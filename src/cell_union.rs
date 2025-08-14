//! S2CellUnion - A region consisting of cells of various sizes
//!
//! This module provides the S2CellUnion implementation, which represents a region
//! consisting of cells of various sizes. It's used to approximate shapes with
//! configurable tradeoffs between accuracy and cell count.
//!
//! ## Architecture
//!
//! S2CellUnion is represented as a vector of sorted, non-overlapping S2CellIds.
//! The vector is normalized by default, meaning groups of 4 child cells are
//! replaced by their parent cell whenever possible.
//!
//! ## Key Operations
//!
//! - **Construction**: From cell lists, ranges, or other unions
//! - **Containment**: Fast O(log n) point and cell containment testing
//! - **Set Operations**: Union, intersection, and difference with other unions
//! - **Expansion**: Add buffer cells around the boundary
//! - **Normalization**: Optimize representation by merging child cells

use crate::error::{S2Error, S2Result};
use crate::{S2CellId, S2Point, S1Angle, S2Cell, S2Cap, S2LatLngRect};
use crate::cell_id::MAX_LEVEL;
use std::cmp::{max, min};
use std::fmt;

/// Maximum number of cells allowed in decode operations (matching C++ default)
const DECODE_MAX_NUM_CELLS: usize = 1_000_000;

/// Current version number for lossless encoding
const ENCODING_VERSION: u8 = 1;

/// S2CellUnion represents a region consisting of cells of various sizes
///
/// ## Design Principles
///
/// - **Normalized representation**: Groups of 4 child cells replaced by parent when possible
/// - **Sorted ordering**: Cell IDs stored in increasing order along Hilbert curve
/// - **Fast operations**: Containment and intersection use binary search
/// - **Memory efficient**: Optional normalization and packing operations
///
/// ## Examples
///
/// ```rust,ignore
/// use s2geometry_rust::{S2CellUnion, S2CellId};
///
/// // Create from cell IDs (automatically normalized)
/// let cells = vec![S2CellId::from_face(0), S2CellId::from_face(1)];
/// let union = S2CellUnion::new(cells);
///
/// // Test containment
/// let point = S2Point::new(1.0, 0.0, 0.0)?;
/// assert!(union.contains_point(&point));
///
/// // Set operations
/// let other = S2CellUnion::new(vec![S2CellId::from_face(2)]);
/// let combined = union.union(&other);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct S2CellUnion {
    /// Vector of sorted, non-overlapping S2CellIds
    cell_ids: Vec<S2CellId>,
}

impl S2CellUnion {
    /// Creates a new S2CellUnion from the given cell IDs
    ///
    /// The input is automatically normalized: sorted, deduplicated, and
    /// groups of 4 child cells are replaced by their parent.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// let cells = vec![S2CellId::from_face(0), S2CellId::from_face(1)];
    /// let union = S2CellUnion::new(cells);
    /// ```
    pub fn new(mut cell_ids: Vec<S2CellId>) -> Self {
        Self::normalize_static(&mut cell_ids);
        Self { cell_ids }
    }

    /// Creates an empty S2CellUnion
    pub fn empty() -> Self {
        Self {
            cell_ids: Vec::new(),
        }
    }

    /// Creates a cell union covering the whole sphere
    ///
    /// Returns a union containing all 6 face cells (level 0).
    pub fn whole_sphere() -> Self {
        Self::new(vec![
            S2CellId::from_face(0),
            S2CellId::from_face(1),
            S2CellId::from_face(2),
            S2CellId::from_face(3),
            S2CellId::from_face(4),
            S2CellId::from_face(5),
        ])
    }

    /// Creates a cell union from normalized cell IDs without validation
    ///
    /// # Safety
    ///
    /// The caller must ensure that `cell_ids` satisfies the requirements of
    /// `is_normalized()`. Use this only when you know the input is already
    /// normalized (e.g., from another S2CellUnion).
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// let normalized_cells = existing_union.release();
    /// let union = S2CellUnion::from_normalized(normalized_cells);
    /// ```
    pub fn from_normalized(cell_ids: Vec<S2CellId>) -> Self {
        debug_assert!(Self::is_valid_normalized(&cell_ids));
        Self { cell_ids }
    }

    /// Creates a cell union from cell IDs without normalization
    ///
    /// The input must be valid (sorted, non-overlapping) but doesn't need
    /// to be normalized (groups of 4 children can be present instead of parent).
    ///
    /// # Safety
    ///
    /// The caller must ensure that `cell_ids` satisfies the requirements of
    /// `is_valid()` but not necessarily `is_normalized()`.
    pub fn from_verbatim(cell_ids: Vec<S2CellId>) -> Self {
        debug_assert!(Self::is_valid_static(&cell_ids));
        Self { cell_ids }
    }

    /// Creates a cell union covering a continuous range of leaf cells
    ///
    /// The output covers all leaf cells from `min_id` to `max_id` inclusive.
    /// Both IDs must be leaf cells, and `min_id <= max_id`.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// let min_leaf = S2CellId::from_face_pos_level(0, 0, MAX_LEVEL)?;
    /// let max_leaf = min_leaf.next();
    /// let union = S2CellUnion::from_min_max(min_leaf, max_leaf)?;
    /// ```
    pub fn from_min_max(min_id: S2CellId, max_id: S2CellId) -> S2Result<Self> {
        if !min_id.is_leaf() || !max_id.is_leaf() {
            return Err(S2Error::invalid_cell_id(min_id.id(), "Both min_id and max_id must be leaf cells"));
        }
        if min_id > max_id {
            return Err(S2Error::invalid_cell_id(min_id.id(), "min_id must be <= max_id"));
        }
        
        Self::from_begin_end(min_id, max_id.next())
    }

    /// Creates a cell union covering a range from begin (inclusive) to end (exclusive)
    ///
    /// Like `from_min_max()` but uses half-open interval semantics.
    /// If `begin == end`, returns an empty union.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// let begin = S2CellId::begin(MAX_LEVEL);
    /// let end = S2CellId::end(MAX_LEVEL);
    /// let whole_sphere = S2CellUnion::from_begin_end(begin, end)?;
    /// ```
    pub fn from_begin_end(begin: S2CellId, end: S2CellId) -> S2Result<Self> {
        if !begin.is_leaf() || !end.is_leaf() {
            return Err(S2Error::invalid_cell_id(begin.id(), "Both begin and end must be leaf cells"));
        }
        if begin > end {
            return Err(S2Error::invalid_cell_id(begin.id(), "begin must be <= end"));
        }

        let mut cell_ids = Vec::new();
        let mut id = begin;
        
        while id != end {
            let next_id = id.maximum_tile(end);
            cell_ids.push(next_id);
            id = next_id.next();
        }
        
        // The output is already normalized by construction
        Ok(Self::from_normalized(cell_ids))
    }

    /// Returns the number of cells in the union
    pub fn num_cells(&self) -> usize {
        self.cell_ids.len()
    }

    /// Returns true if the union is empty
    pub fn is_empty(&self) -> bool {
        self.cell_ids.is_empty()
    }

    /// Returns the cell ID at the given index
    ///
    /// # Panics
    ///
    /// Panics if `index >= num_cells()`.
    pub fn cell_id(&self, index: usize) -> S2CellId {
        self.cell_ids[index]
    }

    /// Returns a reference to the underlying vector of cell IDs
    pub fn cell_ids(&self) -> &[S2CellId] {
        &self.cell_ids
    }

    /// Gives ownership of the cell IDs vector and clears the union
    pub fn release(&mut self) -> Vec<S2CellId> {
        std::mem::take(&mut self.cell_ids)
    }

    /// Clears the union and minimizes memory usage
    pub fn clear(&mut self) {
        self.cell_ids.clear();
        self.cell_ids.shrink_to_fit();
    }

    /// Reduces memory usage if there's excess capacity
    ///
    /// If there are more than `excess` unused elements in the vector,
    /// reallocates to eliminate the excess space.
    pub fn pack(&mut self, excess: usize) {
        if self.cell_ids.capacity() - self.cell_ids.len() > excess {
            self.cell_ids.shrink_to_fit();
        }
    }

    /// Returns true if the cell union is valid
    ///
    /// Valid means the S2CellIds are valid, non-overlapping, and sorted.
    pub fn is_valid(&self) -> bool {
        Self::is_valid_static(&self.cell_ids)
    }

    /// Returns true if the cell union is normalized
    ///
    /// Normalized means valid plus no four cells have a common parent.
    pub fn is_normalized(&self) -> bool {
        Self::is_valid_normalized(&self.cell_ids)
    }

    /// Normalizes the cell union in place
    ///
    /// Discards cells contained by other cells, replaces groups of 4 child
    /// cells by their parent, and sorts in increasing order.
    pub fn normalize(&mut self) {
        Self::normalize_static(&mut self.cell_ids);
    }

    /// Returns true if the union contains the given cell ID
    ///
    /// This is a fast O(log n) operation using binary search.
    pub fn contains_cell_id(&self, id: S2CellId) -> bool {
        if !id.is_valid() {
            return false;
        }

        // Binary search for the first cell that might contain id
        let pos = self.cell_ids.binary_search_by(|cell| {
            if cell.range_max() < id.range_min() {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });
        
        match pos {
            Ok(_) => unreachable!(), // Equal case handled by Greater
            Err(index) => {
                index < self.cell_ids.len() && self.cell_ids[index].contains(&id)
            }
        }
    }

    /// Returns true if the union contains the given point
    ///
    /// This is a fast O(log n) operation.
    pub fn contains_point(&self, point: &S2Point) -> bool {
        self.contains_cell_id(S2CellId::from_point(point))
    }

    /// Returns true if the union contains the given cell
    pub fn contains_cell(&self, cell: &S2Cell) -> bool {
        self.contains_cell_id(cell.id())
    }

    /// Returns true if the union intersects the given cell ID
    ///
    /// This is a fast O(log n) operation using binary search.
    pub fn intersects_cell_id(&self, id: S2CellId) -> bool {
        if !id.is_valid() {
            return false;
        }

        // Binary search for the first cell that might intersect id
        let pos = self.cell_ids.binary_search_by(|cell| {
            if cell.range_max() < id.range_min() {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });
        
        match pos {
            Ok(_) => unreachable!(), // Equal case handled by Greater
            Err(index) => {
                index < self.cell_ids.len() && self.cell_ids[index].intersects(&id)
            }
        }
    }

    /// Returns true if the union intersects the given cell
    pub fn intersects_cell(&self, cell: &S2Cell) -> bool {
        self.intersects_cell_id(cell.id())
    }

    /// Returns true if this union contains the other union
    ///
    /// An empty union is contained by any union (including empty).
    pub fn contains(&self, other: &S2CellUnion) -> bool {
        if other.is_empty() {
            return true;
        }
        if self.is_empty() {
            return false;
        }

        let mut i = 0;
        for &other_id in &other.cell_ids {
            // Advance i to the first cell that might contain other_id
            while i < self.cell_ids.len() && self.cell_ids[i].range_max() < other_id.range_min() {
                i += 1;
            }
            
            if i >= self.cell_ids.len() || !self.cell_ids[i].contains(&other_id) {
                return false;
            }
        }
        true
    }

    /// Returns true if this union intersects the other union
    pub fn intersects(&self, other: &S2CellUnion) -> bool {
        let mut i = 0;
        let mut j = 0;
        
        while i < self.cell_ids.len() && j < other.cell_ids.len() {
            let self_id = self.cell_ids[i];
            let other_id = other.cell_ids[j];
            
            if self_id.range_max() < other_id.range_min() {
                // self_id is entirely before other_id
                i += 1;
            } else if other_id.range_max() < self_id.range_min() {
                // other_id is entirely before self_id  
                j += 1;
            } else {
                // Neither is entirely before the other, so they intersect
                return true;
            }
        }
        false
    }

    /// Returns the union of this and another cell union
    pub fn union(&self, other: &S2CellUnion) -> Self {
        let mut combined = Vec::with_capacity(self.cell_ids.len() + other.cell_ids.len());
        combined.extend_from_slice(&self.cell_ids);
        combined.extend_from_slice(&other.cell_ids);
        Self::new(combined)
    }

    /// Returns the intersection of this union with a cell ID
    pub fn intersection_with_cell_id(&self, id: S2CellId) -> Self {
        if !id.is_valid() {
            return Self::empty();
        }

        if self.contains_cell_id(id) {
            return Self::from_normalized(vec![id]);
        }

        let mut result = Vec::new();
        let id_range_min = id.range_min();
        let id_range_max = id.range_max();
        
        for &cell_id in &self.cell_ids {
            if cell_id.id() >= id_range_min && cell_id.id() <= id_range_max {
                result.push(cell_id);
            } else if cell_id.id() > id_range_max {
                break;
            }
        }
        
        Self::from_normalized(result)
    }

    /// Returns the intersection of this and another cell union
    pub fn intersection(&self, other: &S2CellUnion) -> Self {
        let mut result = Vec::new();
        Self::get_intersection(&self.cell_ids, &other.cell_ids, &mut result);
        Self::from_normalized(result)
    }

    /// Returns the difference between this and another cell union
    pub fn difference(&self, other: &S2CellUnion) -> Self {
        let mut result = Vec::new();
        for &id in &self.cell_ids {
            Self::get_difference_internal(id, other, &mut result);
        }
        Self::from_normalized(result)
    }

    /// Expands the union by adding neighboring cells at the given level
    ///
    /// For each cell in the union, adds all neighboring cells at `expand_level`.
    /// If a cell is smaller than `expand_level`, first promotes it to `expand_level`.
    pub fn expand(&mut self, expand_level: i32) {
        let mut output = Vec::new();
        let level_lsb = S2CellId::lsb_for_level(expand_level);
        
        // Process cells in reverse order to optimize for small regions
        for &id in self.cell_ids.iter().rev() {
            let mut expand_id = id;
            if id.lsb() < level_lsb {
                expand_id = id.parent_at_level(expand_level);
            }
            
            output.push(expand_id);
            expand_id.append_all_neighbors(expand_level, &mut output);
        }
        
        self.cell_ids = output;
        self.normalize();
    }

    /// Expands the union with minimum radius and level constraints
    ///
    /// Expands so all points within `min_radius` are covered, but doesn't
    /// use cells more than `max_level_diff` levels higher than the largest
    /// input cell.
    pub fn expand_with_radius(&mut self, min_radius: S1Angle, max_level_diff: i32) {
        let min_level = self.cell_ids.iter()
            .map(|id| id.level())
            .min()
            .unwrap_or(MAX_LEVEL);
            
        // Find the level where cells are at least min_radius wide
        let radius_level = S2CellId::level_for_min_width(min_radius.radians());
        
        if radius_level == 0 && min_radius.radians() > S2CellId::min_width_at_level(0) {
            // Requested expansion is greater than face cell width
            self.expand(0);
        }
        
        let expand_level = min(min_level + max_level_diff, radius_level);
        self.expand(expand_level);
    }

    /// Returns the number of leaf cells covered by this union
    ///
    /// This will be at most 6 * 2^60 for the whole sphere.
    pub fn leaf_cells_covered(&self) -> u64 {
        self.cell_ids.iter()
            .map(|id| {
                let inverted_level = MAX_LEVEL - id.level();
                1u64 << (inverted_level << 1)
            })
            .sum()
    }

    /// Returns the approximate area covered by the union in steradians
    ///
    /// Uses the average area of each cell level, which may be off by up to
    /// a factor of 1.7 due to cell distortion.
    pub fn average_area(&self) -> f64 {
        S2Cell::average_area_at_level(MAX_LEVEL) * self.leaf_cells_covered() as f64
    }

    /// Returns the approximate area using each cell's approximate area
    pub fn approx_area(&self) -> f64 {
        self.cell_ids.iter()
            .map(|&id| S2Cell::from(id).approx_area())
            .sum()
    }

    /// Returns the exact area using exact geometric calculations
    pub fn exact_area(&self) -> f64 {
        self.cell_ids.iter()
            .map(|&id| S2Cell::from(id).exact_area())
            .sum()
    }

    /// Returns a bounding cap for the union
    pub fn get_cap_bound(&self) -> S2Cap {
        if self.is_empty() {
            return S2Cap::empty();
        }

        // Compute approximate centroid weighted by cell area
        let mut centroid = crate::math::DVec3::ZERO;
        for &id in &self.cell_ids {
            let area = S2Cell::average_area_at_level(id.level());
            centroid += area * id.to_point().coords();
        }
        
        if centroid == crate::math::DVec3::ZERO {
            centroid = crate::math::DVec3::X;
        } else {
            centroid = centroid.normalize();
        }

        // Expand cap to contain all cell bounds
        let center_point = S2Point::from_vec3(centroid).unwrap();
        let mut cap = S2Cap::from_point(center_point);
        
        for &id in &self.cell_ids {
            cap.add_cap(&S2Cell::from(id).get_cap_bound());
        }
        
        cap
    }

    /// Returns a bounding rectangle for the union
    pub fn get_rect_bound(&self) -> S2LatLngRect {
        let mut bound = S2LatLngRect::empty();
        for &id in &self.cell_ids {
            bound = bound.union(&S2Cell::from(id).get_rect_bound());
        }
        bound
    }

    /// Iterator over the cell IDs in the union
    pub fn iter(&self) -> impl Iterator<Item = S2CellId> + '_ {
        self.cell_ids.iter().copied()
    }

    // Internal implementation methods

    /// Static method to check if a vector of cell IDs is valid
    fn is_valid_static(cell_ids: &[S2CellId]) -> bool {
        if !cell_ids.is_empty() && !cell_ids[0].is_valid() {
            return false;
        }
        
        for i in 1..cell_ids.len() {
            if !cell_ids[i].is_valid() {
                return false;
            }
            if cell_ids[i - 1].range_max() >= cell_ids[i].range_min() {
                return false;
            }
        }
        true
    }

    /// Static method to check if a vector of cell IDs is normalized
    fn is_valid_normalized(cell_ids: &[S2CellId]) -> bool {
        if !Self::is_valid_static(cell_ids) {
            return false;
        }
        
        for i in 3..cell_ids.len() {
            if Self::are_siblings(cell_ids[i - 3], cell_ids[i - 2], 
                                 cell_ids[i - 1], cell_ids[i]) {
                return false;
            }
        }
        true
    }

    /// Returns true if the four cells are siblings (have a common parent)
    fn are_siblings(a: S2CellId, b: S2CellId, c: S2CellId, d: S2CellId) -> bool {
        // Quick XOR test (necessary but not sufficient)
        if (a.id() ^ b.id() ^ c.id()) != d.id() {
            return false;
        }

        // Exact test using parent mask
        let mask = d.lsb() << 1;
        let mask = !(mask + (mask << 1));
        let d_masked = d.id() & mask;
        
        (a.id() & mask) == d_masked &&
        (b.id() & mask) == d_masked &&
        (c.id() & mask) == d_masked &&
        !d.is_face()
    }

    /// Static method to normalize a vector of cell IDs in place
    fn normalize_static(cell_ids: &mut Vec<S2CellId>) {
        cell_ids.sort();
        let mut out = 0;
        
        for i in 0..cell_ids.len() {
            let mut id = cell_ids[i];
            
            // Skip cells contained by the previous cell
            if out > 0 && cell_ids[out - 1].contains(&id) {
                continue;
            }
            
            // Remove previous cells contained by this cell
            while out > 0 && id.contains(&cell_ids[out - 1]) {
                out -= 1;
            }
            
            // Collapse four siblings into their parent
            while out >= 3 && Self::are_siblings(cell_ids[out - 3], cell_ids[out - 2],
                                                 cell_ids[out - 1], id) {
                id = id.parent(id.level() - 1).unwrap_or(id);
                out -= 3;
            }
            
            cell_ids[out] = id;
            out += 1;
        }
        
        cell_ids.truncate(out);
    }

    /// Helper for computing intersection of two sorted cell ID vectors
    fn get_intersection(x: &[S2CellId], y: &[S2CellId], out: &mut Vec<S2CellId>) {
        out.clear();
        let mut i = 0;
        let mut j = 0;
        
        while i < x.len() && j < y.len() {
            let x_min = x[i].range_min();
            let y_min = y[j].range_min();
            
            if x_min > y_min {
                if x[i].id() <= y[j].range_max() {
                    out.push(x[i]);
                    i += 1;
                } else {
                    j += 1;
                }
            } else if y_min > x_min {
                if y[j].id() <= x[i].range_max() {
                    out.push(y[j]);
                    j += 1;
                } else {
                    i += 1;
                }
            } else {
                // Same range_min, so one contains the other
                if x[i] < y[j] {
                    out.push(x[i]);
                    i += 1;
                } else {
                    out.push(y[j]);
                    j += 1;
                }
            }
        }
    }

    /// Helper for computing cell difference (recursive)
    fn get_difference_internal(cell: S2CellId, y: &S2CellUnion, result: &mut Vec<S2CellId>) {
        if !y.intersects_cell_id(cell) {
            result.push(cell);
        } else if !y.contains_cell_id(cell) {
            // Recursively subtract from children
            for child in cell.children() {
                Self::get_difference_internal(child, y, result);
            }
        }
    }
}

impl Default for S2CellUnion {
    fn default() -> Self {
        Self::empty()
    }
}

impl fmt::Display for S2CellUnion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const MAX_DISPLAY: usize = 500;
        
        write!(f, "Size:{} S2CellIds:", self.num_cells())?;
        
        let limit = min(MAX_DISPLAY, self.num_cells());
        for (i, &cell_id) in self.cell_ids.iter().take(limit).enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{}", cell_id.to_token())?;
        }
        
        if self.num_cells() > MAX_DISPLAY {
            write!(f, ",...")?;
        }
        
        Ok(())
    }
}

impl FromIterator<S2CellId> for S2CellUnion {
    fn from_iter<T: IntoIterator<Item = S2CellId>>(iter: T) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl IntoIterator for S2CellUnion {
    type Item = S2CellId;
    type IntoIter = std::vec::IntoIter<S2CellId>;

    fn into_iter(self) -> Self::IntoIter {
        self.cell_ids.into_iter()
    }
}

impl<'a> IntoIterator for &'a S2CellUnion {
    type Item = S2CellId;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, S2CellId>>;

    fn into_iter(self) -> Self::IntoIter {
        self.cell_ids.iter().copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_union() {
        let union = S2CellUnion::empty();
        assert!(union.is_empty());
        assert_eq!(union.num_cells(), 0);
        assert!(union.is_valid());
        assert!(union.is_normalized());
    }

    #[test]
    fn test_whole_sphere() {
        let union = S2CellUnion::whole_sphere();
        assert_eq!(union.num_cells(), 6);
        assert!(union.is_valid());
        assert!(union.is_normalized());
        
        // Should contain all face cells
        for face in 0..6 {
            assert!(union.contains_cell_id(S2CellId::from_face(face)));
        }
    }

    #[test]
    fn test_single_cell_construction() {
        let face_id = S2CellId::from_face(1);
        let union = S2CellUnion::new(vec![face_id]);
        
        assert_eq!(union.num_cells(), 1);
        assert_eq!(union.cell_id(0), face_id);
        assert!(union.contains_cell_id(face_id));
    }

    #[test]
    fn test_normalization() {
        let parent_id = S2CellId::from_face(0);
        
        // Create union with all 4 children - should normalize to parent
        let children: Vec<S2CellId> = parent_id.children().collect();
        let union = S2CellUnion::new(children);
        
        assert_eq!(union.num_cells(), 1);
        assert_eq!(union.cell_id(0), parent_id);
    }

    #[test]
    fn test_containment() {
        let parent_id = S2CellId::from_face(0);
        let union = S2CellUnion::new(vec![parent_id]);
        
        // Should contain all descendants
        for child in parent_id.children() {
            assert!(union.contains_cell_id(child));
        }
        
        // Should not contain siblings
        let sibling = S2CellId::from_face(1);
        assert!(!union.contains_cell_id(sibling));
    }

    #[test]
    fn test_union_operation() {
        let union1 = S2CellUnion::new(vec![S2CellId::from_face(0)]);
        let union2 = S2CellUnion::new(vec![S2CellId::from_face(1)]);
        
        let combined = union1.union(&union2);
        
        assert_eq!(combined.num_cells(), 2);
        assert!(combined.contains_cell_id(S2CellId::from_face(0)));
        assert!(combined.contains_cell_id(S2CellId::from_face(1)));
    }

    #[test]
    fn test_intersection() {
        let parent = S2CellUnion::new(vec![S2CellId::from_face(0)]);
        let children: Vec<S2CellId> = S2CellId::from_face(0).children().collect();
        let child_union = S2CellUnion::new(children);
        
        let intersection = parent.intersection(&child_union);
        
        // Intersection should equal the child union
        assert_eq!(intersection, child_union);
    }
}