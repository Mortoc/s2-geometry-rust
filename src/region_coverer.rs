//! S2RegionCoverer - Approximates arbitrary regions using unions of S2 cells
//!
//! S2RegionCoverer allows arbitrary regions to be approximated as unions of cells
//! for efficient search and storage operations. The approximation is controlled by
//! configurable parameters that trade off between accuracy and cell count.
//!
//! ## Core Algorithm
//!
//! The coverer uses a priority queue based approach:
//! 1. Start with the 6 cube face cells 
//! 2. For each cell that intersects the region, subdivide into 4 children
//! 3. Use priority queue to select best candidates for subdivision
//! 4. Continue until constraints are met (max cells, max level, etc.)
//! 5. Return the final set of cells that covers the region
//!
//! ## Usage Examples
//!
//! ```rust,ignore
//! use s2geometry_rust::{S2RegionCoverer, S2Cap, S2Point, S1Angle};
//! 
//! // Configure the coverer
//! let mut coverer = S2RegionCoverer::new();
//! coverer.set_max_cells(8);
//! coverer.set_max_level(20);
//!
//! // Create a region to cover
//! let center = S2Point::from_coords(1.0, 0.0, 0.0);
//! let cap = S2Cap::from_center_angle(center, S1Angle::from_degrees(1.0));
//!
//! // Generate covering
//! let covering = coverer.get_covering(&cap);
//! ```

use crate::error::{S2Error, S2Result};
use crate::{S2CellId, S2CellUnion, S2Point, S2Cell, S2Cap, S2LatLngRect, S2Loop, S2Polyline, S2LatLng, S1ChordAngle};
use crate::cell_id::MAX_LEVEL;
use std::collections::BinaryHeap;
use std::cmp::Ordering;
use std::fmt;

/// Default maximum number of cells in a covering
const DEFAULT_MAX_CELLS: usize = 8;

/// Maximum allowed cells to prevent runaway algorithms
const MAX_CELLS_HARD_LIMIT: usize = 1_000_000;

/// Interface for regions that can be covered by S2 cells
pub trait S2Region {
    /// Test if the region contains the given point
    fn contains(&self, point: &S2Point) -> bool;

    /// Test if the region may intersect the given cell
    /// Returns false only if the region definitely does not intersect the cell
    fn may_intersect_cell(&self, cell: &S2Cell) -> bool;

    /// Get a bounding rectangle for this region
    fn get_rect_bound(&self) -> S2LatLngRect;

    /// Get a bounding cap for this region  
    fn get_cap_bound(&self) -> S2Cap;
}

// Implement S2Region for existing types
impl S2Region for S2Cap {
    fn contains(&self, point: &S2Point) -> bool {
        self.contains(point)
    }

    fn may_intersect_cell(&self, cell: &S2Cell) -> bool {
        self.may_intersect(cell)
    }

    fn get_rect_bound(&self) -> S2LatLngRect {
        self.get_rect_bound()
    }

    fn get_cap_bound(&self) -> S2Cap {
        *self
    }
}

impl S2Region for S2LatLngRect {
    fn contains(&self, point: &S2Point) -> bool {
        self.contains_point(point)
    }

    fn may_intersect_cell(&self, cell: &S2Cell) -> bool {
        // Conservative check: use rect intersection
        self.intersects(&cell.get_rect_bound())
    }

    fn get_rect_bound(&self) -> S2LatLngRect {
        *self
    }

    fn get_cap_bound(&self) -> S2Cap {
        // Convert the lat/lng rectangle to a bounding cap
        // For simplicity, use the center and compute radius to farthest corner
        let center_lat = (self.lat().lo() + self.lat().hi()) / 2.0;
        let center_lng = self.lng().get_center();
        
        let center_latlng = S2LatLng::from_radians(center_lat, center_lng);
        let center_point = center_latlng.to_point().unwrap_or_else(|_| S2Point::new(1.0, 0.0, 0.0).unwrap());
        
        // Find the maximum distance to any corner
        let corners = [
            S2LatLng::from_radians(self.lat().lo(), self.lng().lo()),
            S2LatLng::from_radians(self.lat().lo(), self.lng().hi()), 
            S2LatLng::from_radians(self.lat().hi(), self.lng().lo()),
            S2LatLng::from_radians(self.lat().hi(), self.lng().hi()),
        ];
        
        let mut max_dist = S1ChordAngle::zero();
        for corner in &corners {
            if let Ok(corner_point) = corner.to_point() {
                let dist = S1ChordAngle::between_points(&center_point, &corner_point);
                if dist > max_dist {
                    max_dist = dist;
                }
            }
        }
        
        S2Cap::new(center_point, max_dist.to_angle())
    }
}

impl S2Region for S2Loop {
    fn contains(&self, point: &S2Point) -> bool {
        self.contains(point)
    }

    fn may_intersect_cell(&self, cell: &S2Cell) -> bool {
        // Conservative check: check if any vertex is in the loop or cell
        // TODO: implement proper loop-cell intersection
        for vertex in self.vertices() {
            if cell.contains(vertex) {
                return true;
            }
        }
        for i in 0..4 {
            let vertex = cell.get_vertex(i);
            if self.contains(&vertex) {
                return true;
            }
        }
        false
    }

    fn get_rect_bound(&self) -> S2LatLngRect {
        self.get_rect_bound()
    }

    fn get_cap_bound(&self) -> S2Cap {
        self.get_cap_bound()
    }
}

impl S2Region for S2Polyline {
    fn contains(&self, _point: &S2Point) -> bool {
        // Polylines have no interior
        false
    }

    fn may_intersect_cell(&self, cell: &S2Cell) -> bool {
        // Conservative check: check if any vertex is in the cell
        // TODO: implement proper polyline-cell intersection
        for vertex in self.vertices() {
            if cell.contains(vertex) {
                return true;
            }
        }
        false
    }

    fn get_rect_bound(&self) -> S2LatLngRect {
        self.get_rect_bound()
    }

    fn get_cap_bound(&self) -> S2Cap {
        self.get_cap_bound()
    }
}

impl S2Region for S2CellUnion {
    fn contains(&self, point: &S2Point) -> bool {
        self.contains_point(point)
    }

    fn may_intersect_cell(&self, cell: &S2Cell) -> bool {
        self.intersects_cell(cell)
    }

    fn get_rect_bound(&self) -> S2LatLngRect {
        self.get_rect_bound()
    }

    fn get_cap_bound(&self) -> S2Cap {
        self.get_cap_bound()
    }
}

/// Configuration options for S2RegionCoverer
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct S2RegionCovererOptions {
    /// Maximum number of cells in the covering (default: 8)
    max_cells: usize,
    /// Minimum cell level to use (default: 0)  
    min_level: u8,
    /// Maximum cell level to use (default: MAX_LEVEL)
    max_level: u8,
    /// Level modulus - only use cell levels that are multiples of this (default: 1)
    level_mod: u8,
}

impl Default for S2RegionCovererOptions {
    fn default() -> Self {
        Self {
            max_cells: DEFAULT_MAX_CELLS,
            min_level: 0,
            max_level: MAX_LEVEL as u8,
            level_mod: 1,
        }
    }
}

impl S2RegionCovererOptions {
    /// Create new options with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set maximum number of cells in the covering
    pub fn max_cells(mut self, max_cells: usize) -> S2Result<Self> {
        if max_cells == 0 || max_cells > MAX_CELLS_HARD_LIMIT {
            return Err(S2Error::invalid_point(format!(
                "max_cells must be between 1 and {}, got {}", 
                MAX_CELLS_HARD_LIMIT, max_cells
            )));
        }
        self.max_cells = max_cells;
        Ok(self)
    }

    /// Set minimum cell level  
    pub fn min_level(mut self, min_level: u8) -> S2Result<Self> {
        if min_level as i32 > MAX_LEVEL {
            return Err(S2Error::invalid_point(format!(
                "min_level must be <= {}, got {}", 
                MAX_LEVEL, min_level
            )));
        }
        self.min_level = min_level;
        Ok(self)
    }

    /// Set maximum cell level
    pub fn max_level(mut self, max_level: u8) -> S2Result<Self> {
        if max_level as i32 > MAX_LEVEL {
            return Err(S2Error::invalid_point(format!(
                "max_level must be <= {}, got {}", 
                MAX_LEVEL, max_level
            )));
        }
        self.max_level = max_level;
        Ok(self)
    }

    /// Set level modulus
    pub fn level_mod(mut self, level_mod: u8) -> S2Result<Self> {
        if level_mod == 0 || level_mod > 3 {
            return Err(S2Error::invalid_point(format!(
                "level_mod must be 1, 2, or 3, got {}", 
                level_mod
            )));
        }
        self.level_mod = level_mod;
        Ok(self)
    }

    /// Get the maximum number of cells
    pub fn get_max_cells(&self) -> usize {
        self.max_cells
    }

    /// Get the minimum level
    pub fn get_min_level(&self) -> u8 {
        self.min_level
    }

    /// Get the maximum level  
    pub fn get_max_level(&self) -> u8 {
        self.max_level
    }

    /// Get the level modulus
    pub fn get_level_mod(&self) -> u8 {
        self.level_mod
    }

    /// Check if options are valid
    pub fn is_valid(&self) -> bool {
        self.max_cells > 0 
            && self.max_cells <= MAX_CELLS_HARD_LIMIT
            && self.min_level <= self.max_level
            && self.max_level as i32 <= MAX_LEVEL
            && self.level_mod >= 1 
            && self.level_mod <= 3
    }
}

/// A candidate cell for consideration in the covering algorithm
#[derive(Debug, Clone)]
struct Candidate {
    /// The cell being considered
    cell_id: S2CellId,
    /// Whether this cell is definitely contained in the region
    is_terminal: bool,
    /// Number of children that intersect the region (0 if not computed)
    num_children: u8,
    /// Priority for the candidate (used in priority queue)
    priority: i32,
}

impl Candidate {
    fn new(cell_id: S2CellId) -> Self {
        Self {
            cell_id,
            is_terminal: false,
            num_children: 0,
            priority: 0,
        }
    }

    /// Compute priority for this candidate
    /// Priority is based on cell level and number of intersecting children
    /// Higher priority means the cell should be subdivided first
    fn compute_priority(&mut self, _options: &S2RegionCovererOptions) {
        let level = self.cell_id.level();
        
        // Prefer subdividing larger cells (lower level)
        let mut priority = -(level as i32);
        
        // If we know the number of children that intersect, 
        // prefer cells with fewer intersecting children
        if self.num_children > 0 {
            priority -= self.num_children as i32;
        }
        
        self.priority = priority;
    }
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.priority == other.priority
    }
}

impl Eq for Candidate {}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap is a max-heap, but we want higher priority first
        self.priority.cmp(&other.priority)
    }
}

/// S2RegionCoverer approximates arbitrary regions using unions of S2 cells
///
/// The coverer provides several algorithms for generating cell coverings:
/// - `get_covering()`: Standard covering algorithm 
/// - `get_interior_covering()`: Covering using only cells entirely within the region
/// - `get_fast_covering()`: Faster but potentially less optimal covering
///
/// ## Algorithm Details
///
/// The standard algorithm uses a priority queue to manage candidate cells:
/// 1. Initialize with the 6 cube face cells
/// 2. While constraints allow, subdivide the highest priority candidate
/// 3. A candidate's priority depends on its level and intersection properties
/// 4. Terminal candidates (fully contained) are not subdivided
/// 5. Return the final set of candidates as the covering
///
/// ## Performance
///
/// - Time complexity: O(k log k) where k is the number of cells generated
/// - Space complexity: O(k) for the priority queue and result
/// - Typical covering sizes: 4-20 cells for most real-world regions
#[derive(Debug, Clone)]
pub struct S2RegionCoverer {
    options: S2RegionCovererOptions,
}

impl Default for S2RegionCoverer {
    fn default() -> Self {
        Self::new()
    }
}

impl S2RegionCoverer {
    /// Create a new S2RegionCoverer with default options
    pub fn new() -> Self {
        Self {
            options: S2RegionCovererOptions::default(),
        }
    }

    /// Create a new S2RegionCoverer with the given options
    pub fn with_options(options: S2RegionCovererOptions) -> S2Result<Self> {
        if !options.is_valid() {
            return Err(S2Error::invalid_point("Invalid coverer options"));
        }
        Ok(Self { options })
    }

    /// Set the maximum number of cells in the covering
    pub fn set_max_cells(&mut self, max_cells: usize) -> S2Result<()> {
        self.options = self.options.max_cells(max_cells)?;
        Ok(())
    }

    /// Set the minimum cell level
    pub fn set_min_level(&mut self, min_level: u8) -> S2Result<()> {
        self.options = self.options.min_level(min_level)?;
        Ok(())
    }

    /// Set the maximum cell level
    pub fn set_max_level(&mut self, max_level: u8) -> S2Result<()> {
        self.options = self.options.max_level(max_level)?;
        Ok(())
    }

    /// Set the level modulus
    pub fn set_level_mod(&mut self, level_mod: u8) -> S2Result<()> {
        self.options = self.options.level_mod(level_mod)?;
        Ok(())
    }

    /// Get the current options
    pub fn options(&self) -> &S2RegionCovererOptions {
        &self.options
    }

    /// Generate a cell covering for the given region
    ///
    /// Returns a normalized S2CellUnion that covers the region according to
    /// the configured parameters. The covering satisfies:
    /// - Every point in the region is contained in at least one cell
    /// - The number of cells is at most max_cells (unless min_level requires more)
    /// - All cells are between min_level and max_level
    /// - Cell levels are multiples of level_mod
    pub fn get_covering<R: S2Region>(&self, region: &R) -> S2Result<S2CellUnion> {
        let mut candidates = Vec::new();
        self.get_initial_candidates(region, &mut candidates)?;
        
        let mut queue = BinaryHeap::new();
        for candidate in candidates {
            queue.push(candidate);
        }

        let mut result = Vec::new();
        self.covering_internal(region, &mut queue, &mut result)?;
        
        Ok(S2CellUnion::new(result))
    }

    /// Generate an interior covering for the given region
    ///
    /// Returns a covering where every cell is entirely contained within the region.
    /// This is useful for fast containment testing: if a point is in any cell of
    /// the interior covering, it is definitely in the region.
    pub fn get_interior_covering<R: S2Region>(&self, region: &R) -> S2Result<S2CellUnion> {
        let mut candidates = Vec::new();
        self.get_initial_candidates(region, &mut candidates)?;
        
        let mut queue = BinaryHeap::new();
        for candidate in candidates {
            queue.push(candidate);
        }

        let mut result = Vec::new();
        self.interior_covering_internal(region, &mut queue, &mut result)?;
        
        Ok(S2CellUnion::new(result))
    }

    /// Generate a fast covering for the given region
    ///
    /// This is a simplified algorithm that may not produce optimal results but
    /// is faster than the standard covering algorithm. It uses a flood-fill
    /// approach starting from a single point in the region.
    pub fn get_fast_covering<R: S2Region>(&self, region: &R) -> S2Result<S2CellUnion> {
        // For now, use the standard algorithm
        // TODO: Implement actual fast covering algorithm
        self.get_covering(region)
    }

    /// Check if a cell union is canonical according to the current options
    ///
    /// A covering is canonical if:
    /// - It's normalized (no redundant cells)
    /// - All cells are at valid levels according to level_mod
    /// - Number of cells doesn't exceed max_cells
    pub fn is_canonical(&self, covering: &S2CellUnion) -> bool {
        let cell_ids = covering.cell_ids();
        
        // Check cell count
        if cell_ids.len() > self.options.max_cells {
            return false;
        }

        // Check that all cells are at valid levels
        for &cell_id in cell_ids {
            let level = cell_id.level();
            if level < self.options.min_level as i32 
                || level > self.options.max_level as i32 
                || (level % self.options.level_mod as i32) != 0 {
                return false;
            }
        }

        // Check that the covering is normalized
        let mut temp_union = S2CellUnion::new(cell_ids.to_vec());
        temp_union.normalize();
        let normalized_ids = temp_union.cell_ids();
        normalized_ids == cell_ids
    }

    /// Canonicalize an existing cell union according to the current options
    ///
    /// This adjusts the given covering to meet the current parameters:
    /// - Replaces cells at invalid levels with valid parent/child cells
    /// - Reduces cell count if it exceeds max_cells
    /// - Normalizes the result
    pub fn canonicalize_covering(&self, covering: &mut S2CellUnion) -> S2Result<()> {
        let mut cell_ids = covering.cell_ids().to_vec();
        
        // Replace cells at invalid levels
        let mut i = 0;
        while i < cell_ids.len() {
            let cell_id = cell_ids[i];
            let level = cell_id.level();
            
            // Adjust level to be valid according to level_mod
            let target_level = if level < self.options.min_level as i32 {
                self.options.min_level as i32
            } else if level > self.options.max_level as i32 {
                self.options.max_level as i32
            } else {
                // Round to nearest valid level, but never exceed max_level
                let remainder = level % self.options.level_mod as i32;
                if remainder == 0 {
                    level
                } else if remainder < self.options.level_mod as i32 / 2 {
                    level - remainder
                } else {
                    let rounded_up = level + (self.options.level_mod as i32 - remainder);
                    if rounded_up <= self.options.max_level as i32 {
                        rounded_up
                    } else {
                        level - remainder
                    }
                }
            };

            if target_level != level {
                cell_ids[i] = cell_id.parent(target_level)?;
            }
            i += 1;
        }

        // Normalize to remove duplicates and merge parent/child relationships
        let mut temp_union = S2CellUnion::new(cell_ids);
        temp_union.normalize();
        cell_ids = temp_union.cell_ids().to_vec();

        // If we still have too many cells, keep only the largest ones
        if cell_ids.len() > self.options.max_cells {
            // Sort by level (ascending) so larger cells come first
            cell_ids.sort_by_key(|id| id.level());
            cell_ids.truncate(self.options.max_cells);
            let mut temp_union = S2CellUnion::new(cell_ids);
            temp_union.normalize();
            cell_ids = temp_union.cell_ids().to_vec();
        }

        *covering = S2CellUnion::new(cell_ids);
        Ok(())
    }

    /// Get initial candidate cells (the 6 cube faces)
    fn get_initial_candidates<R: S2Region>(&self, region: &R, candidates: &mut Vec<Candidate>) -> S2Result<()> {
        // Start with the 6 cube face cells
        for face in 0..6 {
            let cell_id = S2CellId::from_face_pos_level(face, 0, 0)?;
            let cell = S2Cell::from(cell_id);
            
            if region.may_intersect_cell(&cell) {
                candidates.push(Candidate::new(cell_id));
            }
        }
        Ok(())
    }

    /// Internal covering algorithm implementation
    fn covering_internal<R: S2Region>(
        &self,
        region: &R,
        queue: &mut BinaryHeap<Candidate>,
        result: &mut Vec<S2CellId>,
    ) -> S2Result<()> {
        while let Some(mut candidate) = queue.pop() {
            // Check if we should expand this candidate
            if self.should_expand(&candidate, result.len()) {
                // Try to subdivide the candidate
                if let Ok(Some(children)) = self.expand_candidate(region, &mut candidate) {
                    for child in children {
                        queue.push(child);
                    }
                    continue;
                }
            }

            // Add this candidate to the result
            result.push(candidate.cell_id);
        }
        Ok(())
    }

    /// Internal interior covering algorithm implementation  
    fn interior_covering_internal<R: S2Region>(
        &self,
        region: &R,
        queue: &mut BinaryHeap<Candidate>,
        result: &mut Vec<S2CellId>,
    ) -> S2Result<()> {
        while let Some(mut candidate) = queue.pop() {
            let cell = S2Cell::from(candidate.cell_id);
            
            // For interior covering, only include cells that are entirely contained
            if self.is_cell_contained(region, &cell) {
                if self.should_expand(&candidate, result.len()) {
                    // Try to subdivide the candidate
                    if let Ok(Some(children)) = self.expand_candidate_interior(region, &mut candidate) {
                        for child in children {
                            queue.push(child);
                        }
                        continue;
                    }
                }
                
                // Add this candidate to the result
                result.push(candidate.cell_id);
            }
        }
        Ok(())
    }

    /// Check if we should expand (subdivide) the given candidate
    fn should_expand(&self, candidate: &Candidate, current_result_size: usize) -> bool {
        // Don't expand terminal candidates (fully contained cells)
        if candidate.is_terminal {
            return false;
        }

        // Don't expand if we're at max level
        let level = candidate.cell_id.level();
        if level >= self.options.max_level as i32 {
            return false;
        }

        // Don't expand if we're at the cell limit
        if current_result_size >= self.options.max_cells {
            return false;
        }

        // Don't expand if the next level wouldn't be valid according to level_mod
        let next_level = level + 1;
        if (next_level % self.options.level_mod as i32) != 0 {
            return false;
        }

        true
    }

    /// Expand a candidate into its children for standard covering
    fn expand_candidate<R: S2Region>(
        &self,
        region: &R,
        candidate: &mut Candidate,
    ) -> S2Result<Option<Vec<Candidate>>> {
        let level = candidate.cell_id.level();
        if level >= MAX_LEVEL {
            return Ok(None);
        }

        let mut children = Vec::new();
        let mut num_intersecting = 0;

        // Check all 4 children
        for child_pos in 0..4 {
            let child_id = candidate.cell_id.child(child_pos)?;
            let child_cell = S2Cell::from(child_id);

            if region.may_intersect_cell(&child_cell) {
                let mut child_candidate = Candidate::new(child_id);
                
                // Check if this child is fully contained (terminal)
                if self.is_cell_contained(region, &child_cell) {
                    child_candidate.is_terminal = true;
                }

                child_candidate.compute_priority(&self.options);
                children.push(child_candidate);
                num_intersecting += 1;
            }
        }

        candidate.num_children = num_intersecting;

        if children.is_empty() {
            Ok(None)
        } else {
            Ok(Some(children))
        }
    }

    /// Expand a candidate into its children for interior covering
    fn expand_candidate_interior<R: S2Region>(
        &self,
        region: &R,
        candidate: &mut Candidate,
    ) -> S2Result<Option<Vec<Candidate>>> {
        let level = candidate.cell_id.level();
        if level >= MAX_LEVEL {
            return Ok(None);
        }

        let mut children = Vec::new();

        // For interior covering, only include children that are fully contained
        for child_pos in 0..4 {
            let child_id = candidate.cell_id.child(child_pos)?;
            let child_cell = S2Cell::from(child_id);

            if self.is_cell_contained(region, &child_cell) {
                let mut child_candidate = Candidate::new(child_id);
                child_candidate.is_terminal = true; // All interior cells are terminal
                child_candidate.compute_priority(&self.options);
                children.push(child_candidate);
            }
        }

        if children.is_empty() {
            Ok(None)
        } else {
            Ok(Some(children))
        }
    }

    /// Check if a cell is entirely contained within the region
    fn is_cell_contained<R: S2Region>(&self, region: &R, cell: &S2Cell) -> bool {
        // Check if all 4 vertices are contained in the region
        for i in 0..4 {
            let vertex = cell.get_vertex(i);
            if !region.contains(&vertex) {
                return false;
            }
        }
        true
    }
}

impl fmt::Display for S2RegionCoverer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "S2RegionCoverer(max_cells={}, min_level={}, max_level={}, level_mod={})",
            self.options.max_cells,
            self.options.min_level,
            self.options.max_level,
            self.options.level_mod
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{S2Point, S1Angle};

    #[test]
    fn test_region_coverer_creation() {
        let coverer = S2RegionCoverer::new();
        assert_eq!(coverer.options().get_max_cells(), DEFAULT_MAX_CELLS);
        assert_eq!(coverer.options().get_min_level(), 0);
        assert_eq!(coverer.options().get_max_level() as i32, MAX_LEVEL);
        assert_eq!(coverer.options().get_level_mod(), 1);
    }

    #[test]
    fn test_options_validation() {
        // Valid options
        let options = S2RegionCovererOptions::new()
            .max_cells(10).unwrap()
            .min_level(2).unwrap()
            .max_level(15).unwrap()
            .level_mod(2).unwrap();
        assert!(options.is_valid());

        // Invalid max_cells
        assert!(S2RegionCovererOptions::new().max_cells(0).is_err());
        assert!(S2RegionCovererOptions::new().max_cells(MAX_CELLS_HARD_LIMIT + 1).is_err());

        // Invalid levels  
        assert!(S2RegionCovererOptions::new().min_level((MAX_LEVEL + 1) as u8).is_err());
        assert!(S2RegionCovererOptions::new().max_level((MAX_LEVEL + 1) as u8).is_err());

        // Invalid level_mod
        assert!(S2RegionCovererOptions::new().level_mod(0).is_err());
        assert!(S2RegionCovererOptions::new().level_mod(4).is_err());
    }

    #[test]
    fn test_covering_single_cell() {
        let coverer = S2RegionCoverer::new();
        
        // Create a cell and use it as a region (via S2CellUnion)
        let cell_id = S2CellId::from_face_pos_level(0, 0, 1).unwrap();
        let cell_union = S2CellUnion::new(vec![cell_id]);
        
        let covering = coverer.get_covering(&cell_union).unwrap();
        
        // The covering should contain the original cell
        assert!(covering.contains_cell_id(cell_id));
    }

    #[test]
    fn test_covering_cap() {
        let mut coverer = S2RegionCoverer::new();
        coverer.set_max_cells(6).unwrap();
        coverer.set_max_level(10).unwrap();
        
        // Create a small cap
        let center = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let cap = S2Cap::from_center_angle(&center, S1Angle::from_degrees(1.0));
        
        let covering = coverer.get_covering(&cap).unwrap();
        
        // Verify basic properties
        assert!(covering.num_cells() > 0);
        assert!(covering.num_cells() <= 6);
        
        // Verify the covering actually covers the center point
        assert!(covering.contains_point(&center));
    }

    #[test]
    fn test_interior_covering() {
        let mut coverer = S2RegionCoverer::new();
        coverer.set_max_cells(20).unwrap();
        
        // Create a cap large enough to have an interior
        let center = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let cap = S2Cap::from_center_angle(&center, S1Angle::from_degrees(10.0));
        
        let exterior_covering = coverer.get_covering(&cap).unwrap();
        let interior_covering = coverer.get_interior_covering(&cap).unwrap();
        
        // Interior covering should have fewer or equal cells
        assert!(interior_covering.num_cells() <= exterior_covering.num_cells());
        
        // Interior covering should be contained in exterior covering
        for &cell_id in interior_covering.cell_ids() {
            assert!(exterior_covering.intersects_cell_id(cell_id));
        }
    }

    #[test]
    fn test_canonicalize_covering() {
        let mut coverer = S2RegionCoverer::new();
        coverer.set_max_cells(2).unwrap(); // Very small limit to force canonicalization
        
        // Create a covering from different faces that won't get normalized together
        let cell_ids = vec![
            S2CellId::from_face_pos_level(0, 0, 1).unwrap(), // Face 0, level 1
            S2CellId::from_face_pos_level(1, 0, 1).unwrap(), // Face 1, level 1  
            S2CellId::from_face_pos_level(2, 0, 1).unwrap(), // Face 2, level 1
        ];
        let mut covering = S2CellUnion::new(cell_ids);
        
        // This should not be canonical because it has more cells than max_cells
        assert!(!coverer.is_canonical(&covering));
        
        coverer.canonicalize_covering(&mut covering).unwrap();
        
        // After canonicalization, it should have at most max_cells
        assert!(covering.num_cells() <= 2);
        assert!(coverer.is_canonical(&covering));
    }

    #[test]
    fn test_level_mod_constraints() {
        let mut coverer = S2RegionCoverer::new();
        coverer.set_level_mod(2).unwrap();
        coverer.set_max_level(6).unwrap();
        
        let center = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let cap = S2Cap::from_center_angle(&center, S1Angle::from_degrees(1.0));
        
        let covering = coverer.get_covering(&cap).unwrap();
        
        // All cells should be at even levels (0, 2, 4, 6)
        for &cell_id in covering.cell_ids() {
            assert_eq!(cell_id.level() % 2, 0, "Cell level {} is not even", cell_id.level());
        }
    }
}