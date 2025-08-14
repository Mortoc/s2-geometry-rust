//! S2ShapeIndex - Spatial index for fast geometric queries
//!
//! This module provides the core S2ShapeIndex trait and related types for
//! indexing polygonal geometry to enable fast spatial queries like intersection
//! testing, containment queries, and distance calculations.

use crate::shape::S2Shape;
use crate::cell_id::S2CellId;
use std::sync::Arc;
use std::fmt;

/// A clipped shape represents the portion of a shape that intersects a specific S2Cell
///
/// This is used internally by the index to store shape fragments that have been
/// clipped to specific cells in the spatial hierarchy.
#[derive(Debug, Clone)]
pub struct S2ClippedShape {
    /// ID of the shape this fragment belongs to
    shape_id: i32,
    /// Whether the center of the containing cell is inside this shape
    /// Only meaningful for polygon shapes
    contains_center: bool,
    /// Edge IDs from the original shape that intersect the cell
    /// Uses inline storage for up to 2 edges for memory efficiency
    edge_ids: EdgeIdStorage,
}

/// Storage for edge IDs with inline optimization for small counts
#[derive(Debug, Clone)]
enum EdgeIdStorage {
    /// No edges
    Empty,
    /// Single edge (common case)
    One(i32),
    /// Two edges (also common)
    Two(i32, i32),
    /// More than two edges
    Many(Vec<i32>),
}

impl EdgeIdStorage {
    fn new() -> Self {
        EdgeIdStorage::Empty
    }

    fn len(&self) -> usize {
        match self {
            EdgeIdStorage::Empty => 0,
            EdgeIdStorage::One(_) => 1,
            EdgeIdStorage::Two(_, _) => 2,
            EdgeIdStorage::Many(vec) => vec.len(),
        }
    }

    fn is_empty(&self) -> bool {
        matches!(self, EdgeIdStorage::Empty)
    }

    fn push(&mut self, edge_id: i32) {
        match self {
            EdgeIdStorage::Empty => {
                *self = EdgeIdStorage::One(edge_id);
            }
            EdgeIdStorage::One(first) => {
                *self = EdgeIdStorage::Two(*first, edge_id);
            }
            EdgeIdStorage::Two(first, second) => {
                *self = EdgeIdStorage::Many(vec![*first, *second, edge_id]);
            }
            EdgeIdStorage::Many(vec) => {
                vec.push(edge_id);
            }
        }
    }

    fn get(&self, index: usize) -> Option<i32> {
        match self {
            EdgeIdStorage::Empty => None,
            EdgeIdStorage::One(id) => if index == 0 { Some(*id) } else { None },
            EdgeIdStorage::Two(first, second) => match index {
                0 => Some(*first),
                1 => Some(*second),
                _ => None,
            },
            EdgeIdStorage::Many(vec) => vec.get(index).copied(),
        }
    }

    fn iter(&self) -> EdgeIdIter {
        EdgeIdIter::new(self)
    }
}

/// Iterator over edge IDs in an EdgeIdStorage
pub struct EdgeIdIter<'a> {
    storage: &'a EdgeIdStorage,
    index: usize,
}

impl<'a> EdgeIdIter<'a> {
    fn new(storage: &'a EdgeIdStorage) -> Self {
        Self { storage, index: 0 }
    }
}

impl<'a> Iterator for EdgeIdIter<'a> {
    type Item = i32;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.storage.get(self.index);
        if result.is_some() {
            self.index += 1;
        }
        result
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.storage.len().saturating_sub(self.index);
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for EdgeIdIter<'a> {}

impl S2ClippedShape {
    /// Create a new clipped shape
    pub fn new(shape_id: i32, contains_center: bool) -> Self {
        Self {
            shape_id,
            contains_center,
            edge_ids: EdgeIdStorage::new(),
        }
    }

    /// Get the shape ID
    pub fn shape_id(&self) -> i32 {
        self.shape_id
    }

    /// Check if the cell center is contained in this shape
    pub fn contains_center(&self) -> bool {
        self.contains_center
    }

    /// Get the number of edges
    pub fn num_edges(&self) -> i32 {
        self.edge_ids.len() as i32
    }

    /// Add an edge ID to this clipped shape
    pub fn add_edge(&mut self, edge_id: i32) {
        self.edge_ids.push(edge_id);
    }

    /// Get an edge ID by index
    pub fn edge_id(&self, index: i32) -> Option<i32> {
        self.edge_ids.get(index as usize)
    }

    /// Iterate over all edge IDs
    pub fn edge_ids(&self) -> EdgeIdIter {
        self.edge_ids.iter()
    }

    /// Check if this clipped shape has any edges
    pub fn is_empty(&self) -> bool {
        self.edge_ids.is_empty()
    }
}

/// Contents of an index cell - all the clipped shapes that intersect this cell
#[derive(Debug, Clone, Default)]
pub struct S2ShapeIndexCell {
    /// Clipped shapes that intersect this cell
    shapes: Vec<S2ClippedShape>,
}

impl S2ShapeIndexCell {
    /// Create a new empty cell
    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
        }
    }

    /// Get the number of shapes in this cell
    pub fn num_clipped(&self) -> i32 {
        self.shapes.len() as i32
    }

    /// Get a clipped shape by index
    pub fn clipped(&self, index: i32) -> Option<&S2ClippedShape> {
        self.shapes.get(index as usize)
    }

    /// Add a clipped shape to this cell
    pub fn add_clipped(&mut self, clipped: S2ClippedShape) {
        self.shapes.push(clipped);
    }

    /// Iterate over all clipped shapes
    pub fn clipped_shapes(&self) -> impl Iterator<Item = &S2ClippedShape> {
        self.shapes.iter()
    }

    /// Clear all clipped shapes
    pub fn clear(&mut self) {
        self.shapes.clear();
    }
}

/// Abstract base class for indexing polygonal geometry
///
/// S2ShapeIndex makes it very fast to answer queries such as finding nearby
/// shapes, measuring distances, testing for intersection and containment.
pub trait S2ShapeIndex: fmt::Debug + Send + Sync {
    /// Iterator type for this index
    type Iterator: S2ShapeIndexIterator;

    /// Returns the number of distinct shape IDs in the index
    ///
    /// Shape IDs are assigned sequentially starting from 0, but may have gaps
    /// if shapes have been removed.
    fn num_shape_ids(&self) -> i32;

    /// Returns the shape with the given ID, or None if no such shape exists
    fn shape(&self, id: i32) -> Option<Arc<dyn S2Shape>>;

    /// Returns an iterator positioned at the beginning of the index
    fn iter(&self) -> Self::Iterator;

    /// Returns memory usage in bytes (approximate)
    fn space_used(&self) -> usize;

    /// Returns an iterator positioned at the given cell ID
    ///
    /// If the cell ID is not present in the index, the iterator is positioned
    /// at the next cell that is present, or at the end if no such cell exists.
    fn iter_at(&self, cell_id: S2CellId) -> Self::Iterator {
        let mut iter = self.iter();
        iter.seek(cell_id);
        iter
    }

    /// Begin an iterator positioned at the first cell
    fn begin(&self) -> Self::Iterator {
        self.iter()
    }

    /// End iterator (positioned past the last cell)
    fn end(&self) -> Self::Iterator {
        let mut iter = self.iter();
        iter.seek_to_end();
        iter
    }
}

/// Iterator for traversing the contents of an S2ShapeIndex
///
/// The iterator visits index cells in increasing order of S2CellId.
pub trait S2ShapeIndexIterator: fmt::Debug + Clone {
    /// Returns the S2CellId of the current index cell
    ///
    /// Requires: !done()
    fn cell_id(&self) -> S2CellId;

    /// Returns the contents of the current index cell
    ///
    /// Requires: !done()
    fn cell(&self) -> &S2ShapeIndexCell;

    /// Returns true if the iterator is positioned past the last index cell
    fn done(&self) -> bool;

    /// Advances the iterator to the next index cell
    fn next(&mut self);

    /// Positions the iterator at the first index cell
    fn begin(&mut self);

    /// Positions the iterator past the last index cell
    fn finish(&mut self);

    /// Positions the iterator at the first cell with ID >= target
    ///
    /// If no such cell exists, the iterator is positioned at the end.
    fn seek(&mut self, target: S2CellId);

    /// Positions the iterator past the last index cell
    fn seek_to_end(&mut self) {
        self.finish();
    }

    /// Returns the previous cell ID, or None if at the beginning
    fn prev(&mut self) -> Option<S2CellId>;
}

/// Trait for mutable shape indexes that allow adding and removing shapes
pub trait MutableS2ShapeIndex: S2ShapeIndex {
    /// Add a shape to the index and return its assigned ID
    ///
    /// The shape will be assigned the next available ID. The index will take
    /// ownership of the shape and may rebuild portions of the index as needed.
    fn add_shape(&mut self, shape: Arc<dyn S2Shape>) -> i32;

    /// Remove a shape from the index
    ///
    /// Returns the removed shape if it existed, or None if the ID was invalid.
    /// May trigger a partial rebuild of the index.
    fn remove_shape(&mut self, shape_id: i32) -> Option<Arc<dyn S2Shape>>;

    /// Remove all shapes from the index
    ///
    /// This efficiently resets the index to an empty state.
    fn clear(&mut self);

    /// Force the index to be rebuilt
    ///
    /// This is normally called automatically as needed, but can be called
    /// explicitly to control when the rebuild happens.
    fn force_build(&mut self);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_id_storage() {
        let mut storage = EdgeIdStorage::new();
        assert!(storage.is_empty());
        assert_eq!(storage.len(), 0);

        storage.push(1);
        assert!(!storage.is_empty());
        assert_eq!(storage.len(), 1);
        assert_eq!(storage.get(0), Some(1));

        storage.push(2);
        assert_eq!(storage.len(), 2);
        assert_eq!(storage.get(0), Some(1));
        assert_eq!(storage.get(1), Some(2));

        storage.push(3);
        assert_eq!(storage.len(), 3);
        assert_eq!(storage.get(0), Some(1));
        assert_eq!(storage.get(1), Some(2));
        assert_eq!(storage.get(2), Some(3));

        // Test iterator
        let ids: Vec<i32> = storage.iter().collect();
        assert_eq!(ids, vec![1, 2, 3]);
    }

    #[test]
    fn test_clipped_shape() {
        let mut clipped = S2ClippedShape::new(5, true);
        assert_eq!(clipped.shape_id(), 5);
        assert!(clipped.contains_center());
        assert_eq!(clipped.num_edges(), 0);
        assert!(clipped.is_empty());

        clipped.add_edge(10);
        clipped.add_edge(20);
        assert_eq!(clipped.num_edges(), 2);
        assert!(!clipped.is_empty());
        assert_eq!(clipped.edge_id(0), Some(10));
        assert_eq!(clipped.edge_id(1), Some(20));

        let edge_ids: Vec<i32> = clipped.edge_ids().collect();
        assert_eq!(edge_ids, vec![10, 20]);
    }

    #[test]
    fn test_shape_index_cell() {
        let mut cell = S2ShapeIndexCell::new();
        assert_eq!(cell.num_clipped(), 0);

        let clipped1 = S2ClippedShape::new(1, false);
        let clipped2 = S2ClippedShape::new(2, true);

        cell.add_clipped(clipped1);
        cell.add_clipped(clipped2);

        assert_eq!(cell.num_clipped(), 2);
        assert_eq!(cell.clipped(0).unwrap().shape_id(), 1);
        assert_eq!(cell.clipped(1).unwrap().shape_id(), 2);
        assert!(!cell.clipped(0).unwrap().contains_center());
        assert!(cell.clipped(1).unwrap().contains_center());

        let shape_ids: Vec<i32> = cell.clipped_shapes()
            .map(|c| c.shape_id())
            .collect();
        assert_eq!(shape_ids, vec![1, 2]);
    }
}