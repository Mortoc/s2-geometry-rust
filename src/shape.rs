//! S2Shape - Abstract interface for geometric objects that can be indexed
//!
//! This module provides the core S2Shape trait and related types for representing
//! polygonal geometry (points, polylines, and polygons) in a flexible way that
//! allows different underlying data representations.

use crate::point::S2Point;
use std::fmt;

/// The dimension of a shape (0=point, 1=polyline, 2=polygon)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShapeDimension {
    /// Point geometry (0-dimensional)
    Point = 0,
    /// Polyline geometry (1-dimensional)
    Polyline = 1,
    /// Polygon geometry (2-dimensional)
    Polygon = 2,
}

impl From<u8> for ShapeDimension {
    fn from(value: u8) -> Self {
        match value {
            0 => ShapeDimension::Point,
            1 => ShapeDimension::Polyline,
            2 => ShapeDimension::Polygon,
            _ => panic!("Invalid shape dimension: {}", value),
        }
    }
}

impl From<ShapeDimension> for u8 {
    fn from(dim: ShapeDimension) -> Self {
        dim as u8
    }
}

/// An edge connecting two vertices
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Edge {
    /// First vertex of the edge
    pub v0: S2Point,
    /// Second vertex of the edge
    pub v1: S2Point,
}

impl Edge {
    /// Create a new edge from two points
    pub fn new(v0: S2Point, v1: S2Point) -> Self {
        Self { v0, v1 }
    }
}

/// A contiguous range of edges that form a connected sequence
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chain {
    /// Starting edge index in the shape
    pub start: i32,
    /// Number of edges in this chain
    pub length: i32,
}

impl Chain {
    /// Create a new chain
    pub fn new(start: i32, length: i32) -> Self {
        Self { start, length }
    }

    /// Get the ending edge index (exclusive)
    pub fn end(&self) -> i32 {
        self.start + self.length
    }

    /// Check if this chain contains the given edge index
    pub fn contains(&self, edge_id: i32) -> bool {
        edge_id >= self.start && edge_id < self.end()
    }
}

/// Position of an edge within a specific chain
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainPosition {
    /// Index of the chain
    pub chain_id: i32,
    /// Offset within the chain
    pub offset: i32,
}

impl ChainPosition {
    /// Create a new chain position
    pub fn new(chain_id: i32, offset: i32) -> Self {
        Self { chain_id, offset }
    }
}

/// A reference point with containment information for polygon shapes
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferencePoint {
    /// The reference point
    pub point: S2Point,
    /// Whether the point is contained within the shape
    /// For non-polygon shapes, this is always false
    pub contained: bool,
}

impl ReferencePoint {
    /// Create a new reference point
    pub fn new(point: S2Point, contained: bool) -> Self {
        Self { point, contained }
    }

    /// Create a reference point that is not contained (for non-polygon shapes)
    pub fn not_contained(point: S2Point) -> Self {
        Self::new(point, false)
    }

    /// Create a reference point that is contained (for polygon shapes)
    pub fn contained(point: S2Point) -> Self {
        Self::new(point, true)
    }
}

/// Abstract interface for geometric shapes that can be indexed by S2ShapeIndex
///
/// S2Shape represents polygonal geometry as a collection of edges, where edges
/// may represent:
/// - Points (0-dimensional): Each "edge" represents a single point
/// - Polylines (1-dimensional): Edges form connected line segments
/// - Polygons (2-dimensional): Edges form closed loops with possible holes
///
/// All edges in a shape must have the same dimension.
pub trait S2Shape: fmt::Debug + Send + Sync {
    /// Returns the number of edges in this shape
    /// 
    /// For point shapes, this equals the number of points.
    /// For polyline shapes, this equals the number of line segments.
    /// For polygon shapes, this equals the total number of edges in all loops.
    fn num_edges(&self) -> i32;

    /// Returns the specified edge of this shape
    ///
    /// # Panics
    /// Panics if edge_id is not in the range [0, num_edges())
    fn edge(&self, edge_id: i32) -> Edge;

    /// Returns the dimension of this shape (0, 1, or 2)
    fn dimension(&self) -> ShapeDimension;

    /// Returns a reference point for this shape with containment information
    ///
    /// For polygon shapes, the point may be contained within the shape.
    /// For non-polygon shapes, the point is never contained.
    fn get_reference_point(&self) -> ReferencePoint;

    /// Returns the number of contiguous edge chains in this shape
    ///
    /// Most shapes have a single chain, but shapes with multiple disconnected
    /// components may have multiple chains.
    fn num_chains(&self) -> i32 {
        if self.num_edges() == 0 { 0 } else { 1 }
    }

    /// Returns the specified edge chain
    ///
    /// # Panics
    /// Panics if chain_id is not in the range [0, num_chains())
    fn chain(&self, chain_id: i32) -> Chain {
        if chain_id != 0 {
            panic!("Invalid chain_id: {}", chain_id);
        }
        Chain::new(0, self.num_edges())
    }

    /// Returns an edge from the specified chain
    ///
    /// This is equivalent to edge(chain(chain_id).start + offset)
    ///
    /// # Panics
    /// Panics if chain_id or offset is invalid
    fn chain_edge(&self, chain_id: i32, offset: i32) -> Edge {
        let chain = self.chain(chain_id);
        if offset < 0 || offset >= chain.length {
            panic!("Invalid offset {} for chain {} with length {}", offset, chain_id, chain.length);
        }
        self.edge(chain.start + offset)
    }

    /// Returns the chain position containing the given edge
    ///
    /// # Panics
    /// Panics if edge_id is not in the range [0, num_edges())
    fn chain_position(&self, edge_id: i32) -> ChainPosition {
        if edge_id < 0 || edge_id >= self.num_edges() {
            panic!("Invalid edge_id: {}", edge_id);
        }
        
        // For shapes with a single chain, this is simple
        if self.num_chains() == 1 {
            return ChainPosition::new(0, edge_id);
        }
        
        // For multiple chains, find which chain contains this edge
        for chain_id in 0..self.num_chains() {
            let chain = self.chain(chain_id);
            if chain.contains(edge_id) {
                return ChainPosition::new(chain_id, edge_id - chain.start);
            }
        }
        
        panic!("Edge {} not found in any chain", edge_id);
    }

    /// Returns an optional type tag for this shape
    ///
    /// This allows applications to associate user-defined metadata with shapes.
    /// The default implementation returns None.
    fn type_tag(&self) -> Option<i32> {
        None
    }

    /// Returns whether this shape has an interior
    ///
    /// This is true for polygons (2D) and false for points and polylines.
    fn has_interior(&self) -> bool {
        self.dimension() == ShapeDimension::Polygon
    }

    /// Returns whether this shape contains the given point
    ///
    /// This is only meaningful for polygon shapes. For other shapes, returns false.
    fn contains(&self, _point: &S2Point) -> bool {
        // Default implementation for non-polygon shapes
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::S2Point;
    use crate::math::DVec3;

    // Simple point shape implementation for testing
    #[derive(Debug)]
    struct TestPointShape {
        points: Vec<S2Point>,
    }

    impl TestPointShape {
        fn new(points: Vec<S2Point>) -> Self {
            Self { points }
        }
    }

    impl S2Shape for TestPointShape {
        fn num_edges(&self) -> i32 {
            self.points.len() as i32
        }

        fn edge(&self, edge_id: i32) -> Edge {
            let point = self.points[edge_id as usize];
            Edge::new(point, point) // For points, both vertices are the same
        }

        fn dimension(&self) -> ShapeDimension {
            ShapeDimension::Point
        }

        fn get_reference_point(&self) -> ReferencePoint {
            if self.points.is_empty() {
                // Use origin as default
                let origin = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
                ReferencePoint::not_contained(origin)
            } else {
                ReferencePoint::not_contained(self.points[0])
            }
        }
    }

    #[test]
    fn test_shape_dimension() {
        assert_eq!(ShapeDimension::Point as u8, 0);
        assert_eq!(ShapeDimension::Polyline as u8, 1);
        assert_eq!(ShapeDimension::Polygon as u8, 2);
    }

    #[test]
    fn test_edge() {
        let p1 = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let p2 = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
        let edge = Edge::new(p1, p2);
        
        assert_eq!(edge.v0, p1);
        assert_eq!(edge.v1, p2);
    }

    #[test]
    fn test_chain() {
        let chain = Chain::new(5, 3);
        assert_eq!(chain.start, 5);
        assert_eq!(chain.length, 3);
        assert_eq!(chain.end(), 8);
        assert!(chain.contains(5));
        assert!(chain.contains(7));
        assert!(!chain.contains(4));
        assert!(!chain.contains(8));
    }

    #[test]
    fn test_chain_position() {
        let pos = ChainPosition::new(2, 3);
        assert_eq!(pos.chain_id, 2);
        assert_eq!(pos.offset, 3);
    }

    #[test]
    fn test_reference_point() {
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let ref_point = ReferencePoint::new(point, true);
        
        assert_eq!(ref_point.point, point);
        assert!(ref_point.contained);
        
        let not_contained = ReferencePoint::not_contained(point);
        assert!(!not_contained.contained);
        
        let contained = ReferencePoint::contained(point);
        assert!(contained.contained);
    }

    #[test]
    fn test_point_shape() {
        let points = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        let shape = TestPointShape::new(points.clone());
        
        assert_eq!(shape.num_edges(), 3);
        assert_eq!(shape.dimension(), ShapeDimension::Point);
        assert_eq!(shape.num_chains(), 1);
        assert!(!shape.has_interior());
        
        let edge = shape.edge(1);
        assert_eq!(edge.v0, points[1]);
        assert_eq!(edge.v1, points[1]);
        
        let chain = shape.chain(0);
        assert_eq!(chain.start, 0);
        assert_eq!(chain.length, 3);
        
        let chain_pos = shape.chain_position(1);
        assert_eq!(chain_pos.chain_id, 0);
        assert_eq!(chain_pos.offset, 1);
        
        let ref_point = shape.get_reference_point();
        assert_eq!(ref_point.point, points[0]);
        assert!(!ref_point.contained);
    }
}