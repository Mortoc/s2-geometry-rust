//! Point shapes - S2Shape implementations for point geometry
//!
//! This module provides shape implementations for point-based geometry,
//! including single points and collections of points.

use crate::shape::{S2Shape, Edge, ShapeDimension, ReferencePoint};
use crate::point::S2Point;
use std::fmt;

/// A shape representing a single point
#[derive(Debug, Clone)]
pub struct S2PointShape {
    point: S2Point,
}

impl S2PointShape {
    /// Create a new point shape
    pub fn new(point: S2Point) -> Self {
        Self { point }
    }

    /// Get the point
    pub fn point(&self) -> S2Point {
        self.point
    }
}

impl S2Shape for S2PointShape {
    fn num_edges(&self) -> i32 {
        1
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id != 0 {
            panic!("Invalid edge_id {} for point shape", edge_id);
        }
        Edge::new(self.point, self.point)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Point
    }

    fn get_reference_point(&self) -> ReferencePoint {
        ReferencePoint::not_contained(self.point)
    }
}

/// A shape representing multiple points
#[derive(Debug, Clone)]
pub struct S2MultiPointShape {
    points: Vec<S2Point>,
}

impl S2MultiPointShape {
    /// Create a new multi-point shape
    pub fn new(points: Vec<S2Point>) -> Self {
        Self { points }
    }

    /// Create an empty multi-point shape
    pub fn empty() -> Self {
        Self { points: Vec::new() }
    }

    /// Add a point to this shape
    pub fn add_point(&mut self, point: S2Point) {
        self.points.push(point);
    }

    /// Get the number of points
    pub fn num_points(&self) -> i32 {
        self.points.len() as i32
    }

    /// Get a point by index
    pub fn point(&self, index: i32) -> S2Point {
        self.points[index as usize]
    }

    /// Get all points
    pub fn points(&self) -> &[S2Point] {
        &self.points
    }

    /// Check if this shape is empty
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

impl S2Shape for S2MultiPointShape {
    fn num_edges(&self) -> i32 {
        self.points.len() as i32
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.num_edges() {
            panic!("Invalid edge_id {} for multi-point shape with {} points", 
                   edge_id, self.num_edges());
        }
        let point = self.points[edge_id as usize];
        Edge::new(point, point)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Point
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.points.is_empty() {
            // Use a default point (arbitrary choice)
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            ReferencePoint::not_contained(self.points[0])
        }
    }
}

/// A shape representing a cloud of points with optional metadata
#[derive(Debug, Clone)]
pub struct S2PointCloudShape<T = ()> {
    points: Vec<S2Point>,
    metadata: Vec<T>,
}

impl<T> S2PointCloudShape<T> {
    /// Create a new point cloud shape with metadata
    pub fn new(points: Vec<S2Point>, metadata: Vec<T>) -> Self {
        assert_eq!(points.len(), metadata.len(), 
                   "Points and metadata must have the same length");
        Self { points, metadata }
    }

    /// Get a point by index
    pub fn point(&self, index: i32) -> S2Point {
        self.points[index as usize]
    }

    /// Get metadata by index
    pub fn metadata(&self, index: i32) -> &T {
        &self.metadata[index as usize]
    }

    /// Get the number of points
    pub fn num_points(&self) -> i32 {
        self.points.len() as i32
    }

    /// Get all points
    pub fn points(&self) -> &[S2Point] {
        &self.points
    }

    /// Get all metadata
    pub fn all_metadata(&self) -> &[T] {
        &self.metadata
    }
}

impl<T: fmt::Debug + Send + Sync> S2Shape for S2PointCloudShape<T> {
    fn num_edges(&self) -> i32 {
        self.points.len() as i32
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.num_edges() {
            panic!("Invalid edge_id {} for point cloud with {} points", 
                   edge_id, self.num_edges());
        }
        let point = self.points[edge_id as usize];
        Edge::new(point, point)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Point
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.points.is_empty() {
            // Use a default point
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            ReferencePoint::not_contained(self.points[0])
        }
    }
}

impl S2PointCloudShape<()> {
    /// Create a point cloud without metadata
    pub fn from_points(points: Vec<S2Point>) -> Self {
        let metadata = vec![(); points.len()];
        Self { points, metadata }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;

    #[test]
    fn test_point_shape() {
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let shape = S2PointShape::new(point);

        assert_eq!(shape.num_edges(), 1);
        assert_eq!(shape.dimension(), ShapeDimension::Point);
        assert!(!shape.has_interior());

        let edge = shape.edge(0);
        assert_eq!(edge.v0, point);
        assert_eq!(edge.v1, point);

        let ref_point = shape.get_reference_point();
        assert_eq!(ref_point.point, point);
        assert!(!ref_point.contained);

        assert_eq!(shape.point(), point);
    }

    #[test]
    #[should_panic(expected = "Invalid edge_id 1 for point shape")]
    fn test_point_shape_invalid_edge() {
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let shape = S2PointShape::new(point);
        shape.edge(1);
    }

    #[test]
    fn test_multi_point_shape() {
        let points = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        let shape = S2MultiPointShape::new(points.clone());

        assert_eq!(shape.num_edges(), 3);
        assert_eq!(shape.num_points(), 3);
        assert_eq!(shape.dimension(), ShapeDimension::Point);
        assert!(!shape.has_interior());
        assert!(!shape.is_empty());

        for i in 0..3 {
            let edge = shape.edge(i);
            assert_eq!(edge.v0, points[i as usize]);
            assert_eq!(edge.v1, points[i as usize]);
            assert_eq!(shape.point(i), points[i as usize]);
        }

        let ref_point = shape.get_reference_point();
        assert_eq!(ref_point.point, points[0]);
        assert!(!ref_point.contained);

        assert_eq!(shape.points(), &points);
    }

    #[test]
    fn test_multi_point_shape_empty() {
        let shape = S2MultiPointShape::empty();
        assert_eq!(shape.num_edges(), 0);
        assert_eq!(shape.num_points(), 0);
        assert!(shape.is_empty());

        let ref_point = shape.get_reference_point();
        assert!(!ref_point.contained);
    }

    #[test]
    fn test_multi_point_shape_add() {
        let mut shape = S2MultiPointShape::empty();
        let point1 = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let point2 = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));

        shape.add_point(point1);
        assert_eq!(shape.num_points(), 1);
        assert_eq!(shape.point(0), point1);

        shape.add_point(point2);
        assert_eq!(shape.num_points(), 2);
        assert_eq!(shape.point(1), point2);
    }

    #[test]
    fn test_point_cloud_shape() {
        let points = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        let metadata = vec!["point1", "point2"];
        let shape = S2PointCloudShape::new(points.clone(), metadata.clone());

        assert_eq!(shape.num_edges(), 2);
        assert_eq!(shape.num_points(), 2);
        assert_eq!(shape.dimension(), ShapeDimension::Point);

        assert_eq!(shape.point(0), points[0]);
        assert_eq!(shape.point(1), points[1]);
        assert_eq!(shape.metadata(0), &"point1");
        assert_eq!(shape.metadata(1), &"point2");

        assert_eq!(shape.points(), &points);
        assert_eq!(shape.all_metadata(), &metadata);
    }

    #[test]
    fn test_point_cloud_shape_from_points() {
        let points = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        let shape = S2PointCloudShape::from_points(points.clone());

        assert_eq!(shape.num_points(), 2);
        assert_eq!(shape.point(0), points[0]);
        assert_eq!(shape.point(1), points[1]);
        assert_eq!(shape.metadata(0), &());
        assert_eq!(shape.metadata(1), &());
    }

    #[test]
    #[should_panic(expected = "Points and metadata must have the same length")]
    fn test_point_cloud_shape_mismatched_lengths() {
        let points = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
        ];
        let metadata = vec!["point1", "point2"];
        S2PointCloudShape::new(points, metadata);
    }
}