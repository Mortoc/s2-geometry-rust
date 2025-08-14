//! Polyline shapes - S2Shape implementations for polyline geometry
//!
//! This module provides shape implementations for polyline-based geometry,
//! including single polylines and collections of polylines.

use crate::shape::{S2Shape, Edge, Chain, ChainPosition, ShapeDimension, ReferencePoint};
use crate::point::S2Point;
use crate::polyline::S2Polyline;
use std::sync::Arc;

/// A shape representing a single polyline
#[derive(Debug, Clone)]
pub struct S2PolylineShape {
    polyline: Arc<S2Polyline>,
}

impl S2PolylineShape {
    /// Create a new polyline shape
    pub fn new(polyline: Arc<S2Polyline>) -> Self {
        Self { polyline }
    }

    /// Get the underlying polyline
    pub fn polyline(&self) -> &S2Polyline {
        &self.polyline
    }
}

impl S2Shape for S2PolylineShape {
    fn num_edges(&self) -> i32 {
        if self.polyline.num_vertices() == 0 {
            0
        } else {
            (self.polyline.num_vertices() - 1) as i32
        }
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.num_edges() {
            panic!("Invalid edge_id {} for polyline with {} edges", 
                   edge_id, self.num_edges());
        }
        
        let v0 = self.polyline.vertex(edge_id as usize).expect("Invalid edge_id");
        let v1 = self.polyline.vertex((edge_id + 1) as usize).expect("Invalid edge_id");
        Edge::new(v0, v1)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polyline
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.polyline.num_vertices() == 0 {
            // Use a default point
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            ReferencePoint::not_contained(self.polyline.vertex(0).expect("Empty polyline"))
        }
    }
}

/// A shape representing multiple disconnected polylines
#[derive(Debug, Clone)]
pub struct S2MultiPolylineShape {
    polylines: Vec<Arc<S2Polyline>>,
    /// Cached information about chains
    chain_starts: Vec<i32>,
    total_edges: i32,
}

impl S2MultiPolylineShape {
    /// Create a new multi-polyline shape
    pub fn new(polylines: Vec<Arc<S2Polyline>>) -> Self {
        let mut chain_starts = Vec::with_capacity(polylines.len() + 1);
        let mut total_edges = 0;
        
        chain_starts.push(0);
        for polyline in &polylines {
            let edges = if polyline.num_vertices() == 0 { 0 } else { polyline.num_vertices() - 1 };
            total_edges += edges as i32;
            chain_starts.push(total_edges);
        }
        
        Self {
            polylines,
            chain_starts,
            total_edges,
        }
    }

    /// Create an empty multi-polyline shape
    pub fn empty() -> Self {
        Self {
            polylines: Vec::new(),
            chain_starts: vec![0],
            total_edges: 0,
        }
    }

    /// Add a polyline to this shape
    pub fn add_polyline(&mut self, polyline: Arc<S2Polyline>) {
        let edges = if polyline.num_vertices() == 0 { 0 } else { polyline.num_vertices() - 1 };
        self.total_edges += edges as i32;
        self.chain_starts.push(self.total_edges);
        self.polylines.push(polyline);
    }

    /// Get the number of polylines
    pub fn num_polylines(&self) -> i32 {
        self.polylines.len() as i32
    }

    /// Get a polyline by index
    pub fn polyline(&self, index: i32) -> &S2Polyline {
        &self.polylines[index as usize]
    }

    /// Get all polylines
    pub fn polylines(&self) -> &[Arc<S2Polyline>] {
        &self.polylines
    }

    /// Check if this shape is empty
    pub fn is_empty(&self) -> bool {
        self.polylines.is_empty()
    }
}

impl S2Shape for S2MultiPolylineShape {
    fn num_edges(&self) -> i32 {
        self.total_edges
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for multi-polyline with {} edges", 
                   edge_id, self.total_edges);
        }
        
        // Find which polyline contains this edge
        let chain_pos = self.chain_position(edge_id);
        let polyline = &self.polylines[chain_pos.chain_id as usize];
        
        let v0 = polyline.vertex(chain_pos.offset as usize).expect("Invalid offset");
        let v1 = polyline.vertex((chain_pos.offset + 1) as usize).expect("Invalid offset");
        Edge::new(v0, v1)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polyline
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.polylines.is_empty() {
            // Use a default point
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            let first_polyline = &self.polylines[0];
            if first_polyline.num_vertices() > 0 {
                ReferencePoint::not_contained(first_polyline.vertex(0).expect("Empty polyline"))
            } else {
                let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
                ReferencePoint::not_contained(default_point)
            }
        }
    }

    fn num_chains(&self) -> i32 {
        self.polylines.len() as i32
    }

    fn chain(&self, chain_id: i32) -> Chain {
        if chain_id < 0 || chain_id >= self.num_chains() {
            panic!("Invalid chain_id {} for multi-polyline with {} chains", 
                   chain_id, self.num_chains());
        }
        
        let start = self.chain_starts[chain_id as usize];
        let end = self.chain_starts[chain_id as usize + 1];
        Chain::new(start, end - start)
    }

    fn chain_position(&self, edge_id: i32) -> ChainPosition {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for multi-polyline with {} edges", 
                   edge_id, self.total_edges);
        }
        
        // Binary search to find the chain containing this edge
        let chain_id = match self.chain_starts.binary_search(&edge_id) {
            Ok(index) => index as i32,
            Err(index) => (index - 1) as i32,
        };
        
        let offset = edge_id - self.chain_starts[chain_id as usize];
        ChainPosition::new(chain_id, offset)
    }
}

/// A wrapper that converts any S2Polyline into an S2Shape
#[derive(Debug, Clone)]
pub struct S2LaxPolylineShape {
    polyline: Arc<S2Polyline>,
    /// Whether to ignore degenerate edges (consecutive duplicate vertices)
    ignore_degeneracies: bool,
}

impl S2LaxPolylineShape {
    /// Create a new lax polyline shape that may contain degenerate edges
    pub fn new(polyline: Arc<S2Polyline>) -> Self {
        Self {
            polyline,
            ignore_degeneracies: true,
        }
    }

    /// Create a new lax polyline shape with explicit degeneracy handling
    pub fn with_degeneracy_handling(polyline: Arc<S2Polyline>, ignore_degeneracies: bool) -> Self {
        Self {
            polyline,
            ignore_degeneracies,
        }
    }

    /// Get the underlying polyline
    pub fn polyline(&self) -> &S2Polyline {
        &self.polyline
    }

    /// Check if degeneracies are ignored
    pub fn ignores_degeneracies(&self) -> bool {
        self.ignore_degeneracies
    }
}

impl S2Shape for S2LaxPolylineShape {
    fn num_edges(&self) -> i32 {
        if self.ignore_degeneracies {
            // Count only non-degenerate edges
            let mut count = 0;
            let num_edges = if self.polyline.num_vertices() == 0 { 0 } else { self.polyline.num_vertices() - 1 };
            for i in 0..num_edges {
                let v0 = self.polyline.vertex(i).expect("Invalid vertex index");
                let v1 = self.polyline.vertex(i + 1).expect("Invalid vertex index");
                if v0 != v1 {
                    count += 1;
                }
            }
            count
        } else {
            if self.polyline.num_vertices() == 0 { 0 } else { (self.polyline.num_vertices() - 1) as i32 }
        }
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if self.ignore_degeneracies {
            // Find the edge_id-th non-degenerate edge
            let mut count = 0;
            let num_edges = if self.polyline.num_vertices() == 0 { 0 } else { self.polyline.num_vertices() - 1 };
            for i in 0..num_edges {
                let v0 = self.polyline.vertex(i).expect("Invalid vertex index");
                let v1 = self.polyline.vertex(i + 1).expect("Invalid vertex index");
                if v0 != v1 {
                    if count == edge_id {
                        return Edge::new(v0, v1);
                    }
                    count += 1;
                }
            }
            panic!("Invalid edge_id {} for lax polyline", edge_id);
        } else {
            let num_edges = if self.polyline.num_vertices() == 0 { 0 } else { (self.polyline.num_vertices() - 1) as i32 };
            if edge_id < 0 || edge_id >= num_edges {
                panic!("Invalid edge_id {} for lax polyline with {} edges", 
                       edge_id, num_edges);
            }
            
            let v0 = self.polyline.vertex(edge_id as usize).expect("Invalid vertex index");
            let v1 = self.polyline.vertex((edge_id + 1) as usize).expect("Invalid vertex index");
            Edge::new(v0, v1)
        }
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polyline
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.polyline.num_vertices() == 0 {
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            ReferencePoint::not_contained(self.polyline.vertex(0).expect("Empty polyline"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;

    #[test]
    fn test_polyline_shape() {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        
        let polyline = Arc::new(S2Polyline::new(vertices.clone()).unwrap());
        let shape = S2PolylineShape::new(polyline.clone());

        assert_eq!(shape.num_edges(), 2);
        assert_eq!(shape.dimension(), ShapeDimension::Polyline);
        assert!(!shape.has_interior());

        let edge0 = shape.edge(0);
        assert_eq!(edge0.v0, vertices[0]);
        assert_eq!(edge0.v1, vertices[1]);

        let edge1 = shape.edge(1);
        assert_eq!(edge1.v0, vertices[1]);
        assert_eq!(edge1.v1, vertices[2]);

        let ref_point = shape.get_reference_point();
        assert_eq!(ref_point.point, vertices[0]);
        assert!(!ref_point.contained);

        assert_eq!(shape.polyline().num_vertices(), 3);
    }

    #[test]
    fn test_multi_polyline_shape() {
        let vertices1 = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        let vertices2 = vec![
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
            S2Point::from_normalized(DVec3::new(1.0, 1.0, 0.0).normalize()),
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 1.0).normalize()),
        ];

        let polyline1 = Arc::new(S2Polyline::new(vertices1.clone()).unwrap());
        let polyline2 = Arc::new(S2Polyline::new(vertices2.clone()).unwrap());
        
        let shape = S2MultiPolylineShape::new(vec![polyline1.clone(), polyline2.clone()]);

        assert_eq!(shape.num_edges(), 3); // 1 + 2 edges
        assert_eq!(shape.num_chains(), 2);
        assert_eq!(shape.num_polylines(), 2);
        assert_eq!(shape.dimension(), ShapeDimension::Polyline);
        assert!(!shape.is_empty());

        // Test first polyline's edge
        let edge0 = shape.edge(0);
        assert_eq!(edge0.v0, vertices1[0]);
        assert_eq!(edge0.v1, vertices1[1]);

        // Test second polyline's first edge
        let edge1 = shape.edge(1);
        assert_eq!(edge1.v0, vertices2[0]);
        assert_eq!(edge1.v1, vertices2[1]);

        // Test chains
        let chain0 = shape.chain(0);
        assert_eq!(chain0.start, 0);
        assert_eq!(chain0.length, 1);

        let chain1 = shape.chain(1);
        assert_eq!(chain1.start, 1);
        assert_eq!(chain1.length, 2);

        // Test chain positions
        let pos0 = shape.chain_position(0);
        assert_eq!(pos0.chain_id, 0);
        assert_eq!(pos0.offset, 0);

        let pos2 = shape.chain_position(2);
        assert_eq!(pos2.chain_id, 1);
        assert_eq!(pos2.offset, 1);

        assert_eq!(shape.polyline(0).num_vertices(), 2);
        assert_eq!(shape.polyline(1).num_vertices(), 3);
    }

    #[test]
    fn test_multi_polyline_shape_empty() {
        let shape = S2MultiPolylineShape::empty();
        assert_eq!(shape.num_edges(), 0);
        assert_eq!(shape.num_chains(), 0);
        assert_eq!(shape.num_polylines(), 0);
        assert!(shape.is_empty());
    }

    #[test]
    fn test_multi_polyline_shape_add() {
        let mut shape = S2MultiPolylineShape::empty();
        
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        let polyline = Arc::new(S2Polyline::new(vertices).unwrap());
        
        shape.add_polyline(polyline.clone());
        
        assert_eq!(shape.num_polylines(), 1);
        assert_eq!(shape.num_edges(), 1);
        assert_eq!(shape.num_chains(), 1);
        assert!(!shape.is_empty());
    }

    #[test]
    fn test_lax_polyline_shape_with_degeneracies() {
        // For this test, we'll create a simple polyline and test the two modes
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];

        let polyline = Arc::new(S2Polyline::new(vertices.clone()).unwrap());
        
        // With degeneracy handling (default)
        let lax_shape = S2LaxPolylineShape::new(polyline.clone());
        assert!(lax_shape.ignores_degeneracies());
        assert_eq!(lax_shape.num_edges(), 2); // Normal edges, no degeneracies to filter
        
        // Without degeneracy handling  
        let strict_shape = S2LaxPolylineShape::with_degeneracy_handling(polyline, false);
        assert!(!strict_shape.ignores_degeneracies());
        assert_eq!(strict_shape.num_edges(), 2); // Same as above since no degeneracies
    }
}