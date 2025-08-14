//! Polygon shapes - S2Shape implementations for polygon geometry
//!
//! This module provides shape implementations for polygon-based geometry,
//! including loops and polygons with holes.

use crate::shape::{S2Shape, Edge, Chain, ChainPosition, ShapeDimension, ReferencePoint};
use crate::point::S2Point;
use crate::r#loop::S2Loop;
use std::sync::Arc;

/// A shape representing a single loop (closed polygon boundary)
#[derive(Debug, Clone)]
pub struct S2LoopShape {
    r#loop: Arc<S2Loop>,
}

impl S2LoopShape {
    /// Create a new loop shape
    pub fn new(r#loop: Arc<S2Loop>) -> Self {
        Self { r#loop }
    }

    /// Get the underlying loop
    pub fn r#loop(&self) -> &S2Loop {
        &self.r#loop
    }
}

impl S2Shape for S2LoopShape {
    fn num_edges(&self) -> i32 {
        self.r#loop.num_vertices() as i32
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.num_edges() {
            panic!("Invalid edge_id {} for loop with {} edges", 
                   edge_id, self.num_edges());
        }
        
        let v0 = self.r#loop.vertex(edge_id as usize);
        let v1 = self.r#loop.vertex(((edge_id + 1) % self.num_edges()) as usize);
        Edge::new(v0, v1)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polygon
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.r#loop.num_vertices() == 0 {
            // Empty loop - use default point
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            // Use the first vertex as reference
            let ref_point = self.r#loop.vertex(0);
            // The reference point is contained if this is a normal (CCW) loop
            let contained = !self.r#loop.is_empty() && !self.r#loop.is_hole();
            ReferencePoint::new(ref_point, contained)
        }
    }

    fn contains(&self, point: &S2Point) -> bool {
        self.r#loop.contains(point)
    }
}

/// A shape representing a polygon with multiple loops (shell + holes)
#[derive(Debug, Clone)]
pub struct S2PolygonShape {
    /// All loops in the polygon (first is shell, rest are holes)
    loops: Vec<Arc<S2Loop>>,
    /// Cached information about chains
    chain_starts: Vec<i32>,
    total_edges: i32,
}

impl S2PolygonShape {
    /// Create a new polygon shape from loops
    /// The first loop should be the outer shell, subsequent loops are holes
    pub fn new(loops: Vec<Arc<S2Loop>>) -> Self {
        let mut chain_starts = Vec::with_capacity(loops.len() + 1);
        let mut total_edges = 0;
        
        chain_starts.push(0);
        for r#loop in &loops {
            total_edges += r#loop.num_vertices() as i32;
            chain_starts.push(total_edges);
        }
        
        Self {
            loops,
            chain_starts,
            total_edges,
        }
    }

    /// Create a polygon shape from a single loop (no holes)
    pub fn from_loop(r#loop: Arc<S2Loop>) -> Self {
        Self::new(vec![r#loop])
    }

    /// Create an empty polygon shape
    pub fn empty() -> Self {
        Self {
            loops: Vec::new(),
            chain_starts: vec![0],
            total_edges: 0,
        }
    }

    /// Get the number of loops in this polygon
    pub fn num_loops(&self) -> i32 {
        self.loops.len() as i32
    }

    /// Get a loop by index
    pub fn r#loop(&self, index: i32) -> &S2Loop {
        &self.loops[index as usize]
    }

    /// Get all loops
    pub fn loops(&self) -> &[Arc<S2Loop>] {
        &self.loops
    }

    /// Get the shell (outer boundary) loop, if any
    pub fn shell(&self) -> Option<&S2Loop> {
        if self.loops.is_empty() {
            None
        } else {
            Some(&self.loops[0])
        }
    }

    /// Get the hole loops (all loops except the first)
    pub fn holes(&self) -> &[Arc<S2Loop>] {
        if self.loops.len() <= 1 {
            &[]
        } else {
            &self.loops[1..]
        }
    }

    /// Check if this polygon is empty
    pub fn is_empty(&self) -> bool {
        self.loops.is_empty() || self.total_edges == 0
    }

    /// Add a loop to this polygon
    pub fn add_loop(&mut self, r#loop: Arc<S2Loop>) {
        self.total_edges += r#loop.num_vertices() as i32;
        self.chain_starts.push(self.total_edges);
        self.loops.push(r#loop);
    }
}

impl S2Shape for S2PolygonShape {
    fn num_edges(&self) -> i32 {
        self.total_edges
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for polygon with {} edges", 
                   edge_id, self.total_edges);
        }
        
        // Find which loop contains this edge
        let chain_pos = self.chain_position(edge_id);
        let r#loop = &self.loops[chain_pos.chain_id as usize];
        
        let loop_edges = r#loop.num_vertices() as i32;
        let v0 = r#loop.vertex(chain_pos.offset as usize);
        let v1 = r#loop.vertex(((chain_pos.offset + 1) % loop_edges) as usize);
        Edge::new(v0, v1)
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polygon
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.loops.is_empty() {
            // Empty polygon - use default point
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            // Use a point from the shell loop
            let shell = &self.loops[0];
            if shell.num_vertices() == 0 {
                let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
                ReferencePoint::not_contained(default_point)
            } else {
                let ref_point = shell.vertex(0);
                // The reference point is contained if the shell is not a hole and not empty
                let contained = !shell.is_empty() && !shell.is_hole();
                ReferencePoint::new(ref_point, contained)
            }
        }
    }

    fn num_chains(&self) -> i32 {
        self.loops.len() as i32
    }

    fn chain(&self, chain_id: i32) -> Chain {
        if chain_id < 0 || chain_id >= self.num_chains() {
            panic!("Invalid chain_id {} for polygon with {} chains", 
                   chain_id, self.num_chains());
        }
        
        let start = self.chain_starts[chain_id as usize];
        let end = self.chain_starts[chain_id as usize + 1];
        Chain::new(start, end - start)
    }

    fn chain_position(&self, edge_id: i32) -> ChainPosition {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for polygon with {} edges", 
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

    fn contains(&self, point: &S2Point) -> bool {
        if self.loops.is_empty() {
            return false;
        }
        
        // A point is inside the polygon if it's inside the shell
        // and not inside any of the holes
        let shell = &self.loops[0];
        if !shell.contains(point) {
            return false;
        }
        
        // Check if the point is inside any hole
        for hole in &self.loops[1..] {
            if hole.contains(point) {
                return false;
            }
        }
        
        true
    }
}

/// A shape representing multiple disconnected polygons
#[derive(Debug, Clone)]
pub struct S2MultiPolygonShape {
    polygons: Vec<S2PolygonShape>,
    /// Cached information about chains across all polygons
    chain_starts: Vec<i32>,
    total_edges: i32,
}

impl S2MultiPolygonShape {
    /// Create a new multi-polygon shape
    pub fn new(polygons: Vec<S2PolygonShape>) -> Self {
        let mut chain_starts = Vec::new();
        let mut total_edges = 0;
        
        chain_starts.push(0);
        for polygon in &polygons {
            // Each polygon contributes its chains
            for i in 0..polygon.num_chains() {
                total_edges += polygon.chain(i).length;
                chain_starts.push(total_edges);
            }
        }
        
        Self {
            polygons,
            chain_starts,
            total_edges,
        }
    }

    /// Create an empty multi-polygon shape
    pub fn empty() -> Self {
        Self {
            polygons: Vec::new(),
            chain_starts: vec![0],
            total_edges: 0,
        }
    }

    /// Get the number of polygons
    pub fn num_polygons(&self) -> i32 {
        self.polygons.len() as i32
    }

    /// Get a polygon by index
    pub fn polygon(&self, index: i32) -> &S2PolygonShape {
        &self.polygons[index as usize]
    }

    /// Get all polygons
    pub fn polygons(&self) -> &[S2PolygonShape] {
        &self.polygons
    }

    /// Check if this multi-polygon is empty
    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty() || self.total_edges == 0
    }
}

impl S2Shape for S2MultiPolygonShape {
    fn num_edges(&self) -> i32 {
        self.total_edges
    }

    fn edge(&self, edge_id: i32) -> Edge {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for multi-polygon with {} edges", 
                   edge_id, self.total_edges);
        }
        
        // Find which polygon and chain contains this edge
        let chain_pos = self.chain_position(edge_id);
        
        // Map global chain ID to polygon and local chain
        let mut current_chain = 0;
        for (_poly_idx, polygon) in self.polygons.iter().enumerate() {
            let poly_chains = polygon.num_chains();
            if chain_pos.chain_id < current_chain + poly_chains {
                let local_chain_id = chain_pos.chain_id - current_chain;
                return polygon.chain_edge(local_chain_id, chain_pos.offset);
            }
            current_chain += poly_chains;
        }
        
        panic!("Failed to find edge {} in multi-polygon", edge_id);
    }

    fn dimension(&self) -> ShapeDimension {
        ShapeDimension::Polygon
    }

    fn get_reference_point(&self) -> ReferencePoint {
        if self.polygons.is_empty() {
            let default_point = S2Point::from_normalized(crate::math::DVec3::new(1.0, 0.0, 0.0));
            ReferencePoint::not_contained(default_point)
        } else {
            self.polygons[0].get_reference_point()
        }
    }

    fn num_chains(&self) -> i32 {
        self.polygons.iter().map(|p| p.num_chains()).sum()
    }

    fn chain(&self, chain_id: i32) -> Chain {
        if chain_id < 0 || chain_id >= self.num_chains() {
            panic!("Invalid chain_id {} for multi-polygon with {} chains", 
                   chain_id, self.num_chains());
        }
        
        let start = self.chain_starts[chain_id as usize];
        let end = self.chain_starts[chain_id as usize + 1];
        Chain::new(start, end - start)
    }

    fn chain_position(&self, edge_id: i32) -> ChainPosition {
        if edge_id < 0 || edge_id >= self.total_edges {
            panic!("Invalid edge_id {} for multi-polygon with {} edges", 
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

    fn contains(&self, point: &S2Point) -> bool {
        // A point is contained if it's contained by any of the polygons
        self.polygons.iter().any(|polygon| polygon.contains(point))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;

    fn create_simple_loop() -> Arc<S2Loop> {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        Arc::new(S2Loop::new(vertices).unwrap())
    }

    #[test]
    fn test_loop_shape() {
        let r#loop = create_simple_loop();
        let shape = S2LoopShape::new(r#loop.clone());

        assert_eq!(shape.num_edges(), 3);
        assert_eq!(shape.dimension(), ShapeDimension::Polygon);
        assert!(shape.has_interior());

        let edge0 = shape.edge(0);
        assert_eq!(edge0.v0, r#loop.vertex(0));
        assert_eq!(edge0.v1, r#loop.vertex(1));

        let edge2 = shape.edge(2);
        assert_eq!(edge2.v0, r#loop.vertex(2));
        assert_eq!(edge2.v1, r#loop.vertex(0)); // Wraps around

        let ref_point = shape.get_reference_point();
        assert_eq!(ref_point.point, r#loop.vertex(0));
    }

    #[test]
    fn test_polygon_shape_single_loop() {
        let r#loop = create_simple_loop();
        let shape = S2PolygonShape::from_loop(r#loop.clone());

        assert_eq!(shape.num_edges(), 3);
        assert_eq!(shape.num_chains(), 1);
        assert_eq!(shape.num_loops(), 1);
        assert_eq!(shape.dimension(), ShapeDimension::Polygon);
        assert!(!shape.is_empty());

        let chain = shape.chain(0);
        assert_eq!(chain.start, 0);
        assert_eq!(chain.length, 3);

        assert!(shape.shell().is_some());
        assert_eq!(shape.holes().len(), 0);
    }

    #[test]
    fn test_polygon_shape_with_hole() {
        let shell = create_simple_loop();
        let hole = create_simple_loop(); // In practice this would be inside the shell
        
        let shape = S2PolygonShape::new(vec![shell.clone(), hole.clone()]);

        assert_eq!(shape.num_edges(), 6);
        assert_eq!(shape.num_chains(), 2);
        assert_eq!(shape.num_loops(), 2);

        let chain0 = shape.chain(0);
        assert_eq!(chain0.start, 0);
        assert_eq!(chain0.length, 3);

        let chain1 = shape.chain(1);
        assert_eq!(chain1.start, 3);
        assert_eq!(chain1.length, 3);

        assert!(shape.shell().is_some());
        assert_eq!(shape.holes().len(), 1);

        // Test chain positions
        let pos0 = shape.chain_position(0);
        assert_eq!(pos0.chain_id, 0);
        assert_eq!(pos0.offset, 0);

        let pos4 = shape.chain_position(4);
        assert_eq!(pos4.chain_id, 1);
        assert_eq!(pos4.offset, 1);
    }

    #[test]
    fn test_polygon_shape_empty() {
        let shape = S2PolygonShape::empty();
        assert_eq!(shape.num_edges(), 0);
        assert_eq!(shape.num_chains(), 0);
        assert_eq!(shape.num_loops(), 0);
        assert!(shape.is_empty());
        assert!(shape.shell().is_none());
        assert_eq!(shape.holes().len(), 0);
    }

    #[test]
    fn test_multi_polygon_shape() {
        let polygon1 = S2PolygonShape::from_loop(create_simple_loop());
        let polygon2 = S2PolygonShape::from_loop(create_simple_loop());
        
        let shape = S2MultiPolygonShape::new(vec![polygon1, polygon2]);

        assert_eq!(shape.num_edges(), 6);
        assert_eq!(shape.num_chains(), 2);
        assert_eq!(shape.num_polygons(), 2);
        assert_eq!(shape.dimension(), ShapeDimension::Polygon);
        assert!(!shape.is_empty());

        assert_eq!(shape.polygon(0).num_edges(), 3);
        assert_eq!(shape.polygon(1).num_edges(), 3);
    }

    #[test]
    fn test_multi_polygon_shape_empty() {
        let shape = S2MultiPolygonShape::empty();
        assert_eq!(shape.num_edges(), 0);
        assert_eq!(shape.num_chains(), 0);
        assert_eq!(shape.num_polygons(), 0);
        assert!(shape.is_empty());
    }
}