//! S2Builder Module - Robust Geometric Construction System
//!
//! This module provides the complete S2Builder system for constructing reliable
//! geometric objects from imperfect input data through snap rounding and
//! topology preservation.

use crate::{S2Point, S2Loop, S2Polyline, S2Error, S2Result, S1Angle};
use std::collections::HashMap;

pub mod snap_functions;
pub mod graph;
pub mod layers;

pub use snap_functions::*;
pub use graph::*;
pub use layers::*;

/// Type-safe identifier for edges in the S2Builder graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeId(pub usize);

/// Type-safe identifier for vertices in the S2Builder graph  
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexId(pub usize);

/// Edge type classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    /// Undirected edge
    Undirected,
    /// Directed edge
    Directed,
}

/// Configuration options for S2Builder
#[derive(Debug, Clone)]
pub struct Options {
    /// Snap function for vertex positioning
    pub snap_function: SnapFunction,
    /// Whether to validate output geometry
    pub validate: bool,
    /// Whether to split self-intersecting edges
    pub split_crossing_edges: bool,
    /// Intersection tolerance for edge intersection calculations
    pub intersection_tolerance: S1Angle,
    /// Whether to simplify edge chains
    pub simplify_edge_chains: bool,
    /// Whether the builder should be idempotent
    pub idempotent: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            snap_function: SnapFunction::Identity(IdentitySnapFunction::with_no_separation()),
            validate: true,
            split_crossing_edges: false,
            intersection_tolerance: S1Angle::zero(),
            simplify_edge_chains: false,
            idempotent: true,
        }
    }
}

impl Options {
    /// Create new options with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the snap function
    pub fn with_snap_function(mut self, snap_function: SnapFunction) -> Self {
        self.snap_function = snap_function;
        self
    }

    /// Set whether to validate output
    pub fn with_validate(mut self, validate: bool) -> Self {
        self.validate = validate;
        self
    }

    /// Set whether to split crossing edges
    pub fn with_split_crossing_edges(mut self, split_crossing_edges: bool) -> Self {
        self.split_crossing_edges = split_crossing_edges;
        self
    }

    /// Set the intersection tolerance
    pub fn with_intersection_tolerance(mut self, intersection_tolerance: S1Angle) -> Self {
        self.intersection_tolerance = intersection_tolerance;
        self
    }

    /// Set whether to simplify edge chains
    pub fn with_simplify_edge_chains(mut self, simplify_edge_chains: bool) -> Self {
        self.simplify_edge_chains = simplify_edge_chains;
        self
    }

    /// Set whether the builder should be idempotent
    pub fn with_idempotent(mut self, idempotent: bool) -> Self {
        self.idempotent = idempotent;
        self
    }

    /// Get whether to split crossing edges
    pub fn split_crossing_edges(&self) -> bool {
        self.split_crossing_edges
    }

    /// Get whether to validate output
    pub fn validate(&self) -> bool {
        self.validate
    }

    /// Get whether to simplify edge chains
    pub fn simplify_edge_chains(&self) -> bool {
        self.simplify_edge_chains
    }

    /// Get whether the builder should be idempotent
    pub fn idempotent(&self) -> bool {
        self.idempotent
    }

    /// Get the intersection tolerance (accounting for split_crossing_edges)
    pub fn intersection_tolerance(&self) -> S1Angle {
        if !self.split_crossing_edges {
            self.intersection_tolerance
        } else {
            // When split_crossing_edges is true, use a minimum of intersection error
            // For now, use a small default value - this should be S2::kIntersectionError
            self.intersection_tolerance.max(S1Angle::from_radians(1e-15))
        }
    }
}

/// Main S2Builder class for robust geometry construction
pub struct S2Builder {
    options: Options,
    input_vertices: Vec<S2Point>,
    input_edges: Vec<(VertexId, VertexId, EdgeType)>,
    layers: Vec<Box<dyn Layer>>,
    built: bool,
}

impl S2Builder {
    /// Create new S2Builder with default options
    pub fn new() -> Self {
        Self::with_options(Options::default())
    }

    /// Create new S2Builder with custom options
    pub fn with_options(options: Options) -> Self {
        Self {
            options,
            input_vertices: Vec::new(),
            input_edges: Vec::new(),
            layers: Vec::new(),
            built: false,
        }
    }

    /// Add a vertex to the builder
    pub fn add_vertex(&mut self, vertex: S2Point) -> S2Result<VertexId> {
        if self.built {
            return Err(S2Error::BuilderError { reason: "Cannot add vertex after build() has been called".to_string() });
        }
        let id = VertexId(self.input_vertices.len());
        self.input_vertices.push(vertex);
        Ok(id)
    }

    /// Add an edge between two points directly
    pub fn add_edge(&mut self, point_a: S2Point, point_b: S2Point) -> S2Result<EdgeId> {
        // Check for antipodal points which form invalid edges
        let dot_product = point_a.coords().dot(point_b.coords());
        if (dot_product + 1.0).abs() < 1e-15 {
            return Err(S2Error::BuilderError { 
                reason: "Cannot add edge between antipodal points".to_string() 
            });
        }
        
        let v1 = self.add_vertex(point_a)?;
        let v2 = self.add_vertex(point_b)?;
        let id = EdgeId(self.input_edges.len());
        self.input_edges.push((v1, v2, EdgeType::Directed));
        Ok(id)
    }

    /// Add an edge between two vertices
    pub fn add_edge_vertices(&mut self, v1: VertexId, v2: VertexId, edge_type: EdgeType) -> EdgeId {
        let id = EdgeId(self.input_edges.len());
        self.input_edges.push((v1, v2, edge_type));
        id
    }

    /// Add a polyline to the builder
    pub fn add_polyline(&mut self, polyline: &S2Polyline, edge_type: EdgeType) -> S2Result<()> {
        if polyline.num_vertices() < 2 {
            return Ok(());
        }

        let mut vertex_ids = Vec::new();
        for i in 0..polyline.num_vertices() {
            let vertex_id = self.add_vertex(polyline.vertex(i).expect("Invalid vertex index"))?;
            vertex_ids.push(vertex_id);
        }

        for i in 0..vertex_ids.len() - 1 {
            self.add_edge_vertices(vertex_ids[i], vertex_ids[i + 1], edge_type);
        }

        Ok(())
    }

    /// Add a loop to the builder
    pub fn add_loop(&mut self, loop_: &S2Loop, edge_type: EdgeType) -> S2Result<()> {
        if loop_.num_vertices() < 3 {
            return Err(S2Error::invalid_loop("Loop must have at least 3 vertices".to_string()));
        }

        let mut vertex_ids = Vec::new();
        for i in 0..loop_.num_vertices() {
            let vertex_id = self.add_vertex(loop_.vertex(i))?;
            vertex_ids.push(vertex_id);
        }

        // Connect all vertices in a loop
        for i in 0..vertex_ids.len() {
            let next_i = (i + 1) % vertex_ids.len();
            self.add_edge_vertices(vertex_ids[i], vertex_ids[next_i], edge_type);
        }

        Ok(())
    }

    /// Add a layer for output processing
    pub fn add_layer(&mut self, layer: Box<dyn Layer>) {
        self.layers.push(layer);
    }

    /// Get number of input vertices
    pub fn num_input_vertices(&self) -> usize {
        self.input_vertices.len()
    }

    /// Get number of input edges
    pub fn num_input_edges(&self) -> usize {
        self.input_edges.len()
    }

    /// Get number of layers
    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }

    /// Check if builder has been built
    pub fn is_built(&self) -> bool {
        self.built
    }

    /// Set build state (for testing)
    #[cfg(test)]
    pub fn set_built_state(&mut self, built: bool) {
        self.built = built;
    }

    /// Build the geometry and process all layers
    pub fn build(&mut self) -> S2Result<()> {
        self.built = true;
        // Create graph from input
        let mut graph = Graph::new();
        
        // Add vertices with snapping
        let mut vertex_map = HashMap::new();
        for (i, vertex) in self.input_vertices.iter().enumerate() {
            let snapped = self.options.snap_function.snap_point(*vertex);
            let graph_vertex_id = graph.add_vertex(snapped);
            vertex_map.insert(VertexId(i), graph_vertex_id);
        }

        // Add edges
        for (v1, v2, edge_type) in &self.input_edges {
            let graph_v1 = vertex_map[v1];
            let graph_v2 = vertex_map[v2];
            graph.add_edge(graph_v1, graph_v2, *edge_type);
        }

        // Process each layer
        for layer in &mut self.layers {
            layer.build(&graph)?;
        }

        Ok(())
    }
}

impl Default for S2Builder {
    fn default() -> Self {
        Self::new()
    }
}