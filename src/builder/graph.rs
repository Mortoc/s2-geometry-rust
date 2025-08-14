//! Internal Graph Representation for S2Builder
//!
//! The Graph module provides the internal representation used by S2Builder
//! to manage snapped edges and vertices during the construction process.
//! This graph maintains topological relationships and provides efficient
//! access to geometric information needed by output layers.
//!
//! # Core Components
//!
//! - [`Graph`]: Main graph structure with vertices and edges
//! - [`Vertex`]: Graph vertex with position and adjacency information
//! - [`Edge`]: Graph edge with endpoints and properties
//! - [`GraphOptions`]: Configuration for graph construction
//!
//! # Key Features
//!
//! - **Topology Preservation**: Maintains proper edge-vertex relationships
//! - **Efficient Queries**: Fast access to vertex neighbors and incident edges
//! - **Label Tracking**: Support for tracking edge sources through labels
//! - **Validation**: Comprehensive checking of graph consistency
//!
//! # Example
//!
//! ```rust,ignore
//! use s2geometry_rust::builder::graph::*;
//!
//! let options = GraphOptions::default();
//! let graph = Graph::from_snapped_edges(snapped_edges, &sites)?;
//! 
//! // Query graph structure
//! for vertex_id in graph.vertices() {
//!     let vertex = graph.vertex(vertex_id);
//!     println!("Vertex {} at {:?}", vertex_id.0, vertex.point());
//! }
//! ```

use crate::error::{S2Error, S2Result};
use crate::point::S2Point;
use crate::builder::{EdgeId, VertexId, EdgeType};
use std::collections::{HashMap, HashSet};

/// Configuration options for graph construction
///
/// Controls how the graph is built from snapped edges and what
/// validation and processing steps are performed.
#[derive(Debug, Clone)]
pub struct GraphOptions {
    /// Whether to remove degenerate edges (zero length)
    remove_degenerate_edges: bool,
    /// Whether to remove duplicate edges  
    remove_duplicate_edges: bool,
    /// Whether to remove sibling pairs (opposite directed edges)
    remove_sibling_pairs: bool,
    /// Whether to validate graph consistency
    validate: bool,
}

impl Default for GraphOptions {
    fn default() -> Self {
        Self {
            remove_degenerate_edges: true,
            remove_duplicate_edges: true,
            remove_sibling_pairs: false, // Usually handled by layers
            validate: true,
        }
    }
}

impl GraphOptions {
    /// Create new graph options with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set whether to remove degenerate edges
    pub fn with_remove_degenerate_edges(mut self, remove: bool) -> Self {
        self.remove_degenerate_edges = remove;
        self
    }

    /// Set whether to remove duplicate edges
    pub fn with_remove_duplicate_edges(mut self, remove: bool) -> Self {
        self.remove_duplicate_edges = remove;
        self
    }

    /// Set whether to remove sibling pairs
    pub fn with_remove_sibling_pairs(mut self, remove: bool) -> Self {
        self.remove_sibling_pairs = remove;
        self
    }

    /// Set whether to validate the graph
    pub fn with_validate(mut self, validate: bool) -> Self {
        self.validate = validate;
        self
    }
}

/// Internal representation of a snapped edge for graph construction
#[derive(Debug, Clone)]
pub struct SnappedEdge {
    /// Source vertex
    pub source: S2Point,
    /// Target vertex
    pub target: S2Point,
    /// Edge type (directed or undirected)
    pub edge_type: EdgeType,
    /// Optional label for tracking edge sources
    pub label: Option<u32>,
}

/// Graph vertex with position and adjacency information
#[derive(Debug, Clone)]
pub struct Vertex {
    /// Position of the vertex on the sphere
    point: S2Point,
    /// IDs of outgoing edges from this vertex
    outgoing_edges: Vec<EdgeId>,
    /// IDs of incoming edges to this vertex
    incoming_edges: Vec<EdgeId>,
}

impl Vertex {
    /// Create a new vertex at the given point
    pub fn new(point: S2Point) -> Self {
        Self {
            point,
            outgoing_edges: Vec::new(),
            incoming_edges: Vec::new(),
        }
    }

    /// Get the vertex position
    pub fn point(&self) -> S2Point {
        self.point
    }

    /// Get outgoing edge IDs
    pub fn outgoing_edges(&self) -> &[EdgeId] {
        &self.outgoing_edges
    }

    /// Get incoming edge IDs
    pub fn incoming_edges(&self) -> &[EdgeId] {
        &self.incoming_edges
    }

    /// Get the degree (total number of incident edges)
    pub fn degree(&self) -> usize {
        self.outgoing_edges.len() + self.incoming_edges.len()
    }

    /// Add an outgoing edge
    fn add_outgoing_edge(&mut self, edge_id: EdgeId) {
        self.outgoing_edges.push(edge_id);
    }

    /// Add an incoming edge
    fn add_incoming_edge(&mut self, edge_id: EdgeId) {
        self.incoming_edges.push(edge_id);
    }
}

/// Graph edge with endpoints and properties
#[derive(Debug, Clone)]
pub struct Edge {
    /// Source vertex ID
    source: VertexId,
    /// Target vertex ID
    target: VertexId,
    /// Edge type (directed or undirected)
    edge_type: EdgeType,
    /// Optional label for tracking edge sources
    label: Option<u32>,
}

impl Edge {
    /// Create a new edge
    pub fn new(source: VertexId, target: VertexId, edge_type: EdgeType, label: Option<u32>) -> Self {
        Self {
            source,
            target,
            edge_type,
            label,
        }
    }

    /// Get the source vertex ID
    pub fn source(&self) -> VertexId {
        self.source
    }

    /// Get the target vertex ID
    pub fn target(&self) -> VertexId {
        self.target
    }

    /// Get the edge type
    pub fn edge_type(&self) -> EdgeType {
        self.edge_type
    }

    /// Get the edge label
    pub fn label(&self) -> Option<u32> {
        self.label
    }

    /// Check if this edge is directed
    pub fn is_directed(&self) -> bool {
        self.edge_type == EdgeType::Directed
    }

    /// Check if this edge is undirected
    pub fn is_undirected(&self) -> bool {
        self.edge_type == EdgeType::Undirected
    }

    /// Get the reversed edge (target becomes source)
    pub fn reversed(&self) -> Self {
        Self {
            source: self.target,
            target: self.source,
            edge_type: self.edge_type,
            label: self.label,
        }
    }
}

/// Main graph structure for S2Builder
///
/// The Graph maintains the topological structure of snapped edges and vertices,
/// providing efficient access to geometric and topological information needed
/// by output layers.
#[derive(Debug)]
pub struct Graph {
    /// All vertices in the graph
    vertices: HashMap<VertexId, Vertex>,
    /// All edges in the graph
    edges: HashMap<EdgeId, Edge>,
    /// Map from vertex positions to vertex IDs (for deduplication)
    position_to_vertex: HashMap<PositionKey, VertexId>,
    /// Next available vertex ID
    next_vertex_id: u32,
    /// Next available edge ID
    next_edge_id: u32,
    /// Graph options used for construction
    options: GraphOptions,
}

/// Hash key for vertex positions (to handle floating-point precision)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PositionKey {
    x_bits: u64,
    y_bits: u64,
    z_bits: u64,
}

impl PositionKey {
    /// Create a position key from an S2Point
    fn from_point(point: S2Point) -> Self {
        let coords = point.coords();
        Self {
            x_bits: coords.x.to_bits(),
            y_bits: coords.y.to_bits(),
            z_bits: coords.z.to_bits(),
        }
    }
}

impl Graph {
    /// Create a new empty graph with default options
    pub fn new() -> Self {
        Self::with_options(GraphOptions::default())
    }

    /// Create a new empty graph with custom options
    pub fn with_options(options: GraphOptions) -> Self {
        Self {
            vertices: HashMap::new(),
            edges: HashMap::new(),
            position_to_vertex: HashMap::new(),
            next_vertex_id: 0,
            next_edge_id: 0,
            options,
        }
    }

    /// Create a graph from snapped edges
    pub fn from_snapped_edges(snapped_edges: Vec<SnappedEdge>, _sites: &[S2Point]) -> S2Result<Self> {
        let options = GraphOptions::default();
        let mut graph = Self::with_options(options);

        // Add all edges to the graph
        for snapped_edge in snapped_edges {
            graph.add_snapped_edge(snapped_edge)?;
        }

        // Apply graph processing options
        graph.process_graph()?;

        if graph.options.validate {
            graph.validate()?;
        }

        Ok(graph)
    }

    /// Get the number of vertices
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get the number of edges
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Get all vertex IDs
    pub fn vertices(&self) -> impl Iterator<Item = VertexId> + '_ {
        self.vertices.keys().copied()
    }

    /// Get all edge IDs
    pub fn edges(&self) -> impl Iterator<Item = EdgeId> + '_ {
        self.edges.keys().copied()
    }

    /// Get a vertex by ID
    pub fn vertex(&self, vertex_id: VertexId) -> Option<&Vertex> {
        self.vertices.get(&vertex_id)
    }

    /// Get an edge by ID
    pub fn edge(&self, edge_id: EdgeId) -> Option<&Edge> {
        self.edges.get(&edge_id)
    }

    /// Add a vertex to the graph and return its ID
    /// 
    /// This is a convenience method that creates a vertex at the given position.
    /// If a vertex already exists at this position, returns the existing vertex ID.
    pub fn add_vertex(&mut self, point: S2Point) -> VertexId {
        self.find_or_create_vertex(point)
    }

    /// Add an edge to the graph between two vertices
    /// 
    /// # Arguments
    /// * `source` - Source vertex ID
    /// * `target` - Target vertex ID  
    /// * `edge_type` - Type of edge (directed or undirected)
    /// 
    /// # Returns
    /// * `EdgeId` of the newly created edge
    /// 
    /// # Panics
    /// Panics if either vertex ID is invalid
    pub fn add_edge(&mut self, source: VertexId, target: VertexId, edge_type: EdgeType) -> EdgeId {
        // Get the actual points for the vertices
        let source_point = self.vertices.get(&source)
            .expect("Invalid source vertex ID")
            .point();
        let target_point = self.vertices.get(&target)
            .expect("Invalid target vertex ID")
            .point();

        // Create a snapped edge and add it
        let snapped_edge = SnappedEdge {
            source: source_point,
            target: target_point,
            edge_type,
            label: None,
        };

        self.add_snapped_edge(snapped_edge)
            .expect("Failed to add edge")
    }

    /// Find or create a vertex at the given position
    pub fn find_or_create_vertex(&mut self, point: S2Point) -> VertexId {
        let key = PositionKey::from_point(point);
        
        if let Some(&vertex_id) = self.position_to_vertex.get(&key) {
            return vertex_id;
        }

        // Create new vertex
        let vertex_id = VertexId(self.next_vertex_id as usize);
        self.next_vertex_id += 1;

        let vertex = Vertex::new(point);
        self.vertices.insert(vertex_id, vertex);
        self.position_to_vertex.insert(key, vertex_id);

        vertex_id
    }

    /// Add a snapped edge to the graph
    fn add_snapped_edge(&mut self, snapped_edge: SnappedEdge) -> S2Result<EdgeId> {
        // Skip degenerate edges if configured
        if self.options.remove_degenerate_edges {
            let distance = snapped_edge.source.angle(&snapped_edge.target);
            if distance < 1e-15 {
                return Err(S2Error::BuilderError {
                    reason: "Degenerate edge removed".to_string(),
                });
            }
        }

        // Find or create vertices
        let source_id = self.find_or_create_vertex(snapped_edge.source);
        let target_id = self.find_or_create_vertex(snapped_edge.target);

        // Create edge
        let edge_id = EdgeId(self.next_edge_id as usize);
        self.next_edge_id += 1;

        let edge = Edge::new(source_id, target_id, snapped_edge.edge_type, snapped_edge.label);

        // Check for duplicates if configured
        if self.options.remove_duplicate_edges {
            for existing_edge in self.edges.values() {
                if existing_edge.source == edge.source 
                   && existing_edge.target == edge.target
                   && existing_edge.edge_type == edge.edge_type {
                    return Err(S2Error::BuilderError {
                        reason: "Duplicate edge removed".to_string(),
                    });
                }
            }
        }

        // Add edge to adjacency lists
        if let Some(source_vertex) = self.vertices.get_mut(&source_id) {
            source_vertex.add_outgoing_edge(edge_id);
        }
        if let Some(target_vertex) = self.vertices.get_mut(&target_id) {
            target_vertex.add_incoming_edge(edge_id);
        }

        // Store edge
        self.edges.insert(edge_id, edge);

        Ok(edge_id)
    }

    /// Apply graph processing options
    fn process_graph(&mut self) -> S2Result<()> {
        if self.options.remove_sibling_pairs {
            self.remove_sibling_pairs()?;
        }

        Ok(())
    }

    /// Remove sibling pairs (opposite directed edges between same vertices)
    fn remove_sibling_pairs(&mut self) -> S2Result<()> {
        let mut edges_to_remove = Vec::new();
        let edge_ids: Vec<EdgeId> = self.edges.keys().copied().collect();

        for &edge_id in &edge_ids {
            if let Some(edge) = self.edges.get(&edge_id) {
                // Look for the reverse edge
                for &other_edge_id in &edge_ids {
                    if edge_id != other_edge_id {
                        if let Some(other_edge) = self.edges.get(&other_edge_id) {
                            if edge.source == other_edge.target 
                               && edge.target == other_edge.source
                               && edge.edge_type == other_edge.edge_type {
                                // Found sibling pair
                                edges_to_remove.push(edge_id);
                                edges_to_remove.push(other_edge_id);
                            }
                        }
                    }
                }
            }
        }

        // Remove duplicate entries and remove edges
        edges_to_remove.sort();
        edges_to_remove.dedup();
        
        for edge_id in edges_to_remove {
            self.remove_edge(edge_id)?;
        }

        Ok(())
    }

    /// Remove an edge from the graph
    fn remove_edge(&mut self, edge_id: EdgeId) -> S2Result<()> {
        if let Some(edge) = self.edges.remove(&edge_id) {
            // Remove from vertex adjacency lists
            if let Some(source_vertex) = self.vertices.get_mut(&edge.source) {
                source_vertex.outgoing_edges.retain(|&id| id != edge_id);
            }
            if let Some(target_vertex) = self.vertices.get_mut(&edge.target) {
                target_vertex.incoming_edges.retain(|&id| id != edge_id);
            }
        }

        Ok(())
    }

    /// Validate graph consistency
    fn validate(&self) -> S2Result<()> {
        // Check that all edge endpoints reference valid vertices
        for (edge_id, edge) in &self.edges {
            if !self.vertices.contains_key(&edge.source) {
                return Err(S2Error::BuilderError {
                    reason: format!("Edge {:?} references invalid source vertex {:?}", edge_id, edge.source),
                });
            }
            if !self.vertices.contains_key(&edge.target) {
                return Err(S2Error::BuilderError {
                    reason: format!("Edge {:?} references invalid target vertex {:?}", edge_id, edge.target),
                });
            }
        }

        // Check that vertex adjacency lists are consistent
        for (vertex_id, vertex) in &self.vertices {
            for &edge_id in &vertex.outgoing_edges {
                if let Some(edge) = self.edges.get(&edge_id) {
                    if edge.source != *vertex_id {
                        return Err(S2Error::BuilderError {
                            reason: format!("Vertex {:?} claims outgoing edge {:?} but edge source is {:?}", 
                                           vertex_id, edge_id, edge.source),
                        });
                    }
                } else {
                    return Err(S2Error::BuilderError {
                        reason: format!("Vertex {:?} references non-existent outgoing edge {:?}", 
                                       vertex_id, edge_id),
                    });
                }
            }

            for &edge_id in &vertex.incoming_edges {
                if let Some(edge) = self.edges.get(&edge_id) {
                    if edge.target != *vertex_id {
                        return Err(S2Error::BuilderError {
                            reason: format!("Vertex {:?} claims incoming edge {:?} but edge target is {:?}", 
                                           vertex_id, edge_id, edge.target),
                        });
                    }
                } else {
                    return Err(S2Error::BuilderError {
                        reason: format!("Vertex {:?} references non-existent incoming edge {:?}", 
                                       vertex_id, edge_id),
                    });
                }
            }
        }

        Ok(())
    }

    /// Get all edges incident to a vertex (both incoming and outgoing)
    pub fn incident_edges(&self, vertex_id: VertexId) -> Vec<EdgeId> {
        if let Some(vertex) = self.vertices.get(&vertex_id) {
            let mut edges = vertex.incoming_edges.clone();
            edges.extend_from_slice(&vertex.outgoing_edges);
            edges
        } else {
            Vec::new()
        }
    }

    /// Get neighboring vertices of a given vertex
    pub fn neighboring_vertices(&self, vertex_id: VertexId) -> Vec<VertexId> {
        let mut neighbors = HashSet::new();
        
        if let Some(vertex) = self.vertices.get(&vertex_id) {
            // Add targets of outgoing edges
            for &edge_id in &vertex.outgoing_edges {
                if let Some(edge) = self.edges.get(&edge_id) {
                    neighbors.insert(edge.target);
                }
            }
            
            // Add sources of incoming edges
            for &edge_id in &vertex.incoming_edges {
                if let Some(edge) = self.edges.get(&edge_id) {
                    neighbors.insert(edge.source);
                }
            }
        }

        neighbors.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;

    #[test]
    fn test_empty_graph() {
        let graph = Graph::new();
        assert_eq!(graph.num_vertices(), 0);
        assert_eq!(graph.num_edges(), 0);
    }

    #[test]
    fn test_add_vertex() {
        let mut graph = Graph::new();
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        
        let vertex_id = graph.find_or_create_vertex(point);
        assert_eq!(graph.num_vertices(), 1);
        
        let vertex = graph.vertex(vertex_id).unwrap();
        assert_eq!(vertex.point(), point);
    }

    #[test]
    fn test_add_edge() {
        let mut graph = Graph::new();
        let point_a = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let point_b = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
        
        let snapped_edge = SnappedEdge {
            source: point_a,
            target: point_b,
            edge_type: EdgeType::Directed,
            label: None,
        };
        
        let edge_id = graph.add_snapped_edge(snapped_edge).unwrap();
        assert_eq!(graph.num_edges(), 1);
        assert_eq!(graph.num_vertices(), 2);
        
        let edge = graph.edge(edge_id).unwrap();
        assert_eq!(edge.edge_type(), EdgeType::Directed);
    }

    #[test]
    fn test_vertex_deduplication() {
        let mut graph = Graph::new();
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        
        let vertex_id1 = graph.find_or_create_vertex(point);
        let vertex_id2 = graph.find_or_create_vertex(point);
        
        assert_eq!(vertex_id1, vertex_id2);
        assert_eq!(graph.num_vertices(), 1);
    }

    #[test]
    fn test_edge_adjacency() {
        let mut graph = Graph::new();
        let point_a = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let point_b = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
        
        let source_id = graph.find_or_create_vertex(point_a);
        let target_id = graph.find_or_create_vertex(point_b);
        
        let snapped_edge = SnappedEdge {
            source: point_a,
            target: point_b,
            edge_type: EdgeType::Directed,
            label: None,
        };
        
        let edge_id = graph.add_snapped_edge(snapped_edge).unwrap();
        
        let source_vertex = graph.vertex(source_id).unwrap();
        let target_vertex = graph.vertex(target_id).unwrap();
        
        assert_eq!(source_vertex.outgoing_edges().len(), 1);
        assert_eq!(source_vertex.incoming_edges().len(), 0);
        assert_eq!(target_vertex.outgoing_edges().len(), 0);
        assert_eq!(target_vertex.incoming_edges().len(), 1);
        
        assert_eq!(source_vertex.outgoing_edges()[0], edge_id);
        assert_eq!(target_vertex.incoming_edges()[0], edge_id);
    }

    #[test]
    fn test_neighboring_vertices() {
        let mut graph = Graph::new();
        let point_a = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        let point_b = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
        let point_c = S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0));
        
        let vertex_a = graph.find_or_create_vertex(point_a);
        let vertex_b = graph.find_or_create_vertex(point_b);
        let vertex_c = graph.find_or_create_vertex(point_c);
        
        // Add edges A->B and A->C
        graph.add_snapped_edge(SnappedEdge {
            source: point_a,
            target: point_b,
            edge_type: EdgeType::Directed,
            label: None,
        }).unwrap();
        
        graph.add_snapped_edge(SnappedEdge {
            source: point_a,
            target: point_c,
            edge_type: EdgeType::Directed,
            label: None,
        }).unwrap();
        
        let neighbors = graph.neighboring_vertices(vertex_a);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&vertex_b));
        assert!(neighbors.contains(&vertex_c));
    }
}