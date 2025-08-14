//! Output Layers for S2Builder
//!
//! Layers define how the snapped graph is converted into specific output
//! geometry types. Each layer is responsible for assembling edges from
//! the internal graph into the desired output format while handling
//! edge types, validation, and error recovery.
//!
//! # Core Concepts
//!
//! - **Layer Trait**: Abstract interface for all output layers
//! - **Graph Processing**: Layers work with the snapped graph representation
//! - **Edge Types**: Layers handle directed/undirected edges appropriately
//! - **Validation**: Layers can validate their output geometry
//! - **Error Handling**: Robust error reporting for invalid configurations
//!
//! # Available Layers
//!
//! - [`S2PolygonLayer`]: Builds S2Polygon from closed loops
//! - [`S2PolylineLayer`]: Builds single S2Polyline from edge chain
//! - [`S2PolylineVectorLayer`]: Builds multiple S2Polylines
//! - [`IndexedS2PolylineVectorLayer`]: Builds polylines with source tracking
//!
//! # Example
//!
//! ```rust,ignore
//! use s2geometry_rust::builder::layers::*;
//! use s2geometry_rust::{S2Polygon, S2Builder};
//!
//! let mut polygon = S2Polygon::empty();
//! let layer = S2PolygonLayer::new(&mut polygon);
//! builder.add_layer(Box::new(layer));
//! builder.build()?;
//! ```

use crate::error::{S2Error, S2Result};
use crate::builder::graph::Graph;
use crate::builder::{EdgeType, EdgeId, VertexId};
use crate::polyline::S2Polyline;
use crate::r#loop::S2Loop;
use crate::point::S2Point;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Debug;

/// Abstract trait for S2Builder output layers
///
/// Layers are responsible for converting the snapped graph into specific
/// output geometry types. Each layer defines its own requirements for
/// edge types, topology, and validation.
pub trait Layer: Debug + Send {
    /// Build output geometry from the snapped graph
    ///
    /// This is the main method where layers process the graph edges
    /// and construct their specific output geometry.
    ///
    /// # Arguments
    /// * `graph` - The snapped graph containing edges and vertices
    ///
    /// # Returns
    /// * `Ok(())` if the layer was built successfully
    /// * `Err(S2Error)` if there was an error during construction
    fn build(&mut self, graph: &Graph) -> S2Result<()>;

    /// Get the graph options required by this layer
    ///
    /// Layers can specify requirements for graph processing such as
    /// whether to remove duplicate edges or sibling pairs.
    fn graph_options(&self) -> LayerGraphOptions {
        LayerGraphOptions::default()
    }
}

/// Options for graph processing by layers
#[derive(Debug, Clone)]
pub struct LayerGraphOptions {
    /// Required edge type for this layer
    pub edge_type: EdgeType,
    /// Whether degenerate edges should be removed
    pub remove_degenerate_edges: bool,
    /// Whether duplicate edges should be removed
    pub remove_duplicate_edges: bool,
    /// Whether sibling pairs should be removed
    pub remove_sibling_pairs: bool,
}

impl Default for LayerGraphOptions {
    fn default() -> Self {
        Self {
            edge_type: EdgeType::Directed,
            remove_degenerate_edges: true,
            remove_duplicate_edges: true,
            remove_sibling_pairs: false,
        }
    }
}

/// Configuration options for polygon layers
#[derive(Debug, Clone)]
pub struct PolygonLayerOptions {
    /// Edge type (directed or undirected)
    edge_type: EdgeType,
    /// Whether to validate the output polygon
    validate: bool,
}

impl Default for PolygonLayerOptions {
    fn default() -> Self {
        Self {
            edge_type: EdgeType::Undirected,
            validate: true,
        }
    }
}

impl PolygonLayerOptions {
    /// Create new options with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the edge type
    pub fn with_edge_type(mut self, edge_type: EdgeType) -> Self {
        self.edge_type = edge_type;
        self
    }

    /// Set whether to validate output
    pub fn with_validate(mut self, validate: bool) -> Self {
        self.validate = validate;
        self
    }

    /// Get the edge type
    pub fn edge_type(&self) -> EdgeType {
        self.edge_type
    }

    /// Get whether validation is enabled
    pub fn validate(&self) -> bool {
        self.validate
    }
}

/// Layer for building an S2Polygon from closed loops
///
/// This layer assembles graph edges into closed loops that form
/// the boundary of a spherical polygon. It handles multiple loops
/// and can determine their nesting relationships.
#[derive(Debug)]
pub struct S2PolygonLayer<'a> {
    /// Output polygon to build into
    polygon: &'a mut Vec<S2Loop>,
    /// Configuration options
    options: PolygonLayerOptions,
    /// Optional label tracking
    label_set_ids: Option<&'a mut Vec<u32>>,
}

impl<'a> S2PolygonLayer<'a> {
    /// Create a new polygon layer
    ///
    /// # Arguments
    /// * `polygon` - Mutable reference to output polygon (as vector of loops)
    /// * `options` - Configuration options
    pub fn new(polygon: &'a mut Vec<S2Loop>, options: PolygonLayerOptions) -> Self {
        Self {
            polygon,
            options,
            label_set_ids: None,
        }
    }

    /// Create a polygon layer with default options
    pub fn with_defaults(polygon: &'a mut Vec<S2Loop>) -> Self {
        Self::new(polygon, PolygonLayerOptions::default())
    }

    /// Create a polygon layer with label tracking
    pub fn with_labels(
        polygon: &'a mut Vec<S2Loop>,
        label_set_ids: &'a mut Vec<u32>,
        options: PolygonLayerOptions,
    ) -> Self {
        Self {
            polygon,
            options,
            label_set_ids: Some(label_set_ids),
        }
    }
}

impl<'a> Layer for S2PolygonLayer<'a> {
    fn build(&mut self, graph: &Graph) -> S2Result<()> {
        // Find all cycles in the graph
        let cycles = self.find_cycles(graph)?;
        
        // Convert cycles to S2Loops
        self.polygon.clear();
        for cycle in cycles {
            let loop_vertices = self.cycle_to_vertices(graph, &cycle)?;
            if loop_vertices.len() >= 3 {
                let s2_loop = S2Loop::new(loop_vertices)?;
                if !self.options.validate || s2_loop.is_valid() {
                    self.polygon.push(s2_loop);
                }
            }
        }

        Ok(())
    }

    fn graph_options(&self) -> LayerGraphOptions {
        LayerGraphOptions {
            edge_type: self.options.edge_type,
            remove_degenerate_edges: true,
            remove_duplicate_edges: true,
            remove_sibling_pairs: self.options.edge_type == EdgeType::Undirected,
        }
    }
}

impl<'a> S2PolygonLayer<'a> {
    /// Find all cycles in the graph
    fn find_cycles(&self, graph: &Graph) -> S2Result<Vec<Vec<EdgeId>>> {
        let mut cycles = Vec::new();
        let mut visited_edges = HashSet::new();

        for edge_id in graph.edges() {
            if visited_edges.contains(&edge_id) {
                continue;
            }

            if let Some(cycle) = self.find_cycle_from_edge(graph, edge_id, &mut visited_edges)? {
                cycles.push(cycle);
            }
        }

        Ok(cycles)
    }

    /// Find a cycle starting from a specific edge
    fn find_cycle_from_edge(
        &self,
        graph: &Graph,
        start_edge: EdgeId,
        visited_edges: &mut HashSet<EdgeId>,
    ) -> S2Result<Option<Vec<EdgeId>>> {
        if visited_edges.contains(&start_edge) {
            return Ok(None);
        }

        let mut cycle = Vec::new();
        let mut current_edge = start_edge;
        let start_vertex = graph.edge(start_edge)
            .ok_or_else(|| S2Error::BuilderError {
                reason: "Invalid edge ID in graph".to_string(),
            })?
            .source();

        loop {
            if visited_edges.contains(&current_edge) {
                break;
            }

            visited_edges.insert(current_edge);
            cycle.push(current_edge);

            // Find next edge in the cycle
            let current_target = graph.edge(current_edge)
                .ok_or_else(|| S2Error::BuilderError {
                    reason: "Invalid edge ID in graph".to_string(),
                })?
                .target();

            // Look for outgoing edge from current target
            let next_edge = self.find_next_edge(graph, current_target, current_edge)?;
            
            match next_edge {
                Some(edge_id) => {
                    current_edge = edge_id;
                    // Check if we've completed a cycle
                    let next_target = graph.edge(edge_id)
                        .ok_or_else(|| S2Error::BuilderError {
                            reason: "Invalid edge ID in graph".to_string(),
                        })?
                        .target();
                    if next_target == start_vertex {
                        visited_edges.insert(edge_id);
                        cycle.push(edge_id);
                        break;
                    }
                }
                None => {
                    // Dead end - not a cycle
                    return Ok(None);
                }
            }
        }

        if cycle.len() >= 3 {
            Ok(Some(cycle))
        } else {
            Ok(None)
        }
    }

    /// Find the next edge in a cycle from a given vertex
    fn find_next_edge(
        &self,
        graph: &Graph,
        vertex_id: VertexId,
        previous_edge: EdgeId,
    ) -> S2Result<Option<EdgeId>> {
        if let Some(vertex) = graph.vertex(vertex_id) {
            // For undirected edges, we need to find an edge that isn't the one we came from
            let mut candidates = Vec::new();
            
            // Add outgoing edges
            for &edge_id in vertex.outgoing_edges() {
                if edge_id != previous_edge {
                    candidates.push(edge_id);
                }
            }

            // For undirected edges, also consider incoming edges as potential next edges
            if self.options.edge_type == EdgeType::Undirected {
                for &edge_id in vertex.incoming_edges() {
                    if edge_id != previous_edge {
                        candidates.push(edge_id);
                    }
                }
            }

            // Choose the first available candidate
            // TODO: For better polygon construction, we should choose based on angles
            Ok(candidates.into_iter().next())
        } else {
            Ok(None)
        }
    }

    /// Convert a cycle of edge IDs to vertex coordinates
    fn cycle_to_vertices(&self, graph: &Graph, cycle: &[EdgeId]) -> S2Result<Vec<S2Point>> {
        let mut vertices = Vec::new();

        for &edge_id in cycle {
            if let Some(edge) = graph.edge(edge_id) {
                if let Some(source_vertex) = graph.vertex(edge.source()) {
                    vertices.push(source_vertex.point());
                }
            }
        }

        Ok(vertices)
    }
}

/// Configuration options for polyline layers
#[derive(Debug, Clone)]
pub struct PolylineLayerOptions {
    /// Edge type (directed or undirected)
    edge_type: EdgeType,
    /// Whether to validate the output polyline
    validate: bool,
}

impl Default for PolylineLayerOptions {
    fn default() -> Self {
        Self {
            edge_type: EdgeType::Directed,
            validate: true,
        }
    }
}

impl PolylineLayerOptions {
    /// Create new options with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the edge type
    pub fn with_edge_type(mut self, edge_type: EdgeType) -> Self {
        self.edge_type = edge_type;
        self
    }

    /// Set whether to validate output
    pub fn with_validate(mut self, validate: bool) -> Self {
        self.validate = validate;
        self
    }

    /// Get the edge type
    pub fn edge_type(&self) -> EdgeType {
        self.edge_type
    }

    /// Get whether validation is enabled
    pub fn validate(&self) -> bool {
        self.validate
    }
}

/// Layer for building a single S2Polyline from a chain of edges
///
/// This layer assembles graph edges into a single connected polyline.
/// It requires that the edges form a single unbroken chain.
#[derive(Debug)]
pub struct S2PolylineLayer<'a> {
    /// Output polyline to build into
    polyline: &'a mut Option<S2Polyline>,
    /// Configuration options
    options: PolylineLayerOptions,
}

impl<'a> S2PolylineLayer<'a> {
    /// Create a new polyline layer
    pub fn new(polyline: &'a mut Option<S2Polyline>, options: PolylineLayerOptions) -> Self {
        Self {
            polyline,
            options,
        }
    }

    /// Create a polyline layer with default options
    pub fn with_defaults(polyline: &'a mut Option<S2Polyline>) -> Self {
        Self::new(polyline, PolylineLayerOptions::default())
    }
}

impl<'a> Layer for S2PolylineLayer<'a> {
    fn build(&mut self, graph: &Graph) -> S2Result<()> {
        // Find a single path through all edges
        let path = self.find_edge_path(graph)?;
        
        if path.is_empty() {
            *self.polyline = None;
            return Ok(());
        }

        // Convert path to vertices
        let vertices = self.path_to_vertices(graph, &path)?;
        
        if vertices.len() >= 2 {
            let s2_polyline = S2Polyline::new(vertices)?;
            if !self.options.validate || s2_polyline.is_valid() {
                *self.polyline = Some(s2_polyline);
            } else {
                *self.polyline = None;
            }
        } else {
            *self.polyline = None;
        }

        Ok(())
    }

    fn graph_options(&self) -> LayerGraphOptions {
        LayerGraphOptions {
            edge_type: self.options.edge_type,
            remove_degenerate_edges: true,
            remove_duplicate_edges: true,
            remove_sibling_pairs: false,
        }
    }
}

impl<'a> S2PolylineLayer<'a> {
    /// Find a path through all edges in the graph
    fn find_edge_path(&self, graph: &Graph) -> S2Result<Vec<EdgeId>> {
        let edge_ids: Vec<EdgeId> = graph.edges().collect();
        
        if edge_ids.is_empty() {
            return Ok(Vec::new());
        }

        // Start from the first edge and try to build a chain
        let mut path = Vec::new();
        let mut remaining_edges: HashSet<EdgeId> = edge_ids.iter().copied().collect();
        
        // Find a starting edge (one with no incoming edges, or just pick first)
        let start_edge = self.find_start_edge(graph, &edge_ids)?;
        
        let mut current_edge = start_edge;
        
        while remaining_edges.contains(&current_edge) {
            remaining_edges.remove(&current_edge);
            path.push(current_edge);
            
            // Find next edge
            if let Some(edge) = graph.edge(current_edge) {
                let next_vertex = edge.target();
                
                // Look for an edge starting from next_vertex
                if let Some(next_edge) = self.find_outgoing_edge(graph, next_vertex, &remaining_edges) {
                    current_edge = next_edge;
                } else {
                    break; // End of chain
                }
            } else {
                break;
            }
        }

        // Check if we used all edges (for a valid single polyline)
        if !remaining_edges.is_empty() {
            return Err(S2Error::BuilderError {
                reason: "Edges do not form a single connected polyline".to_string(),
            });
        }

        Ok(path)
    }

    /// Find a good starting edge for the polyline
    fn find_start_edge(&self, graph: &Graph, edge_ids: &[EdgeId]) -> S2Result<EdgeId> {
        // Look for an edge whose source vertex has no incoming edges
        for &edge_id in edge_ids {
            if let Some(edge) = graph.edge(edge_id) {
                if let Some(vertex) = graph.vertex(edge.source()) {
                    if vertex.incoming_edges().is_empty() {
                        return Ok(edge_id);
                    }
                }
            }
        }

        // If no clear start found, just use the first edge
        edge_ids.first().copied().ok_or_else(|| S2Error::BuilderError {
            reason: "No edges available".to_string(),
        })
    }

    /// Find an outgoing edge from a vertex that's in the remaining set
    fn find_outgoing_edge(&self, graph: &Graph, vertex_id: VertexId, remaining_edges: &HashSet<EdgeId>) -> Option<EdgeId> {
        if let Some(vertex) = graph.vertex(vertex_id) {
            for &edge_id in vertex.outgoing_edges() {
                if remaining_edges.contains(&edge_id) {
                    return Some(edge_id);
                }
            }
        }
        None
    }

    /// Convert an edge path to vertex coordinates
    fn path_to_vertices(&self, graph: &Graph, path: &[EdgeId]) -> S2Result<Vec<S2Point>> {
        let mut vertices = Vec::new();

        for &edge_id in path {
            if let Some(edge) = graph.edge(edge_id) {
                if let Some(source_vertex) = graph.vertex(edge.source()) {
                    vertices.push(source_vertex.point());
                }
            }
        }

        // Add the final target vertex
        if let Some(&last_edge_id) = path.last() {
            if let Some(edge) = graph.edge(last_edge_id) {
                if let Some(target_vertex) = graph.vertex(edge.target()) {
                    vertices.push(target_vertex.point());
                }
            }
        }

        Ok(vertices)
    }
}

/// Layer for building multiple S2Polylines from graph edges
///
/// This layer assembles graph edges into as few polylines as possible,
/// handling disconnected components and multiple edge chains.
#[derive(Debug)]
pub struct S2PolylineVectorLayer<'a> {
    /// Output vector of polylines
    polylines: &'a mut Vec<S2Polyline>,
    /// Configuration options
    options: PolylineLayerOptions,
}

impl<'a> S2PolylineVectorLayer<'a> {
    /// Create a new polyline vector layer
    pub fn new(polylines: &'a mut Vec<S2Polyline>, options: PolylineLayerOptions) -> Self {
        Self {
            polylines,
            options,
        }
    }

    /// Create a polyline vector layer with default options
    pub fn with_defaults(polylines: &'a mut Vec<S2Polyline>) -> Self {
        Self::new(polylines, PolylineLayerOptions::default())
    }
}

impl<'a> Layer for S2PolylineVectorLayer<'a> {
    fn build(&mut self, graph: &Graph) -> S2Result<()> {
        self.polylines.clear();
        
        let mut remaining_edges: HashSet<EdgeId> = graph.edges().collect();
        
        while !remaining_edges.is_empty() {
            // Find a connected component
            let component = self.extract_connected_component(graph, &mut remaining_edges)?;
            
            if !component.is_empty() {
                let vertices = self.path_to_vertices(graph, &component)?;
                if vertices.len() >= 2 {
                    let polyline = S2Polyline::new(vertices)?;
                    if !self.options.validate || polyline.is_valid() {
                        self.polylines.push(polyline);
                    }
                }
            }
        }

        Ok(())
    }

    fn graph_options(&self) -> LayerGraphOptions {
        LayerGraphOptions {
            edge_type: self.options.edge_type,
            remove_degenerate_edges: true,
            remove_duplicate_edges: true,
            remove_sibling_pairs: false,
        }
    }
}

impl<'a> S2PolylineVectorLayer<'a> {
    /// Extract a connected component from the remaining edges
    fn extract_connected_component(&self, graph: &Graph, remaining_edges: &mut HashSet<EdgeId>) -> S2Result<Vec<EdgeId>> {
        if remaining_edges.is_empty() {
            return Ok(Vec::new());
        }

        // Start with any remaining edge
        let start_edge = *remaining_edges.iter().next().unwrap();
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();

        queue.push_back(start_edge);

        while let Some(edge_id) = queue.pop_front() {
            if visited.contains(&edge_id) || !remaining_edges.contains(&edge_id) {
                continue;
            }

            visited.insert(edge_id);
            remaining_edges.remove(&edge_id);
            component.push(edge_id);

            // Add connected edges to queue
            if let Some(edge) = graph.edge(edge_id) {
                // Add edges from source vertex
                if let Some(vertex) = graph.vertex(edge.source()) {
                    for &connected_edge in vertex.outgoing_edges() {
                        if remaining_edges.contains(&connected_edge) {
                            queue.push_back(connected_edge);
                        }
                    }
                    for &connected_edge in vertex.incoming_edges() {
                        if remaining_edges.contains(&connected_edge) {
                            queue.push_back(connected_edge);
                        }
                    }
                }

                // Add edges from target vertex
                if let Some(vertex) = graph.vertex(edge.target()) {
                    for &connected_edge in vertex.outgoing_edges() {
                        if remaining_edges.contains(&connected_edge) {
                            queue.push_back(connected_edge);
                        }
                    }
                    for &connected_edge in vertex.incoming_edges() {
                        if remaining_edges.contains(&connected_edge) {
                            queue.push_back(connected_edge);
                        }
                    }
                }
            }
        }

        // Sort component to create a reasonable path
        Ok(self.sort_edges_into_path(graph, component)?)
    }

    /// Sort edges into a reasonable path order
    fn sort_edges_into_path(&self, graph: &Graph, mut edges: Vec<EdgeId>) -> S2Result<Vec<EdgeId>> {
        if edges.is_empty() {
            return Ok(edges);
        }

        let mut path = Vec::new();
        let mut remaining: HashSet<EdgeId> = edges.iter().copied().collect();

        // Find a starting edge (prefer one with no incoming connections)
        let start_edge = self.find_path_start(graph, &edges)?;
        let mut current_edge = start_edge;

        while remaining.contains(&current_edge) {
            remaining.remove(&current_edge);
            path.push(current_edge);

            // Find next connected edge
            if let Some(next_edge) = self.find_next_connected_edge(graph, current_edge, &remaining) {
                current_edge = next_edge;
            } else {
                // No more connected edges, pick any remaining edge if available
                if let Some(&next_edge) = remaining.iter().next() {
                    current_edge = next_edge;
                } else {
                    break;
                }
            }
        }

        Ok(path)
    }

    /// Find a good starting edge for a path
    fn find_path_start(&self, graph: &Graph, edges: &[EdgeId]) -> S2Result<EdgeId> {
        // Look for an edge whose source has degree 1 (end of a path)
        for &edge_id in edges {
            if let Some(edge) = graph.edge(edge_id) {
                if let Some(vertex) = graph.vertex(edge.source()) {
                    if vertex.degree() == 1 {
                        return Ok(edge_id);
                    }
                }
            }
        }

        // Just use the first edge
        edges.first().copied().ok_or_else(|| S2Error::BuilderError {
            reason: "No edges available for path".to_string(),
        })
    }

    /// Find the next edge connected to the current edge
    fn find_next_connected_edge(&self, graph: &Graph, current_edge: EdgeId, remaining: &HashSet<EdgeId>) -> Option<EdgeId> {
        if let Some(edge) = graph.edge(current_edge) {
            let target_vertex = edge.target();
            
            // Look for an edge starting from the target vertex
            if let Some(vertex) = graph.vertex(target_vertex) {
                for &candidate_edge in vertex.outgoing_edges() {
                    if remaining.contains(&candidate_edge) {
                        return Some(candidate_edge);
                    }
                }
            }
        }

        None
    }

    /// Convert an edge path to vertices (same as S2PolylineLayer)
    fn path_to_vertices(&self, graph: &Graph, path: &[EdgeId]) -> S2Result<Vec<S2Point>> {
        let mut vertices = Vec::new();

        for &edge_id in path {
            if let Some(edge) = graph.edge(edge_id) {
                if let Some(source_vertex) = graph.vertex(edge.source()) {
                    vertices.push(source_vertex.point());
                }
            }
        }

        // Add the final target vertex
        if let Some(&last_edge_id) = path.last() {
            if let Some(edge) = graph.edge(last_edge_id) {
                if let Some(target_vertex) = graph.vertex(edge.target()) {
                    vertices.push(target_vertex.point());
                }
            }
        }

        Ok(vertices)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::graph::*;
    use crate::math::DVec3;

    #[test]
    fn test_polygon_layer_creation() {
        let mut polygon = Vec::new();
        let options = PolygonLayerOptions::default();
        let layer = S2PolygonLayer::new(&mut polygon, options);
        
        assert_eq!(layer.graph_options().edge_type, EdgeType::Undirected);
    }

    #[test]
    fn test_polyline_layer_creation() {
        let mut polyline = None;
        let options = PolylineLayerOptions::default();
        let layer = S2PolylineLayer::new(&mut polyline, options);
        
        assert_eq!(layer.graph_options().edge_type, EdgeType::Directed);
    }

    #[test]
    fn test_polyline_vector_layer_creation() {
        let mut polylines = Vec::new();
        let options = PolylineLayerOptions::default();
        let layer = S2PolylineVectorLayer::new(&mut polylines, options);
        
        assert_eq!(layer.graph_options().edge_type, EdgeType::Directed);
    }

    // TODO: Add integration tests with actual graph building
}