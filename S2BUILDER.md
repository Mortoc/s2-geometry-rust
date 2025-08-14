# S2Builder Implementation Guide

This document describes the S2Builder implementation in the Rust S2 geometry library, providing a comprehensive toolkit for robust geometric construction from imperfect input data.

## Overview

S2Builder is a sophisticated geometry construction and repair toolkit that handles topology preservation, snap rounding, and building reliable geometric objects from imperfect input data. It serves as the foundation for complex geometric operations like polygon boolean operations, simplification, and format conversion.

## Architecture

The S2Builder system follows a multi-stage construction process:

1. **Input Collection**: Edges and vertices are collected and normalized
2. **Site Selection**: Choose snap sites using the configured snap function  
3. **Edge Snapping**: Snap input edges to the selected sites
4. **Graph Construction**: Build internal graph representation with topology
5. **Layer Processing**: Each layer builds its output from the snapped graph
6. **Validation**: Verify geometric consistency and report errors

## Core Components

### S2Builder

The main builder class that orchestrates the entire construction process.

```rust
use s2geometry_rust::builder::*;

// Create builder with default options
let mut builder = S2Builder::default();

// Add input geometry
builder.add_edge(point_a, point_b)?;
builder.add_loop(&vertices, EdgeType::Undirected)?;

// Add output layers
let mut polygon = Vec::new();
let layer = S2PolygonLayer::with_defaults(&mut polygon);
builder.add_layer(Box::new(layer))?;

// Build all outputs
builder.build()?;
```

### Options

Configuration for S2Builder operation, controlling snap rounding, topology, and validation.

```rust
let options = Options::new()
    .with_snap_function(S2CellIdSnapFunction::new(10))
    .with_split_crossing_edges(true)
    .with_intersection_tolerance(1e-15)
    .with_validate_output(true);

let builder = S2Builder::new(options);
```

### Snap Functions

Functions that determine how input vertices are mapped to discrete snap sites.

#### Available Snap Functions

- **IdentitySnapFunction**: No snapping (identity transformation)
- **IntLatLngSnapFunction**: Snap to integer lat/lng coordinates
- **S2CellIdSnapFunction**: Snap to S2CellId centers at a given level
- **MinEdgeLengthSnapFunction**: Enforces minimum edge length constraints

```rust
// Snap to level 10 S2Cell centers (roughly 600m resolution)
let snap_fn = S2CellIdSnapFunction::new(10);

// Snap to microdegree precision
let snap_fn = IntLatLngSnapFunction::microdegrees();

// No snapping with minimum separation
let snap_fn = IdentitySnapFunction::new(1e-15);
```

### Graph

Internal representation managing snapped edges and vertices with topological relationships.

```rust
// Graph is built automatically by S2Builder
let graph = builder.graph().unwrap(); // After build()

// Query graph structure
for vertex_id in graph.vertices() {
    let vertex = graph.vertex(vertex_id).unwrap();
    println!("Vertex at {:?}", vertex.point());
}
```

### Layers

Output layers define how the snapped graph is converted into specific geometry types.

#### S2PolygonLayer

Builds S2Polygon from closed loops.

```rust
let mut loops = Vec::new();
let options = PolygonLayerOptions::default()
    .with_edge_type(EdgeType::Undirected)
    .with_validate(true);
let layer = S2PolygonLayer::new(&mut loops, options);
```

#### S2PolylineLayer

Builds a single S2Polyline from edge chain.

```rust
let mut polyline = None;
let options = PolylineLayerOptions::default()
    .with_edge_type(EdgeType::Directed);
let layer = S2PolylineLayer::new(&mut polyline, options);
```

#### S2PolylineVectorLayer

Builds multiple S2Polylines from disconnected components.

```rust
let mut polylines = Vec::new();
let layer = S2PolylineVectorLayer::with_defaults(&mut polylines);
```

## Usage Examples

### Basic Triangle Construction

```rust
use s2geometry_rust::builder::*;
use s2geometry_rust::{S2Builder, S2Loop, EdgeType};

let mut builder = S2Builder::default();

// Define triangle vertices
let vertices = vec![
    S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
    S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
    S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
];

// Add as undirected loop
builder.add_loop(&vertices, EdgeType::Undirected)?;

// Set up polygon output
let mut loops = Vec::new();
let layer = S2PolygonLayer::with_defaults(&mut loops);
builder.add_layer(Box::new(layer))?;

// Build the polygon
builder.build()?;
```

### Robust Polyline Construction

```rust
// Use S2CellId snapping for numerical robustness
let snap_function = S2CellIdSnapFunction::new(15);
let options = Options::default()
    .with_snap_function(snap_function)
    .with_validate_input(true);

let mut builder = S2Builder::new(options);

// Add path vertices
let path = vec![
    S2LatLng::from_degrees(0.0, 0.0).to_point()?,
    S2LatLng::from_degrees(0.0, 30.0).to_point()?,
    S2LatLng::from_degrees(0.0, 60.0).to_point()?,
];

builder.add_polyline(&path, EdgeType::Directed)?;

// Build polyline with validation
let mut polyline = None;
let layer = S2PolylineLayer::with_defaults(&mut polyline);
builder.add_layer(Box::new(layer))?;
builder.build()?;
```

### Multiple Output Layers

```rust
let mut builder = S2Builder::default();

// Add complex input geometry
builder.add_loop(&square_vertices, EdgeType::Undirected)?;
builder.add_polyline(&diagonal_path, EdgeType::Directed)?;

// Multiple output layers
let mut polygon_loops = Vec::new();
let mut all_polylines = Vec::new();

builder.add_layer(Box::new(S2PolygonLayer::with_defaults(&mut polygon_loops)))?;
builder.add_layer(Box::new(S2PolylineVectorLayer::with_defaults(&mut all_polylines)))?;

builder.build()?;
// Now have both polygon and polyline outputs
```

## Key Features

### Numerical Robustness

- **Snap Rounding**: Ensures vertices are placed at discrete, well-separated locations
- **Error Bounds**: Conservative analysis of geometric precision limits
- **Exact Arithmetic**: Fallback for critical geometric predicates

### Topology Preservation

- **Edge Relationships**: Maintains proper vertex-edge adjacency
- **Loop Closure**: Ensures polygonal loops are properly closed
- **Intersection Handling**: Configurable edge intersection processing

### Flexible Configuration

- **Multiple Snap Functions**: Choose appropriate discretization strategy
- **Edge Types**: Support for directed and undirected edges
- **Validation Levels**: Configurable input/output validation
- **Performance Tuning**: Options for speed vs. robustness trade-offs

### Error Handling

- **Comprehensive Validation**: Input geometry validation with detailed error messages
- **Graceful Degradation**: Robust handling of degenerate cases
- **Error Recovery**: Continue processing when possible, report issues clearly

## Implementation Status

### Completed Components

- ✅ Core S2Builder architecture
- ✅ Snap function implementations (Identity, IntLatLng, S2CellId, MinEdgeLength)
- ✅ Graph construction and management
- ✅ Basic layer implementations (Polygon, Polyline, PolylineVector)
- ✅ Configuration options and error handling
- ✅ Comprehensive test suite
- ✅ Documentation and examples

### Future Enhancements

- 🔄 Advanced snap functions (custom grids, geographic projections)
- 🔄 Edge intersection splitting algorithms
- 🔄 Edge chain simplification
- 🔄 Label tracking and source edge mapping
- 🔄 Performance optimizations for large datasets
- 🔄 Additional layer types (indexed polylines, shape collections)

## Design Philosophy

The S2Builder implementation follows these principles:

1. **Safety First**: Leverage Rust's ownership system for memory safety
2. **Robustness**: Numerical precision issues are handled explicitly
3. **Flexibility**: Support diverse input types and output requirements
4. **Performance**: Efficient algorithms with careful memory management
5. **Usability**: Clear APIs with comprehensive error messages
6. **Compatibility**: Match Google's C++ S2 library behavior where possible

## Testing

The implementation includes comprehensive tests covering:

- Unit tests for all components
- Integration tests for complete workflows
- Error condition testing
- Performance benchmarks
- Compatibility tests with reference implementations

Run tests with:
```bash
cargo test test_s2builder
```

Run the demo example:
```bash
cargo run --example s2builder_demo
```

## Performance Considerations

- **Memory Usage**: Efficient graph representation with minimal overhead
- **Algorithmic Complexity**: O(n log n) vertex processing, O(m) edge processing
- **Cache Efficiency**: Spatial locality through S2CellId ordering
- **Parallelization**: Thread-safe design enables parallel layer processing

## Integration with S2 Ecosystem

S2Builder integrates seamlessly with other S2 geometry components:

- **Input Sources**: S2Point, S2LatLng, S2Polyline, S2Loop
- **Output Targets**: S2Polygon, S2Polyline, custom geometry types
- **Spatial Indexing**: Compatible with S2ShapeIndex for query operations
- **Boolean Operations**: Foundation for union, intersection, difference operations

## Conclusion

The S2Builder implementation provides a robust, flexible foundation for geometric construction in Rust. It handles the complex challenges of numerical precision, topology preservation, and error recovery while maintaining high performance and usability.

This implementation enables reliable processing of real-world geographic data with the safety and performance benefits of Rust, making it suitable for applications ranging from GIS systems to computational geometry research.