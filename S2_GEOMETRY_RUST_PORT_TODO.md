# S2 Geometry Library - Rust Port TODO (TDD-Based)

This document outlines a comprehensive plan for porting the Google S2 Geometry library from C++ to Rust using strict Test-Driven Development (TDD) principles. Each major component follows the Red-Green-Refactor cycle.

## Project Structure Overview

Based on analysis of the C++ codebase at https://github.com/google/s2geometry, the library is organized into:

- **Core geometric types**: S2Point, S2CellId, S1Angle, S1ChordAngle
- **Fundamental utilities**: Vector math, exact arithmetic, bit operations
- **Regional types**: S2Region, S2Cap, S2Cell, S2Loop, S2Polygon 
- **Indexing structures**: S2CellIndex, S2ShapeIndex, MutableS2ShapeIndex
- **Query operations**: Distance queries, containment, intersection
- **Builder utilities**: S2Builder and various layer implementations
- **Advanced operations**: Boolean operations, buffering, validation

## Phase 1: Foundation and Core Utilities (Level 0 Dependencies)

### 1.1 Base System Types and Error Handling

#### 1.1.1 Error Handling Framework
- **RED**: Write failing tests for S2Error types covering:
  - Data loss errors, invalid input errors, computation failures
  - Error chaining and context preservation
  - Conversion from/to standard Rust error types
- **GREEN**: Implement minimal S2Error enum with thiserror derive
- **REFACTOR**: Add comprehensive error variants and documentation

#### 1.1.2 Basic Type Definitions and Constants
- **RED**: Write tests for fundamental constants and type aliases
  - Precision constants, geometric tolerances
  - Type aliases for common numeric types
- **GREEN**: Define basic constants and type system
- **REFACTOR**: Organize into appropriate modules with clear visibility

### 1.2 Bit Manipulation Utilities (`s2/util/bits/`)

#### 1.2.1 Core Bit Operations
- **RED**: Write tests for bit manipulation operations:
  - Bit counting (population count, leading/trailing zeros)
  - Bit interleaving operations for spatial indexing
  - Endianness handling
- **GREEN**: Implement using Rust's built-in bit operations and external crates
- **REFACTOR**: Optimize for performance, add SIMD where beneficial

### 1.3 Mathematical Utilities (`s2/util/math/`)

#### 1.3.1 Vector Math Foundation
- **RED**: Write comprehensive tests for Vector2, Vector3, Vector4:
  - Construction, component access, basic arithmetic
  - Dot product, cross product, normalization
  - Comparison operations with floating-point tolerance
  - Conversion between different numeric types
- **GREEN**: Implement generic Vector types with operator overloading via traits
- **REFACTOR**: Optimize for performance, ensure zero-cost abstractions

#### 1.3.2 ExactFloat Arithmetic
- **RED**: Write tests for exact arithmetic operations:
  - Construction from integers, rationals, floating-point
  - Addition, subtraction, multiplication (no division)
  - Comparison operations, sign determination
  - Conversion to/from standard floating-point types
- **GREEN**: Implement using BigInt crates (num-bigint) or similar
- **REFACTOR**: Optimize for geometric predicate use cases

#### 1.3.3 Mathematical Utilities
- **RED**: Write tests for mathutil functions:
  - Floating-point comparisons with tolerance
  - Angle normalization and wrapping
  - Min/max operations, clamping
- **GREEN**: Implement utility functions
- **REFACTOR**: Ensure numerical stability and correctness

### 1.4 Encoding/Decoding Framework (`s2/util/coding/`)

#### 1.4.1 Binary Encoding Support
- **RED**: Write tests for encoding/decoding primitives:
  - Variable-length integer encoding (varint)
  - Fixed-width type encoding/decoding
  - Buffer management and error handling
- **GREEN**: Implement Encoder/Decoder traits and basic implementations  
- **REFACTOR**: Optimize for performance, ensure zero-copy where possible

## Phase 2: Core Geometric Types (Level 1 Dependencies)

### 2.1 Fundamental Angle Types

#### 2.1.1 S1Angle Implementation
- **RED**: Write comprehensive tests for S1Angle:
  - Construction from radians, degrees, E5/E6/E7 representations
  - Conversion between representations with precision guarantees
  - Arithmetic operations (addition, subtraction, multiplication by scalar)
  - Comparison operations, normalization
  - Trigonometric operations where exact
- **GREEN**: Implement S1Angle as wrapper around f64 with conversion methods
- **REFACTOR**: Optimize conversions, add comprehensive documentation

#### 2.1.2 S1ChordAngle Implementation
- **RED**: Write tests for S1ChordAngle:
  - Construction from S1Angle, directly from squared chord length
  - Conversion to/from S1Angle with appropriate precision
  - Comparison operations, arithmetic where applicable
  - Special cases (zero, straight, infinity angles)
- **GREEN**: Implement S1ChordAngle with internal squared-chord representation
- **REFACTOR**: Optimize for distance calculations, ensure numerical stability

### 2.2 Core Point Representation

#### 2.2.1 S2Point Implementation
- **RED**: Write comprehensive tests for S2Point:
  - Construction from coordinates, conversion from Vector3
  - Normalization operations, magnitude calculations
  - Arithmetic operations (addition, subtraction, scalar multiplication)
  - Dot product, cross product operations
  - Comparison with floating-point tolerance
  - Serialization/deserialization
- **GREEN**: Implement S2Point as wrapper around Vector3<f64>
- **REFACTOR**: Ensure zero-cost conversion between S2Point and Vector3

### 2.3 Coordinate Systems and Projections

#### 2.3.1 S2Coords System
- **RED**: Write tests for coordinate transformations:
  - Cube face projections (6 faces of unit cube)
  - UV coordinate system on each face
  - ST coordinate normalization
  - Conversion between (face,u,v) and S2Point
  - Jacobian calculations for area computations
- **GREEN**: Implement coordinate transformation functions
- **REFACTOR**: Optimize for performance, ensure numerical accuracy

## Phase 3: Cell System and Hierarchical Indexing (Level 2 Dependencies)

### 3.1 S2CellId - Core Indexing Structure

#### 3.1.1 S2CellId Basic Operations
- **RED**: Write comprehensive tests for S2CellId:
  - Construction from S2Point, from face/level/position
  - Hilbert curve position encoding/decoding
  - Level operations (parent, children, neighbors)
  - Range queries and containment checks
  - Token string conversion for debugging
  - Serialization with compact representations
- **GREEN**: Implement S2CellId as wrapper around u64 with bit manipulation
- **REFACTOR**: Optimize bit operations, add extensive documentation

#### 3.1.2 S2CellId Advanced Features
- **RED**: Write tests for advanced S2CellId operations:
  - Range iterator functionality
  - Coverage operations for geometric shapes
  - Distance calculations between cells
  - Union and intersection of cell ranges
- **GREEN**: Implement advanced algorithms
- **REFACTOR**: Optimize for large-scale indexing operations

### 3.2 S2Cell - Geometric Cell Representation

#### 3.2.1 S2Cell Implementation
- **RED**: Write tests for S2Cell geometry:
  - Construction from S2CellId
  - Vertex computation and edge geometry
  - Center point and area calculations
  - Bounding rectangle computation
  - Containment and intersection tests
- **GREEN**: Implement S2Cell with lazy computation of geometric properties
- **REFACTOR**: Optimize geometric calculations, cache expensive operations

### 3.3 Interval Types

#### 3.3.1 R1Interval (1D Intervals)
- **RED**: Write tests for 1D interval operations:
  - Construction, bounds checking, emptiness tests
  - Union, intersection, containment operations
  - Expansion operations with margins
- **GREEN**: Implement R1Interval for angle/coordinate ranges
- **REFACTOR**: Ensure robust handling of edge cases

#### 3.3.2 S1Interval (Circular Intervals)
- **RED**: Write tests for circular interval operations:
  - Handling of wraparound at 2π
  - Union and intersection with wraparound logic
  - Complement operations
  - Distance calculations
- **GREEN**: Implement S1Interval with wraparound logic
- **REFACTOR**: Optimize for common longitude/latitude operations

## Phase 4: Regional Types and Geometric Shapes (Level 3 Dependencies)

### 4.1 S2Region Abstract Interface

#### 4.1.1 S2Region Trait Definition
- **RED**: Write tests for S2Region trait:
  - Bounding cap and rectangle computation
  - Containment tests for points and cells
  - Intersection tests with other regions
  - Area and centroid calculations
- **GREEN**: Define S2Region trait with required methods
- **REFACTOR**: Ensure trait is ergonomic for implementors

### 4.2 S2Cap - Spherical Cap Implementation

#### 4.2.1 S2Cap Geometric Operations
- **RED**: Write comprehensive tests for S2Cap:
  - Construction from center point and height/angle
  - Containment tests for points and other caps
  - Intersection and union operations
  - Bounding rectangle computation
  - Area and circumference calculations
- **GREEN**: Implement S2Cap with center point and height
- **REFACTOR**: Optimize for common query patterns

### 4.3 S2Loop - Simple Polygon Implementation

#### 4.3.1 S2Loop Basic Structure
- **RED**: Write tests for S2Loop construction and validation:
  - Vertex ordering and orientation validation
  - Self-intersection detection
  - Hole detection and nesting validation
  - Area and perimeter calculations
- **GREEN**: Implement basic S2Loop structure
- **REFACTOR**: Optimize validation algorithms

#### 4.3.2 S2Loop Geometric Queries
- **RED**: Write tests for S2Loop queries:
  - Point containment using winding number
  - Edge intersection detection
  - Distance calculations to points and edges
  - Bounding rectangle computation
- **GREEN**: Implement geometric query algorithms
- **REFACTOR**: Optimize for performance-critical paths

### 4.4 S2Polygon - Complex Polygon Implementation

#### 4.4.1 S2Polygon Structure and Validation
- **RED**: Write tests for S2Polygon:
  - Multi-loop polygon construction
  - Hole nesting validation
  - Orientation consistency checking
  - Area calculation with holes
- **GREEN**: Implement S2Polygon as collection of S2Loops
- **REFACTOR**: Optimize multi-loop operations

### 4.5 S2Polyline - Path Representation

#### 4.5.1 S2Polyline Operations
- **RED**: Write tests for S2Polyline:
  - Construction from vertex sequence
  - Length and interpolation calculations
  - Closest point queries
  - Simplification algorithms
- **GREEN**: Implement S2Polyline with vertex storage
- **REFACTOR**: Optimize for path-based queries

## Phase 5: Indexing and Spatial Queries (Level 4 Dependencies)

### 5.1 S2Shape Abstract Interface

#### 5.1.1 S2Shape Trait Definition
- **RED**: Write tests for S2Shape trait:
  - Edge iteration and access
  - Reference point computation
  - Dimension and type identification
- **GREEN**: Define S2Shape trait for different geometry types
- **REFACTOR**: Ensure efficient implementation for all shape types

### 5.2 S2ShapeIndex - Spatial Index

#### 5.2.1 Basic Index Operations
- **RED**: Write tests for S2ShapeIndex:
  - Shape insertion and removal
  - Cell-based spatial partitioning
  - Iterator functionality over index contents
  - Memory management and updates
- **GREEN**: Implement basic index structure
- **REFACTOR**: Optimize for insert/query performance balance

#### 5.2.2 Advanced Index Queries
- **RED**: Write tests for complex queries:
  - Range queries over cell regions
  - Intersection detection between shapes
  - Nearest neighbor searches
  - Batch query operations
- **GREEN**: Implement query algorithms
- **REFACTOR**: Optimize query performance, add query statistics

### 5.3 Distance Query Operations

#### 5.3.1 S2ClosestEdgeQuery
- **RED**: Write tests for closest edge queries:
  - Distance calculations to points, edges, shapes
  - K-nearest neighbor searches
  - Distance bounds and early termination
  - Query result ranking and filtering
- **GREEN**: Implement distance query algorithms
- **REFACTOR**: Optimize for large-scale datasets

#### 5.3.2 S2ClosestPointQuery and S2ClosestCellQuery
- **RED**: Write tests for point and cell queries
- **GREEN**: Implement specialized query types
- **REFACTOR**: Share common query infrastructure

## Phase 6: Builder Framework and Construction (Level 5 Dependencies)

### 6.1 S2Builder - Robust Geometry Construction

#### 6.1.1 S2Builder Core Framework
- **RED**: Write tests for S2Builder:
  - Input geometry processing and validation
  - Snapping and simplification operations
  - Error handling and recovery
  - Layer-based output generation
- **GREEN**: Implement core S2Builder framework
- **REFACTOR**: Optimize for large input datasets

#### 6.1.2 S2Builder Layer Implementations
- **RED**: Write tests for various layer types:
  - S2PolygonLayer, S2PolylineLayer
  - Point and shape collection layers
  - Custom layer implementations
- **GREEN**: Implement standard layer types
- **REFACTOR**: Ensure consistent API across layer types

## Phase 7: Advanced Geometric Operations (Level 6 Dependencies)

### 7.1 Boolean Operations

#### 7.1.1 S2BooleanOperation Implementation
- **RED**: Write comprehensive tests for boolean operations:
  - Union, intersection, difference, symmetric difference
  - Multi-polygon operations
  - Precision and robustness testing
  - Edge case handling (touching edges, overlaps)
- **GREEN**: Implement boolean operation algorithms
- **REFACTOR**: Optimize for complex polygon operations

### 7.2 Buffer Operations

#### 7.2.1 S2BufferOperation Implementation
- **RED**: Write tests for buffer operations:
  - Positive and negative buffering
  - Join styles and end cap handling
  - Buffer of points, lines, and polygons
- **GREEN**: Implement buffering algorithms
- **REFACTOR**: Optimize buffer generation

### 7.3 Geometric Predicates and Validation

#### 7.3.1 S2Predicates Implementation
- **RED**: Write tests for robust geometric predicates:
  - Orientation tests with exact arithmetic
  - Edge intersection tests
  - Point-in-polygon tests
  - Degeneracy handling
- **GREEN**: Implement predicates using ExactFloat when needed
- **REFACTOR**: Optimize common cases while maintaining robustness

## Phase 8: High-Level APIs and Integration (Level 7 Dependencies)

### 8.1 Text Format and Debugging Support

#### 8.1.1 S2TextFormat Implementation
- **RED**: Write tests for text format parsing and generation:
  - WKT-style format support
  - Debug string generation
  - Parsing error handling and recovery
- **GREEN**: Implement text format support
- **REFACTOR**: Ensure compatibility with existing formats

### 8.2 Performance and Benchmarking

#### 8.2.1 Comprehensive Benchmarking Suite
- **RED**: Write performance tests comparing to C++ implementation:
  - Core geometric operations benchmarks
  - Index construction and query benchmarks
  - Memory usage profiling
  - Scalability testing with large datasets
- **GREEN**: Implement benchmarking framework
- **REFACTOR**: Optimize identified performance bottlenecks

### 8.3 Python Bindings (Optional)

#### 8.3.1 PyO3-based Python Interface
- **RED**: Write tests for Python API compatibility:
  - Core type exposure to Python
  - Memory management across language boundary
  - Error handling integration
- **GREEN**: Implement basic Python bindings
- **REFACTOR**: Optimize for Python performance patterns

## Testing Strategy Throughout Development

### Unit Testing Requirements
- Each component must achieve >95% code coverage
- Property-based testing for geometric operations using quickcheck
- Comparison testing against C++ reference implementation
- Numerical stability testing with edge cases

### Integration Testing
- End-to-end workflows using real geographic datasets
- Performance regression testing
- Memory safety verification under stress
- Thread safety validation for concurrent operations

### Documentation Standards
- Comprehensive rustdoc for all public APIs
- Mathematical background and algorithm explanations
- Migration guide from C++ to Rust API
- Performance characteristics documentation

## Cargo Workspace Structure

```
s2geometry-rust/
├── Cargo.toml (workspace)
├── s2-core/           # Core geometric types
├── s2-index/          # Indexing structures  
├── s2-query/          # Query operations
├── s2-builder/        # Geometry construction
├── s2-ops/            # Advanced operations
├── s2-text/           # Text format support
├── s2-python/         # Python bindings (optional)
├── benchmarks/        # Performance testing
└── integration-tests/ # End-to-end testing
```

## External Dependencies

### Required Crates
- `num-bigint` - For ExactFloat arithmetic
- `serde` - Serialization support
- `thiserror` - Error handling
- `approx` - Floating-point comparisons
- `geo-types` - Integration with Rust geo ecosystem

### Development Dependencies
- `criterion` - Benchmarking
- `proptest` - Property-based testing
- `quickcheck` - Additional property testing
- `pyo3` - Python bindings (optional)

## Success Criteria

The port is considered successful when:

1. **Correctness**: All geometric operations produce results equivalent to C++ implementation within numerical precision limits
2. **Performance**: Rust implementation performs within 10% of C++ performance for core operations
3. **Safety**: Zero unsafe code in public APIs, comprehensive error handling
4. **Usability**: APIs feel natural to Rust developers while preserving essential functionality
5. **Documentation**: Complete documentation with examples and migration guide
6. **Testing**: Comprehensive test suite with high coverage and property-based testing
7. **Integration**: Compatible with Rust geo ecosystem and provides Python bindings

## Implementation Timeline Estimate

- **Phase 1-2** (Foundation + Core Types): 8-10 weeks
- **Phase 3-4** (Cells + Regions): 10-12 weeks  
- **Phase 5** (Indexing): 8-10 weeks
- **Phase 6** (Builder): 6-8 weeks
- **Phase 7** (Advanced Ops): 10-12 weeks
- **Phase 8** (Integration): 4-6 weeks

**Total Estimated Time**: 46-58 weeks (approximately 1 year)

This timeline assumes dedicated development effort and may be reduced with parallel development of independent components.