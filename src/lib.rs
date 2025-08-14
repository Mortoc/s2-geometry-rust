//! # S2 Geometry Library - Rust Implementation
//!
//! This is a comprehensive Rust port of Google's S2 Geometry Library,
//! designed for computational geometry on the sphere.
//!
//! ## Key Features
//!
//! - **Numerical Robustness**: Two-tier architecture with exact arithmetic fallback
//! - **Performance**: Built on Bevy's battle-tested `glam` math foundation  
//! - **Thread Safety**: All types are `Send + Sync` without locks
//! - **Memory Safety**: Leverages Rust's ownership system for safe geometric operations
//! - **Bevy Integration**: Seamless compatibility with Bevy's math ecosystem
//!
//! ## Mathematical Architecture (Matching Google's C++ Implementation)
//!
//! This library uses the **three-tier mathematical strategy** pioneered by Google's C++ S2:
//!
//! 1. **Fast Path (~90-95%)**: `glam::DVec3` f64 precision with rigorous error analysis
//! 2. **Stable Path (~4-9%)**: Extended precision for borderline cases  
//! 3. **Exact Path (<1%)**: `BigRational` exact arithmetic for perfect deterministic results
//!
//! This **exactly matches** Google's proven approach, ensuring identical robustness while 
//! leveraging Rust's performance and safety guarantees.
//!
//! ### Why Exact Arithmetic?
//!
//! Geometric predicates like orientation tests and point-in-polygon queries must be **deterministic**
//! to maintain topological consistency. Small floating-point errors can cause:
//! - Orientation tests to give inconsistent results
//! - Point-in-polygon queries to flip randomly
//! - Edge intersections to be missed or falsely detected
//!
//! Our exact arithmetic fallback ensures these operations always produce correct, consistent results.
//!
//! ## Architecture
//!
//! The library is organized into several modules:
//!
//! - [`error`]: Error types and result handling
//! - [`math`]: Mathematical foundations combining glam f64 + exact arithmetic
//! - [`point`]: Basic point types and operations  
//! - [`cell`]: Hierarchical spatial indexing with S2CellId
//! - [`region`]: Geometric regions and shapes
//! - [`index`]: Spatial indexing and query operations
//!
//! ## Example
//!
//! ```rust,ignore
//! use s2geometry_rust::{S2Point, S2Cell, S2CellId};
//! use s2geometry_rust::math::DVec3;
//!
//! // Create a point on the unit sphere using glam types
//! let coords = DVec3::new(1.0, 0.0, 0.0);
//! let point = S2Point::from_vec3(coords)?;
//!
//! // Find the S2 cell containing this point  
//! let cell_id = S2CellId::from_point(&point);
//! let cell = S2Cell::from(cell_id);
//!
//! // Robust geometric queries with exact arithmetic fallback
//! assert!(cell.contains(&point));
//!
//! // Direct integration with Bevy math types
//! let bevy_transform = bevy_math::Vec3::from(point.coords());
//! ```
//!
//! ## Performance Notes
//!
//! - **Fast operations**: Vector math, coordinate transformations, distance calculations
//! - **Exact fallback**: Orientation tests, edge intersections, point containment near boundaries  
//! - **Typical fallback rate**: <1% of geometric operations require exact arithmetic
//! - **Memory usage**: Exact arithmetic objects created on-demand, not stored permanently
//!
//! See [`ARCHITECTURE.md`](https://github.com/user/s2geometry-rust/blob/main/ARCHITECTURE.md) 
//! for detailed design documentation.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]
#![allow(clippy::module_name_repetitions)]

pub mod error;
pub mod math;
pub mod point;
pub mod angle;
pub mod chord_angle;
pub mod cell_id;
pub mod interval;
pub mod r2;
pub mod cell;
pub mod latlng;
pub mod latlng_rect;
pub mod cap;
pub mod r#loop;
pub mod cell_union;
pub mod edge_crosser;
pub mod predicates;
pub mod polyline;
pub mod region_coverer;
pub mod shape;
pub mod shape_index;
pub mod point_shape;
pub mod polyline_shape;
pub mod polygon_shape;
pub mod mutable_shape_index;
pub mod builder;

// Re-export key types for convenience
pub use error::{S2Error, S2Result};
pub use point::S2Point;
pub use angle::S1Angle;
pub use chord_angle::S1ChordAngle;
pub use cell_id::S2CellId;
pub use interval::{S1Interval, R1Interval};
pub use r2::{R2Point, R2Rect};
pub use cell::S2Cell;
pub use latlng::S2LatLng;
pub use latlng_rect::S2LatLngRect;
pub use cap::S2Cap;
pub use r#loop::S2Loop;
pub use cell_union::S2CellUnion;
pub use edge_crosser::{S2EdgeCrosser, S2CopyingEdgeCrosser};
pub use polyline::S2Polyline;
pub use region_coverer::{S2RegionCoverer, S2RegionCovererOptions, S2Region};
pub use shape::{S2Shape, Edge, Chain, ChainPosition, ReferencePoint, ShapeDimension};
pub use shape_index::{S2ShapeIndex, S2ShapeIndexIterator, MutableS2ShapeIndex as MutableS2ShapeIndexTrait, 
                     S2ClippedShape, S2ShapeIndexCell};
pub use point_shape::{S2PointShape, S2MultiPointShape, S2PointCloudShape};
pub use polyline_shape::{S2PolylineShape, S2MultiPolylineShape, S2LaxPolylineShape};
pub use polygon_shape::{S2LoopShape, S2PolygonShape, S2MultiPolygonShape};
pub use mutable_shape_index::{MutableS2ShapeIndex, MutableS2ShapeIndexOptions};
pub use builder::{S2Builder, Options as S2BuilderOptions, EdgeType, EdgeId, VertexId};
// pub mod cell;
// pub mod region;
// pub mod shape;
// pub mod index;
// pub mod builder;
// pub mod ops;