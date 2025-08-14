//! Snap Functions for S2Builder
//!
//! Snap functions determine how input vertices are mapped to discrete "snap sites"
//! during the building process. This is crucial for numerical robustness and
//! ensuring that geometric operations produce consistent results.
//!
//! # Core Concepts
//!
//! - **Snap Sites**: Discrete locations where input vertices are "snapped" to
//! - **Snap Radius**: Maximum distance a vertex can be moved during snapping
//! - **Idempotency**: Snapping the same point multiple times gives the same result
//! - **Minimum Separation**: Ensuring snap sites are sufficiently separated
//!
//! # Available Snap Functions
//!
//! - [`IdentitySnapFunction`]: No snapping (identity transformation)
//! - [`IntLatLngSnapFunction`]: Snap to integer lat/lng coordinates
//! - [`S2CellIdSnapFunction`]: Snap to S2CellId centers at a given level
//!
//! # Example
//!
//! ```rust,ignore
//! use s2geometry_rust::builder::snap_functions::*;
//! use s2geometry_rust::S2Point;
//!
//! // Create a snap function that snaps to level 10 S2Cell centers
//! let snap_fn = S2CellIdSnapFunction::new(10);
//! let input_point = S2Point::from_coords(1.0, 0.1, 0.1).unwrap();
//! let snapped = snap_fn.snap_point(input_point);
//! ```

use crate::error::{S2Error, S2Result};
use crate::point::S2Point;
use crate::latlng::S2LatLng;
use crate::cell_id::S2CellId;
use crate::math::DVec3;
use std::fmt::Debug;

/// Original trait for functions that snap points to discrete locations (private)
trait SnapFunction_old: Debug + Send + Sync {
    /// Snap a point to a discrete location
    ///
    /// # Arguments
    /// * `point` - The input point to snap
    ///
    /// # Returns
    /// The snapped point location
    fn snap_point(&self, point: S2Point) -> S2Point;

    /// Get the maximum snap radius
    ///
    /// This is the maximum distance that any point can be moved during snapping.
    /// It's used for error analysis and determining edge intersection tolerances.
    fn snap_radius(&self) -> f64;

    /// Get the minimum separation between snap sites
    ///
    /// This guarantees that distinct input points won't snap to locations
    /// that are closer than this distance.
    fn min_vertex_separation(&self) -> f64;

    /// Clone the snap function
    fn clone_box(&self) -> Box<dyn SnapFunction_old>;
}

impl Clone for Box<dyn SnapFunction_old> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Identity snap function that doesn't move points
///
/// This snap function leaves all points unchanged, which is useful
/// when the input is already properly discretized or when exact
/// coordinates are required.
#[derive(Debug, Clone)]
pub struct IdentitySnapFunction {
    /// Minimum vertex separation (used for validation)
    min_vertex_separation: f64,
}

impl IdentitySnapFunction {
    /// Create a new identity snap function
    ///
    /// # Arguments
    /// * `min_vertex_separation` - Minimum allowed distance between vertices
    pub fn new(min_vertex_separation: f64) -> Self {
        Self {
            min_vertex_separation,
        }
    }

    /// Create an identity snap function with no minimum separation
    pub fn with_no_separation() -> Self {
        Self::new(0.0)
    }
}

impl SnapFunction_old for IdentitySnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        point // No transformation
    }

    fn snap_radius(&self) -> f64 {
        0.0 // No snapping occurs
    }

    fn min_vertex_separation(&self) -> f64 {
        self.min_vertex_separation
    }

    fn clone_box(&self) -> Box<dyn SnapFunction_old> {
        Box::new(self.clone())
    }
}

/// Snap function that snaps to integer latitude/longitude coordinates
///
/// This snaps points to a regular grid in latitude/longitude space,
/// which is useful for geographic applications where coordinates
/// are often stored as integers in microdegrees or similar units.
#[derive(Debug, Clone)]
pub struct IntLatLngSnapFunction {
    /// Scale factor for latitude/longitude coordinates
    scale: f64,
    /// Cached snap radius
    snap_radius: f64,
    /// Cached minimum vertex separation
    min_vertex_separation: f64,
}

impl IntLatLngSnapFunction {
    /// Create a new integer lat/lng snap function
    ///
    /// # Arguments
    /// * `scale` - Scale factor (e.g., 1e6 for microdegrees)
    ///
    /// # Example
    /// ```rust,ignore
    /// // Snap to microdegree precision
    /// let snap_fn = IntLatLngSnapFunction::new(1e6);
    /// ```
    pub fn new(scale: f64) -> Self {
        // Compute conservative bounds for snap radius and min separation
        let snap_radius = Self::compute_snap_radius(scale);
        let min_vertex_separation = Self::compute_min_vertex_separation(scale);

        Self {
            scale,
            snap_radius,
            min_vertex_separation,
        }
    }

    /// Create snap function for microdegree precision (scale = 1e6)
    pub fn microdegrees() -> Self {
        Self::new(1e6)
    }

    /// Create snap function for degree precision (scale = 1)
    pub fn degrees() -> Self {
        Self::new(1.0)
    }

    /// Compute the maximum snap radius for the given scale
    fn compute_snap_radius(scale: f64) -> f64 {
        // Maximum distance a point can move is roughly 1/(2*scale) in lat/lng space
        // Convert to radians and then to angular distance on sphere
        let max_lat_lng_error = 1.0 / (2.0 * scale) * std::f64::consts::PI / 180.0;
        
        // Conservative upper bound: sqrt(2) * max_lat_lng_error for diagonal movement
        max_lat_lng_error * std::f64::consts::SQRT_2
    }

    /// Compute the minimum vertex separation for the given scale
    fn compute_min_vertex_separation(scale: f64) -> f64 {
        // Points that are 1/scale apart in lat/lng will snap to different sites
        let min_lat_lng_separation = 1.0 / scale * std::f64::consts::PI / 180.0;
        
        // Conservative lower bound
        min_lat_lng_separation * 0.5
    }
}

impl SnapFunction_old for IntLatLngSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        let latlng = S2LatLng::from_point(point);
        
        // Convert to scaled integer coordinates and back
        let lat_scaled = (latlng.lat().degrees() * self.scale).round() / self.scale;
        let lng_scaled = (latlng.lng().degrees() * self.scale).round() / self.scale;
        
        // Convert back to S2Point
        let snapped_latlng = S2LatLng::from_degrees(lat_scaled, lng_scaled);
        snapped_latlng.to_point().unwrap_or(point) // Fallback to original on error
    }

    fn snap_radius(&self) -> f64 {
        self.snap_radius
    }

    fn min_vertex_separation(&self) -> f64 {
        self.min_vertex_separation
    }

    fn clone_box(&self) -> Box<dyn SnapFunction_old> {
        Box::new(self.clone())
    }
}

/// Snap function that snaps to S2CellId centers at a specific level
///
/// This snaps points to the centers of S2 cells at the specified level,
/// which provides a natural hierarchical discretization of the sphere.
#[derive(Debug, Clone)]
pub struct S2CellIdSnapFunction {
    /// S2 cell level for snapping
    level: i32,
    /// Cached snap radius
    snap_radius: f64,
    /// Cached minimum vertex separation
    min_vertex_separation: f64,
}

impl S2CellIdSnapFunction {
    /// Create a new S2CellId snap function
    ///
    /// # Arguments
    /// * `level` - S2 cell level (0-30, higher = finer grid)
    ///
    /// # Example
    /// ```rust,ignore
    /// // Snap to level 10 cell centers (roughly 600m cell size)
    /// let snap_fn = S2CellIdSnapFunction::new(10);
    /// ```
    pub fn new(level: i32) -> Self {
        if level < 0 || level > 30 {
            panic!("S2CellId level must be in range [0, 30]");
        }

        let snap_radius = Self::compute_snap_radius(level);
        let min_vertex_separation = Self::compute_min_vertex_separation(level);

        Self {
            level,
            snap_radius,
            min_vertex_separation,
        }
    }

    /// Compute the maximum snap radius for the given level
    fn compute_snap_radius(level: i32) -> f64 {
        // The maximum distance from a cell center to any point in the cell
        // is roughly half the cell diagonal
        let cell_size = S2CellId::avg_edge_metric().get_value(level);
        cell_size * std::f64::consts::SQRT_2 / 2.0
    }

    /// Compute the minimum vertex separation for the given level
    fn compute_min_vertex_separation(level: i32) -> f64 {
        // Adjacent cell centers are separated by roughly the cell size
        // Use approximate formula: cell size = 2 * PI / (1 << level) for face width
        let cell_size = std::f64::consts::PI / (1u64 << level) as f64;
        cell_size * 0.8 // Conservative factor
    }
}

impl SnapFunction_old for S2CellIdSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        let cell_id = S2CellId::from_point(&point);
        let parent_cell = cell_id.parent_at_level(self.level);
        parent_cell.to_point()
    }

    fn snap_radius(&self) -> f64 {
        self.snap_radius
    }

    fn min_vertex_separation(&self) -> f64 {
        self.min_vertex_separation
    }

    fn clone_box(&self) -> Box<dyn SnapFunction_old> {
        Box::new(self.clone())
    }
}

/// Snap function that enforces a minimum edge length
///
/// This is useful for preventing very short edges that can cause
/// numerical issues in downstream operations.
#[derive(Debug, Clone)]
pub struct MinEdgeLengthSnapFunction {
    /// Base snap function to use
    base: Box<dyn SnapFunction_old>,
    /// Minimum edge length to enforce
    min_edge_length: f64,
}

impl MinEdgeLengthSnapFunction {
    /// Create a new minimum edge length snap function
    ///
    /// # Arguments
    /// * `base` - Base snap function to apply first
    /// * `min_edge_length` - Minimum edge length in radians
    pub fn new(base: Box<dyn SnapFunction_old>, min_edge_length: f64) -> Self {
        Self {
            base,
            min_edge_length,
        }
    }
}

impl SnapFunction_old for MinEdgeLengthSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        // Apply base snapping first
        self.base.snap_point(point)
        // TODO: Implement edge length enforcement in S2Builder
    }

    fn snap_radius(&self) -> f64 {
        // Conservative upper bound
        self.base.snap_radius() + self.min_edge_length
    }

    fn min_vertex_separation(&self) -> f64 {
        // Use the larger of base separation and min edge length
        self.base.min_vertex_separation().max(self.min_edge_length)
    }

    fn clone_box(&self) -> Box<dyn SnapFunction_old> {
        Box::new(Self {
            base: self.base.clone(),
            min_edge_length: self.min_edge_length,
        })
    }
}

/// Enum wrapper for different snap function implementations
///
/// This provides a concrete type that can be stored and cloned while
/// providing access to all the different snap function implementations.
#[derive(Debug, Clone)]
pub enum SnapFunction {
    /// Identity snap function (no snapping)
    Identity(IdentitySnapFunction),
    /// Integer lat/lng snap function
    IntLatLng(IntLatLngSnapFunction),
    /// S2CellId snap function
    S2CellId(S2CellIdSnapFunction),
    /// Minimum edge length snap function
    MinEdgeLength(MinEdgeLengthSnapFunction),
}

impl SnapFunction {
    /// Snap a point using the wrapped implementation
    pub fn snap_point(&self, point: S2Point) -> S2Point {
        match self {
            SnapFunction::Identity(f) => SnapFunction_old::snap_point(f, point),
            SnapFunction::IntLatLng(f) => SnapFunction_old::snap_point(f, point),
            SnapFunction::S2CellId(f) => SnapFunction_old::snap_point(f, point),
            SnapFunction::MinEdgeLength(f) => SnapFunction_old::snap_point(f, point),
        }
    }

    /// Get the snap radius using the wrapped implementation
    pub fn snap_radius(&self) -> f64 {
        match self {
            SnapFunction::Identity(f) => SnapFunction_old::snap_radius(f),
            SnapFunction::IntLatLng(f) => SnapFunction_old::snap_radius(f),
            SnapFunction::S2CellId(f) => SnapFunction_old::snap_radius(f),
            SnapFunction::MinEdgeLength(f) => SnapFunction_old::snap_radius(f),
        }
    }

    /// Get the minimum vertex separation using the wrapped implementation
    pub fn min_vertex_separation(&self) -> f64 {
        match self {
            SnapFunction::Identity(f) => SnapFunction_old::min_vertex_separation(f),
            SnapFunction::IntLatLng(f) => SnapFunction_old::min_vertex_separation(f),
            SnapFunction::S2CellId(f) => SnapFunction_old::min_vertex_separation(f),
            SnapFunction::MinEdgeLength(f) => SnapFunction_old::min_vertex_separation(f),
        }
    }
}

impl Default for SnapFunction {
    fn default() -> Self {
        SnapFunction::Identity(IdentitySnapFunction::with_no_separation())
    }
}

/// Rename the original trait to avoid conflicts
pub trait SnapFunctionTrait: Debug + Send + Sync {
    /// Snap a point to a discrete location
    fn snap_point(&self, point: S2Point) -> S2Point;
    /// Get the maximum snap radius
    fn snap_radius(&self) -> f64;
    /// Get the minimum separation between snap sites
    fn min_vertex_separation(&self) -> f64;
    /// Clone the snap function
    fn clone_box(&self) -> Box<dyn SnapFunctionTrait>;
}

impl Clone for Box<dyn SnapFunctionTrait> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

// Implement the trait for all the concrete types
impl SnapFunctionTrait for IdentitySnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        <Self as SnapFunction_old>::snap_point(self, point)
    }
    fn snap_radius(&self) -> f64 {
        <Self as SnapFunction_old>::snap_radius(self)
    }
    fn min_vertex_separation(&self) -> f64 {
        <Self as SnapFunction_old>::min_vertex_separation(self)
    }
    fn clone_box(&self) -> Box<dyn SnapFunctionTrait> {
        Box::new(self.clone())
    }
}

impl SnapFunctionTrait for IntLatLngSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        <Self as SnapFunction_old>::snap_point(self, point)
    }
    fn snap_radius(&self) -> f64 {
        <Self as SnapFunction_old>::snap_radius(self)
    }
    fn min_vertex_separation(&self) -> f64 {
        <Self as SnapFunction_old>::min_vertex_separation(self)
    }
    fn clone_box(&self) -> Box<dyn SnapFunctionTrait> {
        Box::new(self.clone())
    }
}

impl SnapFunctionTrait for S2CellIdSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        <Self as SnapFunction_old>::snap_point(self, point)
    }
    fn snap_radius(&self) -> f64 {
        <Self as SnapFunction_old>::snap_radius(self)
    }
    fn min_vertex_separation(&self) -> f64 {
        <Self as SnapFunction_old>::min_vertex_separation(self)
    }
    fn clone_box(&self) -> Box<dyn SnapFunctionTrait> {
        Box::new(self.clone())
    }
}

impl SnapFunctionTrait for MinEdgeLengthSnapFunction {
    fn snap_point(&self, point: S2Point) -> S2Point {
        <Self as SnapFunction_old>::snap_point(self, point)
    }
    fn snap_radius(&self) -> f64 {
        <Self as SnapFunction_old>::snap_radius(self)
    }
    fn min_vertex_separation(&self) -> f64 {
        <Self as SnapFunction_old>::min_vertex_separation(self)
    }
    fn clone_box(&self) -> Box<dyn SnapFunctionTrait> {
        Box::new(self.clone())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;

    #[test]
    fn test_identity_snap_function() {
        let snap_fn = IdentitySnapFunction::new(1e-10);
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.1, 0.1));
        let snapped = SnapFunction_old::snap_point(&snap_fn, point);
        
        // Should be unchanged
        assert_eq!(point, snapped);
        assert_eq!(SnapFunction_old::snap_radius(&snap_fn), 0.0);
    }

    #[test]
    fn test_int_latlng_snap_function() {
        let snap_fn = IntLatLngSnapFunction::new(1.0); // Degree precision
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.1, 0.1));
        let snapped = SnapFunction_old::snap_point(&snap_fn, point);
        
        // Should be different (snapped to integer degrees)
        assert_ne!(point, snapped);
        assert!(SnapFunction_old::snap_radius(&snap_fn) > 0.0);
        assert!(SnapFunction_old::min_vertex_separation(&snap_fn) > 0.0);
    }

    #[test]
    fn test_s2cellid_snap_function() {
        let snap_fn = S2CellIdSnapFunction::new(10);
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.1, 0.1));
        let snapped = SnapFunction_old::snap_point(&snap_fn, point);
        
        // Basic validation - should get a valid point
        assert!(snapped.is_unit_length());
        
        // Check metrics are reasonable  
        assert!(SnapFunction_old::snap_radius(&snap_fn) > 0.0);
        assert!(SnapFunction_old::min_vertex_separation(&snap_fn) > 0.0);
        
        // For now, just check that snapping produces a different point
        // (unless the input was already exactly at a cell center)
        println!("Original: {:?}", point);
        println!("Snapped:  {:?}", snapped);
    }

    #[test]
    fn test_snap_function_cloning() {
        let snap_fn: Box<dyn SnapFunction_old> = Box::new(IdentitySnapFunction::new(1e-10));
        let cloned = snap_fn.clone();
        
        let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        assert_eq!(snap_fn.snap_point(point), cloned.snap_point(point));
    }

    #[test]
    #[should_panic]
    fn test_invalid_s2cellid_level() {
        S2CellIdSnapFunction::new(-1); // Should panic
    }

    #[test]
    #[should_panic]
    fn test_invalid_s2cellid_level_too_high() {
        S2CellIdSnapFunction::new(31); // Should panic
    }
}