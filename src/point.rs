//! S2Point - Points on the unit sphere
//!
//! S2Point represents a point on the unit sphere as a 3D unit vector.
//! This is the fundamental type for all S2 geometry operations.

use crate::math::{DVec3, float_utils::*};
use crate::error::{S2Error, S2Result};
use std::fmt;

/// A point on the unit sphere, represented as a 3D unit vector.
///
/// S2Point is the fundamental geometric primitive in S2 geometry.
/// It wraps a DVec3 and provides geometric operations.
///
/// # Mathematical Properties
/// - Can represent both normalized and unnormalized vectors (for intermediate calculations)
/// - Provides normalization when needed for sphere operations
/// - Supports exact arithmetic operations for robustness
///
/// # Example
/// ```rust,ignore
/// use s2geometry_rust::S2Point;
/// use s2geometry_rust::math::DVec3;
///
/// let coords = DVec3::new(1.0, 0.0, 0.0);
/// let point = S2Point::from_vec3(coords)?;
/// assert_eq!(point.x(), 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct S2Point {
    coords: DVec3,
}

impl S2Point {
    /// Create an S2Point from normalized coordinates
    /// 
    /// This is used when you know the coordinates are already normalized.
    /// Does not perform normalization checks in release mode.
    #[inline]
    pub fn from_normalized(coords: DVec3) -> Self {
        Self { coords }
    }

    /// Create an S2Point from coordinates with automatic normalization
    /// 
    /// Returns an error if the input vector has zero length.
    pub fn from_vec3(coords: DVec3) -> S2Result<Self> {
        let length_sq = coords.length_squared();
        if length_sq < f64::MIN_POSITIVE {
            return Err(S2Error::invalid_point("Zero-length vector"));
        }
        
        let normalized = coords / length_sq.sqrt();
        Ok(Self { coords: normalized })
    }

    /// Create an S2Point from individual coordinates
    pub fn new(x: f64, y: f64, z: f64) -> S2Result<Self> {
        Self::from_vec3(DVec3::new(x, y, z))
    }

    /// Create an S2Point directly from coordinates without normalization
    /// Used for testing and intermediate calculations
    #[inline]
    pub fn from_coords_raw(x: f64, y: f64, z: f64) -> Self {
        Self { coords: DVec3::new(x, y, z) }
    }

    /// Get the underlying coordinate vector
    #[inline]
    pub fn coords(&self) -> DVec3 {
        self.coords
    }

    /// Get X coordinate  
    #[inline]
    pub fn x(&self) -> f64 {
        self.coords.x
    }

    /// Get Y coordinate
    #[inline] 
    pub fn y(&self) -> f64 {
        self.coords.y
    }

    /// Get Z coordinate
    #[inline]
    pub fn z(&self) -> f64 {
        self.coords.z
    }

    /// Normalize this point to unit length
    /// Returns a new normalized S2Point
    pub fn normalize(&self) -> Self {
        let length = self.coords.length();
        if approx_eq(length, 1.0) {
            *self
        } else if length > f64::MIN_POSITIVE {
            Self { coords: self.coords / length }
        } else {
            // Return original point for zero-length vectors
            *self
        }
    }

    /// Compute dot product with another point
    #[inline]
    pub fn dot(&self, other: &S2Point) -> f64 {
        self.coords.dot(other.coords)
    }

    /// Compute cross product with another point  
    #[inline]
    pub fn cross(&self, other: &S2Point) -> DVec3 {
        self.coords.cross(other.coords)
    }

    /// Compute angle between two points (assumes both are normalized)
    /// Uses atan2 for better numerical precision, especially for small angles
    pub fn angle(&self, other: &S2Point) -> f64 {
        let cross_prod = self.coords().cross(other.coords());
        let dot_prod = self.coords().dot(other.coords());
        cross_prod.length().atan2(dot_prod)
    }

    /// Compute Euclidean distance between two points
    #[inline]
    pub fn distance(&self, other: &S2Point) -> f64 {
        (self.coords - other.coords).length()
    }

    /// Compute squared Euclidean distance between two points
    #[inline]
    pub fn distance_squared(&self, other: &S2Point) -> f64 {
        (self.coords - other.coords).length_squared()
    }

    /// Check if this point is unit length (normalized)
    #[inline]
    pub fn is_unit_length(&self) -> bool {
        const TOLERANCE: f64 = 1e-15;
        (self.coords.length_squared() - 1.0).abs() <= TOLERANCE
    }

    /// Interpolate between two points with parameter t ∈ [0,1]
    /// Uses spherical linear interpolation (slerp) for normalized points
    pub fn interpolate(&self, other: &S2Point, t: f64) -> S2Point {
        if approx_eq(t, 0.0) {
            return *self;
        }
        if approx_eq(t, 1.0) {
            return *other;
        }

        let dot_prod = self.dot(other).clamp(-1.0, 1.0);
        let angle = dot_prod.acos();
        
        if approx_zero(angle) {
            // Points are identical or nearly so
            return *self;
        }
        
        let sin_angle = angle.sin();
        if approx_zero(sin_angle) {
            // Points are antipodal, use linear interpolation 
            let result = self.coords * (1.0 - t) + other.coords * t;
            S2Point::from_vec3(result).unwrap_or(*self)
        } else {
            // Standard slerp
            let factor1 = ((1.0 - t) * angle).sin() / sin_angle;
            let factor2 = (t * angle).sin() / sin_angle;
            let result = self.coords * factor1 + other.coords * factor2;
            S2Point::from_normalized(result)
        }
    }

    /// Component-wise division (matching C++ Vector3::DivComponents)
    #[inline]
    pub fn div_components(&self, other: &S2Point) -> S2Point {
        let result = DVec3::new(
            self.x() / other.x(),
            self.y() / other.y(), 
            self.z() / other.z()
        );
        S2Point { coords: result }
    }

    /// Component-wise square root (matching C++ Vector3::Sqrt) 
    #[inline]
    pub fn sqrt(&self) -> S2Point {
        let result = DVec3::new(
            self.x().sqrt(),
            self.y().sqrt(),
            self.z().sqrt()
        );
        // Don't normalize - this is component-wise operation
        S2Point { coords: result }
    }

    /// Component-wise floor (matching C++ Vector3::Floor)
    #[inline]
    pub fn floor(&self) -> S2Point {
        let result = DVec3::new(
            self.x().floor(),
            self.y().floor(),
            self.z().floor()
        );
        S2Point { coords: result }
    }

    /// Component-wise ceil (matching C++ Vector3::Ceil)  
    #[inline]
    pub fn ceil(&self) -> S2Point {
        let result = DVec3::new(
            self.x().ceil(),
            self.y().ceil(),
            self.z().ceil()
        );
        S2Point { coords: result }
    }

    /// Component-wise round (matching C++ Vector3::FRound)
    #[inline] 
    pub fn fround(&self) -> S2Point {
        let result = DVec3::new(
            self.x().round(),
            self.y().round(),
            self.z().round()
        );
        S2Point { coords: result }
    }
}

// Arithmetic operations (matching C++ Vector3 behavior exactly)
impl std::ops::Add for S2Point {
    type Output = S2Point;
    
    #[inline]
    fn add(self, other: S2Point) -> S2Point {
        S2Point { coords: self.coords + other.coords }
    }
}

impl std::ops::Sub for S2Point {
    type Output = S2Point;
    
    #[inline]
    fn sub(self, other: S2Point) -> S2Point {
        S2Point { coords: self.coords - other.coords }
    }
}

impl std::ops::SubAssign for S2Point {
    #[inline]
    fn sub_assign(&mut self, other: S2Point) {
        self.coords -= other.coords;
    }
}

impl std::ops::Mul<f64> for S2Point {
    type Output = S2Point;
    
    #[inline]
    fn mul(self, scalar: f64) -> S2Point {
        S2Point { coords: self.coords * scalar }
    }
}

impl std::ops::Mul<S2Point> for f64 {
    type Output = S2Point;
    
    #[inline]
    fn mul(self, point: S2Point) -> S2Point {
        S2Point { coords: self * point.coords }
    }
}

impl std::ops::Div<f64> for S2Point {
    type Output = S2Point;
    
    #[inline]
    fn div(self, scalar: f64) -> S2Point {
        S2Point { coords: self.coords / scalar }
    }
}

impl std::ops::Neg for S2Point {
    type Output = S2Point;
    
    #[inline]
    fn neg(self) -> S2Point {
        S2Point { coords: -self.coords }
    }
}

// Hash implementation
impl std::hash::Hash for S2Point {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Use bit representation for exact hash consistency
        self.coords.x.to_bits().hash(state);
        self.coords.y.to_bits().hash(state);
        self.coords.z.to_bits().hash(state);
    }
}

// Display implementation
impl fmt::Display for S2Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.10}, {:.10}, {:.10})", self.x(), self.y(), self.z())
    }
}

// Conversion from DVec3
impl From<DVec3> for S2Point {
    fn from(coords: DVec3) -> Self {
        Self { coords }
    }
}

// Conversion to DVec3  
impl From<S2Point> for DVec3 {
    fn from(point: S2Point) -> Self {
        point.coords
    }
}

/// Hash function for S2Point coordinates
pub struct S2PointHash;

impl S2PointHash {
    /// Creates a new S2PointHash instance
    pub fn new() -> Self {
        Self
    }
    
    /// Computes a hash value for the given S2Point
    pub fn hash_point(&self, point: &S2Point) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let mut hasher = DefaultHasher::new();
        point.hash(&mut hasher);
        hasher.finish()
    }
}

impl Default for S2PointHash {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::constants::PI_2;

    #[test]
    fn test_point_creation() {
        let p = S2Point::new(1.0, 0.0, 0.0).unwrap();
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 0.0);
        assert_eq!(p.z(), 0.0);
    }

    #[test]
    fn test_point_normalization() {
        let p = S2Point::new(2.0, 0.0, 0.0).unwrap();
        assert!(approx_eq(p.x(), 1.0));
        assert!(approx_eq(p.coords().length(), 1.0));
    }

    #[test]
    fn test_point_operations() {
        let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let p2 = S2Point::new(0.0, 1.0, 0.0).unwrap();
        
        assert!(approx_zero(p1.dot(&p2)));
        assert!(approx_eq(p1.angle(&p2), PI_2));
        
        let cross = p1.cross(&p2);
        assert!(approx_eq(cross.z, 1.0));
    }

    #[test]
    fn test_raw_operations() {
        // Test that component-wise operations work without normalization
        let p1 = S2Point::from_coords_raw(4.0, 8.0, 16.0);
        let p2 = S2Point::from_coords_raw(2.0, 2.0, 2.0);
        let result = p1.div_components(&p2);
        
        assert_eq!(result.x(), 2.0);
        assert_eq!(result.y(), 4.0); 
        assert_eq!(result.z(), 8.0);
    }
}