//! S2LatLng - Latitude/longitude coordinates on the unit sphere
//!
//! S2LatLng represents a point on the unit sphere as a pair of latitude-longitude
//! coordinates. This provides the bridge between traditional geographic coordinates
//! and S2's spherical geometry system.

use crate::math::{constants::*, float_utils::*};
use crate::{S1Angle, S2Point, R2Point, S2Error, S2Result};
use std::fmt;
use std::ops::{Add, Sub, Mul};
use std::hash::{Hash, Hasher};

/// A point on the unit sphere represented as latitude-longitude coordinates.
///
/// S2LatLng stores coordinates internally as radians in an R2Point, following
/// the C++ implementation exactly. The latitude is restricted to [-π/2, π/2]
/// and longitude to [-π, π] for valid coordinates.
///
/// # Example
/// ```rust,ignore
/// use s2geometry_rust::S2LatLng;
///
/// let seattle = S2LatLng::from_degrees(47.6062, -122.3321);
/// assert!(seattle.is_valid());
/// assert_eq!(seattle.lat().degrees(), 47.6062);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct S2LatLng {
    coords: R2Point, // [lat_radians, lng_radians]
}

impl S2LatLng {
    /// Create S2LatLng from S1Angle latitude and longitude
    #[inline]
    pub fn new(lat: S1Angle, lng: S1Angle) -> Self {
        Self {
            coords: R2Point::new(lat.radians(), lng.radians()),
        }
    }

    /// Default constructor - creates (0, 0) lat/lng
    #[inline]
    pub fn default() -> Self {
        Self {
            coords: R2Point::new(0.0, 0.0),
        }
    }

    /// Create S2LatLng from radians
    #[inline]
    pub fn from_radians(lat_radians: f64, lng_radians: f64) -> Self {
        Self {
            coords: R2Point::new(lat_radians, lng_radians),
        }
    }

    /// Create S2LatLng from degrees
    #[inline]
    pub fn from_degrees(lat_degrees: f64, lng_degrees: f64) -> Self {
        Self::new(
            S1Angle::from_degrees(lat_degrees),
            S1Angle::from_degrees(lng_degrees),
        )
    }

    /// Create S2LatLng from E5 representation (degrees * 1e5)
    #[inline]
    pub fn from_e5(lat_e5: i32, lng_e5: i32) -> Self {
        Self::new(S1Angle::from_e5(lat_e5), S1Angle::from_e5(lng_e5))
    }

    /// Create S2LatLng from E6 representation (degrees * 1e6)
    #[inline]
    pub fn from_e6(lat_e6: i32, lng_e6: i32) -> Self {
        Self::new(S1Angle::from_e6(lat_e6), S1Angle::from_e6(lng_e6))
    }

    /// Create S2LatLng from E7 representation (degrees * 1e7)
    #[inline]
    pub fn from_e7(lat_e7: i32, lng_e7: i32) -> Self {
        Self::new(S1Angle::from_e7(lat_e7), S1Angle::from_e7(lng_e7))
    }

    /// Create S2LatLng from unsigned E6 representation
    #[inline]
    pub fn from_unsigned_e6(lat_e6: u32, lng_e6: u32) -> Self {
        Self::new(
            S1Angle::from_unsigned_e6(lat_e6),
            S1Angle::from_unsigned_e6(lng_e6),
        )
    }

    /// Create S2LatLng from unsigned E7 representation
    #[inline]
    pub fn from_unsigned_e7(lat_e7: u32, lng_e7: u32) -> Self {
        Self::new(
            S1Angle::from_unsigned_e7(lat_e7),
            S1Angle::from_unsigned_e7(lng_e7),
        )
    }

    /// Create S2LatLng from S2Point
    pub fn from_point(p: S2Point) -> Self {
        Self {
            coords: R2Point::new(Self::latitude_radians(p), Self::longitude_radians(p)),
        }
    }

    /// Returns an invalid S2LatLng (coordinates outside valid range)
    #[inline]
    pub fn invalid() -> Self {
        // Use coordinates outside the valid range (matching C++ implementation)
        Self::from_radians(PI, 2.0 * PI)
    }

    /// Compute latitude from S2Point (static method matching C++ Latitude)
    pub fn latitude(p: S2Point) -> S1Angle {
        S1Angle::from_radians(Self::latitude_radians(p))
    }

    /// Compute longitude from S2Point (static method matching C++ Longitude)
    pub fn longitude(p: S2Point) -> S1Angle {
        S1Angle::from_radians(Self::longitude_radians(p))
    }

    /// Internal method to compute latitude in radians from S2Point
    fn latitude_radians(p: S2Point) -> f64 {
        // Use atan2 for better accuracy near poles (matching C++ implementation)
        // The "+ 0.0" ensures -0.0 becomes +0.0 for consistent formatting
        let coords = p.coords();
        (coords.z + 0.0).atan2((coords.x * coords.x + coords.y * coords.y).sqrt())
    }

    /// Internal method to compute longitude in radians from S2Point
    fn longitude_radians(p: S2Point) -> f64 {
        // The "+ 0.0" ensures -0.0 becomes +0.0 for consistent formatting
        // atan2(0, 0) is defined to be 0
        let coords = p.coords();
        (coords.y + 0.0).atan2(coords.x + 0.0)
    }

    /// Get latitude as S1Angle
    #[inline]
    pub fn lat(&self) -> S1Angle {
        S1Angle::from_radians(self.coords.x())
    }

    /// Get longitude as S1Angle
    #[inline]
    pub fn lng(&self) -> S1Angle {
        S1Angle::from_radians(self.coords.y())
    }

    /// Get internal coordinates as R2Point
    #[inline]
    pub fn coords(&self) -> R2Point {
        self.coords
    }

    /// Check if coordinates are valid
    /// 
    /// Valid means latitude in [-90°, 90°] and longitude in [-180°, 180°]
    pub fn is_valid(&self) -> bool {
        let lat_rad = self.coords.x().abs();
        let lng_rad = self.coords.y().abs();
        lat_rad <= PI_2 && lng_rad <= PI && 
            self.coords.x().is_finite() && self.coords.y().is_finite()
    }

    /// Normalize coordinates to valid range
    /// 
    /// Clamps latitude to [-90°, 90°] and wraps longitude to [-180°, 180°]
    /// Returns Invalid() if coordinates contain NaN or infinity
    pub fn normalized(&self) -> Self {
        let lat_rad = self.coords.x();
        let lng_rad = self.coords.y();

        // Check for non-finite values
        if !lat_rad.is_finite() || !lng_rad.is_finite() {
            return Self::invalid();
        }

        // Normalize latitude to [-π/2, π/2]
        let mut norm_lat = lat_rad;
        let mut norm_lng = lng_rad;

        // Handle latitude clamping (matching C++ implementation)
        if norm_lat > PI_2 {
            norm_lat = PI_2;  // Clamp to 90 degrees
        } else if norm_lat < -PI_2 {
            norm_lat = -PI_2; // Clamp to -90 degrees
        }

        // Clamp latitude to exact bounds
        norm_lat = norm_lat.clamp(-PI_2, PI_2);

        // Normalize longitude to [-π, π]
        norm_lng = norm_lng.rem_euclid(2.0 * PI);
        if norm_lng > PI {
            norm_lng -= 2.0 * PI;
        }

        Self::from_radians(norm_lat, norm_lng)
    }

    /// Convert to S2Point
    pub fn to_point(&self) -> S2Result<S2Point> {
        let lat_rad = self.coords.x();
        let lng_rad = self.coords.y();

        // Check for finite values
        if !lat_rad.is_finite() || !lng_rad.is_finite() {
            return Err(S2Error::invalid_point("Non-finite lat/lng coordinates"));
        }

        // Convert to Cartesian coordinates
        let cos_lat = lat_rad.cos();
        let sin_lat = lat_rad.sin();
        let cos_lng = lng_rad.cos();
        let sin_lng = lng_rad.sin();

        let x = cos_lat * cos_lng;
        let y = cos_lat * sin_lng;
        let z = sin_lat;

        S2Point::new(x, y, z)
    }

    /// Get distance to another S2LatLng using Haversine formula
    ///
    /// This implementation matches the C++ GetDistance method and uses
    /// the Haversine formula for great-circle distance calculation.
    /// Both S2LatLng points should be normalized for best accuracy.
    pub fn get_distance(&self, other: &Self) -> S1Angle {
        let lat1 = self.coords.x();
        let lat2 = other.coords.x();
        let lng1 = self.coords.y();
        let lng2 = other.coords.y();

        let dlat = lat2 - lat1;
        let dlng = lng2 - lng1;

        let a = (dlat * 0.5).sin().powi(2) + 
                lat1.cos() * lat2.cos() * (dlng * 0.5).sin().powi(2);
        
        // Use atan2 for better numerical stability
        let distance = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
        
        S1Angle::from_radians(distance)
    }

    /// Check if approximately equal to another S2LatLng
    /// 
    /// This matches the C++ ApproxEquals behavior exactly
    pub fn approx_equals(&self, other: &Self, max_error: S1Angle) -> bool {
        self.coords.aequal(other.coords, max_error.radians())
    }

    /// Export coordinates as degrees separated by comma
    /// Format: "lat,lng" e.g. "47.6062,-122.3321"
    pub fn to_string_in_degrees(&self) -> String {
        // Use normalized coordinates for string output
        let normalized = if self.is_valid() { *self } else { self.normalized() };
        format!("{:.6},{:.6}", 
                normalized.lat().degrees(), 
                normalized.lng().degrees())
    }
}

impl Default for S2LatLng {
    fn default() -> Self {
        Self::default()
    }
}

// Arithmetic operations
impl Add for S2LatLng {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            coords: self.coords + other.coords,
        }
    }
}

impl Sub for S2LatLng {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            coords: self.coords - other.coords,
        }
    }
}

impl Mul<f64> for S2LatLng {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self {
        Self {
            coords: self.coords * scalar,
        }
    }
}

impl Mul<S2LatLng> for f64 {
    type Output = S2LatLng;

    fn mul(self, latlng: S2LatLng) -> S2LatLng {
        S2LatLng {
            coords: latlng.coords * self,
        }
    }
}

// Comparison operations
impl PartialOrd for S2LatLng {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // Use lexicographic ordering: first by latitude, then by longitude
        match self.coords.x().partial_cmp(&other.coords.x()) {
            Some(std::cmp::Ordering::Equal) => self.coords.y().partial_cmp(&other.coords.y()),
            other => other,
        }
    }
}

// Hash implementation for use in hash maps
impl Hash for S2LatLng {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash the bit representation for consistent results
        self.coords.x().to_bits().hash(state);
        self.coords.y().to_bits().hash(state);
    }
}

impl Eq for S2LatLng {
    // Eq implementation that is consistent with PartialEq
    // This is safe because PartialEq for S2LatLng uses exact equality
}

// Display implementation
impl fmt::Display for S2LatLng {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string_in_degrees())
    }
}

/// Legacy hash function for backwards compatibility (matches C++ S2LatLngHash)
#[derive(Default)]
pub struct S2LatLngHash;

impl S2LatLngHash {
    pub fn hash(&self, latlng: &S2LatLng) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        let mut hasher = DefaultHasher::new();
        latlng.hash(&mut hasher);
        hasher.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_basic_construction() {
        let ll = S2LatLng::from_degrees(45.0, 90.0);
        assert_relative_eq!(ll.lat().degrees(), 45.0, epsilon = 1e-15);
        assert_relative_eq!(ll.lng().degrees(), 90.0, epsilon = 1e-15);
        assert!(ll.is_valid());
    }

    #[test]
    fn test_invalid_construction() {
        let invalid = S2LatLng::from_degrees(-91.0, 0.0);
        assert!(!invalid.is_valid());
        
        let invalid2 = S2LatLng::from_degrees(0.0, 181.0);
        assert!(!invalid2.is_valid());
    }

    #[test]
    fn test_normalization() {
        let bad = S2LatLng::from_degrees(120.0, 200.0);
        assert!(!bad.is_valid());
        
        let good = bad.normalized();
        assert!(good.is_valid());
        assert_relative_eq!(good.lat().degrees(), 90.0, epsilon = 1e-13);
        assert_relative_eq!(good.lng().degrees(), -160.0, epsilon = 1e-13);
    }

    #[test]
    fn test_arithmetic() {
        let a = S2LatLng::from_degrees(10.0, 20.0);
        let b = S2LatLng::from_degrees(20.0, 30.0);
        
        let sum = a + b;
        assert_relative_eq!((sum).lat().degrees(), 30.0, epsilon = 1e-15);
        assert_relative_eq!((sum).lng().degrees(), 50.0, epsilon = 1e-15);
        
        let diff = a - b;
        assert_relative_eq!(diff.lat().degrees(), -10.0, epsilon = 1e-15);
        assert_relative_eq!(diff.lng().degrees(), -10.0, epsilon = 1e-15);
        
        let scaled = 0.5 * a;
        assert_relative_eq!(scaled.lat().degrees(), 5.0, epsilon = 1e-15);
        assert_relative_eq!(scaled.lng().degrees(), 10.0, epsilon = 1e-15);
    }

    #[test]
    fn test_point_conversion() {
        // Test conversion from S2Point to S2LatLng and back
        let original_point = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let ll = S2LatLng::from_point(original_point);
        let converted_back = ll.to_point().unwrap();
        
        assert_relative_eq!(
            original_point.coords().x, 
            converted_back.coords().x, 
            epsilon = 1e-15
        );
    }

    #[test] 
    fn test_distance() {
        let ll1 = S2LatLng::from_degrees(90.0, 0.0);  // North pole
        let ll2 = S2LatLng::from_degrees(90.0, 0.0);  // Same point
        
        let distance = ll1.get_distance(&ll2);
        assert_relative_eq!(distance.radians(), 0.0, epsilon = 1e-15);
        
        // Test a known distance
        let seattle = S2LatLng::from_degrees(47.6062, -122.3321);
        let vancouver = S2LatLng::from_degrees(49.2827, -123.1207);
        let dist = seattle.get_distance(&vancouver);
        
        // Distance should be around 1.9 degrees (rough verification)
        assert!(dist.degrees() > 1.0 && dist.degrees() < 3.0);
    }

    #[test]
    fn test_string_formatting() {
        let ll = S2LatLng::from_degrees(47.6062, -122.3321);
        let formatted = ll.to_string_in_degrees();
        
        // Should contain both coordinates separated by comma
        assert!(formatted.contains("47.606"));
        assert!(formatted.contains("-122.332"));
        assert!(formatted.contains(","));
    }

    #[test]
    fn test_special_cases() {
        // Test default constructor
        let default_ll = S2LatLng::default();
        assert!(default_ll.is_valid());
        assert_eq!(default_ll.lat().degrees(), 0.0);
        assert_eq!(default_ll.lng().degrees(), 0.0);

        // Test invalid constant
        let invalid = S2LatLng::invalid();
        assert!(!invalid.is_valid());

        // Test infinity/NaN handling
        let inf_ll = S2LatLng::from_degrees(f64::INFINITY, 0.0);
        assert!(!inf_ll.is_valid());
        assert!(!inf_ll.normalized().is_valid());

        let nan_ll = S2LatLng::from_degrees(f64::NAN, 0.0);
        assert!(!nan_ll.is_valid());
        assert!(!nan_ll.normalized().is_valid());
    }
}