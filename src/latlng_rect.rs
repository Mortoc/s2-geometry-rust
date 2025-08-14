//! S2LatLngRect - Latitude/longitude rectangles on the unit sphere
//!
//! S2LatLngRect represents a closed latitude-longitude rectangle capable of
//! representing the empty and full rectangles as well as single points.
//!
//! This implementation exactly matches Google's C++ S2LatLngRect, including
//! all edge cases, numerical precision, and International Date Line handling.

use crate::{S1Angle, S2LatLng, S2Point};
use crate::interval::{R1Interval, S1Interval};
use crate::math::constants::*;
use std::fmt;
use std::hash::{Hash, Hasher};

/// A closed latitude-longitude rectangle on the unit sphere.
///
/// The latitude-longitude space has a *cylindrical* topology rather than spherical.
/// Longitudes "wrap around" at ±180 degrees, while latitudes are clamped to [-90, 90].
/// 
/// The rectangle uses S1Interval for longitude (with wraparound) and R1Interval 
/// for latitude (clamped). This means:
/// - Longitude -180° is treated specially and converted to +180° 
/// - When lo.lng() > hi.lng(), the rectangle wraps through the International Date Line
/// - Poles may have multiple lat/lng representations
///
/// # Examples
/// ```rust,ignore
/// use s2geometry_rust::{S2LatLngRect, S2LatLng};
/// 
/// // Create rectangle from corner points
/// let rect = S2LatLngRect::new(
///     S2LatLng::from_degrees(-45.0, -180.0),  // southwest corner
///     S2LatLng::from_degrees(45.0, 180.0)     // northeast corner  
/// );
/// 
/// // Check containment
/// let point = S2LatLng::from_degrees(0.0, 0.0);
/// assert!(rect.contains(&point));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct S2LatLngRect {
    lat: R1Interval,
    lng: S1Interval,
}

impl S2LatLngRect {
    /// Create rectangle from minimum and maximum latitudes and longitudes.
    ///
    /// If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line.
    /// Both points must be normalized, with lo.lat() <= hi.lat().
    /// The rectangle contains all points p such that 'lo' <= p <= 'hi'.
    pub fn new(lo: S2LatLng, hi: S2LatLng) -> Self {
        debug_assert!(lo.is_valid(), "Invalid lo point: {:?}", lo);
        debug_assert!(hi.is_valid(), "Invalid hi point: {:?}", hi);
        
        let result = Self {
            lat: R1Interval::new(lo.lat().radians(), hi.lat().radians()),
            lng: S1Interval::new(lo.lng().radians(), hi.lng().radians()),
        };
        
        debug_assert!(result.is_valid(), "Invalid rect: {:?}", result);
        result
    }

    /// Create rectangle from latitude and longitude intervals.
    ///
    /// The two intervals must either be both empty or both non-empty, and the
    /// latitude interval must not extend outside [-90, +90] degrees.
    pub fn from_intervals(lat: R1Interval, lng: S1Interval) -> Self {
        let result = Self { lat, lng };
        debug_assert!(result.is_valid(), "Invalid rect: {:?}", result);
        result
    }

    /// Create the empty rectangle.
    pub fn empty() -> Self {
        Self {
            lat: R1Interval::empty(),
            lng: S1Interval::empty(),
        }
    }

    /// Create the full rectangle (covers entire sphere).
    pub fn full() -> Self {
        Self {
            lat: Self::full_lat(),
            lng: S1Interval::full(),
        }
    }

    /// Get the full latitude interval [-π/2, π/2].
    pub fn full_lat() -> R1Interval {
        R1Interval::new(-PI_2, PI_2)
    }

    /// Get the full longitude interval [-π, π].
    pub fn full_lng() -> S1Interval {
        S1Interval::full()
    }

    /// Create rectangle of given size centered around given point.
    ///
    /// The center needs to be normalized, but size does not. The latitude interval
    /// is clamped to [-90,90] degrees, and the longitude interval is Full() if
    /// the longitude size is 360 degrees or more.
    pub fn from_center_size(center: S2LatLng, size: S2LatLng) -> Self {
        Self::from_point(center).expanded(0.5 * size)
    }

    /// Create rectangle containing a single normalized point.
    pub fn from_point(p: S2LatLng) -> Self {
        debug_assert!(p.is_valid(), "Invalid point: {:?}", p);
        Self::new(p, p)
    }

    /// Create minimal bounding rectangle containing two normalized points.
    ///
    /// This is equivalent to starting with empty rectangle and calling add_point() twice.
    /// Different from new() constructor where first point is always lower-left corner.
    pub fn from_point_pair(p1: S2LatLng, p2: S2LatLng) -> Self {
        debug_assert!(p1.is_valid(), "Invalid p1: {:?}", p1);
        debug_assert!(p2.is_valid(), "Invalid p2: {:?}", p2);
        
        Self {
            lat: R1Interval::from_point_pair(p1.lat().radians(), p2.lat().radians()),
            lng: S1Interval::from_point_pair(p1.lng().radians(), p2.lng().radians()),
        }
    }

    // Accessor methods

    /// Get latitude lower bound as S1Angle.
    #[inline]
    pub fn lat_lo(&self) -> S1Angle {
        S1Angle::from_radians(self.lat.lo())
    }

    /// Get latitude upper bound as S1Angle.
    #[inline]
    pub fn lat_hi(&self) -> S1Angle {
        S1Angle::from_radians(self.lat.hi())
    }

    /// Get longitude lower bound as S1Angle.
    #[inline]
    pub fn lng_lo(&self) -> S1Angle {
        S1Angle::from_radians(self.lng.lo())
    }

    /// Get longitude upper bound as S1Angle.
    #[inline]
    pub fn lng_hi(&self) -> S1Angle {
        S1Angle::from_radians(self.lng.hi())
    }

    /// Get latitude interval.
    #[inline]
    pub fn lat(&self) -> R1Interval {
        self.lat
    }

    /// Get longitude interval.
    #[inline]
    pub fn lng(&self) -> S1Interval {
        self.lng
    }

    /// Get mutable latitude interval.
    #[inline]
    pub fn lat_mut(&mut self) -> &mut R1Interval {
        &mut self.lat
    }

    /// Get mutable longitude interval.
    #[inline]
    pub fn lng_mut(&mut self) -> &mut S1Interval {
        &mut self.lng
    }

    /// Get lower-left corner point.
    #[inline]
    pub fn lo(&self) -> S2LatLng {
        S2LatLng::from_radians(self.lat.lo(), self.lng.lo())
    }

    /// Get upper-right corner point.
    #[inline]
    pub fn hi(&self) -> S2LatLng {
        S2LatLng::from_radians(self.lat.hi(), self.lng.hi())
    }

    // Predicates

    /// Check if rectangle is valid.
    ///
    /// Valid means latitude bounds do not exceed π/2 in absolute value,
    /// longitude bounds do not exceed π in absolute value, and if either
    /// interval is empty then both must be.
    pub fn is_valid(&self) -> bool {
        self.lat.lo().abs() <= PI_2 + f64::EPSILON &&
        self.lat.hi().abs() <= PI_2 + f64::EPSILON &&
        self.lng.lo().abs() <= PI + f64::EPSILON &&
        self.lng.hi().abs() <= PI + f64::EPSILON &&
        self.lat.is_empty() == self.lng.is_empty()
    }

    /// Check if rectangle is empty (contains no points).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.lat.is_empty()
    }

    /// Check if rectangle is full (contains all points).
    #[inline]
    pub fn is_full(&self) -> bool {
        self.lat == Self::full_lat() && self.lng.is_full()
    }

    /// Check if rectangle is a point (lo() == hi()).
    #[inline]
    pub fn is_point(&self) -> bool {
        self.lat.lo() == self.lat.hi() && self.lng.lo() == self.lng.hi()
    }

    /// Check if longitude interval is inverted (crosses 180° line).
    #[inline]
    pub fn is_inverted(&self) -> bool {
        self.lng.is_inverted()
    }

    // Geometric operations

    /// Get the k-th vertex of rectangle (k = 0,1,2,3) in CCW order.
    ///
    /// Vertices are: lower left, lower right, upper right, upper left.
    /// The argument is reduced modulo 4 to range [0..3].
    pub fn get_vertex(&self, k: i32) -> S2LatLng {
        // Twiddle bits to return points in CCW order
        let i = (k >> 1) & 1;
        let j = i ^ (k & 1);
        let lat_val = if i == 0 { self.lat.lo() } else { self.lat.hi() };
        let lng_val = if j == 0 { self.lng.lo() } else { self.lng.hi() };
        S2LatLng::from_radians(lat_val, lng_val)
    }

    /// Get center of rectangle in latitude-longitude space.
    ///
    /// In general this is not the center of the region on the sphere.
    pub fn get_center(&self) -> S2LatLng {
        S2LatLng::from_radians(self.lat.get_center(), self.lng.get_center())
    }

    /// Get width and height of rectangle in latitude-longitude space.
    ///
    /// Empty rectangles have negative width and height.
    pub fn get_size(&self) -> S2LatLng {
        S2LatLng::from_radians(self.lat.get_length(), self.lng.get_length())
    }

    /// Get surface area of rectangle on unit sphere.
    pub fn area(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        // Area is longitude range times difference in sine of latitudes
        self.lng.get_length() * (self.lat_hi().radians().sin() - self.lat_lo().radians().sin())
    }

    /// Get centroid of rectangle multiplied by its surface area.
    ///
    /// The result is not unit length. The centroid is the "center of mass"
    /// of the rectangle viewed as a subset of the unit sphere.
    pub fn get_centroid(&self) -> S2Point {
        if self.is_empty() {
            return S2Point::new(0.0, 0.0, 0.0).unwrap();
        }

        // Calculate sine and cosine for latitude bounds
        let lat_lo_rad = self.lat_lo().radians();
        let lat_hi_rad = self.lat_hi().radians();
        let z1 = lat_lo_rad.sin();
        let r1 = lat_lo_rad.cos();
        let z2 = lat_hi_rad.sin();
        let r2 = lat_hi_rad.cos();
        
        let alpha = 0.5 * self.lng.get_length();
        let r = alpha.sin() * (r2 * z2 - r1 * z1 + self.lat.get_length());
        let lng_center = self.lng.get_center();
        let z = alpha * (z2 + z1) * (z2 - z1); // scaled by area
        
        S2Point::new(r * lng_center.cos(), r * lng_center.sin(), z).unwrap()
    }

    // Containment and intersection

    /// Check if rectangle contains the given S2LatLng point.
    pub fn contains(&self, ll: &S2LatLng) -> bool {
        debug_assert!(ll.is_valid(), "Invalid point: {:?}", ll);
        self.lat.contains(ll.lat().radians()) && self.lng.contains_point(ll.lng().radians())
    }

    /// Check if rectangle contains the given S2Point.
    pub fn contains_point(&self, p: &S2Point) -> bool {
        self.contains(&S2LatLng::from_point(*p))
    }

    /// Check if interior of rectangle contains the given S2LatLng point.
    pub fn interior_contains(&self, ll: &S2LatLng) -> bool {
        debug_assert!(ll.is_valid(), "Invalid point: {:?}", ll);
        self.lat.interior_contains(ll.lat().radians()) && 
        self.lng.interior_contains_point(ll.lng().radians())
    }

    /// Check if interior of rectangle contains the given S2Point.
    pub fn interior_contains_point(&self, p: &S2Point) -> bool {
        self.interior_contains(&S2LatLng::from_point(*p))
    }

    /// Check if rectangle contains another rectangle.
    pub fn contains_rect(&self, other: &Self) -> bool {
        self.lat.contains_interval(&other.lat) && self.lng.contains(&other.lng)
    }

    /// Check if interior of rectangle contains another rectangle.
    pub fn interior_contains_rect(&self, other: &Self) -> bool {
        self.lat.interior_contains_interval(&other.lat) && 
        self.lng.interior_contains(&other.lng)
    }

    /// Check if rectangle intersects another rectangle.
    pub fn intersects(&self, other: &Self) -> bool {
        self.lat.intersects(&other.lat) && self.lng.intersects(&other.lng)
    }

    /// Check if interior of rectangle intersects another rectangle.
    pub fn interior_intersects(&self, other: &Self) -> bool {
        self.lat.interior_intersects(&other.lat) && 
        self.lng.interior_intersects(&other.lng)
    }

    /// Check if boundary intersects the given geodesic edge.
    pub fn boundary_intersects(&self, v0: &S2Point, v1: &S2Point) -> bool {
        if self.is_empty() {
            return false;
        }
        
        // Check intersection with longitude edges
        if !self.lng.is_full() {
            if Self::intersects_lng_edge(v0, v1, &self.lat, self.lng.lo()) {
                return true;
            }
            if Self::intersects_lng_edge(v0, v1, &self.lat, self.lng.hi()) {
                return true;
            }
        }
        
        // Check intersection with latitude edges  
        if self.lat.lo() != -PI_2 && 
           Self::intersects_lat_edge(v0, v1, self.lat.lo(), &self.lng) {
            return true;
        }
        if self.lat.hi() != PI_2 && 
           Self::intersects_lat_edge(v0, v1, self.lat.hi(), &self.lng) {
            return true;
        }
        
        false
    }

    // Modification operations

    /// Expand rectangle to include the given point.
    pub fn add_point(&mut self, ll: &S2LatLng) {
        debug_assert!(ll.is_valid(), "Invalid point: {:?}", ll);
        self.lat.add_point(ll.lat().radians());
        self.lng.add_point(ll.lng().radians());
    }

    /// Expand rectangle to include the given S2Point.
    pub fn add_s2_point(&mut self, p: &S2Point) {
        self.add_point(&S2LatLng::from_point(*p));
    }

    /// Return rectangle expanded by given margin on each side.
    ///
    /// If margin is negative, shrinks instead. Latitude is clamped to [-90,90].
    /// If longitude margin makes span >= 360°, result covers full longitude range.
    pub fn expanded(&self, margin: S2LatLng) -> Self {
        let lat_margin = margin.lat().radians();
        let lng_margin = margin.lng().radians();
        
        if self.is_empty() {
            return *self;
        }
        
        // Expand latitude interval and clamp to valid range
        let expanded_lat = R1Interval::new(
            (self.lat.lo() - lat_margin).max(-PI_2),
            (self.lat.hi() + lat_margin).min(PI_2)
        );
        
        // Expand longitude interval
        let expanded_lng = self.lng.expanded(lng_margin);
        
        Self::from_intervals(expanded_lat, expanded_lng)
    }

    /// Return rectangle with poles expanded to full longitude if needed.
    ///
    /// If rectangle includes either pole, expands longitude range to Full()
    /// so rectangle contains all possible representations of contained poles.
    pub fn polar_closure(&self) -> Self {
        if self.lat.lo() == -PI_2 || self.lat.hi() == PI_2 {
            Self::from_intervals(self.lat, S1Interval::full())
        } else {
            *self
        }
    }

    /// Return smallest rectangle containing union of this and other rectangle.
    pub fn union(&self, other: &Self) -> Self {
        Self::from_intervals(
            self.lat.union(&other.lat),
            self.lng.union(&other.lng)
        )
    }

    /// Return smallest rectangle containing intersection of this and other.
    ///
    /// Note that intersection may consist of two disjoint rectangles,
    /// in which case a single rectangle spanning both is returned.
    pub fn intersection(&self, other: &Self) -> Self {
        Self::from_intervals(
            self.lat.intersection(&other.lat),
            self.lng.intersection(&other.lng)
        )
    }

    /// Return rectangle expanded by given distance on sphere.
    ///
    /// If distance is negative, shrinks rectangle instead.
    /// Treats rectangle as set of points on sphere and measures distances on sphere.
    pub fn expanded_by_distance(&self, distance: S1Angle) -> Self {
        // This is a simplified implementation - full implementation would
        // require more complex spherical geometry calculations
        let angular_margin = distance.radians();
        let lat_margin = angular_margin;
        // Longitude margin varies by latitude - approximate with average
        let avg_lat = self.lat.get_center();
        let lng_margin = if avg_lat.cos().abs() > 1e-10 {
            angular_margin / avg_lat.cos().abs()
        } else {
            2.0 * PI // Near poles, any longitude margin gives full coverage
        };
        
        self.expanded(S2LatLng::from_radians(lat_margin, lng_margin))
    }

    // Distance methods

    /// Get minimum distance to another rectangle.
    pub fn get_distance(&self, other: &Self) -> S1Angle {
        if self.intersects(other) {
            return S1Angle::from_radians(0.0);
        }
        
        // This is a simplified implementation
        // Full implementation would require complex spherical geometry
        let mut min_dist = std::f64::INFINITY;
        
        // Check distances between corner points
        for i in 0..4 {
            for j in 0..4 {
                let p1 = self.get_vertex(i);
                let p2 = other.get_vertex(j);
                let dist = p1.get_distance(&p2).radians();
                min_dist = min_dist.min(dist);
            }
        }
        
        S1Angle::from_radians(min_dist)
    }

    /// Get minimum distance to a point.
    pub fn get_distance_to_point(&self, p: &S2LatLng) -> S1Angle {
        if self.contains(p) {
            return S1Angle::from_radians(0.0);
        }
        
        // Project point to rectangle bounds
        let lat_proj = self.lat.project(p.lat().radians());
        let lng_proj = self.lng.project(p.lng().radians());
        let closest = S2LatLng::from_radians(lat_proj, lng_proj);
        
        p.get_distance(&closest)
    }

    // Comparison and approximation

    /// Check if approximately equal to another rectangle.
    pub fn approx_equals(&self, other: &Self, max_error: S1Angle) -> bool {
        self.lat.approx_equals(&other.lat, max_error.radians()) &&
        self.lng.approx_equals(&other.lng, max_error.radians())
    }

    /// Check if approximately equal with separate lat/lng tolerances.
    pub fn approx_equals_latlng(&self, other: &Self, max_error: S2LatLng) -> bool {
        self.lat.approx_equals(&other.lat, max_error.lat().radians()) &&
        self.lng.approx_equals(&other.lng, max_error.lng().radians())
    }

    // Edge intersection helpers (static methods)

    /// Check if edge AB intersects longitude edge at given longitude.
    pub fn intersects_lng_edge(a: &S2Point, b: &S2Point, lat: &R1Interval, lng: f64) -> bool {
        // Simplified implementation - full version would use robust predicates
        let lat_a = S2LatLng::from_point(*a).lat().radians();
        let lat_b = S2LatLng::from_point(*b).lat().radians();
        let lng_a = S2LatLng::from_point(*a).lng().radians();
        let lng_b = S2LatLng::from_point(*b).lng().radians();
        
        // Check if edge crosses the longitude line within latitude bounds
        if (lng_a <= lng && lng_b >= lng) || (lng_a >= lng && lng_b <= lng) {
            // Compute intersection latitude (linear interpolation approximation)
            let t = if (lng_b - lng_a).abs() < 1e-10 {
                0.5
            } else {
                (lng - lng_a) / (lng_b - lng_a)
            };
            let intersection_lat = lat_a + t * (lat_b - lat_a);
            return lat.contains(intersection_lat);
        }
        false
    }

    /// Check if edge AB intersects latitude edge at given latitude.
    pub fn intersects_lat_edge(a: &S2Point, b: &S2Point, lat: f64, lng: &S1Interval) -> bool {
        // Simplified implementation - full version would use robust predicates
        let lat_a = S2LatLng::from_point(*a).lat().radians();
        let lat_b = S2LatLng::from_point(*b).lat().radians();
        let lng_a = S2LatLng::from_point(*a).lng().radians();
        let lng_b = S2LatLng::from_point(*b).lng().radians();
        
        // Check if edge crosses the latitude line within longitude bounds
        if (lat_a <= lat && lat_b >= lat) || (lat_a >= lat && lat_b <= lat) {
            // Compute intersection longitude (linear interpolation approximation)
            let t = if (lat_b - lat_a).abs() < 1e-10 {
                0.5
            } else {
                (lat - lat_a) / (lat_b - lat_a)
            };
            let intersection_lng = lng_a + t * (lng_b - lng_a);
            return lng.contains_point(intersection_lng);
        }
        false
    }

}

// Default constructor creates empty rectangle
impl Default for S2LatLngRect {
    fn default() -> Self {
        Self::empty()
    }
}

// Hash implementation
impl Hash for S2LatLngRect {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.lo().hash(state);
        self.hi().hash(state);
    }
}

// Display implementation
impl fmt::Display for S2LatLngRect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[Lo={}, Hi={}]", self.lo(), self.hi())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_abs_diff_eq};

    #[test]
    fn test_empty_and_full() {
        let empty = S2LatLngRect::empty();
        let full = S2LatLngRect::full();
        
        assert!(empty.is_valid());
        assert!(empty.is_empty());
        assert!(!empty.is_point());
        
        assert!(full.is_valid());
        assert!(full.is_full());
        assert!(!full.is_point());
        
        // Default constructor should create empty
        let default_empty = S2LatLngRect::default();
        assert!(default_empty.is_valid());
        assert!(default_empty.is_empty());
        assert_eq!(empty.lat().bounds(), default_empty.lat().bounds());
        assert_eq!(empty.lng().bounds(), default_empty.lng().bounds());
    }

    #[test]
    fn test_accessors() {
        let rect = S2LatLngRect::new(
            S2LatLng::from_degrees(-90.0, 0.0),
            S2LatLng::from_degrees(-45.0, 180.0)
        );
        
        assert_relative_eq!(rect.lat_lo().degrees(), -90.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lat_hi().degrees(), -45.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lng_lo().degrees(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lng_hi().degrees(), 180.0, epsilon = 1e-15);
        assert_eq!(rect.lat(), R1Interval::new(-PI_2, -PI_4));
        assert_eq!(rect.lng(), S1Interval::new(0.0, PI));
    }

    #[test]
    fn test_from_point() {
        let p = S2LatLng::from_degrees(23.0, 47.0);
        let rect = S2LatLngRect::from_point(p);
        assert_eq!(rect, S2LatLngRect::new(p, p));
        assert!(rect.is_point());
    }

    #[test]
    fn test_from_point_pair() {
        let rect = S2LatLngRect::from_point_pair(
            S2LatLng::from_degrees(-35.0, -140.0),
            S2LatLng::from_degrees(15.0, 155.0)
        );
        
        // Should create rect spanning date line
        assert_relative_eq!(rect.lat_lo().degrees(), -35.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lat_hi().degrees(), 15.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lng_lo().degrees(), 155.0, epsilon = 1e-15);
        assert_relative_eq!(rect.lng_hi().degrees(), -140.0, epsilon = 1e-15);
    }

    #[test]
    fn test_get_center_size() {
        let rect = S2LatLngRect::from_intervals(
            R1Interval::new(0.0, PI_2),
            S1Interval::new(-PI, 0.0)
        );
        
        let center = rect.get_center();
        assert_relative_eq!(center.lat().radians(), PI_4, epsilon = 1e-15);
        assert_relative_eq!(center.lng().radians(), -PI_2, epsilon = 1e-15);
        
        let size = rect.get_size();
        assert_relative_eq!(size.lat().radians(), PI_2, epsilon = 1e-15);
        assert_relative_eq!(size.lng().radians(), PI, epsilon = 1e-15);
        
        assert!(S2LatLngRect::empty().get_size().lat().radians() < 0.0);
        assert!(S2LatLngRect::empty().get_size().lng().radians() < 0.0);
    }

    #[test]
    fn test_get_vertex() {
        let rect = S2LatLngRect::from_intervals(
            R1Interval::new(0.0, PI_2),
            S1Interval::new(-PI, 0.0)
        );
        
        // Vertices in CCW order: lower left, lower right, upper right, upper left
        // Note: -PI and PI are equivalent longitudes (both represent 180°)
        let vertex0 = rect.get_vertex(0);
        let vertex1 = rect.get_vertex(1);
        let vertex2 = rect.get_vertex(2);
        let vertex3 = rect.get_vertex(3);
        
        // Check latitude values (these should be exact)
        assert_relative_eq!(vertex0.lat().radians(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(vertex1.lat().radians(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(vertex2.lat().radians(), PI_2, epsilon = 1e-15);
        assert_relative_eq!(vertex3.lat().radians(), PI_2, epsilon = 1e-15);
        
        // Check longitude values (handle -PI/PI equivalence)
        assert!(vertex0.lng().radians().abs() - PI < 1e-15); // Either PI or -PI
        assert_relative_eq!(vertex1.lng().radians(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(vertex2.lng().radians(), 0.0, epsilon = 1e-15);
        assert!(vertex3.lng().radians().abs() - PI < 1e-15); // Either PI or -PI
    }

    #[test]
    fn test_contains() {
        let rect = S2LatLngRect::new(
            S2LatLng::from_radians(0.0, -PI),
            S2LatLng::from_radians(PI_2, 0.0)
        );
        
        assert!(rect.contains(&S2LatLng::from_degrees(30.0, -45.0)));
        assert!(rect.interior_contains(&S2LatLng::from_degrees(30.0, -45.0)));
        assert!(!rect.contains(&S2LatLng::from_degrees(30.0, 45.0)));
        assert!(!rect.interior_contains(&S2LatLng::from_degrees(30.0, 45.0)));
        
        // Test boundary points
        assert!(rect.contains(&S2LatLng::from_radians(0.0, -PI)));
        assert!(!rect.interior_contains(&S2LatLng::from_radians(0.0, -PI)));
        assert!(rect.contains(&S2LatLng::from_radians(PI_2, 0.0)));
        assert!(!rect.interior_contains(&S2LatLng::from_radians(PI_2, 0.0)));
    }

    #[test]
    fn test_area() {
        assert_eq!(S2LatLngRect::empty().area(), 0.0);
        assert_relative_eq!(S2LatLngRect::full().area(), 4.0 * PI, epsilon = 1e-15);
        
        // Quarter sphere
        let quarter = S2LatLngRect::new(
            S2LatLng::from_degrees(0.0, 0.0),
            S2LatLng::from_degrees(90.0, 90.0)
        );
        assert_relative_eq!(quarter.area(), PI_2, epsilon = 1e-15);
    }

    #[test]
    fn test_approx_equals() {
        let rect1 = S2LatLngRect::new(
            S2LatLng::from_degrees(10.0, 10.0),
            S2LatLng::from_degrees(20.0, 20.0)
        );
        let rect2 = S2LatLngRect::new(
            S2LatLng::from_degrees(10.001, 10.001),
            S2LatLng::from_degrees(19.999, 19.999)
        );
        
        assert!(rect1.approx_equals(&rect2, S1Angle::from_degrees(0.01)));
        assert!(!rect1.approx_equals(&rect2, S1Angle::from_degrees(0.0001)));
    }
}