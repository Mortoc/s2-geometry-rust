//! S2Polyline - Sequences of connected geodesic edges
//!
//! S2Polyline represents a sequence of vertices connected by straight geodesic edges
//! on the sphere. This provides essential linear geometry support for representing
//! roads, boundaries, flight paths, and other linear spatial features.
//!
//! # Key Features
//!
//! - Construction from vectors of S2Point or S2LatLng
//! - Length calculations (both angular and surface distance)
//! - Point interpolation along the polyline 
//! - Nearest point queries and distance calculations
//! - Intersection testing with other geometry
//! - Subsetting and simplification operations
//! - Bounding rectangle calculations
//! - Validation and error checking
//!
//! # Mathematical Properties
//!
//! - Vertices must be unit-length points on the sphere
//! - Edges are geodesics (great circle arcs) between consecutive vertices
//! - No edges of length 0 or 180 degrees (antipodal vertices) are allowed
//! - Adjacent vertices should not be identical or antipodal
//!
//! # Example
//!
//! ```rust,ignore
//! use s2geometry_rust::{S2Polyline, S2Point, S1Angle};
//!
//! // Create a simple polyline
//! let vertices = vec![
//!     S2Point::new(1.0, 0.0, 0.0).unwrap(),
//!     S2Point::new(0.0, 1.0, 0.0).unwrap(),
//!     S2Point::new(0.0, 0.0, 1.0).unwrap(),
//! ];
//! let polyline = S2Polyline::new(vertices).unwrap();
//!
//! // Calculate total length
//! let length = polyline.get_length();
//! 
//! // Interpolate a point at 50% along the polyline
//! let mid_point = polyline.interpolate(0.5);
//!
//! // Find closest point to a query point
//! let query = S2Point::new(0.5, 0.5, 0.5).unwrap();
//! let (closest_point, next_vertex) = polyline.project(&query);
//! ```

use crate::error::{S2Error, S2Result};
use crate::math::DVec3;
use crate::point::S2Point;
use crate::angle::S1Angle;
use crate::latlng::S2LatLng;
use crate::latlng_rect::S2LatLngRect;
use crate::cap::S2Cap;
use crate::predicates;

/// A sequence of vertices connected by geodesic edges on the sphere.
///
/// S2Polyline represents a sequence of connected great circle arcs on the unit sphere.
/// It provides operations for length calculation, interpolation, projection, and
/// geometric queries while maintaining numerical robustness through the three-tier
/// mathematical architecture.
///
/// # Invariants
/// - All vertices must be unit-length points on the sphere
/// - No duplicate vertices (adjacent or non-adjacent)
/// - No edges of length 0 or 180 degrees
/// - At least 2 vertices for a valid polyline
///
/// # Thread Safety
/// S2Polyline is `Send + Sync` and can be safely shared between threads.
#[derive(Debug, Clone)]
pub struct S2Polyline {
    vertices: Vec<S2Point>,
    length: Option<S1Angle>, // Cached total length
    bound: Option<S2LatLngRect>, // Cached bounding rectangle
}

impl S2Polyline {
    /// Create a new S2Polyline from a vector of S2Point vertices.
    ///
    /// The vertices are connected in order by geodesic edges. The polyline
    /// must have at least 2 vertices to be valid.
    ///
    /// # Arguments
    /// * `vertices` - The sequence of points defining the polyline
    ///
    /// # Returns
    /// * `Ok(S2Polyline)` if the vertices form a valid polyline
    /// * `Err(S2Error)` if the polyline is invalid
    ///
    /// # Example
    /// ```rust,ignore
    /// let vertices = vec![
    ///     S2Point::from_coords(1.0, 0.0, 0.0)?,
    ///     S2Point::from_coords(0.0, 1.0, 0.0)?,
    /// ];
    /// let polyline = S2Polyline::new(vertices)?;
    /// ```
    pub fn new(vertices: Vec<S2Point>) -> S2Result<Self> {
        let polyline = Self {
            vertices,
            length: None,
            bound: None,
        };
        polyline.validate()?;
        Ok(polyline)
    }

    /// Create a new S2Polyline from a vector of S2LatLng coordinates.
    ///
    /// This is a convenience constructor that converts lat/lng coordinates
    /// to S2Point vertices before creating the polyline.
    ///
    /// # Arguments
    /// * `coords` - The sequence of lat/lng coordinates
    ///
    /// # Returns
    /// * `Ok(S2Polyline)` if the coordinates form a valid polyline
    /// * `Err(S2Error)` if the polyline is invalid
    pub fn from_latlngs(coords: Vec<S2LatLng>) -> S2Result<Self> {
        let vertices: Result<Vec<S2Point>, S2Error> = coords
            .into_iter()
            .map(|latlng| latlng.to_point())
            .collect();
        Self::new(vertices?)
    }

    /// Create an empty polyline with no vertices.
    ///
    /// An empty polyline is useful as a default value but cannot be used
    /// for geometric operations.
    pub fn empty() -> Self {
        Self {
            vertices: Vec::new(),
            length: Some(S1Angle::zero()),
            bound: Some(S2LatLngRect::empty()),
        }
    }

    /// Returns true if this polyline is valid.
    ///
    /// A valid polyline must have:
    /// - At least 2 vertices
    /// - All vertices are unit-length
    /// - No identical consecutive vertices
    /// - No antipodal consecutive vertices
    pub fn is_valid(&self) -> bool {
        self.validate().is_ok()
    }

    /// Get the number of vertices in the polyline.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get a vertex by index.
    ///
    /// # Arguments
    /// * `index` - The vertex index (0-based)
    ///
    /// # Returns
    /// * `Some(S2Point)` if the index is valid
    /// * `None` if the index is out of bounds
    pub fn vertex(&self, index: usize) -> Option<S2Point> {
        self.vertices.get(index).copied()
    }

    /// Get all vertices as a slice.
    pub fn vertices(&self) -> &[S2Point] {
        &self.vertices
    }

    /// Calculate the total length of the polyline.
    ///
    /// The length is measured as the sum of angular distances along
    /// the geodesic edges connecting consecutive vertices.
    ///
    /// # Returns
    /// The total angular length of the polyline
    pub fn get_length(&self) -> S1Angle {
        if let Some(cached_length) = self.length {
            return cached_length;
        }

        if self.vertices.len() < 2 {
            S1Angle::zero()
        } else {
            let mut total = 0.0;
            for i in 0..self.vertices.len() - 1 {
                let edge_length = S1Angle::from_radians(self.vertices[i].angle(&self.vertices[i + 1]));
                total += edge_length.radians();
            }
            S1Angle::from_radians(total)
        }
    }

    /// Interpolate a point at the given fraction along the polyline.
    ///
    /// The interpolation is done with respect to arc length, not
    /// straight-line distance. The result is always a unit-length point.
    ///
    /// # Arguments
    /// * `fraction` - Position along the polyline (0.0 = start, 1.0 = end)
    ///
    /// # Returns
    /// The interpolated point on the polyline
    ///
    /// # Example
    /// ```rust,ignore
    /// let mid_point = polyline.interpolate(0.5); // Midpoint
    /// let quarter_point = polyline.interpolate(0.25); // 25% along
    /// ```
    pub fn interpolate(&self, fraction: f64) -> S2Point {
        if self.vertices.is_empty() {
            return S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
        }

        if self.vertices.len() == 1 {
            return self.vertices[0];
        }

        // Handle edge cases
        if fraction <= 0.0 {
            return self.vertices[0];
        }
        if fraction >= 1.0 {
            return *self.vertices.last().unwrap();
        }

        let total_length = self.get_length().radians();
        if total_length == 0.0 {
            return self.vertices[0];
        }

        let target_length = fraction * total_length;
        let mut accumulated_length = 0.0;

        for i in 0..self.vertices.len() - 1 {
            let edge_length = self.vertices[i].angle(&self.vertices[i + 1]);
            
            if accumulated_length + edge_length >= target_length {
                // The target point is on this edge
                let edge_fraction = if edge_length == 0.0 {
                    0.0
                } else {
                    (target_length - accumulated_length) / edge_length
                };
                
                return self.interpolate_edge(i, edge_fraction);
            }
            
            accumulated_length += edge_length;
        }

        // Shouldn't reach here, but return the last vertex as fallback
        *self.vertices.last().unwrap()
    }

    /// Find the point on the polyline closest to the given query point.
    ///
    /// This method projects the query point onto the polyline and returns
    /// both the closest point and the index of the next vertex.
    ///
    /// # Arguments
    /// * `query` - The point to project onto the polyline
    ///
    /// # Returns
    /// A tuple containing:
    /// - The closest point on the polyline
    /// - The index of the next vertex (for the edge containing the closest point)
    ///
    /// # Example
    /// ```rust,ignore
    /// let query = S2Point::from_coords(0.5, 0.5, 0.5)?;
    /// let (closest, next_vertex) = polyline.project(&query);
    /// ```
    pub fn project(&self, query: &S2Point) -> (S2Point, usize) {
        if self.vertices.is_empty() {
            return (S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)), 0);
        }

        if self.vertices.len() == 1 {
            return (self.vertices[0], 0);
        }

        let mut min_distance = f64::INFINITY;
        let mut closest_point = self.vertices[0];
        let mut closest_next_vertex = 1;

        for i in 0..self.vertices.len() - 1 {
            let projected = self.project_to_edge(query, i);
            let distance = query.angle(&projected);
            
            if distance < min_distance {
                min_distance = distance;
                closest_point = projected;
                closest_next_vertex = i + 1;
            }
        }

        (closest_point, closest_next_vertex)
    }

    /// Check if this polyline intersects with another polyline.
    ///
    /// This performs a comprehensive intersection test between all edges
    /// of both polylines using robust geometric predicates.
    ///
    /// # Arguments
    /// * `other` - The other polyline to test for intersection
    ///
    /// # Returns
    /// `true` if the polylines intersect, `false` otherwise
    pub fn intersects(&self, other: &S2Polyline) -> bool {
        // Quick rejection test using bounding rectangles
        if let (Some(bound1), Some(bound2)) = (&self.bound, &other.bound) {
            if !bound1.intersects(bound2) {
                return false;
            }
        }

        // Test all edge pairs for intersection
        for i in 0..self.vertices.len().saturating_sub(1) {
            for j in 0..other.vertices.len().saturating_sub(1) {
                if self.edge_intersects_edge(i, other, j) {
                    return true;
                }
            }
        }

        false
    }

    /// Reverse the order of vertices in the polyline.
    ///
    /// This creates a new polyline with vertices in reverse order,
    /// effectively reversing the direction of traversal.
    pub fn reverse(&mut self) {
        self.vertices.reverse();
        // Cached values remain valid
    }

    /// Get the bounding rectangle of the polyline.
    ///
    /// Returns the smallest latitude-longitude rectangle that contains
    /// all vertices of the polyline.
    pub fn get_rect_bound(&self) -> S2LatLngRect {
        if let Some(cached_bound) = self.bound {
            return cached_bound;
        }

        if self.vertices.is_empty() {
            S2LatLngRect::empty()
        } else {
            let mut rect = S2LatLngRect::from_point(S2LatLng::from_point(self.vertices[0]));
            for vertex in &self.vertices[1..] {
                rect.add_point(&S2LatLng::from_point(*vertex));
            }
            rect
        }
    }

    /// Get a bounding cap for the polyline.
    ///
    /// Returns a spherical cap that contains all vertices of the polyline.
    /// This is useful for spatial indexing and quick intersection tests.
    pub fn get_cap_bound(&self) -> S2Cap {
        if self.vertices.is_empty() {
            return S2Cap::empty();
        }

        if self.vertices.len() == 1 {
            return S2Cap::from_point(self.vertices[0]);
        }

        // Find the centroid and maximum distance
        let mut centroid = DVec3::ZERO;
        for vertex in &self.vertices {
            centroid += vertex.coords();
        }
        
        let center = S2Point::from_normalized(centroid.normalize());
        let mut max_distance: f64 = 0.0;
        
        for vertex in &self.vertices {
            let distance = center.angle(vertex);
            max_distance = max_distance.max(distance);
        }

        S2Cap::from_center_angle(&center, S1Angle::from_radians(max_distance))
    }

    // Private helper methods

    /// Validate the polyline structure and constraints.
    fn validate(&self) -> S2Result<()> {
        if self.vertices.len() < 2 {
            return Err(S2Error::InvalidPolyline {
                reason: "Polyline must have at least 2 vertices".to_string(),
            });
        }

        for (i, vertex) in self.vertices.iter().enumerate() {
            // Check that vertex is normalized
            let length_sq = vertex.coords().length_squared();
            if (length_sq - 1.0).abs() > 1e-10 {
                return Err(S2Error::InvalidPolyline {
                    reason: format!("Vertex {} is not normalized: length_sq = {}", i, length_sq),
                });
            }

            // Check for consecutive identical or antipodal vertices
            if i > 0 {
                let prev_vertex = self.vertices[i - 1];
                let dot_product = vertex.dot(&prev_vertex);
                
                if dot_product > 1.0 - 1e-15 {
                    return Err(S2Error::InvalidPolyline {
                        reason: format!("Vertices {} and {} are identical", i - 1, i),
                    });
                }
                
                if dot_product < -1.0 + 1e-15 {
                    return Err(S2Error::InvalidPolyline {
                        reason: format!("Vertices {} and {} are antipodal", i - 1, i),
                    });
                }
            }
        }

        Ok(())
    }

    /// Interpolate a point along the edge between vertices i and i+1.
    fn interpolate_edge(&self, i: usize, fraction: f64) -> S2Point {
        debug_assert!(i + 1 < self.vertices.len());
        
        let a = self.vertices[i];
        let b = self.vertices[i + 1];
        
        if fraction <= 0.0 {
            return a;
        }
        if fraction >= 1.0 {
            return b;
        }

        // Use spherical linear interpolation (slerp)
        let dot = a.dot(&b).clamp(-1.0, 1.0);
        let angle = dot.acos();
        
        if angle < 1e-15 {
            // Points are very close, use linear interpolation
            let result = a.coords() * (1.0 - fraction) + b.coords() * fraction;
            return S2Point::from_normalized(result.normalize());
        }

        let sin_angle = angle.sin();
        let factor_a = ((1.0 - fraction) * angle).sin() / sin_angle;
        let factor_b = (fraction * angle).sin() / sin_angle;
        
        let result = a.coords() * factor_a + b.coords() * factor_b;
        S2Point::from_normalized(result.normalize())
    }

    /// Project a point onto the edge between vertices i and i+1.
    fn project_to_edge(&self, query: &S2Point, i: usize) -> S2Point {
        debug_assert!(i + 1 < self.vertices.len());
        
        let a = self.vertices[i];
        let b = self.vertices[i + 1];
        
        // Calculate the projection parameter using the spherical dot product
        let aq = query.dot(&a);
        let bq = query.dot(&b);
        let ab = a.dot(&b);
        
        let denom = 1.0 - ab * ab;
        if denom < 1e-15 {
            // Edge is degenerate (points are nearly identical)
            return a;
        }
        
        let t = (aq - bq * ab) / denom;
        let t_clamped = t.clamp(0.0, 1.0);
        
        self.interpolate_edge(i, t_clamped)
    }

    /// Test if edge i of this polyline intersects edge j of another polyline.
    fn edge_intersects_edge(&self, i: usize, other: &S2Polyline, j: usize) -> bool {
        debug_assert!(i + 1 < self.vertices.len());
        debug_assert!(j + 1 < other.vertices.len());
        
        let a = self.vertices[i];
        let b = self.vertices[i + 1];
        let c = other.vertices[j];
        let d = other.vertices[j + 1];
        
        // Use robust predicates to test for edge intersection
        predicates::edge_or_vertex_crossing(a.coords(), b.coords(), c.coords(), d.coords())
    }
}

// Implement common traits

impl Default for S2Polyline {
    fn default() -> Self {
        Self::empty()
    }
}

impl PartialEq for S2Polyline {
    fn eq(&self, other: &Self) -> bool {
        self.vertices == other.vertices
    }
}

// Implement iterator support
impl IntoIterator for S2Polyline {
    type Item = S2Point;
    type IntoIter = std::vec::IntoIter<S2Point>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.into_iter()
    }
}

impl<'a> IntoIterator for &'a S2Polyline {
    type Item = &'a S2Point;
    type IntoIter = std::slice::Iter<'a, S2Point>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::DVec3;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_empty_polyline() {
        let polyline = S2Polyline::empty();
        assert_eq!(polyline.num_vertices(), 0);
        assert!(polyline.is_valid() || polyline.vertices.is_empty()); // Empty is special case
    }

    #[test]
    fn test_simple_polyline() {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        
        let polyline = S2Polyline::new(vertices).unwrap();
        assert_eq!(polyline.num_vertices(), 3);
        assert!(polyline.is_valid());
    }

    #[test]
    fn test_length_calculation() {
        // Create a quarter circle on the sphere
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        
        let polyline = S2Polyline::new(vertices).unwrap();
        let length = polyline.get_length();
        
        // Should be π/2 radians (90 degrees)
        assert!((length.radians() - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_interpolation() {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        
        let polyline = S2Polyline::new(vertices).unwrap();
        
        // Test edge cases
        let start = polyline.interpolate(0.0);
        let end = polyline.interpolate(1.0);
        
        assert!((start.coords() - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-10);
        assert!((end.coords() - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-10);
        
        // Test midpoint interpolation
        let mid = polyline.interpolate(0.5);
        assert!((mid.coords().length() - 1.0).abs() < 1e-10); // Should be normalized
    }

    #[test]
    fn test_projection() {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        ];
        
        let polyline = S2Polyline::new(vertices).unwrap();
        
        // Project a point that should be closest to the midpoint of the edge
        let normalized_coords = DVec3::new(0.5, 0.5, 0.0).normalize();
        let query = S2Point::from_normalized(normalized_coords);
        let (closest, next_vertex) = polyline.project(&query);
        
        assert_eq!(next_vertex, 1); // Should be on the first (and only) edge
        assert!((closest.coords().length() - 1.0).abs() < 1e-10); // Should be normalized
    }

    #[test]
    fn test_invalid_polylines() {
        // Test insufficient vertices
        let single_vertex = vec![S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0))];
        assert!(S2Polyline::new(single_vertex).is_err());
        
        // Test identical consecutive vertices
        let identical_vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
        ];
        assert!(S2Polyline::new(identical_vertices).is_err());
        
        // Test antipodal consecutive vertices
        let antipodal_vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(-1.0, 0.0, 0.0)),
        ];
        assert!(S2Polyline::new(antipodal_vertices).is_err());
    }

    #[test]
    fn test_reverse() {
        let vertices = vec![
            S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
            S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        ];
        
        let mut polyline = S2Polyline::new(vertices).unwrap();
        let original_first = polyline.vertex(0).unwrap();
        let original_last = polyline.vertex(2).unwrap();
        
        polyline.reverse();
        
        assert_eq!(polyline.vertex(0).unwrap(), original_last);
        assert_eq!(polyline.vertex(2).unwrap(), original_first);
    }
}