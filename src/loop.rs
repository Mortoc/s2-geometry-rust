//! S2Loop - Simple spherical polygons
//!
//! S2Loop represents a simple closed polygon on the sphere. It consists of a 
//! single chain of vertices where the first vertex is implicitly connected to 
//! the last. All loops are defined to have a CCW orientation, i.e. the interior 
//! of the loop is on the left side of the edges.

use crate::error::{S2Error, S2Result};
use crate::math::{DVec3, predicates::*};
use crate::point::S2Point;
use crate::angle::S1Angle;
use crate::chord_angle::S1ChordAngle;
use crate::latlng_rect::S2LatLngRect;
use crate::cap::S2Cap;
use crate::cell::S2Cell;
use std::f64::consts::PI;

/// A simple spherical polygon represented as a chain of vertices.
///
/// S2Loop represents a simple closed polygon on the unit sphere. The interior
/// is defined to be on the left side of the edges (CCW orientation). Loops 
/// must have at least 3 vertices except for special empty and full loops.
///
/// # Special Cases
/// - Empty loop: Contains no points (represented with a single vertex)
/// - Full loop: Contains all points (represented with a single vertex)
///
/// # Invariants
/// - No duplicate vertices (adjacent or non-adjacent)
/// - No self-intersecting edges
/// - No edges of length 180 degrees (antipodal vertices)
/// - At least 3 vertices (except empty/full)
///
/// # Example
/// ```rust,ignore
/// use s2geometry_rust::{S2Loop, S2Point};
/// 
/// let vertices = vec![
///     S2Point::new(1.0, 0.0, 0.0).unwrap(),
///     S2Point::new(0.0, 1.0, 0.0).unwrap(), 
///     S2Point::new(0.0, 0.0, 1.0).unwrap(),
/// ];
/// let loop_poly = S2Loop::new(vertices).unwrap();
/// assert!(loop_poly.is_valid());
/// ```
#[derive(Debug, Clone)]
pub struct S2Loop {
    vertices: Vec<S2Point>,
    depth: i32,
    origin_inside: bool,
    bound: S2LatLngRect,
}

impl S2Loop {
    /// Create a new S2Loop from a vector of vertices.
    ///
    /// The vertices should be in counter-clockwise order around the interior.
    /// The last vertex is implicitly connected to the first.
    pub fn new(vertices: Vec<S2Point>) -> S2Result<Self> {
        let mut loop_obj = Self {
            vertices,
            depth: 0,
            origin_inside: false,
            bound: S2LatLngRect::empty(),
        };
        loop_obj.init()?;
        Ok(loop_obj)
    }

    /// Create an empty loop (contains no points).
    pub fn empty() -> Self {
        Self {
            vertices: vec![Self::empty_vertex()],
            depth: 0,
            origin_inside: false,
            bound: S2LatLngRect::empty(),
        }
    }

    /// Create a full loop (contains all points).
    pub fn full() -> Self {
        Self {
            vertices: vec![Self::full_vertex()],
            depth: 0,
            origin_inside: true,
            bound: S2LatLngRect::full(),
        }
    }

    /// Create an S2Loop from an S2Cell.
    pub fn from_cell(cell: &S2Cell) -> Self {
        let vertices = vec![
            cell.get_vertex(0),
            cell.get_vertex(1), 
            cell.get_vertex(2),
            cell.get_vertex(3),
        ];
        Self::new(vertices).expect("Cell vertices should form valid loop")
    }

    /// Special vertex for empty loop (northern hemisphere).
    fn empty_vertex() -> S2Point {
        S2Point::from_coords_raw(0.0, 0.0, 1.0)
    }

    /// Special vertex for full loop (southern hemisphere).
    fn full_vertex() -> S2Point {
        S2Point::from_coords_raw(0.0, 0.0, -1.0)
    }

    /// Initialize the loop after setting vertices.
    fn init(&mut self) -> S2Result<()> {
        // Validate the loop
        self.validate()?;
        
        // Initialize origin containment and bounds
        self.init_origin_and_bound();
        
        Ok(())
    }

    /// Validate that this is a well-formed loop.
    fn validate(&self) -> S2Result<()> {
        // Handle special cases
        if self.is_empty_or_full() {
            return Ok(());
        }

        // Must have at least 3 vertices
        if self.vertices.len() < 3 {
            return Err(S2Error::invalid_loop("Loop must have at least 3 vertices"));
        }

        // Check for duplicate vertices
        for i in 0..self.vertices.len() {
            for j in (i + 1)..self.vertices.len() {
                if self.vertices[i].coords().abs_diff_eq(self.vertices[j].coords(), f64::EPSILON) {
                    return Err(S2Error::invalid_loop("Loop has duplicate vertices"));
                }
            }
        }

        // Check for degenerate edges
        for i in 0..self.vertices.len() {
            let next = (i + 1) % self.vertices.len();
            if self.vertices[i].coords().abs_diff_eq(self.vertices[next].coords(), f64::EPSILON) {
                return Err(S2Error::invalid_loop("Loop has degenerate edge"));
            }
        }

        // Check for self-intersections (simplified)
        // Full implementation would use robust edge crossing tests
        for i in 0..self.vertices.len() {
            let next_i = (i + 1) % self.vertices.len();
            for j in (i + 2)..self.vertices.len() {
                if j == self.vertices.len() - 1 && i == 0 {
                    continue; // Skip adjacent edges
                }
                let next_j = (j + 1) % self.vertices.len();
                
                let crossing = crossing_sign(
                    self.vertices[i].coords(),
                    self.vertices[next_i].coords(),
                    self.vertices[j].coords(), 
                    self.vertices[next_j].coords()
                );
                
                if crossing > 0 {
                    return Err(S2Error::invalid_loop("Loop edges cross"));
                }
            }
        }

        Ok(())
    }

    /// Initialize origin containment and bounding rectangle.
    fn init_origin_and_bound(&mut self) {
        if self.is_empty_or_full() {
            self.origin_inside = self.vertices[0].z() < 0.0;
            self.bound = if self.is_empty() {
                S2LatLngRect::empty()
            } else {
                S2LatLngRect::full()
            };
            return;
        }

        // Compute origin containment using robust predicates
        self.origin_inside = self.contains_origin();
        
        // Compute bounding rectangle
        self.bound = self.compute_rect_bound();
    }

    /// Check if the loop contains the origin using point-in-polygon test.
    fn contains_origin(&self) -> bool {
        if self.is_empty_or_full() {
            return self.vertices[0].z() < 0.0;
        }

        // Use winding number / sign accumulation
        let origin = DVec3::new(0.0, 0.0, 1.0);
        let mut sign_sum = 0;
        
        for i in 0..self.vertices.len() {
            let next = (i + 1) % self.vertices.len();
            sign_sum += robust_sign(
                origin,
                self.vertices[i].coords(),
                self.vertices[next].coords()
            );
        }
        
        sign_sum != 0
    }

    /// Compute the bounding rectangle for this loop.
    fn compute_rect_bound(&self) -> S2LatLngRect {
        if self.is_empty() {
            return S2LatLngRect::empty();
        }
        if self.is_full() {
            return S2LatLngRect::full();
        }

        let mut rect = S2LatLngRect::empty();
        for vertex in &self.vertices {
            let latlng = crate::S2LatLng::from_point(*vertex);
            rect.add_point(&latlng);
        }
        
        // Expand for edges that may extend beyond vertices
        // Simplified implementation - full version would check edge crossings
        let margin = crate::S2LatLng::from_radians(1e-15, 1e-15);
        rect.expanded(margin)
    }

    /// Number of vertices in the loop.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get vertex at index i. Index wrapping is supported.
    pub fn vertex(&self, i: usize) -> S2Point {
        self.vertices[i % self.vertices.len()]
    }

    /// Get all vertices as a slice.
    pub fn vertices(&self) -> &[S2Point] {
        &self.vertices
    }

    /// Check if this is the empty loop.
    pub fn is_empty(&self) -> bool {
        self.is_empty_or_full() && !self.origin_inside
    }

    /// Check if this is the full loop.
    pub fn is_full(&self) -> bool {
        self.is_empty_or_full() && self.origin_inside
    }

    /// Check if this is either empty or full.
    pub fn is_empty_or_full(&self) -> bool {
        self.vertices.len() == 1
    }

    /// Get/set the nesting depth (used by S2Polygon).
    pub fn depth(&self) -> i32 {
        self.depth
    }

    pub fn set_depth(&mut self, depth: i32) {
        self.depth = depth;
    }

    /// Check if this loop represents a hole.
    pub fn is_hole(&self) -> bool {
        (self.depth & 1) != 0
    }

    /// Get the sign of the loop (+1 for shell, -1 for hole).
    pub fn sign(&self) -> i32 {
        if self.is_hole() { -1 } else { 1 }
    }

    /// Check if the loop is normalized (area <= 2π).
    pub fn is_normalized(&self) -> bool {
        self.get_area() <= 2.0 * PI
    }

    /// Normalize the loop so its area is at most 2π.
    pub fn normalize(&mut self) {
        if !self.is_normalized() {
            self.invert();
        }
    }

    /// Invert the loop (complement the region).
    pub fn invert(&mut self) {
        if self.is_empty_or_full() {
            self.vertices[0] = if self.is_empty() {
                Self::full_vertex()
            } else {
                Self::empty_vertex()
            };
            self.origin_inside = !self.origin_inside;
            self.bound = if self.is_empty() {
                S2LatLngRect::empty()
            } else {
                S2LatLngRect::full()
            };
        } else {
            self.vertices.reverse();
            self.origin_inside = !self.origin_inside;
            // Bound stays the same for normal loops
        }
    }

    /// Get the area of the loop interior.
    pub fn get_area(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        if self.is_full() {
            return 4.0 * PI;
        }

        // Use the standard spherical polygon area formula
        let mut area = 0.0;
        for i in 0..self.vertices.len() {
            let next = (i + 1) % self.vertices.len();
            area += robust_sign(
                DVec3::ZERO,
                self.vertices[i].coords(),
                self.vertices[next].coords()
            ) as f64 * self.vertices[i].coords().dot(self.vertices[next].coords()).acos();
        }
        
        (area.abs() - (self.vertices.len() as f64 - 2.0) * PI).abs()
    }

    /// Get the centroid of the loop (weighted by area).
    pub fn get_centroid(&self) -> S2Point {
        if self.is_empty() {
            return S2Point::from_coords_raw(0.0, 0.0, 0.0);
        }
        if self.is_full() {
            return S2Point::from_coords_raw(0.0, 0.0, 0.0);
        }

        // Simplified centroid calculation
        let mut centroid = DVec3::ZERO;
        for vertex in &self.vertices {
            centroid += vertex.coords();
        }
        
        if centroid.length_squared() > 0.0 {
            S2Point::from_normalized(centroid.normalize())
        } else {
            S2Point::from_coords_raw(0.0, 0.0, 0.0)
        }
    }

    /// Get the geodesic curvature of the loop.
    pub fn get_curvature(&self) -> f64 {
        2.0 * PI - self.get_area()
    }

    /// Check if the loop contains the given point.
    pub fn contains(&self, point: &S2Point) -> bool {
        if self.is_empty() {
            return false;
        }
        if self.is_full() {
            return true;
        }

        // Point-in-polygon test using winding number
        let mut sign_sum = 0;
        let p = point.coords();
        
        for i in 0..self.vertices.len() {
            let next = (i + 1) % self.vertices.len();
            sign_sum += robust_sign(
                p,
                self.vertices[i].coords(),
                self.vertices[next].coords()
            );
        }
        
        sign_sum != 0
    }

    /// Check if this loop contains another loop.
    pub fn contains_loop(&self, other: &S2Loop) -> bool {
        // Handle special cases
        if self.is_full() || other.is_empty() {
            return true;
        }
        if self.is_empty() || other.is_full() {
            return false;
        }

        // Check if all vertices of other are contained
        for vertex in &other.vertices {
            if !self.contains(vertex) {
                return false;
            }
        }

        // TODO: Also need to check for edge intersections
        true
    }

    /// Check if this loop intersects another loop.
    pub fn intersects(&self, other: &S2Loop) -> bool {
        // Handle special cases
        if self.is_empty() || other.is_empty() {
            return false;
        }
        if self.is_full() || other.is_full() {
            return true;
        }

        // Check if any vertex of either loop is contained by the other
        for vertex in &self.vertices {
            if other.contains(vertex) {
                return true;
            }
        }
        for vertex in &other.vertices {
            if self.contains(vertex) {
                return true;
            }
        }

        // TODO: Check for edge intersections
        false
    }

    /// Check if two loops have the same boundary.
    pub fn boundary_equals(&self, other: &S2Loop) -> bool {
        if self.num_vertices() != other.num_vertices() {
            return false;
        }

        // Handle special cases
        if self.is_empty_or_full() && other.is_empty_or_full() {
            return self.is_empty() == other.is_empty();
        }

        // Find the first matching vertex
        for offset in 0..self.num_vertices() {
            let mut matches = true;
            for i in 0..self.num_vertices() {
                let self_vertex = self.vertex(i);
                let other_vertex = other.vertex((i + offset) % other.num_vertices());
                if !self_vertex.coords().abs_diff_eq(other_vertex.coords(), 1e-15) {
                    matches = false;
                    break;
                }
            }
            if matches {
                return true;
            }
        }

        false
    }

    /// Get the bounding rectangle.
    pub fn get_rect_bound(&self) -> S2LatLngRect {
        self.bound.clone()
    }

    /// Get a bounding cap.
    pub fn get_cap_bound(&self) -> S2Cap {
        if self.is_empty() {
            return S2Cap::empty();
        }
        if self.is_full() {
            return S2Cap::full();
        }

        // Simplified: Use the centroid and maximum distance
        let center = self.get_centroid();
        let mut max_dist_sq: f64 = 0.0;
        
        for vertex in &self.vertices {
            let dist_sq = (center.coords() - vertex.coords()).length_squared();
            max_dist_sq = max_dist_sq.max(dist_sq);
        }
        
        let radius = S1ChordAngle::from_length_squared(max_dist_sq);
        S2Cap::from_center_chord_angle(center, radius)
    }

    /// Check if the loop is valid.
    pub fn is_valid(&self) -> bool {
        self.validate().is_ok()
    }

    /// Get maximum error for GetCurvature computation.
    pub fn get_curvature_max_error(&self) -> f64 {
        // Simplified error bound - full implementation would be more complex
        1e-14 * self.num_vertices() as f64
    }

    /// Compute distance from point to loop interior.
    pub fn get_distance(&self, point: &S2Point) -> S1Angle {
        if self.is_empty() {
            return S1Angle::infinity();
        }
        if self.contains(point) {
            return S1Angle::zero();
        }
        self.get_distance_to_boundary(point)
    }

    /// Compute distance from point to loop boundary.
    pub fn get_distance_to_boundary(&self, point: &S2Point) -> S1Angle {
        if self.is_empty_or_full() {
            return S1Angle::infinity();
        }

        // Find minimum distance to any edge
        let mut min_distance = S1Angle::infinity();
        for i in 0..self.vertices.len() {
            let next = (i + 1) % self.vertices.len();
            let edge_distance = self.distance_to_edge(point, &self.vertices[i], &self.vertices[next]);
            if edge_distance < min_distance {
                min_distance = edge_distance;
            }
        }
        min_distance
    }

    /// Helper to compute distance from point to an edge.
    fn distance_to_edge(&self, point: &S2Point, a: &S2Point, b: &S2Point) -> S1Angle {
        // Simplified distance calculation - use great circle distance to edge
        // Full implementation would use more sophisticated edge distance calculation
        let to_a = S1Angle::from_radians(point.coords().dot(a.coords()).acos());
        let to_b = S1Angle::from_radians(point.coords().dot(b.coords()).acos());
        to_a.min(to_b)
    }

    /// Project point to loop interior if contained, otherwise to boundary.
    pub fn project(&self, point: &S2Point) -> S2Point {
        if self.contains(point) {
            *point
        } else {
            self.project_to_boundary(point)
        }
    }

    /// Project point to loop boundary.
    pub fn project_to_boundary(&self, point: &S2Point) -> S2Point {
        if self.is_empty_or_full() {
            return *point;
        }

        // Find closest point on boundary - simplified implementation
        let mut closest_point = self.vertices[0];
        let mut min_distance_sq = (point.coords() - closest_point.coords()).length_squared();

        for vertex in &self.vertices {
            let distance_sq = (point.coords() - vertex.coords()).length_squared();
            if distance_sq < min_distance_sq {
                min_distance_sq = distance_sq;
                closest_point = *vertex;
            }
        }

        closest_point
    }

    /// Create a regular polygon loop.
    pub fn make_regular_loop(center: S2Point, radius: S1Angle, num_vertices: i32) -> S2Result<Self> {
        if num_vertices < 3 {
            return Err(S2Error::invalid_loop("Regular loop needs at least 3 vertices"));
        }

        let mut vertices = Vec::with_capacity(num_vertices as usize);
        
        // Create an orthonormal frame with center as the z-axis
        let z = center.coords();
        let x = if z.z.abs() < 0.9 {
            DVec3::new(0.0, 0.0, 1.0).cross(z).normalize()
        } else {
            DVec3::new(1.0, 0.0, 0.0).cross(z).normalize()
        };
        let y = z.cross(x);

        let radius_rad = radius.radians();
        let angle_step = 2.0 * PI / num_vertices as f64;

        for i in 0..num_vertices {
            let angle = i as f64 * angle_step;
            let cos_angle = angle.cos();
            let sin_angle = angle.sin();
            
            // Point on the circle in the tangent plane
            let local_point = cos_angle * x + sin_angle * y;
            
            // Project onto the sphere
            let point_on_sphere = radius_rad.cos() * z + radius_rad.sin() * local_point;
            vertices.push(S2Point::from_normalized(point_on_sphere));
        }

        Self::new(vertices)
    }
}

impl PartialEq for S2Loop {
    fn eq(&self, other: &Self) -> bool {
        self.boundary_equals(other)
    }
}

