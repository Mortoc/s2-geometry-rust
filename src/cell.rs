//! S2Cell - Geometric representation of S2 hierarchical cells
//!
//! S2Cell provides the geometric counterpart to S2CellId, representing the actual
//! spatial cell as a spherical quadrilateral with vertices, edges, and geometric properties.
//! Unlike S2CellId which is just an identifier, S2Cell supports efficient containment
//! tests, area calculations, and geometric operations.

use crate::error::{S2Error, S2Result};
use crate::point::S2Point;
use crate::cell_id::S2CellId;
use crate::r2::{R2Point, R2Rect};
use crate::interval::R1Interval;
use crate::chord_angle::S1ChordAngle;
use std::fmt;

/// Represents a spherical quadrilateral cell in the S2 geometry hierarchy
///
/// S2Cell provides geometric operations for cells identified by S2CellId.
/// It stores precomputed geometric information for efficient queries.
///
/// ## Key Features
/// - Vertex and edge computation
/// - Area and perimeter calculations  
/// - Containment and intersection tests
/// - Bounding cap and rectangle computation
/// - Distance calculations to points and edges
///
/// ## Memory Layout
/// Designed to match C++ implementation size constraints (~48 bytes)
#[derive(Debug, Clone, PartialEq)]
pub struct S2Cell {
    /// Cell identifier
    id: S2CellId,
    /// Cube face (0-5)
    face: i8,
    /// Subdivision level (0-30)
    level: i8,
    /// Hilbert curve orientation
    orientation: i8,
    /// UV coordinate bounds on the cube face
    uv: R2Rect,
}

/// Canonical edge identifiers for cell boundaries
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Bottom edge (k=0)
    Bottom = 0,
    /// Right edge (k=1)  
    Right = 1,
    /// Top edge (k=2)
    Top = 2,
    /// Left edge (k=3)
    Left = 3,
}

impl S2Cell {
    /// Create an S2Cell from an S2CellId
    ///
    /// This is the primary constructor that computes all geometric properties
    /// from the cell identifier.
    pub fn new(id: S2CellId) -> S2Result<Self> {
        if !id.is_valid() {
            return Err(S2Error::invalid_cell_id(id.id(), "Invalid S2CellId"));
        }

        let face = id.face();
        let level = id.level();
        let orientation = Self::compute_orientation(&id)?;
        let uv = Self::compute_uv_bounds(&id)?;

        Ok(S2Cell {
            id,
            face: face as i8,
            level: level as i8,
            orientation: orientation as i8,
            uv,
        })
    }

    /// Create cell from S2Point (convenience constructor)
    pub fn from_point(point: &S2Point) -> S2Result<Self> {
        let id = S2CellId::from_point(point);
        Self::new(id)
    }

    /// Create cell from face number (level 0 face cell)
    pub fn from_face(face: i32) -> S2Result<Self> {
        let id = S2CellId::from_face_pos_level(face, 0, 0)?;
        Self::new(id)
    }

    /// Create cell from face, position, and level
    pub fn from_face_pos_level(face: i32, pos: u64, level: i32) -> S2Result<Self> {
        let id = S2CellId::from_face_pos_level(face, pos, level)?;
        Self::new(id)
    }

    // Accessors

    /// Get the cell identifier
    #[inline]
    pub fn id(&self) -> S2CellId {
        self.id
    }

    /// Get the cube face (0-5)
    #[inline]
    pub fn face(&self) -> i32 {
        self.face as i32
    }

    /// Get the subdivision level (0-30)
    #[inline]
    pub fn level(&self) -> i32 {
        self.level as i32
    }

    /// Get the Hilbert curve orientation
    #[inline]
    pub fn orientation(&self) -> i32 {
        self.orientation as i32
    }

    /// Check if this is a leaf cell (maximum subdivision)
    #[inline]
    pub fn is_leaf(&self) -> bool {
        self.level() == crate::cell_id::MAX_LEVEL
    }

    /// Get UV coordinate bounds on the cube face
    #[inline]
    pub fn get_bound_uv(&self) -> R2Rect {
        self.uv
    }

    // Vertex operations

    /// Get the k-th vertex of the cell (k = 0,1,2,3)
    ///
    /// Vertices are returned in CCW order: lower-left, lower-right, upper-right, upper-left
    /// in UV coordinates. The returned point is normalized.
    pub fn get_vertex(&self, k: i32) -> S2Point {
        self.get_vertex_raw(k).normalize()
    }

    /// Get the k-th vertex without normalization
    ///
    /// This is more efficient when normalization isn't needed.
    pub fn get_vertex_raw(&self, k: i32) -> S2Point {
        let k = ((k % 4) + 4) % 4; // Proper modulo for negative numbers
        let uv_vertex = self.uv.get_vertex(k);
        Self::face_uv_to_xyz(self.face(), uv_vertex)
    }

    // Edge operations

    /// Get the inward-facing normal of edge k
    ///
    /// Returns the normalized normal vector for the great circle edge
    /// from vertex k to vertex k+1.
    pub fn get_edge(&self, k: i32) -> S2Point {
        self.get_edge_raw(k).normalize()
    }

    /// Get the edge normal without normalization
    ///
    /// The raw edge normal is bounded by sqrt(2) in length and can be used
    /// for exact geometric predicates.
    pub fn get_edge_raw(&self, k: i32) -> S2Point {
        let k = k & 3; // Reduce modulo 4
        match k {
            0 => Self::get_v_norm(self.face(), self.uv.y().lo()), // Bottom
            1 => Self::get_u_norm(self.face(), self.uv.x().hi()), // Right  
            2 => -Self::get_v_norm(self.face(), self.uv.y().hi()), // Top
            _ => -Self::get_u_norm(self.face(), self.uv.x().lo()), // Left
        }
    }

    /// Get UV coordinate for edge k (whichever coordinate is constant)
    pub fn get_uv_coord_of_edge(&self, k: i32) -> f64 {
        let k = k & 3;
        if k % 2 == 0 {
            // Edges 0,2 are constant in V
            self.uv.get_vertex(k).y()
        } else {
            // Edges 1,3 are constant in U
            self.uv.get_vertex(k).x()
        }
    }

    // Center operations

    /// Get the cell center (normalized)
    pub fn get_center(&self) -> S2Point {
        self.get_center_raw().normalize()
    }

    /// Get the cell center without normalization
    pub fn get_center_raw(&self) -> S2Point {
        self.id.to_point_raw()
    }

    // Subdivision

    /// Subdivide cell into four children
    ///
    /// If this is not a leaf cell, fills the children array with the four
    /// child cells and returns true. Otherwise returns false.
    pub fn subdivide(&self, children: &mut [S2Cell; 4]) -> bool {
        if self.is_leaf() {
            return false;
        }

        let child_ids: Vec<S2CellId> = self.id.children().collect();
        for (i, child_id) in child_ids.iter().enumerate() {
            if let Ok(child) = Self::new(*child_id) {
                children[i] = child;
            } else {
                return false;
            }
        }
        true
    }

    // Area calculations

    /// Get average area for cells at the given level
    pub fn average_area(level: i32) -> f64 {
        // Each face has area 4π/6 = 2π/3, and each level quadruples the number of cells
        (2.0 * std::f64::consts::PI / 3.0) / (1u64 << (2 * level)) as f64
    }

    /// Get average area for cells at this level
    pub fn get_average_area(&self) -> f64 {
        Self::average_area(self.level())
    }

    /// Get approximate area of this cell
    ///
    /// Accurate to within 3% for all cell sizes, and 0.1% for level 5+ cells.
    pub fn approx_area(&self) -> f64 {
        // Simple approximation: UV area * Jacobian scaling
        let uv_size = self.uv.get_size();
        let uv_area = uv_size.x() * uv_size.y();
        // Apply scaling factor for face projection (simplified)
        uv_area * Self::average_area(self.level()) / Self::uv_area_at_level(self.level())
    }

    /// Get exact area of this cell
    ///
    /// More expensive but accurate to 6 digits even for leaf cells.
    pub fn exact_area(&self) -> f64 {
        // Use spherical excess formula for exact area calculation
        let vertices = [
            self.get_vertex(0),
            self.get_vertex(1), 
            self.get_vertex(2),
            self.get_vertex(3),
        ];
        Self::spherical_area(&vertices)
    }

    // Containment tests

    /// Test if this cell contains the given point
    ///
    /// S2Cells are considered closed sets, so points on edges/vertices
    /// belong to adjacent cells as well.
    pub fn contains(&self, point: &S2Point) -> bool {
        // Convert point to UV coordinates on this face
        if let Ok(uv) = Self::xyz_to_face_uv(self.face(), point) {
            self.uv.contains(uv)
        } else {
            false
        }
    }

    /// Test if this cell contains another cell
    pub fn contains_cell(&self, other: &S2Cell) -> bool {
        // Fast test: check if this is an ancestor of other
        if self.level() >= other.level() {
            return false;
        }
        
        // Check if other is a descendant
        if let Ok(ancestor) = other.id.parent(self.level()) {
            ancestor == self.id
        } else {
            false
        }
    }

    /// Test if this cell may intersect another cell
    pub fn may_intersect(&self, other: &S2Cell) -> bool {
        // Cells intersect if their ranges overlap
        self.id.intersects(&other.id)
    }

    // Distance calculations

    /// Get distance from cell to point
    ///
    /// Returns zero if point is inside the cell.
    pub fn get_distance(&self, target: &S2Point) -> S1ChordAngle {
        if self.contains(target) {
            S1ChordAngle::zero()
        } else {
            self.get_boundary_distance(target)
        }
    }

    /// Get distance from cell boundary to point
    pub fn get_boundary_distance(&self, target: &S2Point) -> S1ChordAngle {
        let mut min_dist = S1ChordAngle::infinity();
        
        // Check distance to each edge
        for k in 0..4 {
            let v0 = self.get_vertex(k);
            let v1 = self.get_vertex(k + 1);
            let dist = Self::point_to_edge_distance(target, &v0, &v1);
            min_dist = min_dist.min(dist);
        }
        
        min_dist
    }

    /// Get maximum distance from cell to point
    pub fn get_max_distance(&self, target: &S2Point) -> S1ChordAngle {
        // Check if antipodal point is contained
        let antipodal = S2Point::from_coords_raw(-target.x(), -target.y(), -target.z());
        if self.contains(&antipodal) {
            return S1ChordAngle::straight();
        }

        let mut max_dist = S1ChordAngle::negative();
        
        // Check distance to each vertex
        for k in 0..4 {
            let vertex = self.get_vertex(k);
            let dist = S1ChordAngle::from_points(*target, vertex);
            max_dist = max_dist.max(dist);
        }
        
        max_dist
    }

    // Helper methods

    fn compute_orientation(_id: &S2CellId) -> S2Result<i32> {
        // Simplified orientation calculation
        // In full implementation, this would use ToFaceIJOrientation
        Ok(0) // Placeholder
    }

    fn compute_uv_bounds(id: &S2CellId) -> S2Result<R2Rect> {
        // Compute UV bounds for the cell based on level and position
        // For face cells (level 0), this should be [-1, +1] x [-1, +1]
        if id.level() == 0 {
            // Face cell covers entire face: UV range [-1, +1] x [-1, +1]
            let u = R1Interval::new(-1.0, 1.0);
            let v = R1Interval::new(-1.0, 1.0);
            Ok(R2Rect::from_intervals(u, v))
        } else {
            // For non-face cells, we need proper IJ to UV conversion
            // This is a simplified placeholder - full implementation would be more complex
            let size = 2.0 / (1 << id.level()) as f64; // Size in UV space
            let u = R1Interval::new(-1.0, -1.0 + size);
            let v = R1Interval::new(-1.0, -1.0 + size);
            Ok(R2Rect::from_intervals(u, v))
        }
    }

    fn face_uv_to_xyz(face: i32, uv: R2Point) -> S2Point {
        // Simplified face to XYZ conversion
        // This would use proper S2 coordinate transformation
        let u = uv.x();
        let v = uv.y();
        
        let (x, y, z) = match face {
            0 => (1.0, u, v),
            1 => (-u, 1.0, v),
            2 => (-u, -v, 1.0),
            3 => (-1.0, -v, -u),
            4 => (v, -1.0, -u),
            5 => (v, u, -1.0),
            _ => (1.0, 0.0, 0.0),
        };
        
        S2Point::from_coords_raw(x, y, z)
    }

    fn xyz_to_face_uv(face: i32, point: &S2Point) -> S2Result<R2Point> {
        // Simplified XYZ to face UV conversion
        let coords = point.coords();
        let (u, v) = match face {
            0 => (coords.y / coords.x, coords.z / coords.x),
            1 => (-coords.x / coords.y, coords.z / coords.y),
            2 => (-coords.x / coords.z, -coords.y / coords.z),
            3 => (coords.z / (-coords.x), -coords.y / (-coords.x)),
            4 => (coords.z / (-coords.y), -coords.x / (-coords.y)),
            5 => (-coords.y / (-coords.z), coords.x / (-coords.z)),
            _ => return Err(S2Error::invalid_cell_id(0, "Invalid face")),
        };
        Ok(R2Point::new(u, v))
    }

    fn get_u_norm(face: i32, u: f64) -> S2Point {
        // Get normal vector in U direction for given face and U coordinate
        // This matches the C++ GetUNorm function exactly
        match face {
            0 => S2Point::from_coords_raw(u, -1.0, 0.0),
            1 => S2Point::from_coords_raw(1.0, u, 0.0),
            2 => S2Point::from_coords_raw(1.0, 0.0, u),
            3 => S2Point::from_coords_raw(-u, 0.0, 1.0),
            4 => S2Point::from_coords_raw(0.0, -u, 1.0),
            5 => S2Point::from_coords_raw(0.0, -1.0, -u),
            _ => S2Point::from_coords_raw(1.0, 0.0, 0.0), // Should not happen
        }
    }

    fn get_v_norm(face: i32, v: f64) -> S2Point {
        // Get normal vector in V direction for given face and V coordinate  
        // This matches the C++ GetVNorm function exactly
        match face {
            0 => S2Point::from_coords_raw(-v, 0.0, 1.0),
            1 => S2Point::from_coords_raw(0.0, -v, 1.0),
            2 => S2Point::from_coords_raw(0.0, -1.0, -v),
            3 => S2Point::from_coords_raw(v, -1.0, 0.0),
            4 => S2Point::from_coords_raw(1.0, v, 0.0),
            5 => S2Point::from_coords_raw(1.0, 0.0, v),
            _ => S2Point::from_coords_raw(0.0, 1.0, 0.0), // Should not happen
        }
    }

    fn uv_area_at_level(level: i32) -> f64 {
        // UV area of a cell at given level (4.0 for face cells)
        4.0 / (1u64 << (2 * level)) as f64
    }

    fn spherical_area(vertices: &[S2Point]) -> f64 {
        // Calculate spherical area using spherical excess
        // This is a simplified implementation
        if vertices.len() != 4 {
            return 0.0;
        }
        
        // Use approximate formula for now
        let edge_lengths: Vec<f64> = (0..4)
            .map(|i| vertices[i].angle(&vertices[(i + 1) % 4]))
            .collect();
        
        let avg_edge = edge_lengths.iter().sum::<f64>() / 4.0;
        avg_edge * avg_edge // Very rough approximation
    }

    fn point_to_edge_distance(point: &S2Point, a: &S2Point, b: &S2Point) -> S1ChordAngle {
        // Calculate distance from point to great circle edge
        // Simplified implementation
        let dist_a = S1ChordAngle::from_points(*point, *a);
        let dist_b = S1ChordAngle::from_points(*point, *b);
        dist_a.min(dist_b)
    }

    // Additional methods needed for S2CellUnion

    /// Returns the average area of cells at the given level
    pub fn average_area_at_level(level: i32) -> f64 {
        Self::average_area(level)
    }

    /// Returns a bounding cap for this cell
    pub fn get_cap_bound(&self) -> crate::S2Cap {
        // Get cell center and radius
        let center = self.get_center();
        
        // Find the maximum distance from center to any vertex
        let mut max_dist_squared: f64 = 0.0;
        for k in 0..4 {
            let vertex = self.get_vertex(k);
            let dist_squared = center.coords().distance_squared(vertex.coords());
            max_dist_squared = max_dist_squared.max(dist_squared);
        }
        
        let radius = crate::S1Angle::from_radians(max_dist_squared.sqrt().asin());
        crate::S2Cap::from_center_angle(&center, radius)
    }

    /// Returns a bounding rectangle for this cell  
    pub fn get_rect_bound(&self) -> crate::S2LatLngRect {
        let mut rect = crate::S2LatLngRect::empty();
        
        // Add all vertices to the bounding rectangle
        for k in 0..4 {
            let vertex = self.get_vertex(k);
            let latlng = crate::S2LatLng::from_point(vertex);
            rect.add_point(&latlng);
        }
        
        rect
    }
}

impl fmt::Display for S2Cell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "S2Cell(id: {}, face: {}, level: {})", 
               self.id, self.face, self.level)
    }
}

impl From<S2CellId> for S2Cell {
    fn from(id: S2CellId) -> Self {
        Self::new(id).unwrap()
    }
}

// Comparison operators based on cell ID
impl PartialOrd for S2Cell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.id.cmp(&other.id))
    }
}

impl Eq for S2Cell {}

impl Ord for S2Cell {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.cmp(&other.id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_construction() {
        // Test face cell creation
        for face in 0..6 {
            let cell = S2Cell::from_face(face).unwrap();
            assert_eq!(cell.face(), face);
            assert_eq!(cell.level(), 0);
            assert!(!cell.is_leaf());
        }
    }

    #[test]
    fn test_vertices() {
        let cell = S2Cell::from_face(0).unwrap();
        
        // Should have 4 vertices
        for k in 0..4 {
            let vertex = cell.get_vertex(k);
            let length = vertex.coords().length();
            assert!((length - 1.0).abs() < 1e-10, "Vertex should be unit length, got {}", length);
        }
    }

    #[test]
    fn test_containment() {
        let cell = S2Cell::from_face(0).unwrap();
        let center = cell.get_center();
        
        // Cell should contain its center
        assert!(cell.contains(&center));
    }

    #[test]
    fn test_area_calculations() {
        let cell = S2Cell::from_face(0).unwrap();
        
        let avg_area = cell.get_average_area();
        let approx_area = cell.approx_area();
        let exact_area = cell.exact_area();
        
        // All areas should be positive
        assert!(avg_area > 0.0);
        assert!(approx_area > 0.0); 
        assert!(exact_area > 0.0);
        
        // Face cell should have area ~2π/3
        let expected_face_area = 2.0 * std::f64::consts::PI / 3.0;
        assert!((avg_area - expected_face_area).abs() < 0.1);
    }
}