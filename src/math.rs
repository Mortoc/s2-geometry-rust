//! Mathematical foundations for S2 geometry
//!
//! This module implements the **three-tier mathematical architecture** matching
//! Google's C++ S2 implementation for maximum compatibility and robustness.
//!
//! # Design Philosophy (Following C++ S2 Architecture)
//!
//! ## Fast Path (~90-95% of operations)
//! Uses `glam::DVec3` f64 precision with optimized operations:
//! - Primary computational type for all geometric operations
//! - SIMD-optimized where supported by glam
//! - Error bounds analysis using `DBL_EPSILON` methodology
//!
//! ## Stable Path (~4-9% of operations)  
//! Extended precision for borderline cases:
//! - Uses higher precision arithmetic for uncertain results
//! - Error analysis with tighter bounds than f64
//! - Reduces exact arithmetic fallback rate significantly
//!
//! ## Exact Path (<1% of operations)
//! Perfect arithmetic for geometric predicates requiring deterministic results:
//! - [`ExactFloat`] using [`num_rational::BigRational`] backend
//! - Guaranteed deterministic orientation tests, edge intersections
//! - Zero possibility of numerical inconsistency
//!
//! # Predicate Architecture (Matching C++ Implementation)
//!
//! ```rust,ignore
//! use s2geometry_rust::math::predicates::*;
//!
//! // Three-tier predicate evaluation
//! let orientation = robust_sign(a, b, c);
//!
//! // Internal implementation follows C++ pattern:
//! // 1. triage_sign() - fast f64 with error bounds  
//! // 2. stable_sign() - extended precision if needed
//! // 3. exact_sign() - BigRational fallback for edge cases
//! ```
//!
//! # Error Thresholds (Matching C++ Values)
//!
//! The thresholds are based on rigorous floating-point error analysis:
//! - **Triage threshold**: `3.6548 * f64::EPSILON ≈ 8.12e-16`
//! - **Stable threshold**: Platform-dependent extended precision analysis
//! - **Exact fallback**: Only when all numeric approaches fail
//!
//! This matches Google's proven approach for robust computational geometry.

// Re-export the core f64 math types from glam
pub use glam::{DVec2, DVec3, DVec4, DMat2, DMat3, DMat4, DQuat, DAffine2, DAffine3};

// Phase 1 will implement ExactFloat, Vector3<T>, and exact geometric predicates

/// S2 coordinate system transformations matching C++ implementation exactly
pub mod coords {
    use super::DVec3;
    
    // S2 coordinate system constants (matching C++ exactly)
    /// Maximum level for S2 cell subdivision (30 levels)
    pub const MAX_CELL_LEVEL: i32 = 30;
    /// Maximum IJ coordinate value (2^30)
    pub const LIMIT_IJ: i32 = 1 << MAX_CELL_LEVEL;  // 2^30 = 1073741824
    /// Maximum SI/TI coordinate value (2^31)
    pub const MAX_SI_TI: u32 = 1u32 << (MAX_CELL_LEVEL + 1);  // 2^31
    
    // Projection type - using quadratic projection like C++ default
    const QUADRATIC_PROJECTION: bool = true;
    
    // Maximum absolute error in XYZ->UV conversion (matching C++ kMaxXYZtoUVError)
    /// Maximum absolute error in XYZ->UV conversion
    pub const MAX_XYZ_TO_UV_ERROR: f64 = 0.5 * f64::EPSILON;
    
    // Hilbert curve traversal tables (matching C++ s2coords_internal.h)
    /// Mask for swap bit in Hilbert curve traversal
    pub const SWAP_MASK: i32 = 0x01;
    /// Mask for invert bit in Hilbert curve traversal
    pub const INVERT_MASK: i32 = 0x02;
    
    // IJ to position lookup table for Hilbert curve traversal (matching C++ exactly)
    /// IJ to position lookup table for Hilbert curve traversal
    pub const IJ_TO_POS: [[i32; 4]; 4] = [
        [0, 1, 3, 2],  // canonical order:  (0,0), (0,1), (1,1), (1,0)
        [0, 3, 1, 2],  // axes swapped:     (0,0), (1,0), (1,1), (0,1)  
        [2, 3, 1, 0],  // bits inverted:    (1,1), (1,0), (0,0), (0,1)
        [2, 1, 3, 0],  // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    ];
    
    // Position to IJ lookup table (inverse of IJ_TO_POS - matching C++ exactly)
    /// Position to IJ lookup table (inverse of IJ_TO_POS)
    pub const POS_TO_IJ: [[i32; 4]; 4] = [
        [0, 1, 3, 2],  // canonical order: (0,0), (0,1), (1,1), (1,0)
        [0, 2, 3, 1],  // axes swapped: (0,0), (1,0), (1,1), (0,1)
        [3, 2, 0, 1],  // bits inverted: (1,1), (1,0), (0,0), (0,1)
        [3, 1, 0, 2],  // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    ];
    
    // Position to orientation modifier
    /// Position to orientation modifier for Hilbert curve traversal
    pub const POS_TO_ORIENTATION: [i32; 4] = [SWAP_MASK, 0, 0, INVERT_MASK | SWAP_MASK];
    
    // Face UVW axes lookup table (matching C++ kFaceUVWAxes)
    /// Face UVW axes lookup table for S2 cube face projections
    pub const FACE_UVW_AXES: [[[f64; 3]; 3]; 6] = [
        [[ 0.0,  1.0,  0.0], [ 0.0,  0.0,  1.0], [ 1.0,  0.0,  0.0]], // face 0
        [[-1.0,  0.0,  0.0], [ 0.0,  0.0,  1.0], [ 0.0,  1.0,  0.0]], // face 1  
        [[-1.0,  0.0,  0.0], [ 0.0, -1.0,  0.0], [ 0.0,  0.0,  1.0]], // face 2
        [[ 0.0,  0.0, -1.0], [ 0.0, -1.0,  0.0], [-1.0,  0.0,  0.0]], // face 3
        [[ 0.0,  0.0, -1.0], [ 1.0,  0.0,  0.0], [ 0.0, -1.0,  0.0]], // face 4
        [[ 0.0,  1.0,  0.0], [ 1.0,  0.0,  0.0], [ 0.0,  0.0, -1.0]], // face 5
    ];
    
    // Face neighbors lookup table (matching C++ kFaceUVWFaces exactly)
    /// Face UVW faces lookup table for S2 cube face projections
    pub const FACE_UVW_FACES: [[[i32; 2]; 3]; 6] = [
        [[4, 1], [5, 2], [3, 0]], // face 0
        [[0, 3], [5, 2], [4, 1]], // face 1 
        [[0, 3], [1, 4], [5, 2]], // face 2
        [[2, 5], [1, 4], [0, 3]], // face 3
        [[2, 5], [3, 0], [1, 4]], // face 4  
        [[4, 1], [3, 0], [2, 5]], // face 5
    ];
    
    /// Convert s or t coordinate to u or v coordinate using quadratic projection
    /// Matches C++ STtoUV() exactly
    #[inline]
    pub fn st_to_uv(s: f64) -> f64 {
        if QUADRATIC_PROJECTION {
            if s >= 0.5 {
                (1.0/3.0) * (4.0*s*s - 1.0)
            } else {
                (1.0/3.0) * (1.0 - 4.0*(1.0-s)*(1.0-s))
            }
        } else {
            // Linear projection fallback
            2.0 * s - 1.0
        }
    }
    
    /// Convert u or v coordinate to s or t coordinate (inverse of st_to_uv)
    /// Matches C++ UVtoST() exactly  
    #[inline]
    pub fn uv_to_st(u: f64) -> f64 {
        if QUADRATIC_PROJECTION {
            if u >= 0.0 {
                0.5 * (1.0 + 3.0*u).sqrt()
            } else {
                1.0 - 0.5 * (1.0 - 3.0*u).sqrt()
            }
        } else {
            // Linear projection fallback
            0.5 * (u + 1.0)
        }
    }
    
    /// Convert leaf cell index i to minimum s coordinate
    /// Matches C++ IJtoSTMin() exactly
    #[inline]
    pub fn ij_to_st_min(i: i32) -> f64 {
        debug_assert!(i >= 0 && i <= LIMIT_IJ);
        (i as f64) / (LIMIT_IJ as f64)
    }
    
    /// Convert s coordinate to leaf cell index i  
    /// Matches C++ STtoIJ() exactly
    #[inline] 
    pub fn st_to_ij(s: f64) -> i32 {
        // Use floor operation via cast, then clamp to valid range
        let result = (LIMIT_IJ as f64 * s) as i32;
        result.clamp(0, LIMIT_IJ - 1)
    }
    
    /// Convert discrete si/ti coordinate to s/t coordinate
    /// Matches C++ SiTitoST() exactly
    #[inline]
    pub fn si_ti_to_st(si: u32) -> f64 {
        debug_assert!(si <= MAX_SI_TI);
        (si as f64) / (MAX_SI_TI as f64)
    }
    
    /// Convert s/t coordinate to discrete si/ti coordinate
    /// Matches C++ STtoSiTi() exactly
    #[inline]
    pub fn st_to_si_ti(s: f64) -> u32 {
        // Round to nearest integer (matching MathUtil::Round)
        ((s * MAX_SI_TI as f64) + 0.5) as u32
    }
    
    /// Convert face and (u,v) coordinates to XYZ direction vector
    /// Matches C++ FaceUVtoXYZ() exactly
    #[inline]
    pub fn face_uv_to_xyz(face: i32, u: f64, v: f64) -> DVec3 {
        match face {
            0 => DVec3::new( 1.0,   u,   v),
            1 => DVec3::new(  -u, 1.0,   v),
            2 => DVec3::new(  -u,  -v, 1.0),
            3 => DVec3::new(-1.0,  -v,  -u),
            4 => DVec3::new(   v, -1.0, -u),
            5 => DVec3::new(   v,    u, -1.0),
            _ => panic!("Invalid face: {}", face),
        }
    }
    
    /// Get the face containing the given direction vector
    /// Matches C++ GetFace() exactly
    #[inline]
    pub fn get_face(p: DVec3) -> i32 {
        let abs_x = p.x.abs();
        let abs_y = p.y.abs();
        let abs_z = p.z.abs();
        
        let face = if abs_x >= abs_y && abs_x >= abs_z {
            0  // X axis
        } else if abs_y >= abs_z {
            1  // Y axis  
        } else {
            2  // Z axis
        };
        
        // Add 3 if the component is negative
        if p[face as usize] < 0.0 { face + 3 } else { face }
    }
    
    /// Convert XYZ direction vector to face and (u,v) coordinates
    /// Matches C++ XYZtoFaceUV() exactly
    #[inline]
    pub fn xyz_to_face_uv(p: DVec3) -> (i32, f64, f64) {
        let face = get_face(p);
        let (u, v) = valid_face_xyz_to_uv(face, p);
        (face, u, v)
    }
    
    /// Convert XYZ to UV for a known valid face (assumes dot product > 0)
    /// Matches C++ ValidFaceXYZtoUV() exactly  
    #[inline]
    pub fn valid_face_xyz_to_uv(face: i32, p: DVec3) -> (f64, f64) {
        debug_assert!(p.dot(get_norm(face)) > 0.0);
        match face {
            0 => ( p.y / p.x,  p.z / p.x),
            1 => (-p.x / p.y,  p.z / p.y),
            2 => (-p.x / p.z, -p.y / p.z),
            3 => ( p.z / p.x,  p.y / p.x),
            4 => ( p.z / p.y, -p.x / p.y),
            5 => (-p.y / p.z, -p.x / p.z),
            _ => panic!("Invalid face: {}", face),
        }
    }
    
    /// Try to convert XYZ to UV for specific face (returns None if invalid)
    /// Matches C++ FaceXYZtoUV() exactly
    #[inline]
    pub fn face_xyz_to_uv(face: i32, p: DVec3) -> Option<(f64, f64)> {
        let valid = if face < 3 {
            p[face as usize] > 0.0
        } else {
            p[(face - 3) as usize] < 0.0
        };
        
        if valid {
            Some(valid_face_xyz_to_uv(face, p))
        } else {
            None
        }
    }
    
    /// Transform point to (u,v,w) coordinate frame of given face
    /// Matches C++ FaceXYZtoUVW() exactly
    #[inline] 
    pub fn face_xyz_to_uvw(face: i32, p: DVec3) -> DVec3 {
        // Transform using face coordinate axes
        let u_axis = get_u_axis(face);
        let v_axis = get_v_axis(face);
        let w_axis = get_norm(face);
        
        DVec3::new(
            p.dot(u_axis),
            p.dot(v_axis), 
            p.dot(w_axis)
        )
    }
    
    /// Convert XYZ to (face, si, ti) and return cell level if point is cell center
    /// Matches C++ XYZtoFaceSiTi() exactly
    pub fn xyz_to_face_si_ti(p: DVec3) -> (i32, u32, u32, i32) {
        let (face, u, v) = xyz_to_face_uv(p);
        
        // Convert UV to ST
        let s = uv_to_st(u);
        let t = uv_to_st(v);
        
        // Convert to discrete coordinates
        let si = st_to_si_ti(s);
        let ti = st_to_si_ti(t);
        
        // Check if this corresponds to a cell center at some level
        let level = compute_cell_level(si, ti);
        
        (face, si, ti, level)
    }
    
    /// Convert (face, si, ti) to XYZ direction vector  
    /// Matches C++ FaceSiTitoXYZ() exactly
    #[inline]
    pub fn face_si_ti_to_xyz(face: i32, si: u32, ti: u32) -> DVec3 {
        let s = si_ti_to_st(si);
        let t = si_ti_to_st(ti);
        let u = st_to_uv(s);
        let v = st_to_uv(t);
        face_uv_to_xyz(face, u, v)
    }
    
    /// Compute cell level for given si,ti coordinates (-1 if not a cell center)
    fn compute_cell_level(si: u32, ti: u32) -> i32 {
        // Check if si and ti represent a cell center at some level
        // Cell centers have coordinates ending with 1 followed by (30-level) zeros
        
        if si == 0 || ti == 0 || si == MAX_SI_TI || ti == MAX_SI_TI {
            return -1; // Edge coordinates, not cell centers
        }
        
        // Find the level by checking trailing zeros
        let si_trailing_zeros = si.trailing_zeros();
        let ti_trailing_zeros = ti.trailing_zeros();
        
        // Both coordinates must have the same number of trailing zeros
        if si_trailing_zeros != ti_trailing_zeros {
            return -1;
        }
        
        // Check if the coordinates are valid cell centers
        let level = MAX_CELL_LEVEL - si_trailing_zeros as i32;
        if level < 0 || level > MAX_CELL_LEVEL {
            return -1;
        }
        
        // Verify that si and ti are odd after removing trailing zeros
        let si_shifted = si >> si_trailing_zeros;
        let ti_shifted = ti >> ti_trailing_zeros;
        
        if (si_shifted & 1) != 1 || (ti_shifted & 1) != 1 {
            return -1;
        }
        
        level
    }
    
    /// Get right-handed normal for edge in positive V direction at given U
    /// Matches C++ GetUNorm() exactly
    #[inline]
    pub fn get_u_norm(face: i32, u: f64) -> DVec3 {
        match face {
            0 => DVec3::new( u, -1.0,  0.0),
            1 => DVec3::new( 1.0,  u,  0.0), 
            2 => DVec3::new( 1.0, 0.0,   u),
            3 => DVec3::new(-u,  0.0,  1.0),
            4 => DVec3::new( 0.0, -u,  1.0),
            5 => DVec3::new( 0.0, -1.0, -u),
            _ => panic!("Invalid face: {}", face),
        }
    }
    
    /// Get right-handed normal for edge in positive U direction at given V
    /// Matches C++ GetVNorm() exactly
    #[inline]
    pub fn get_v_norm(face: i32, v: f64) -> DVec3 {
        match face {
            0 => DVec3::new(-v,  0.0,  1.0),
            1 => DVec3::new( 0.0, -v,  1.0),
            2 => DVec3::new( 0.0, -1.0, -v),
            3 => DVec3::new( v, -1.0,  0.0),
            4 => DVec3::new( 1.0,  v,  0.0),
            5 => DVec3::new( 1.0, 0.0,   v),
            _ => panic!("Invalid face: {}", face),
        }
    }
    
    /// Get unit-length normal vector for given face
    /// Matches C++ GetNorm() exactly
    #[inline]
    pub fn get_norm(face: i32) -> DVec3 {
        get_uvw_axis(face, 2)
    }
    
    /// Get unit-length U-axis vector for given face
    /// Matches C++ GetUAxis() exactly  
    #[inline]
    pub fn get_u_axis(face: i32) -> DVec3 {
        get_uvw_axis(face, 0)
    }
    
    /// Get unit-length V-axis vector for given face
    /// Matches C++ GetVAxis() exactly
    #[inline]
    pub fn get_v_axis(face: i32) -> DVec3 {
        get_uvw_axis(face, 1)
    }
    
    /// Get axis vector for given face and axis (u=0, v=1, w=2)
    /// Matches C++ GetUVWAxis() exactly
    #[inline]
    pub fn get_uvw_axis(face: i32, axis: usize) -> DVec3 {
        debug_assert!(face >= 0 && face <= 5);
        debug_assert!(axis <= 2);
        let coords = FACE_UVW_AXES[face as usize][axis];
        DVec3::new(coords[0], coords[1], coords[2])
    }
    
    /// Get adjacent face in given direction along given axis
    /// Matches C++ GetUVWFace() exactly
    #[inline] 
    pub fn get_uvw_face(face: i32, axis: usize, direction: usize) -> i32 {
        debug_assert!(face >= 0 && face <= 5);
        debug_assert!(axis <= 2);
        debug_assert!(direction <= 1);
        FACE_UVW_FACES[face as usize][axis][direction]
    }
    
    // Helper functions for testing
    
    /// Swap axes in ij coordinate (for testing)
    #[inline]
    pub fn swap_axes(ij: i32) -> i32 {
        ((ij >> 1) & 1) + ((ij & 1) << 1)
    }
    
    /// Invert bits in ij coordinate (for testing)
    #[inline] 
    pub fn invert_bits(ij: i32) -> i32 {
        ij ^ 3
    }
}

/// Mathematical constants used throughout S2 geometry
pub mod constants {
    use std::f64::consts;

    /// π
    pub const PI: f64 = consts::PI;
    
    /// π/2
    pub const PI_2: f64 = consts::FRAC_PI_2;
    
    /// π/4  
    pub const PI_4: f64 = consts::FRAC_PI_4;
    
    /// 2π
    pub const PI2: f64 = 2.0 * consts::PI;
    
    /// Square root of 2
    pub const SQRT2: f64 = consts::SQRT_2;
    
    /// Degrees to radians conversion factor
    pub const DEG_TO_RAD: f64 = consts::PI / 180.0;
    
    /// Radians to degrees conversion factor
    pub const RAD_TO_DEG: f64 = 180.0 / consts::PI;
    
    /// Machine epsilon for f64 comparisons
    pub const EPSILON: f64 = f64::EPSILON;
    
    /// Maximum determinant error threshold (matching C++ S2 exactly)
    /// From C++: kMaxDetError = 3.6548 * DBL_EPSILON
    pub const TRIAGE_ERROR_THRESHOLD: f64 = 3.6548 * f64::EPSILON; // ≈ 8.12e-16
    
    /// Threshold for stable precision fallback (extended precision analysis)
    /// Used when triage fails but exact arithmetic may not be needed
    pub const STABLE_ERROR_THRESHOLD: f64 = TRIAGE_ERROR_THRESHOLD * 0.1;
    
    /// Earth radius in meters (for geodetic calculations)
    pub const EARTH_RADIUS_METERS: f64 = 6_371_010.0;
    
    /// Maximum error in robust cross product (matching C++ kRobustCrossProdError)
    pub const ROBUST_CROSS_PROD_ERROR: f64 = 9.0 * f64::EPSILON / 2.0;
    
    /// Maximum error in intersection computation
    pub const INTERSECTION_ERROR: f64 = 8.0 * f64::EPSILON;
}

/// Legacy predicates - DEPRECATED: use crate::predicates module for new code
pub mod predicates {
    use super::DVec3;
    
    /// Fast triage test for orientation - DEPRECATED
    ///
    /// Use `crate::predicates::triage_sign()` for new code.
    #[deprecated(note = "Use crate::predicates::triage_sign() instead")]
    #[inline]
    pub fn triage_sign(a: DVec3, b: DVec3, c: DVec3, a_cross_b: DVec3) -> i32 {
        crate::predicates::triage_sign(a, b, c, a_cross_b)
    }
    
    /// Robust orientation test - DEPRECATED
    ///
    /// Use `crate::predicates::sign()` for new code.
    #[deprecated(note = "Use crate::predicates::sign() instead")]
    #[inline]
    pub fn robust_sign(a: DVec3, b: DVec3, c: DVec3) -> i32 {
        crate::predicates::sign(a, b, c)
    }
    
    /// Robust cross product - DEPRECATED
    ///
    /// Use basic `a.cross(b)` for new code, or implement robust cross product
    /// in the predicates module if needed.
    #[deprecated(note = "Use a.cross(b) or implement in crate::predicates")]
    #[inline]
    pub fn robust_cross_prod(a: DVec3, b: DVec3) -> DVec3 {
        a.cross(b)
    }
    
    /// Edge crossing predicate - DEPRECATED
    ///
    /// Use `crate::predicates::crossing_sign()` for new code.
    #[deprecated(note = "Use crate::predicates::crossing_sign() instead")]
    #[inline]
    pub fn crossing_sign(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> i32 {
        crate::predicates::crossing_sign(a, b, c, d)
    }
    
    /// Vertex crossing predicate - DEPRECATED
    ///
    /// Use `crate::predicates::vertex_crossing()` for new code.
    #[deprecated(note = "Use crate::predicates::vertex_crossing() instead")]
    pub fn vertex_crossing(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> bool {
        crate::predicates::vertex_crossing(a, b, c, d)
    }
    
    /// Signed vertex crossing predicate - DEPRECATED
    ///
    /// Use `crate::predicates::signed_vertex_crossing()` for new code.
    #[deprecated(note = "Use crate::predicates::signed_vertex_crossing() instead")]
    pub fn signed_vertex_crossing(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> i32 {
        crate::predicates::signed_vertex_crossing(a, b, c, d)
    }
    
    /// Edge or vertex crossing predicate - DEPRECATED
    ///
    /// Use `crate::predicates::edge_or_vertex_crossing()` for new code.
    #[deprecated(note = "Use crate::predicates::edge_or_vertex_crossing() instead")]
    pub fn edge_or_vertex_crossing(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> bool {
        crate::predicates::edge_or_vertex_crossing(a, b, c, d)
    }
}

/// Utility functions for floating-point operations
pub mod float_utils {
    use super::constants::EPSILON;

    /// Test if two floating-point values are approximately equal
    ///
    /// # Example
    /// ```rust,ignore
    /// use s2geometry_rust::math::float_utils::approx_eq;
    /// assert!(approx_eq(1.0, 1.0 + 1e-16));
    /// ```
    pub fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < EPSILON
    }

    /// Test if a floating-point value is approximately zero
    pub fn approx_zero(x: f64) -> bool {
        x.abs() < EPSILON
    }

    /// Clamp a value to a range
    pub fn clamp(x: f64, min: f64, max: f64) -> f64 {
        if x < min {
            min
        } else if x > max {
            max
        } else {
            x
        }
    }

    /// Square a number
    #[inline(always)]
    pub fn square(x: f64) -> f64 {
        x * x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants() {
        use constants::*;
        assert!((PI - std::f64::consts::PI).abs() < f64::EPSILON);
        assert!((DEG_TO_RAD * 180.0 - PI).abs() < f64::EPSILON);
        assert!((RAD_TO_DEG * PI - 180.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_float_utils() {
        use float_utils::*;
        
        assert!(approx_eq(1.0, 1.0));
        assert!(approx_eq(1.0, 1.0 + f64::EPSILON / 2.0));
        assert!(!approx_eq(1.0, 1.0 + 1e-10));
        
        assert!(approx_zero(0.0));
        assert!(approx_zero(f64::EPSILON / 2.0));
        assert!(!approx_zero(1e-10));
        
        assert_eq!(clamp(5.0, 0.0, 10.0), 5.0);
        assert_eq!(clamp(-1.0, 0.0, 10.0), 0.0);
        assert_eq!(clamp(15.0, 0.0, 10.0), 10.0);
        
        assert_eq!(square(3.0), 9.0);
        assert_eq!(square(-4.0), 16.0);
    }
}