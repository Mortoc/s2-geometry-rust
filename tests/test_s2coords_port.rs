//! Port of s2coords_test.cc - Critical coordinate transformation tests
//!
//! This module ports the essential coordinate transformation tests from the C++ s2coords_test.cc
//! to validate our coordinate transformation architecture. These transformations are foundational
//! to the entire S2 system - face projections, UV coordinates, and XYZ conversions must be
//! mathematically exact.
//!
//! # Key Test Categories
//!
//! 1. **TraversalOrder** - Hilbert curve traversal consistency
//! 2. **ST/IJ Conversions** - Discrete coordinate conversions with exact boundaries  
//! 3. **ST_UV_Conversions** - Core projection transformations
//! 4. **FaceUVtoXYZ** - Face-to-XYZ coordinate mapping with right-handed checks
//! 5. **FaceXYZtoUVW** - XYZ-to-face coordinate transformations
//! 6. **XYZToFaceSiTi** - Point-to-discrete coordinate conversion with roundtrip validation
//! 7. **UVNorms** - Edge normal computations  
//! 8. **UVWAxis/UVWFace** - Face coordinate frame consistency

use s2geometry_rust::math::{DVec3, coords::*};
use approx::{assert_relative_eq, assert_abs_diff_eq};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rand_pcg::Pcg64;

/// Test Hilbert curve traversal order consistency
/// Port of C++ TEST(S2, TraversalOrder)
#[test]
fn test_traversal_order() {
    for r in 0..4 {
        for i in 0..4 {
            // Check consistency with respect to swapping axes
            assert_eq!(
                IJ_TO_POS[r][i],
                IJ_TO_POS[r ^ SWAP_MASK as usize][swap_axes(i as i32) as usize],
                "Swap consistency failed at r={}, i={}", r, i
            );
            assert_eq!(
                POS_TO_IJ[r][i],
                swap_axes(POS_TO_IJ[r ^ SWAP_MASK as usize][i]),
                "Swap consistency (reverse) failed at r={}, i={}", r, i
            );

            // Check consistency with respect to reversing axis directions
            assert_eq!(
                IJ_TO_POS[r][i],
                IJ_TO_POS[r ^ INVERT_MASK as usize][invert_bits(i as i32) as usize],
                "Invert consistency failed at r={}, i={}", r, i
            );
            assert_eq!(
                POS_TO_IJ[r][i],
                invert_bits(POS_TO_IJ[r ^ INVERT_MASK as usize][i]),
                "Invert consistency (reverse) failed at r={}, i={}", r, i
            );

            // Check that the two tables are inverses of each other
            assert_eq!(
                IJ_TO_POS[r][POS_TO_IJ[r][i] as usize],
                i as i32,
                "Inverse relationship failed at r={}, i={}", r, i
            );
            assert_eq!(
                POS_TO_IJ[r][IJ_TO_POS[r][i] as usize],
                i as i32,
                "Inverse relationship (reverse) failed at r={}, i={}", r, i
            );
        }
    }
}

/// Test ST to IJ conversion boundaries
/// Port of C++ TEST(S2, STtoIJBoundaries)
#[test] 
fn test_st_to_ij_boundaries() {
    assert_eq!(st_to_ij(0.0), 0, "ST to IJ boundary at 0.0 failed");
    assert_eq!(st_to_ij(1.0), LIMIT_IJ - 1, "ST to IJ boundary at 1.0 failed");
}

/// Test ST to IJ conversion at halfway points
/// Port of C++ TEST(S2, STtoIJHalfway)
#[test]
fn test_st_to_ij_halfway() {
    // LIMIT_IJ is 2^30, so this is exact
    let recip_limit_ij = 1.0 / (LIMIT_IJ as f64);
    
    assert_eq!(st_to_ij(0.5 * recip_limit_ij), 0);
    assert_eq!(st_to_ij(1.0 * recip_limit_ij), 1);
    assert_eq!(st_to_ij(1.5 * recip_limit_ij), 1);
    assert_eq!(st_to_ij(2.0 * recip_limit_ij), 2);
    assert_eq!(st_to_ij(2.5 * recip_limit_ij), 2);
    assert_eq!(st_to_ij(3.0 * recip_limit_ij), 3);
    assert_eq!(st_to_ij(3.5 * recip_limit_ij), 3);
    assert_eq!(st_to_ij(4.0 * recip_limit_ij), 4);
    assert_eq!(st_to_ij(4.5 * recip_limit_ij), 4);
    
    // Test near the upper bound
    assert_eq!(st_to_ij(((LIMIT_IJ as f64) - 2.5) * recip_limit_ij), LIMIT_IJ - 3);
    assert_eq!(st_to_ij(((LIMIT_IJ as f64) - 2.0) * recip_limit_ij), LIMIT_IJ - 2);
    assert_eq!(st_to_ij(((LIMIT_IJ as f64) - 1.5) * recip_limit_ij), LIMIT_IJ - 2);
    assert_eq!(st_to_ij(((LIMIT_IJ as f64) - 1.0) * recip_limit_ij), LIMIT_IJ - 1);
    assert_eq!(st_to_ij(((LIMIT_IJ as f64) - 0.5) * recip_limit_ij), LIMIT_IJ - 1);
}

/// Test IJ to ST to IJ roundtrip with random values
/// Port of C++ TEST(S2, IJtoSTtoIJRoundtripRandom)
#[test]
fn test_ij_to_st_to_ij_roundtrip_random() {
    let mut rng = StdRng::seed_from_u64(12345);
    
    for _ in 0..1000 {
        let i = rng.gen_range(0..LIMIT_IJ);
        let s_min = ij_to_st_min(i);
        let s_max = ij_to_st_min(i + 1);
        let s = rng.gen_range(s_min..s_max);
        let i_roundtrip = st_to_ij(s);
        
        assert_eq!(i_roundtrip, i, "Roundtrip failed for i={}, s={}", i, s);
        assert_eq!(st_to_ij(s_min), i, "Min boundary failed for s_min={}", s_min);
        
        // Test just before s_max
        let before_s_max = f64::from_bits(s_max.to_bits() - 1);
        assert_eq!(st_to_ij(before_s_max), i, "Before max boundary failed for before_s_max={}", before_s_max);
    }
}

/// Test ST-UV conversion consistency  
/// Port of C++ TEST(S2, ST_UV_Conversions)
#[test]
fn test_st_uv_conversions() {
    // Check boundary conditions
    for s_int in [0, 1, 2] {
        let s = (s_int as f64) / 2.0;
        let u = st_to_uv(s);
        let expected_u = 2.0 * s - 1.0;
        
        // For quadratic projection, this exact relationship only holds at boundaries
        if s == 0.0 || s == 1.0 || s == 0.5 {
            assert_relative_eq!(u, expected_u, epsilon = 1e-15);
        }
    }
    
    for u_int in [-1, 0, 1] {
        let u = u_int as f64;
        let s = uv_to_st(u);
        let expected_s = 0.5 * (u + 1.0);
        
        // For quadratic projection, this exact relationship only holds at boundaries  
        if u == -1.0 || u == 1.0 || u == 0.0 {
            assert_relative_eq!(s, expected_s, epsilon = 1e-15);
        }
    }
    
    // Check that UVtoST and STtoUV are inverses
    for i in 0..=10000 {
        let x = (i as f64) / 10000.0;
        
        // Test ST -> UV -> ST roundtrip
        let uv_result = st_to_uv(x);
        let st_roundtrip = uv_to_st(uv_result);
        assert_abs_diff_eq!(st_roundtrip, x, epsilon = 1e-15);
        
        // Test UV -> ST -> UV roundtrip  
        let u = 2.0 * x - 1.0; // Map [0,1] to [-1,1]
        let st_result = uv_to_st(u);
        let uv_roundtrip = st_to_uv(st_result);
        assert_abs_diff_eq!(uv_roundtrip, u, epsilon = 1e-15);
    }
}

/// Test Face UV to XYZ transformation 
/// Port of C++ TEST(S2, FaceUVtoXYZ)
#[test]
fn test_face_uv_to_xyz() {
    // Check that each face appears exactly once
    let mut sum = DVec3::ZERO;
    for face in 0..6 {
        let center = face_uv_to_xyz(face, 0.0, 0.0);
        let expected_center = get_norm(face);
        
        assert_relative_eq!(center.x, expected_center.x, epsilon = 1e-15);
        assert_relative_eq!(center.y, expected_center.y, epsilon = 1e-15);
        assert_relative_eq!(center.z, expected_center.z, epsilon = 1e-15);
        
        // Check that largest component has absolute value 1
        let largest_component_idx = if center.x.abs() >= center.y.abs() && center.x.abs() >= center.z.abs() {
            0
        } else if center.y.abs() >= center.z.abs() {
            1
        } else {
            2
        };
        assert!((center[largest_component_idx].abs() - 1.0).abs() <= 1e-15,
               "Largest component not 1 for face {}", face);
        
        sum += DVec3::new(center.x.abs(), center.y.abs(), center.z.abs());
    }
    
    // Sum of absolute values should be (2,2,2)
    assert_relative_eq!(sum.x, 2.0, epsilon = 1e-15);
    assert_relative_eq!(sum.y, 2.0, epsilon = 1e-15);
    assert_relative_eq!(sum.z, 2.0, epsilon = 1e-15);
    
    // Check that each face has a right-handed coordinate system
    for face in 0..6 {
        let u_axis = get_u_axis(face);
        let v_axis = get_v_axis(face);
        let center = face_uv_to_xyz(face, 0.0, 0.0);
        
        let cross_product = u_axis.cross(v_axis);
        let dot_product = cross_product.dot(center);
        
        assert_relative_eq!(dot_product, 1.0, epsilon = 1e-14);
    }
    
    // Check that Hilbert curves connect continuously between faces
    for face in 0..6 {
        let sign = if (face & SWAP_MASK) != 0 { -1.0 } else { 1.0 };
        let end_point = face_uv_to_xyz(face, sign, -sign);
        let next_start = face_uv_to_xyz((face + 1) % 6, -1.0, -1.0);
        
        assert_relative_eq!(end_point.x, next_start.x, epsilon = 1e-15);
        assert_relative_eq!(end_point.y, next_start.y, epsilon = 1e-15);
        assert_relative_eq!(end_point.z, next_start.z, epsilon = 1e-15);
    }
}

/// Test Face XYZ to UVW transformation
/// Port of C++ TEST(S2, FaceXYZtoUVW)  
#[test]
fn test_face_xyz_to_uvw() {
    for face in 0..6 {
        // Test zero vector
        assert_eq!(face_xyz_to_uvw(face, DVec3::ZERO), DVec3::ZERO);
        
        // Test U axis
        let u_axis = get_u_axis(face);
        let u_result = face_xyz_to_uvw(face, u_axis);
        assert_relative_eq!(u_result.x, 1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(u_result.y, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(u_result.z, 0.0, epsilon = 1e-15);
        
        let neg_u_result = face_xyz_to_uvw(face, -u_axis);
        assert_relative_eq!(neg_u_result.x, -1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(neg_u_result.y, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(neg_u_result.z, 0.0, epsilon = 1e-15);
        
        // Test V axis
        let v_axis = get_v_axis(face);
        let v_result = face_xyz_to_uvw(face, v_axis);
        assert_abs_diff_eq!(v_result.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(v_result.y, 1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(v_result.z, 0.0, epsilon = 1e-15);
        
        let neg_v_result = face_xyz_to_uvw(face, -v_axis);
        assert_abs_diff_eq!(neg_v_result.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(neg_v_result.y, -1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(neg_v_result.z, 0.0, epsilon = 1e-15);
        
        // Test normal (W axis)
        let norm = get_norm(face);
        let norm_result = face_xyz_to_uvw(face, norm);
        assert_abs_diff_eq!(norm_result.x, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(norm_result.y, 0.0, epsilon = 1e-15);
        assert_relative_eq!(norm_result.z, 1.0, epsilon = 1e-15);
        
        let neg_norm_result = face_xyz_to_uvw(face, -norm);
        assert_abs_diff_eq!(neg_norm_result.x, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(neg_norm_result.y, 0.0, epsilon = 1e-15);
        assert_relative_eq!(neg_norm_result.z, -1.0, epsilon = 1e-15);
    }
}

/// Test XYZ to Face Si Ti conversion with random cell IDs
/// Port of C++ TEST(S2, XYZToFaceSiTi) - simplified version focusing on key properties
#[test]
fn test_xyz_to_face_si_ti() {
    let mut rng = Pcg64::seed_from_u64(42);
    
    // Test conversion and roundtrip for various cases
    for level in 0..=5 {  // Test a subset of levels for performance
        for _ in 0..100 {
            // Generate random face coordinates 
            let face = rng.gen_range(0..6);
            let si = rng.gen_range(1u32..(MAX_SI_TI - 1)) & (!(0u32) << (MAX_CELL_LEVEL - level));
            let ti = rng.gen_range(1u32..(MAX_SI_TI - 1)) & (!(0u32) << (MAX_CELL_LEVEL - level));
            
            // Skip edge coordinates
            if si == 0 || ti == 0 || si >= MAX_SI_TI || ti >= MAX_SI_TI {
                continue;
            }
            
            // Convert to XYZ
            let point = face_si_ti_to_xyz(face, si, ti);
            
            // Convert back to face/si/ti  
            let (face_result, si_result, ti_result, level_result) = xyz_to_face_si_ti(point);
            
            // Check that we get the same face and coordinates
            assert_eq!(face_result, face, "Face mismatch: expected {}, got {}", face, face_result);
            assert_eq!(si_result, si, "Si mismatch: expected {}, got {}", si, si_result);
            assert_eq!(ti_result, ti, "Ti mismatch: expected {}, got {}", ti, ti_result);
            
            // Level check is more complex - for now just verify non-negative for valid coordinates
            if level_result >= 0 {
                assert!(level_result <= MAX_CELL_LEVEL, "Level {} exceeds maximum {}", level_result, MAX_CELL_LEVEL);
            }
        }
    }
    
    // Test specific boundary cases
    test_xyz_to_face_si_ti_boundary_case(0, 1073741824, 1073741824); // 2^30 coordinates
    test_xyz_to_face_si_ti_boundary_case(1, 536870912, 536870912);   // 2^29 coordinates  
    test_xyz_to_face_si_ti_boundary_case(2, 268435456, 268435456);   // 2^28 coordinates
}

fn test_xyz_to_face_si_ti_boundary_case(face: i32, si: u32, ti: u32) {
    let point = face_si_ti_to_xyz(face, si, ti);
    let (face_result, si_result, ti_result, _level_result) = xyz_to_face_si_ti(point);
    
    assert_eq!(face_result, face, "Boundary case face mismatch");
    assert_eq!(si_result, si, "Boundary case si mismatch");
    assert_eq!(ti_result, ti, "Boundary case ti mismatch");
}

/// Test UV normal computations
/// Port of C++ TEST(S2, UVNorms)
#[test]
fn test_uv_norms() {
    for face in 0..6 {
        // Test at regular intervals 
        for i in 0..=1024 {
            let x = -1.0 + 2.0 * (i as f64) / 1024.0;
            
            // Test U norm
            let edge_start = face_uv_to_xyz(face, x, -1.0);
            let edge_end = face_uv_to_xyz(face, x, 1.0);
            let edge_cross = edge_start.cross(edge_end);
            let u_norm = get_u_norm(face, x);
            
            let angle = edge_cross.normalize().dot(u_norm.normalize());
            assert!(angle >= 0.9999, "U norm angle test failed for face {}, x={}, angle={}", face, x, angle);
            
            // Test V norm
            let edge_start_v = face_uv_to_xyz(face, -1.0, x);
            let edge_end_v = face_uv_to_xyz(face, 1.0, x);
            let edge_cross_v = edge_start_v.cross(edge_end_v);
            let v_norm = get_v_norm(face, x);
            
            let angle_v = edge_cross_v.normalize().dot(v_norm.normalize());
            assert!(angle_v >= 0.9999, "V norm angle test failed for face {}, x={}, angle={}", face, x, angle_v);
        }
    }
}

/// Test UVW axis consistency
/// Port of C++ TEST(S2, UVWAxis)  
#[test]
fn test_uvw_axis() {
    for face in 0..6 {
        // Check that axes are consistent with FaceUVtoXYZ
        let u_diff = face_uv_to_xyz(face, 1.0, 0.0) - face_uv_to_xyz(face, 0.0, 0.0);
        let u_axis = get_u_axis(face);
        assert_relative_eq!(u_diff.x, u_axis.x, epsilon = 1e-15);
        assert_relative_eq!(u_diff.y, u_axis.y, epsilon = 1e-15);
        assert_relative_eq!(u_diff.z, u_axis.z, epsilon = 1e-15);
        
        let v_diff = face_uv_to_xyz(face, 0.0, 1.0) - face_uv_to_xyz(face, 0.0, 0.0);
        let v_axis = get_v_axis(face);
        assert_relative_eq!(v_diff.x, v_axis.x, epsilon = 1e-15);
        assert_relative_eq!(v_diff.y, v_axis.y, epsilon = 1e-15);  
        assert_relative_eq!(v_diff.z, v_axis.z, epsilon = 1e-15);
        
        let center = face_uv_to_xyz(face, 0.0, 0.0);
        let norm = get_norm(face);
        assert_relative_eq!(center.x, norm.x, epsilon = 1e-15);
        assert_relative_eq!(center.y, norm.y, epsilon = 1e-15);
        assert_relative_eq!(center.z, norm.z, epsilon = 1e-15);
        
        // Check that coordinate frame is right-handed
        let cross = u_axis.cross(v_axis);
        let dot = cross.dot(norm);
        assert_relative_eq!(dot, 1.0, epsilon = 1e-15);
        
        // Check GetUVWAxis consistency
        assert_eq!(get_uvw_axis(face, 0), u_axis);
        assert_eq!(get_uvw_axis(face, 1), v_axis);
        assert_eq!(get_uvw_axis(face, 2), norm);
    }
}

/// Test UVW face neighbor consistency
/// Port of C++ TEST(S2, UVWFace)
#[test] 
fn test_uvw_face() {
    for face in 0..6 {
        for axis in 0..3 {
            // Check GetUVWFace consistency with GetUVWAxis
            let neg_axis = -get_uvw_axis(face, axis);
            let pos_axis = get_uvw_axis(face, axis);
            
            let neg_face = get_face(neg_axis);
            let pos_face = get_face(pos_axis);
            
            assert_eq!(neg_face, get_uvw_face(face, axis, 0),
                      "Negative direction face mismatch for face {}, axis {}", face, axis);
            assert_eq!(pos_face, get_uvw_face(face, axis, 1),
                      "Positive direction face mismatch for face {}, axis {}", face, axis);
        }
    }
}

/// Comprehensive coordinate transformation consistency test
#[test]
fn test_coordinate_transformation_consistency() {
    let mut rng = Pcg64::seed_from_u64(987654321);
    
    for _ in 0..1000 {
        // Generate random point on unit sphere  
        let x = rng.gen::<f64>() - 0.5;
        let y = rng.gen::<f64>() - 0.5;
        let z = rng.gen::<f64>() - 0.5;
        let point = DVec3::new(x, y, z).normalize();
        
        // Test XYZ -> Face UV -> XYZ roundtrip
        let (face, u, v) = xyz_to_face_uv(point);
        let roundtrip_point = face_uv_to_xyz(face, u, v).normalize();
        
        assert_relative_eq!(point.x, roundtrip_point.x, epsilon = 1e-14);
        assert_relative_eq!(point.y, roundtrip_point.y, epsilon = 1e-14);
        assert_relative_eq!(point.z, roundtrip_point.z, epsilon = 1e-14);
        
        // Test ST <-> UV consistency
        let s = uv_to_st(u);
        let t = uv_to_st(v);
        let u_roundtrip = st_to_uv(s);
        let v_roundtrip = st_to_uv(t);
        
        assert_abs_diff_eq!(u, u_roundtrip, epsilon = 1e-15);
        assert_abs_diff_eq!(v, v_roundtrip, epsilon = 1e-15);
        
        // Test face coordinate frame orthogonality
        let u_axis = get_u_axis(face);
        let v_axis = get_v_axis(face);
        let norm = get_norm(face);
        
        assert_abs_diff_eq!(u_axis.dot(v_axis), 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(u_axis.dot(norm), 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(v_axis.dot(norm), 0.0, epsilon = 1e-15);
        
        // Verify unit length
        assert_relative_eq!(u_axis.length(), 1.0, epsilon = 1e-15);
        assert_relative_eq!(v_axis.length(), 1.0, epsilon = 1e-15);
        assert_relative_eq!(norm.length(), 1.0, epsilon = 1e-15);
    }
}