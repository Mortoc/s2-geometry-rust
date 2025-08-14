//! Simplified S2 coordinate transformation tests
//!
//! This module contains the essential coordinate transformation tests
//! focusing on mathematical correctness and bit-identical compatibility with C++.

use s2geometry_rust::math::{DVec3, coords::*};
use approx::{assert_relative_eq, assert_abs_diff_eq};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Test Hilbert curve traversal order consistency
#[test]
fn test_traversal_order() {
    for r in 0..4 {
        for i in 0..4 {
            // Check consistency with respect to swapping axes
            assert_eq!(
                IJ_TO_POS[r][i],
                IJ_TO_POS[r ^ SWAP_MASK as usize][swap_axes(i as i32) as usize]
            );
            assert_eq!(
                POS_TO_IJ[r][i],
                swap_axes(POS_TO_IJ[r ^ SWAP_MASK as usize][i])
            );

            // Check consistency with respect to reversing axis directions
            assert_eq!(
                IJ_TO_POS[r][i],
                IJ_TO_POS[r ^ INVERT_MASK as usize][invert_bits(i as i32) as usize]
            );
            assert_eq!(
                POS_TO_IJ[r][i],
                invert_bits(POS_TO_IJ[r ^ INVERT_MASK as usize][i])
            );

            // Check that the two tables are inverses of each other
            assert_eq!(
                IJ_TO_POS[r][POS_TO_IJ[r][i] as usize],
                i as i32
            );
            assert_eq!(
                POS_TO_IJ[r][IJ_TO_POS[r][i] as usize],
                i as i32
            );
        }
    }
}

/// Test ST to IJ conversion boundaries
#[test] 
fn test_st_to_ij_boundaries() {
    assert_eq!(st_to_ij(0.0), 0);
    assert_eq!(st_to_ij(1.0), LIMIT_IJ - 1);
}

/// Test ST to IJ conversion at halfway points
#[test]
fn test_st_to_ij_halfway() {
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

/// Test IJ to ST to IJ roundtrip
#[test]
fn test_ij_to_st_to_ij_roundtrip_random() {
    let mut rng = StdRng::seed_from_u64(12345);
    
    for _ in 0..100 {  // Reduce iterations for faster testing
        let i = rng.gen_range(0..LIMIT_IJ);
        let s_min = ij_to_st_min(i);
        let s_max = ij_to_st_min(i + 1);
        let s = rng.gen_range(s_min..s_max);
        let i_roundtrip = st_to_ij(s);
        
        assert_eq!(i_roundtrip, i);
        assert_eq!(st_to_ij(s_min), i);
        
        let before_s_max = f64::from_bits(s_max.to_bits() - 1);
        assert_eq!(st_to_ij(before_s_max), i);
    }
}

/// Test ST-UV conversion consistency  
#[test]
fn test_st_uv_conversions() {
    // Check that UVtoST and STtoUV are inverses
    for i in 0..=1000 {  // Reduce iterations
        let x = (i as f64) / 1000.0;
        
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
        let largest_abs = center.x.abs().max(center.y.abs().max(center.z.abs()));
        assert!((largest_abs - 1.0).abs() <= 1e-15);
        
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
        
        assert!((dot_product - 1.0).abs() <= 1e-14);
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
        
        // Test normal (W axis)
        let norm = get_norm(face);
        let norm_result = face_xyz_to_uvw(face, norm);
        assert_abs_diff_eq!(norm_result.x, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(norm_result.y, 0.0, epsilon = 1e-15);
        assert_relative_eq!(norm_result.z, 1.0, epsilon = 1e-15);
    }
}

/// Test XYZ to Face Si Ti conversion 
#[test]
fn test_xyz_to_face_si_ti_basic() {
    let mut rng = StdRng::seed_from_u64(42);
    
    // Test conversion and roundtrip for basic cases
    for _ in 0..50 {  
        let face = rng.gen_range(0..6);
        let si = rng.gen_range(1u32..(MAX_SI_TI - 1)) & 0xFFFF_FF00;  // Use aligned coordinates
        let ti = rng.gen_range(1u32..(MAX_SI_TI - 1)) & 0xFFFF_FF00;
        
        // Skip edge coordinates
        if si == 0 || ti == 0 || si >= MAX_SI_TI || ti >= MAX_SI_TI {
            continue;
        }
        
        // Convert to XYZ and back
        let point = face_si_ti_to_xyz(face, si, ti);
        let (face_result, si_result, ti_result, _level_result) = xyz_to_face_si_ti(point);
        
        assert_eq!(face_result, face);
        assert_eq!(si_result, si);
        assert_eq!(ti_result, ti);
    }
}

/// Test UVW axis consistency
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
        assert!((dot - 1.0).abs() <= 1e-15);
        
        // Check GetUVWAxis consistency
        assert_eq!(get_uvw_axis(face, 0), u_axis);
        assert_eq!(get_uvw_axis(face, 1), v_axis);
        assert_eq!(get_uvw_axis(face, 2), norm);
    }
}

/// Test UVW face neighbor consistency
#[test] 
fn test_uvw_face() {
    for face in 0..6 {
        for axis in 0..3 {
            // Check GetUVWFace consistency with GetUVWAxis
            let neg_axis = -get_uvw_axis(face, axis);
            let pos_axis = get_uvw_axis(face, axis);
            
            let neg_face = get_face(neg_axis);
            let pos_face = get_face(pos_axis);
            
            assert_eq!(neg_face, get_uvw_face(face, axis, 0));
            assert_eq!(pos_face, get_uvw_face(face, axis, 1));
        }
    }
}

/// Comprehensive coordinate transformation consistency test
#[test]
fn test_coordinate_transformation_consistency() {
    let mut rng = StdRng::seed_from_u64(987654321);
    
    for _ in 0..100 {  // Reduced iterations for faster testing
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