//! General S2Point C++ Compatibility Tests
//!
//! Tests that verify functional equivalence between Rust and C++ S2Point implementations.

use s2geometry_rust::{S2Point, math::DVec3};
use s2geometry_cpp_compat_tests::*;

const TOLERANCE: f64 = 1e-15;

#[test]
fn test_s2point_normalize_equivalence() {
    // Test that point normalization produces identical results
    let test_vectors = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0),
        DVec3::new(-1.0, -1.0, -1.0),
        DVec3::new(2.0, 3.0, 4.0),
        DVec3::new(-5.0, 2.0, -1.0),
        DVec3::new(0.1, 0.2, 0.3),
        DVec3::new(100.0, 200.0, 300.0),
        DVec3::new(1e-10, 1e-10, 1e-10),
    ];
    
    for (vec_idx, vec) in test_vectors.iter().enumerate() {
        // Create unnormalized point
        let unnormalized_point = S2Point::from_coords(*vec);
        
        // Test Rust implementation
        let rust_normalized = unnormalized_point.normalize();
        
        // Test C++ implementation
        let cpp_normalized = cpp_s2point_normalize(unnormalized_point.into());
        
        let rust_coords = rust_normalized.coords();
        assert!(points_equal_approx(
            S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
            cpp_normalized,
            TOLERANCE
        ), "Normalize mismatch for vector {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
           vec_idx,
           rust_coords.x, rust_coords.y, rust_coords.z,
           cpp_normalized.x, cpp_normalized.y, cpp_normalized.z);
        
        // Verify both results are actually normalized
        let rust_norm = (rust_coords.x * rust_coords.x + rust_coords.y * rust_coords.y + rust_coords.z * rust_coords.z).sqrt();
        let cpp_norm = (cpp_normalized.x * cpp_normalized.x + cpp_normalized.y * cpp_normalized.y + cpp_normalized.z * cpp_normalized.z).sqrt();
        
        assert!(angles_equal_approx(rust_norm, 1.0, TOLERANCE),
            "Rust normalized point {} is not unit length: {}", vec_idx, rust_norm);
        assert!(angles_equal_approx(cpp_norm, 1.0, TOLERANCE),
            "C++ normalized point {} is not unit length: {}", vec_idx, cpp_norm);
    }
}

#[test]
fn test_s2point_angle_equivalence() {
    // Test that angle calculations between points match
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(-1.0, 0.0, 0.0),
        DVec3::new(0.0, -1.0, 0.0),
        DVec3::new(0.0, 0.0, -1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
        DVec3::new(-1.0, -1.0, -1.0).normalize(),
        DVec3::new(1.0, -1.0, 0.0).normalize(),
        DVec3::new(0.707, 0.707, 0.0),
    ];
    
    // Test all pairs of points
    for (i, vec_a) in test_points.iter().enumerate() {
        for (j, vec_b) in test_points.iter().enumerate() {
            if i >= j { continue; } // Avoid duplicate pairs and self-angle
            
            let point_a = S2Point::from_normalized(*vec_a);
            let point_b = S2Point::from_normalized(*vec_b);
            
            // Test Rust implementation
            let rust_angle = point_a.angle(&point_b).radians();
            
            // Test C++ implementation
            let cpp_angle = cpp_s2point_angle(point_a.into(), point_b.into());
            
            assert!(angles_equal_approx(rust_angle, cpp_angle, TOLERANCE),
                "Angle mismatch between points {} and {}: Rust={}, C++={}, diff={}", 
                i, j, rust_angle, cpp_angle, (rust_angle - cpp_angle).abs());
        }
    }
}

#[test]
fn test_s2point_cross_prod_equivalence() {
    // Test that cross product calculations match
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 0.0).normalize(),
        DVec3::new(0.0, 1.0, 1.0).normalize(),
        DVec3::new(1.0, 0.0, 1.0).normalize(),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
        DVec3::new(-1.0, 1.0, 0.0).normalize(),
    ];
    
    // Test selected pairs of points for cross product
    let test_pairs = vec![
        (0, 1), // X × Y = Z
        (1, 2), // Y × Z = X  
        (2, 0), // Z × X = Y
        (0, 2), // X × Z = -Y
        (1, 0), // Y × X = -Z
        (2, 1), // Z × Y = -X
        (3, 4), // Arbitrary normalized vectors
        (5, 6),
        (6, 7),
    ];
    
    for (pair_idx, (i, j)) in test_pairs.iter().enumerate() {
        let point_a = S2Point::from_normalized(test_points[*i]);
        let point_b = S2Point::from_normalized(test_points[*j]);
        
        // Test Rust implementation
        let rust_cross = point_a.cross_prod(&point_b);
        
        // Test C++ implementation
        let cpp_cross = cpp_s2point_cross_prod(point_a.into(), point_b.into());
        
        let rust_coords = rust_cross.coords();
        assert!(points_equal_approx(
            S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
            cpp_cross,
            TOLERANCE
        ), "Cross product mismatch for pair {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
           pair_idx,
           rust_coords.x, rust_coords.y, rust_coords.z,
           cpp_cross.x, cpp_cross.y, cpp_cross.z);
    }
}

#[test]
fn test_s2point_known_cross_products() {
    // Test specific known cross product results
    let x_axis = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    let y_axis = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
    let z_axis = S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0));
    
    // X × Y should equal Z
    let rust_x_cross_y = x_axis.cross_prod(&y_axis);
    let cpp_x_cross_y = cpp_s2point_cross_prod(x_axis.into(), y_axis.into());
    
    let rust_coords = rust_x_cross_y.coords();
    assert!(points_equal_approx(
        S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
        cpp_x_cross_y,
        TOLERANCE
    ), "X × Y mismatch: Rust=({}, {}, {}), C++=({}, {}, {})", 
       rust_coords.x, rust_coords.y, rust_coords.z,
       cpp_x_cross_y.x, cpp_x_cross_y.y, cpp_x_cross_y.z);
    
    // Result should be close to Z axis
    assert!(angles_equal_approx(rust_coords.z, 1.0, TOLERANCE),
        "X × Y should produce Z axis, got z component: {}", rust_coords.z);
    assert!(angles_equal_approx(cpp_x_cross_y.z, 1.0, TOLERANCE),
        "C++ X × Y should produce Z axis, got z component: {}", cpp_x_cross_y.z);
}

#[test]
fn test_s2point_special_cases() {
    // Test special cases and edge conditions
    
    // Test with very small vectors
    let tiny_vec = DVec3::new(1e-15, 1e-15, 1e-15);
    let tiny_point = S2Point::from_coords(tiny_vec);
    
    let rust_tiny_normalized = tiny_point.normalize();
    let cpp_tiny_normalized = cpp_s2point_normalize(tiny_point.into());
    
    let rust_coords = rust_tiny_normalized.coords();
    assert!(points_equal_approx(
        S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
        cpp_tiny_normalized,
        TOLERANCE
    ), "Tiny vector normalize mismatch: Rust=({}, {}, {}), C++=({}, {}, {})", 
       rust_coords.x, rust_coords.y, rust_coords.z,
       cpp_tiny_normalized.x, cpp_tiny_normalized.y, cpp_tiny_normalized.z);
    
    // Test with very large vectors
    let huge_vec = DVec3::new(1e15, 1e15, 1e15);
    let huge_point = S2Point::from_coords(huge_vec);
    
    let rust_huge_normalized = huge_point.normalize();
    let cpp_huge_normalized = cpp_s2point_normalize(huge_point.into());
    
    let rust_coords = rust_huge_normalized.coords();
    assert!(points_equal_approx(
        S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
        cpp_huge_normalized,
        TOLERANCE
    ), "Huge vector normalize mismatch: Rust=({}, {}, {}), C++=({}, {}, {})", 
       rust_coords.x, rust_coords.y, rust_coords.z,
       cpp_huge_normalized.x, cpp_huge_normalized.y, cpp_huge_normalized.z);
}

#[test]
fn test_s2point_angle_special_cases() {
    // Test special angle cases
    let point = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    
    // Angle with itself should be 0
    let rust_self_angle = point.angle(&point).radians();
    let cpp_self_angle = cpp_s2point_angle(point.into(), point.into());
    
    assert!(angles_equal_approx(rust_self_angle, 0.0, TOLERANCE),
        "Self angle should be 0, got Rust: {}", rust_self_angle);
    assert!(angles_equal_approx(cpp_self_angle, 0.0, TOLERANCE),
        "Self angle should be 0, got C++: {}", cpp_self_angle);
    
    // Angle with antipodal point should be π
    let antipodal = S2Point::from_normalized(DVec3::new(-1.0, 0.0, 0.0));
    let rust_antipodal_angle = point.angle(&antipodal).radians();
    let cpp_antipodal_angle = cpp_s2point_angle(point.into(), antipodal.into());
    
    assert!(angles_equal_approx(rust_antipodal_angle, std::f64::consts::PI, TOLERANCE),
        "Antipodal angle should be π, got Rust: {}", rust_antipodal_angle);
    assert!(angles_equal_approx(cpp_antipodal_angle, std::f64::consts::PI, TOLERANCE),
        "Antipodal angle should be π, got C++: {}", cpp_antipodal_angle);
}