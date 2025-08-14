//! Conditional C++ Compatibility Tests
//!
//! These tests check at runtime if the C++ submodule is available and skip if not.

use s2geometry_rust::{S2Point, math::DVec3};
use s2geometry_cpp_compat_tests::*;

const TOLERANCE: f64 = 1e-15;

#[test]
fn test_s2point_normalize_equivalence() {
    if !is_cpp_submodule_available() {
        println!("Skipping test_s2point_normalize_equivalence - C++ submodule not available");
        return;
    }
    
    #[cfg(feature = "cpp-compat")]
    {
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
            let cpp_normalized = simple_s2point_normalize(unnormalized_point.into());
            
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
}

#[test]
fn test_s2point_angle_equivalence() {
    if !is_cpp_submodule_available() {
        println!("Skipping test_s2point_angle_equivalence - C++ submodule not available");
        return;
    }
    
    #[cfg(feature = "cpp-compat")]
    {
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
                let cpp_angle = simple_s2point_angle(point_a.into(), point_b.into());
                
                assert!(angles_equal_approx(rust_angle, cpp_angle, TOLERANCE),
                    "Angle mismatch between points {} and {}: Rust={}, C++={}, diff={}", 
                    i, j, rust_angle, cpp_angle, (rust_angle - cpp_angle).abs());
            }
        }
    }
}

#[test]
fn test_s2point_cross_prod_equivalence() {
    if !is_cpp_submodule_available() {
        println!("Skipping test_s2point_cross_prod_equivalence - C++ submodule not available");
        return;
    }
    
    #[cfg(feature = "cpp-compat")]
    {
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
            let cpp_cross = simple_s2point_cross_prod(point_a.into(), point_b.into());
            
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
}

#[test]
fn test_submodule_status_reporting() {
    // This test always runs to show the submodule status
    print_submodule_status();
}