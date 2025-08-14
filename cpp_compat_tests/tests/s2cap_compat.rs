//! S2Cap C++ Compatibility Tests
//!
//! Tests that verify functional equivalence between Rust and C++ S2Cap implementations.
//! This ensures that the same cap with the same arguments produces identical results.

use s2geometry_rust::{S2Point, S2Cap, S1Angle, math::DVec3};
use s2geometry_cpp_compat_tests::*;

const TOLERANCE: f64 = 1e-15;

#[test]
fn test_s2cap_cell_covering_equivalence() {
    // Test that S2Cap cell covering produces identical results between Rust and C++
    let center = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    let radius = S1Angle::from_radians(0.1);
    
    // Get covering from Rust implementation
    let rust_cap = S2Cap::from_center_angle(center, radius);
    let mut rust_covering = Vec::new();
    let mut coverer = s2geometry_rust::S2RegionCoverer::new();
    coverer.set_max_level(10);
    coverer.set_max_cells(100);
    coverer.get_covering(&rust_cap, &mut rust_covering);
    
    // Get covering from C++ implementation
    let cpp_covering = cpp_s2cap_from_point_angle(
        center.into(),
        radius.radians()
    );
    
    // Compare results
    assert_eq!(rust_covering.len(), cpp_covering.len(),
        "Cell covering count mismatch: Rust={}, C++={}", 
        rust_covering.len(), cpp_covering.len());
    
    // Sort both vectors for comparison (order might differ)
    let mut rust_ids: Vec<u64> = rust_covering.iter().map(|c| c.id()).collect();
    let mut cpp_ids: Vec<u64> = cpp_covering.iter().map(|c| c.id).collect();
    rust_ids.sort();
    cpp_ids.sort();
    
    for (i, (rust_id, cpp_id)) in rust_ids.iter().zip(cpp_ids.iter()).enumerate() {
        assert_eq!(rust_id, cpp_id, 
            "Cell ID mismatch at index {}: Rust={}, C++={}", i, rust_id, cpp_id);
    }
}

#[test]
fn test_s2cap_contains_point_equivalence() {
    // Test multiple cap centers and points to ensure containment matches
    let test_cases = vec![
        // (center, radius, test_point)
        (DVec3::new(1.0, 0.0, 0.0), 0.1, DVec3::new(0.95, 0.1, 0.0)),
        (DVec3::new(0.0, 1.0, 0.0), 0.5, DVec3::new(0.2, 0.9, 0.1)),
        (DVec3::new(0.0, 0.0, 1.0), 0.01, DVec3::new(0.0, 0.001, 0.999)),
        (DVec3::new(-1.0, 0.0, 0.0), 1.0, DVec3::new(-0.5, 0.5, 0.5)),
    ];
    
    for (i, (center_vec, radius_rad, point_vec)) in test_cases.iter().enumerate() {
        let center = S2Point::from_normalized(*center_vec);
        let point = S2Point::from_normalized(*point_vec);
        let radius = S1Angle::from_radians(*radius_rad);
        
        // Test Rust implementation
        let rust_cap = S2Cap::from_center_angle(center, radius);
        let rust_contains = rust_cap.contains(&point);
        
        // Test C++ implementation
        let cpp_contains = cpp_s2cap_contains_point(
            center.into(),
            *radius_rad,
            point.into()
        );
        
        assert_eq!(rust_contains, cpp_contains,
            "Contains point mismatch for test case {}: Rust={}, C++={}", 
            i, rust_contains, cpp_contains);
    }
}

#[test]
fn test_s2cap_area_equivalence() {
    // Test cap area calculation matches between implementations
    let test_cases = vec![
        (DVec3::new(1.0, 0.0, 0.0), 0.01),
        (DVec3::new(0.0, 1.0, 0.0), 0.1),
        (DVec3::new(0.0, 0.0, 1.0), 0.5),
        (DVec3::new(1.0, 1.0, 1.0).normalize(), 1.0),
        (DVec3::new(-1.0, -1.0, -1.0).normalize(), 2.0),
    ];
    
    for (i, (center_vec, radius_rad)) in test_cases.iter().enumerate() {
        let center = S2Point::from_normalized(*center_vec);
        let radius = S1Angle::from_radians(*radius_rad);
        
        // Test Rust implementation
        let rust_cap = S2Cap::from_center_angle(center, radius);
        let rust_area = rust_cap.area();
        
        // Test C++ implementation
        let cpp_area = cpp_s2cap_area(center.into(), *radius_rad);
        
        assert!(angles_equal_approx(rust_area, cpp_area, TOLERANCE),
            "Area mismatch for test case {}: Rust={}, C++={}, diff={}", 
            i, rust_area, cpp_area, (rust_area - cpp_area).abs());
    }
}

#[test]
fn test_s2cap_edge_cases() {
    // Test edge cases like empty caps, full caps, etc.
    
    // Empty cap
    let empty_cap = S2Cap::empty();
    let empty_center = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    let empty_radius = S1Angle::from_radians(-1.0); // Negative radius = empty
    
    let test_point = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
    
    let rust_empty_contains = empty_cap.contains(&test_point);
    let cpp_empty_contains = cpp_s2cap_contains_point(
        empty_center.into(),
        -1.0,
        test_point.into()
    );
    
    assert_eq!(rust_empty_contains, cpp_empty_contains,
        "Empty cap containment mismatch: Rust={}, C++={}", 
        rust_empty_contains, cpp_empty_contains);
    
    // Full cap (radius >= π)
    let full_center = S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0));
    let full_radius = S1Angle::from_radians(std::f64::consts::PI);
    let full_cap = S2Cap::from_center_angle(full_center, full_radius);
    
    let rust_full_contains = full_cap.contains(&test_point);
    let cpp_full_contains = cpp_s2cap_contains_point(
        full_center.into(),
        std::f64::consts::PI,
        test_point.into()
    );
    
    assert_eq!(rust_full_contains, cpp_full_contains,
        "Full cap containment mismatch: Rust={}, C++={}", 
        rust_full_contains, cpp_full_contains);
}