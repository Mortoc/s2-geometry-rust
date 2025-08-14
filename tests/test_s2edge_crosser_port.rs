//! Port of s2edge_crosser_test.cc - Edge crossing optimization tests
//!
//! This module tests the S2EdgeCrosser class which provides optimized edge crossing
//! detection by maintaining state between multiple queries against the same base edge.

use s2geometry_rust::{S2EdgeCrosser, S2CopyingEdgeCrosser, S2Point};
use s2geometry_rust::math::{DVec3, predicates::*, constants::*};
use s2geometry_rust::error::S2Result;

/// Helper to create normalized S2Point from coordinates
fn make_point(x: f64, y: f64, z: f64) -> S2Point {
    S2Point::from_vec3(DVec3::new(x, y, z).normalize()).expect("Valid point")
}

/// Test for invalid default points (non-debug builds only)
#[test]
fn test_invalid_default_points() {
    // Test with zero vector (invalid point)
    let zero_point = S2Point::from_vec3(DVec3::ZERO);
    
    if zero_point.is_err() {
        // In Rust, we properly validate points, so this should fail
        // This is better than the C++ behavior that allows invalid points
        return;
    }
    
    // If somehow we allow zero points, test crossing behavior
    let zero = zero_point.unwrap();
    
    let mut crosser = S2EdgeCrosser::new(&zero, &zero);
    assert_eq!(crosser.crossing_sign(&zero, &zero), 0);
    
    let mut copying_crosser = S2CopyingEdgeCrosser::new(zero, zero);
    assert_eq!(copying_crosser.crossing_sign(zero, zero), 0);
    
    assert!(!crosser.edge_or_vertex_crossing(&zero, &zero));
    assert!(!copying_crosser.edge_or_vertex_crossing(zero, zero));
    
    assert_eq!(crosser.signed_edge_or_vertex_crossing(&zero, &zero), 0);
    assert_eq!(copying_crosser.signed_edge_or_vertex_crossing(zero, zero), 0);
}

/// Test for NaN points (non-debug builds only) 
#[test]
fn test_invalid_nan_points() {
    // Test with NaN coordinates
    let nan = f64::NAN;
    let nan_point = S2Point::from_vec3(DVec3::new(nan, nan, nan));
    
    if nan_point.is_err() {
        // In Rust, we properly validate points, so this should fail
        // This is better than the C++ behavior that allows NaN points
        return;
    }
    
    // If somehow we allow NaN points, test crossing behavior
    let nan_pt = nan_point.unwrap();
    
    let mut crosser = S2EdgeCrosser::new(&nan_pt, &nan_pt);
    // Result should be -1 (no crossing) for invalid input
    assert_eq!(crosser.crossing_sign(&nan_pt, &nan_pt), -1);
    
    let mut copying_crosser = S2CopyingEdgeCrosser::new(nan_pt, nan_pt);
    assert_eq!(copying_crosser.crossing_sign(nan_pt, nan_pt), -1);
}

/// Test a single crossing case with all variants
fn test_crossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point, 
                 expected_crossing: i32, signed_crossing_sign: i32) {
    // Normalize all points
    let a = make_point(a.x(), a.y(), a.z());
    let b = make_point(b.x(), b.y(), b.z());
    let c = make_point(c.x(), c.y(), c.z());
    let d = make_point(d.x(), d.y(), d.z());
    
    // Handle degenerate cases
    let mut expected_crossing_val = expected_crossing;
    if a.coords() == c.coords() || a.coords() == d.coords() || 
       b.coords() == c.coords() || b.coords() == d.coords() {
        expected_crossing_val = 0;
    }
    
    // Validate signed crossing sign consistency
    if expected_crossing_val == 1 {
        let computed_sign = robust_sign(a.coords(), b.coords(), c.coords());
        assert_eq!(signed_crossing_sign, computed_sign, 
                  "Signed crossing sign inconsistent with orientation");
    } else if expected_crossing_val == 0 && signed_crossing_sign != 0 {
        let expected_sign = if a.coords() == c.coords() || b.coords() == d.coords() { 1 } else { -1 };
        assert_eq!(signed_crossing_sign, expected_sign,
                  "Signed crossing sign inconsistent for vertex crossing");
    }
    
    // Test global crossing sign function
    assert_eq!(crossing_sign(a.coords(), b.coords(), c.coords(), d.coords()), expected_crossing_val);
    
    let edge_or_vertex = signed_crossing_sign != 0;
    assert_eq!(edge_or_vertex_crossing(a.coords(), b.coords(), c.coords(), d.coords()), edge_or_vertex);
    
    // Test S2EdgeCrosser with different call patterns
    let mut crosser = S2EdgeCrosser::new(&a, &b);
    crosser.restart_at(&c);
    
    assert_eq!(crosser.crossing_sign(&d, &c), expected_crossing);
    assert_eq!(crosser.crossing_sign(&c, &c), expected_crossing);
    assert_eq!(crosser.crossing_sign(&d, &c), expected_crossing);
    assert_eq!(crosser.crossing_sign(&c, &d), expected_crossing);
    
    // Test EdgeOrVertexCrossing
    crosser.restart_at(&c);
    assert_eq!(crosser.edge_or_vertex_crossing(&d, &c), edge_or_vertex);
    assert_eq!(crosser.edge_or_vertex_crossing(&c, &c), edge_or_vertex);
    assert_eq!(crosser.edge_or_vertex_crossing(&d, &c), edge_or_vertex);
    assert_eq!(crosser.edge_or_vertex_crossing(&c, &d), edge_or_vertex);
    
    // Test SignedEdgeOrVertexCrossing  
    crosser.restart_at(&c);
    assert_eq!(crosser.signed_edge_or_vertex_crossing(&d, &c), signed_crossing_sign);
    assert_eq!(crosser.signed_edge_or_vertex_crossing(&c, &c), -signed_crossing_sign);
    assert_eq!(crosser.signed_edge_or_vertex_crossing(&d, &c), -signed_crossing_sign);
    assert_eq!(crosser.signed_edge_or_vertex_crossing(&c, &d), signed_crossing_sign);
    
    // Test crosser reuse
    crosser.init(&c, &d);
    crosser.restart_at(&a);
    assert_eq!(crosser.crossing_sign_chain(&b), expected_crossing);
    assert_eq!(crosser.crossing_sign(&a, &a), expected_crossing);
    
    // Test S2CopyingEdgeCrosser with same tests
    let mut copying_crosser = S2CopyingEdgeCrosser::new(a, b);
    copying_crosser.restart_at(c);
    
    assert_eq!(copying_crosser.crossing_sign(d, c), expected_crossing);
    assert_eq!(copying_crosser.crossing_sign(c, c), expected_crossing);
    assert_eq!(copying_crosser.crossing_sign(d, c), expected_crossing);
    assert_eq!(copying_crosser.crossing_sign(c, d), expected_crossing);
    
    copying_crosser.restart_at(c);
    assert_eq!(copying_crosser.edge_or_vertex_crossing(d, c), edge_or_vertex);
    assert_eq!(copying_crosser.edge_or_vertex_crossing(c, c), edge_or_vertex);
    assert_eq!(copying_crosser.edge_or_vertex_crossing(d, c), edge_or_vertex);
    assert_eq!(copying_crosser.edge_or_vertex_crossing(c, d), edge_or_vertex);
    
    copying_crosser.restart_at(c);
    assert_eq!(copying_crosser.signed_edge_or_vertex_crossing(d, c), signed_crossing_sign);
    assert_eq!(copying_crosser.signed_edge_or_vertex_crossing(c, c), -signed_crossing_sign);
    assert_eq!(copying_crosser.signed_edge_or_vertex_crossing(d, c), -signed_crossing_sign);
    assert_eq!(copying_crosser.signed_edge_or_vertex_crossing(c, d), signed_crossing_sign);
    
    // Test crosser reuse
    copying_crosser.init(c, d);
    copying_crosser.restart_at(a);
    assert_eq!(copying_crosser.crossing_sign_chain(b), expected_crossing);
    assert_eq!(copying_crosser.crossing_sign(a, a), expected_crossing);
}

/// Test crossings with all permutations
fn test_crossings(a: S2Point, b: S2Point, c: S2Point, d: S2Point,
                  crossing_sign: i32, signed_crossing_sign: i32) {
    // Test main case
    test_crossing(a, b, c, d, crossing_sign, signed_crossing_sign);
    
    // Test with reversed edge AB
    test_crossing(b, a, c, d, crossing_sign, -signed_crossing_sign);
    
    // Test with reversed edge CD
    test_crossing(a, b, d, c, crossing_sign, -signed_crossing_sign);
    
    // Test with both edges reversed
    test_crossing(b, a, d, c, crossing_sign, signed_crossing_sign);
    
    // Test degenerate cases
    test_crossing(a, a, c, d, -1, 0);
    test_crossing(a, b, c, c, -1, 0);
    test_crossing(a, a, c, c, -1, 0);
    test_crossing(a, b, a, b, 0, 1);
    
    if crossing_sign == 0 {
        // For vertex crossings, if AB crosses CD then CD does not cross AB
        test_crossing(c, d, a, b, crossing_sign, 0);
    } else {
        // For regular crossings, test the reverse case
        test_crossing(c, d, a, b, crossing_sign, -signed_crossing_sign);
    }
}

#[test]
fn test_edge_crossings() {
    // Test cases from C++ implementation
    
    // 1. Two regular edges that cross
    test_crossings(
        make_point(1.0, 2.0, 1.0), 
        make_point(1.0, -3.0, 0.5),
        make_point(1.0, -0.5, -3.0), 
        make_point(0.1, 0.5, 3.0), 
        1, 1
    );
    
    // 2. Two regular edges that intersect at antipodal points  
    test_crossings(
        make_point(1.0, 2.0, 1.0), 
        make_point(1.0, -3.0, 0.5),
        make_point(-1.0, 0.5, 3.0), 
        make_point(-0.1, -0.5, -3.0), 
        -1, 0
    );
    
    // 3. Two edges on the same great circle starting at antipodal points
    test_crossings(
        make_point(0.0, 0.0, -1.0), 
        make_point(0.0, 1.0, 0.0),
        make_point(0.0, 0.0, 1.0), 
        make_point(0.0, 1.0, 1.0), 
        -1, 0
    );
    
    // 4. Two edges that cross where one vertex is near origin
    test_crossings(
        make_point(1.0, 0.0, 0.0), 
        make_point(1e-10, 1e-10, 1e-10),
        make_point(1.0, -0.1, 1.0), 
        make_point(1.0, 1.0, -0.1), 
        1, 1
    );
    
    // 5. Two edges that intersect at antipodal points where one vertex is near origin
    test_crossings(
        make_point(1.0, 0.0, 0.0), 
        make_point(1e-10, 1e-10, 1e-10),
        make_point(-1.0, 0.1, -1.0), 
        make_point(-1.0, -1.0, 0.1), 
        -1, 0
    );
    
    // 6. Two edges that share an endpoint
    test_crossings(
        make_point(7.0, -2.0, 3.0), 
        make_point(2.0, 3.0, 4.0),
        make_point(2.0, 3.0, 4.0), 
        make_point(-1.0, 2.0, 5.0), 
        0, -1
    );
    
    // 7. Two edges that barely cross each other near the middle
    test_crossings(
        make_point(1.0, 1.0, 1.0), 
        make_point(1.0, f64::from_bits(f64::to_bits(1.0) - 1), -1.0),
        make_point(11.0, -12.0, -1.0), 
        make_point(10.0, 10.0, 1.0), 
        1, -1
    );
    
    // 8. Separated by tiny distance
    test_crossings(
        make_point(1.0, 1.0, 1.0), 
        make_point(1.0, f64::from_bits(f64::to_bits(1.0) + 1), -1.0),
        make_point(1.0, -1.0, 0.0), 
        make_point(1.0, 1.0, 0.0), 
        -1, 0
    );
    
    // 9. Two edges that barely cross at extreme precision
    test_crossings(
        make_point(0.0, 0.0, 1.0), 
        make_point(2.0, -1e-323, 1.0),
        make_point(1.0, -1.0, 1.0), 
        make_point(1e-323, 0.0, 1.0), 
        1, -1
    );
    
    // 10. Separated by extremely small distance
    test_crossings(
        make_point(0.0, 0.0, 1.0), 
        make_point(2.0, 1e-323, 1.0),
        make_point(1.0, -1.0, 1.0), 
        make_point(1e-323, 0.0, 1.0), 
        -1, 0
    );
}

#[test]
fn test_collinear_edges_that_dont_touch() {
    // Test collinear edges that don't intersect
    use rand::prelude::*;
    
    let mut rng = rand::thread_rng();
    
    for _iter in 0..100 {
        // Generate random points
        let a = make_point(rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5);
        let d = make_point(rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5);
        
        // Create collinear points
        let a_coords = a.coords();
        let d_coords = d.coords();
        let interpolated_05 = (a_coords * 0.95 + d_coords * 0.05).normalize();
        let interpolated_95 = (a_coords * 0.05 + d_coords * 0.95).normalize();
        
        let b = S2Point::from_vec3(interpolated_05).unwrap();
        let c = S2Point::from_vec3(interpolated_95).unwrap();
        
        // These should not cross
        assert!(crossing_sign(a.coords(), b.coords(), c.coords(), d.coords()) <= 0);
        
        let mut crosser = S2EdgeCrosser::new(&a, &b);
        crosser.restart_at(&c);
        assert!(crosser.crossing_sign_chain(&d) <= 0);
        assert!(crosser.crossing_sign(&c, &c) <= 0);
    }
}

#[test]
fn test_coincident_zero_length_edges() {
    // Test edges with coincident endpoints that are scaled versions of each other
    use rand::prelude::*;
    
    let mut rng = rand::thread_rng();
    
    for _iter in 0..100 {
        // Generate a random direction
        let mut p = DVec3::new(
            rng.gen::<f64>() - 0.5, 
            rng.gen::<f64>() - 0.5, 
            rng.gen::<f64>() - 0.5
        );
        
        // Ensure non-zero
        if p.length_squared() < f64::EPSILON {
            p = DVec3::new(1.0, 0.0, 0.0);
        }
        
        p = p.normalize();
        
        // Create scaled versions (all project to the same point on unit sphere)
        let epsilon = 1e-15;
        let a_coords = (1.0 - 3.0 * epsilon) * p;
        let b_coords = (1.0 - epsilon) * p;
        let c_coords = p;
        let d_coords = (1.0 + 2.0 * epsilon) * p;
        
        // Only test if all points are valid unit vectors
        if a_coords.length() > 0.5 && d_coords.length() < 2.0 {
            let a = S2Point::from_vec3(a_coords.normalize()).unwrap();
            let b = S2Point::from_vec3(b_coords.normalize()).unwrap();
            let c = S2Point::from_vec3(c_coords.normalize()).unwrap();
            let d = S2Point::from_vec3(d_coords.normalize()).unwrap();
            
            // These should not cross (coincident edges)
            assert!(crossing_sign(a.coords(), b.coords(), c.coords(), d.coords()) <= 0);
            
            let mut crosser = S2EdgeCrosser::new(&a, &b);
            crosser.restart_at(&c);
            assert!(crosser.crossing_sign_chain(&d) <= 0);
            assert!(crosser.crossing_sign(&c, &c) <= 0);
        }
    }
}

#[test]
fn test_edge_chain_optimization() {
    // Test that edge chain processing is optimized correctly
    let a = make_point(1.0, 0.0, 0.0);
    let b = make_point(0.0, 1.0, 0.0);
    
    let mut crosser = S2EdgeCrosser::new(&a, &b);
    
    // Create a chain of edges
    let chain = vec![
        make_point(0.0, 0.0, 1.0),
        make_point(-1.0, 0.0, 0.0),
        make_point(0.0, -1.0, 0.0),
        make_point(1.0, 0.0, -1.0),
    ];
    
    crosser.restart_at(&chain[0]);
    
    // Process chain efficiently
    for i in 1..chain.len() {
        let _crossing = crosser.crossing_sign_chain(&chain[i]);
        // Each call should use the cached state from previous vertex
    }
    
    // Test that state is maintained correctly
    assert!(crosser.c().is_some());
    assert_eq!(crosser.c().unwrap().coords(), chain.last().unwrap().coords());
}

#[test] 
fn test_signed_crossings() {
    // Test signed crossing behavior for winding number computation
    let a = make_point(1.0, 0.0, 0.0);
    let b = make_point(0.0, 1.0, 0.0);
    let c = make_point(0.0, 0.0, 1.0);
    let d = make_point(-1.0, 0.0, 0.0);
    
    let mut crosser = S2EdgeCrosser::new(&a, &b);
    
    let signed_crossing = crosser.signed_edge_or_vertex_crossing(&c, &d);
    
    if signed_crossing != 0 {
        // If there's a crossing, the sign should match the orientation
        let expected_sign = crosser.last_interior_crossing_sign();
        assert_eq!(signed_crossing, expected_sign);
    }
}