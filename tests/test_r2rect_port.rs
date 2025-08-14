//! Port of C++ r2rect_test.cc to Rust
//!
//! This test file ports all the R2Rect tests from the original C++ implementation
//! to verify complete compatibility with the Google S2 library.

use s2geometry_rust::r2::{R2Point, R2Rect};
use s2geometry_rust::interval::R1Interval;

/// Test helper function that mirrors the C++ TestIntervalOps function
///
/// Tests all interval operations on a pair of rectangles:
/// - Contains(), InteriorContains(), Intersects(), InteriorIntersects()
/// - Union() and Intersection() operations
/// - AddRect() and AddPoint() operations
fn test_interval_ops(
    x: &R2Rect,
    y: &R2Rect,
    expected_rexion: &str,
    expected_union: &R2Rect,
    expected_intersection: &R2Rect,
) {
    // Test all of the interval operations on the given pair of intervals.
    // "expected_rexion" is a sequence of "T" and "F" characters corresponding
    // to the expected results of Contains(), InteriorContains(), Intersects(),
    // and InteriorIntersects() respectively.

    let chars: Vec<char> = expected_rexion.chars().collect();
    assert_eq!(chars.len(), 4, "expected_rexion must have exactly 4 characters");

    assert_eq!(chars[0] == 'T', x.contains_rect(y), "Contains() test failed");
    assert_eq!(chars[1] == 'T', x.interior_contains_rect(y), "InteriorContains() test failed");
    assert_eq!(chars[2] == 'T', x.intersects(y), "Intersects() test failed");
    assert_eq!(chars[3] == 'T', x.interior_intersects(y), "InteriorIntersects() test failed");

    assert_eq!(x.union(y) == *x, x.contains_rect(y), "Union contains check failed");
    assert_eq!(!x.intersection(y).is_empty(), x.intersects(y), "Intersection empty check failed");

    assert_eq!(*expected_union, x.union(y), "Union result mismatch");
    assert_eq!(*expected_intersection, x.intersection(y), "Intersection result mismatch");

    let mut r = *x;
    r.add_rect(y);
    assert_eq!(*expected_union, r, "AddRect result mismatch");

    if y.get_size() == R2Point::new(0.0, 0.0) {
        let mut r = *x;
        r.add_point(y.lo());
        assert_eq!(*expected_union, r, "AddPoint result mismatch");
    }
}

#[test]
fn test_empty_rectangles() {
    // Test basic properties of empty rectangles.
    let empty = R2Rect::empty();
    assert!(empty.is_valid());
    assert!(empty.is_empty());
    assert_eq!(empty, empty);
}

#[test]
fn test_constructors_and_accessors() {
    // Check various constructors and accessor methods.
    let r = R2Rect::new(R2Point::new(0.1, 0.0), R2Point::new(0.25, 1.0));
    assert_eq!(0.1, r.x().lo());
    assert_eq!(0.25, r.x().hi());
    assert_eq!(0.0, r.y().lo());
    assert_eq!(1.0, r.y().hi());

    assert_eq!(0.1, r[0][0]);
    assert_eq!(0.25, r[0][1]);
    assert_eq!(0.0, r[1][0]);
    assert_eq!(1.0, r[1][1]);

    assert_eq!(R1Interval::new(0.1, 0.25), r.x());
    assert_eq!(R1Interval::new(0.0, 1.0), r.y());

    assert_eq!(R1Interval::new(0.1, 0.25), r[0]);
    assert_eq!(R1Interval::new(0.0, 1.0), r[1]);

    let mut r = r;
    r[0] = R1Interval::new(3.0, 4.0);
    *r[1].get_mut(0) = 5.0;
    *r[1].get_mut(1) = 6.0;
    assert_eq!(R1Interval::new(3.0, 4.0), r[0]);
    assert_eq!(R1Interval::new(5.0, 6.0), r[1]);

    assert_eq!(r, r);
    assert_ne!(r, R2Rect::empty());

    let r2 = R2Rect::default();
    assert!(r2.is_empty());
    assert_eq!(r2, R2Rect::empty());
}

#[test]
fn test_from_center_size() {
    // FromCenterSize()
    assert!(R2Rect::from_center_size(R2Point::new(0.3, 0.5), R2Point::new(0.2, 0.4))
        .approx_equals(&R2Rect::new(R2Point::new(0.2, 0.3), R2Point::new(0.4, 0.7)), 1e-15));
    
    assert!(R2Rect::from_center_size(R2Point::new(1.0, 0.1), R2Point::new(0.0, 2.0))
        .approx_equals(&R2Rect::new(R2Point::new(1.0, -0.9), R2Point::new(1.0, 1.1)), 1e-15));
}

#[test]
fn test_from_point() {
    // FromPoint(), FromPointPair()
    let d1 = R2Rect::new(R2Point::new(0.1, 0.0), R2Point::new(0.25, 1.0));
    assert_eq!(R2Rect::new(d1.lo(), d1.lo()), R2Rect::from_point(d1.lo()));
    
    assert_eq!(
        R2Rect::new(R2Point::new(0.15, 0.3), R2Point::new(0.35, 0.9)),
        R2Rect::from_point_pair(R2Point::new(0.15, 0.9), R2Point::new(0.35, 0.3))
    );
    
    assert_eq!(
        R2Rect::new(R2Point::new(0.12, 0.0), R2Point::new(0.83, 0.5)),
        R2Rect::from_point_pair(R2Point::new(0.83, 0.0), R2Point::new(0.12, 0.5))
    );
}

#[test]
fn test_simple_predicates() {
    // GetCenter(), GetVertex(), Contains(R2Point), InteriorContains(R2Point).
    let sw1 = R2Point::new(0.0, 0.25);
    let ne1 = R2Point::new(0.5, 0.75);
    let r1 = R2Rect::new(sw1, ne1);

    assert_eq!(R2Point::new(0.25, 0.5), r1.get_center());
    assert_eq!(R2Point::new(0.0, 0.25), r1.get_vertex(0));
    assert_eq!(R2Point::new(0.5, 0.25), r1.get_vertex(1));
    assert_eq!(R2Point::new(0.5, 0.75), r1.get_vertex(2));
    assert_eq!(R2Point::new(0.0, 0.75), r1.get_vertex(3));
    
    assert!(r1.contains(R2Point::new(0.2, 0.4)));
    assert!(!r1.contains(R2Point::new(0.2, 0.8)));
    assert!(!r1.contains(R2Point::new(-0.1, 0.4)));
    assert!(!r1.contains(R2Point::new(0.6, 0.1)));
    assert!(r1.contains(sw1));
    assert!(r1.contains(ne1));
    assert!(!r1.interior_contains(sw1));
    assert!(!r1.interior_contains(ne1));

    // Make sure that GetVertex() returns vertices in CCW order.
    for k in 0..4 {
        let a = r1.get_vertex((k - 1) & 3);
        let b = r1.get_vertex(k);
        let c = r1.get_vertex((k + 1) & 3);
        // Check that (b - a) × (c - a) > 0 (counterclockwise)
        assert!((b - a).cross_prod(&(c - a)) > 0.0, "Vertices not in CCW order at k={}", k);
    }
}

#[test]
fn test_interval_operations() {
    // Contains(R2Rect), InteriorContains(R2Rect),
    // Intersects(), InteriorIntersects(), Union(), Intersection().
    //
    // Much more testing of these methods is done in s1interval_test
    // and r1interval_test.

    let empty = R2Rect::empty();
    let sw1 = R2Point::new(0.0, 0.25);
    let ne1 = R2Point::new(0.5, 0.75);
    let r1 = R2Rect::new(sw1, ne1);
    let r1_mid = R2Rect::new(R2Point::new(0.25, 0.5), R2Point::new(0.25, 0.5));
    let r_sw1 = R2Rect::new(sw1, sw1);
    let r_ne1 = R2Rect::new(ne1, ne1);

    test_interval_ops(&r1, &r1_mid, "TTTT", &r1, &r1_mid);
    test_interval_ops(&r1, &r_sw1, "TFTF", &r1, &r_sw1);
    test_interval_ops(&r1, &r_ne1, "TFTF", &r1, &r_ne1);

    assert_eq!(R2Rect::new(R2Point::new(0.0, 0.25), R2Point::new(0.5, 0.75)), r1);
    
    test_interval_ops(
        &r1,
        &R2Rect::new(R2Point::new(0.45, 0.1), R2Point::new(0.75, 0.3)),
        "FFTT",
        &R2Rect::new(R2Point::new(0.0, 0.1), R2Point::new(0.75, 0.75)),
        &R2Rect::new(R2Point::new(0.45, 0.25), R2Point::new(0.5, 0.3))
    );
    
    test_interval_ops(
        &r1,
        &R2Rect::new(R2Point::new(0.5, 0.1), R2Point::new(0.7, 0.3)),
        "FFTF",
        &R2Rect::new(R2Point::new(0.0, 0.1), R2Point::new(0.7, 0.75)),
        &R2Rect::new(R2Point::new(0.5, 0.25), R2Point::new(0.5, 0.3))
    );
    
    test_interval_ops(
        &r1,
        &R2Rect::new(R2Point::new(0.45, 0.1), R2Point::new(0.7, 0.25)),
        "FFTF",
        &R2Rect::new(R2Point::new(0.0, 0.1), R2Point::new(0.7, 0.75)),
        &R2Rect::new(R2Point::new(0.45, 0.25), R2Point::new(0.5, 0.25))
    );

    test_interval_ops(
        &R2Rect::new(R2Point::new(0.1, 0.2), R2Point::new(0.1, 0.3)),
        &R2Rect::new(R2Point::new(0.15, 0.7), R2Point::new(0.2, 0.8)),
        "FFFF",
        &R2Rect::new(R2Point::new(0.1, 0.2), R2Point::new(0.2, 0.8)),
        &empty
    );

    // Check that the intersection of two rectangles that overlap in x but not y
    // is valid, and vice versa.
    test_interval_ops(
        &R2Rect::new(R2Point::new(0.1, 0.2), R2Point::new(0.4, 0.5)),
        &R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(0.2, 0.1)),
        "FFFF",
        &R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(0.4, 0.5)),
        &empty
    );
    
    test_interval_ops(
        &R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(0.1, 0.3)),
        &R2Rect::new(R2Point::new(0.2, 0.1), R2Point::new(0.3, 0.4)),
        "FFFF",
        &R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(0.3, 0.4)),
        &empty
    );
}

#[test]
fn test_add_point() {
    // AddPoint()
    let sw1 = R2Point::new(0.0, 0.25);
    let ne1 = R2Point::new(0.5, 0.75);
    let r1 = R2Rect::new(sw1, ne1);

    let mut r2 = R2Rect::empty();
    r2.add_point(R2Point::new(0.0, 0.25));
    r2.add_point(R2Point::new(0.5, 0.25));
    r2.add_point(R2Point::new(0.0, 0.75));
    r2.add_point(R2Point::new(0.1, 0.4));
    assert_eq!(r1, r2);
}

#[test]
fn test_project() {
    let r1 = R2Rect::from_intervals(R1Interval::new(0.0, 0.5), R1Interval::new(0.25, 0.75));

    assert_eq!(R2Point::new(0.0, 0.25), r1.project(R2Point::new(-0.01, 0.24)));
    assert_eq!(R2Point::new(0.0, 0.48), r1.project(R2Point::new(-5.0, 0.48)));
    assert_eq!(R2Point::new(0.0, 0.75), r1.project(R2Point::new(-5.0, 2.48)));
    assert_eq!(R2Point::new(0.19, 0.75), r1.project(R2Point::new(0.19, 2.48)));
    assert_eq!(R2Point::new(0.5, 0.75), r1.project(R2Point::new(6.19, 2.48)));
    assert_eq!(R2Point::new(0.5, 0.53), r1.project(R2Point::new(6.19, 0.53)));
    assert_eq!(R2Point::new(0.5, 0.25), r1.project(R2Point::new(6.19, -2.53)));
    assert_eq!(R2Point::new(0.33, 0.25), r1.project(R2Point::new(0.33, -2.53)));
    assert_eq!(R2Point::new(0.33, 0.37), r1.project(R2Point::new(0.33, 0.37)));
}

#[test]
fn test_expanded() {
    // Expanded()
    assert!(R2Rect::empty().expanded(R2Point::new(0.1, 0.3)).is_empty());
    assert!(R2Rect::empty().expanded(R2Point::new(-0.1, -0.3)).is_empty());
    
    assert!(R2Rect::new(R2Point::new(0.2, 0.4), R2Point::new(0.3, 0.7))
        .expanded(R2Point::new(0.1, 0.3))
        .approx_equals(&R2Rect::new(R2Point::new(0.1, 0.1), R2Point::new(0.4, 1.0)), 1e-15));
    
    assert!(R2Rect::new(R2Point::new(0.2, 0.4), R2Point::new(0.3, 0.7))
        .expanded(R2Point::new(-0.1, 0.3))
        .is_empty());
    
    assert!(R2Rect::new(R2Point::new(0.2, 0.4), R2Point::new(0.3, 0.7))
        .expanded(R2Point::new(0.1, -0.2))
        .is_empty());
    
    assert!(R2Rect::new(R2Point::new(0.2, 0.4), R2Point::new(0.3, 0.7))
        .expanded(R2Point::new(0.1, -0.1))
        .approx_equals(&R2Rect::new(R2Point::new(0.1, 0.5), R2Point::new(0.4, 0.6)), 1e-15));
    
    assert!(R2Rect::new(R2Point::new(0.2, 0.4), R2Point::new(0.3, 0.7))
        .expanded_by_margin(0.1)
        .approx_equals(&R2Rect::new(R2Point::new(0.1, 0.3), R2Point::new(0.4, 0.8)), 1e-15));
}

#[test]
fn test_vertex_ordering() {
    // Test that vertices are returned in counter-clockwise order
    let rect = R2Rect::new(R2Point::new(1.0, 2.0), R2Point::new(3.0, 4.0));
    
    // Vertices should be: (1,2), (3,2), (3,4), (1,4) in CCW order
    assert_eq!(rect.get_vertex(0), R2Point::new(1.0, 2.0)); // lower-left
    assert_eq!(rect.get_vertex(1), R2Point::new(3.0, 2.0)); // lower-right  
    assert_eq!(rect.get_vertex(2), R2Point::new(3.0, 4.0)); // upper-right
    assert_eq!(rect.get_vertex(3), R2Point::new(1.0, 4.0)); // upper-left
    
    // Test wrapping
    assert_eq!(rect.get_vertex(4), rect.get_vertex(0));
    assert_eq!(rect.get_vertex(-1), rect.get_vertex(3));
}

#[test]
fn test_vertex_ij() {
    let rect = R2Rect::new(R2Point::new(1.0, 2.0), R2Point::new(3.0, 4.0));
    
    assert_eq!(rect.get_vertex_ij(0, 0), R2Point::new(1.0, 2.0)); // left, bottom
    assert_eq!(rect.get_vertex_ij(1, 0), R2Point::new(3.0, 2.0)); // right, bottom
    assert_eq!(rect.get_vertex_ij(1, 1), R2Point::new(3.0, 4.0)); // right, top
    assert_eq!(rect.get_vertex_ij(0, 1), R2Point::new(1.0, 4.0)); // left, top
}

#[test]
fn test_r2point_operations() {
    let p1 = R2Point::new(3.0, 4.0);
    let p2 = R2Point::new(1.0, 2.0);
    
    // Test basic arithmetic
    assert_eq!(p1 + p2, R2Point::new(4.0, 6.0));
    assert_eq!(p1 - p2, R2Point::new(2.0, 2.0));
    assert_eq!(p1 * 2.0, R2Point::new(6.0, 8.0));
    
    // Test dot product
    assert_eq!(p1.dot_prod(&p2), 11.0); // 3*1 + 4*2 = 11
    
    // Test cross product  
    assert_eq!(p1.cross_prod(&p2), 2.0); // 3*2 - 4*1 = 2
    
    // Test orthogonal vector
    assert_eq!(p1.ortho(), R2Point::new(-4.0, 3.0));
    
    // Test distance
    assert_eq!(p1.distance_squared(&p2), 8.0); // (3-1)^2 + (4-2)^2 = 8
    assert_eq!(p1.distance(&p2), 8.0_f64.sqrt());
    
    // Test norm
    assert_eq!(p1.norm_squared(), 25.0); // 3^2 + 4^2 = 25
    assert_eq!(p1.norm(), 5.0); // sqrt(25) = 5
}

#[test]
fn test_display_formatting() {
    let point = R2Point::new(1.5, -2.3);
    assert_eq!(format!("{}", point), "(1.5, -2.3)");
    
    let rect = R2Rect::new(R2Point::new(0.0, 1.0), R2Point::new(2.0, 3.0));
    assert_eq!(format!("{}", rect), "[0, 2] x [1, 3]");
    
    let empty = R2Rect::empty();
    assert_eq!(format!("{}", empty), "EMPTY");
}

#[test]
fn test_comprehensive_operations() {
    // Create some test rectangles
    let r1 = R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(2.0, 2.0));
    let r2 = R2Rect::new(R2Point::new(1.0, 1.0), R2Point::new(3.0, 3.0));
    let r3 = R2Rect::new(R2Point::new(4.0, 4.0), R2Point::new(5.0, 5.0));
    
    // Test intersections
    assert!(r1.intersects(&r2));
    assert!(!r1.intersects(&r3));
    
    // Test unions
    let union_12 = r1.union(&r2);
    assert_eq!(union_12, R2Rect::new(R2Point::new(0.0, 0.0), R2Point::new(3.0, 3.0)));
    
    // Test intersections
    let intersect_12 = r1.intersection(&r2);
    assert_eq!(intersect_12, R2Rect::new(R2Point::new(1.0, 1.0), R2Point::new(2.0, 2.0)));
    
    let intersect_13 = r1.intersection(&r3);
    assert!(intersect_13.is_empty());
    
    // Test containment
    let small = R2Rect::new(R2Point::new(0.5, 0.5), R2Point::new(1.5, 1.5));
    assert!(r1.contains_rect(&small));
    assert!(!small.contains_rect(&r1));
}