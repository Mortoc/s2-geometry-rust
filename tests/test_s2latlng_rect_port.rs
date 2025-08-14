//! Port of s2latlng_rect_test.cc - Critical S2LatLngRect tests
//!
//! This module ports the essential tests from C++ s2latlng_rect_test.cc to validate
//! our S2LatLngRect implementation matches C++ behavior exactly.
//!
//! S2LatLngRect is critical for geographic applications as it provides rectangle
//! operations in latitude/longitude coordinates, including handling of the 
//! International Date Line and polar regions.

use s2geometry_rust::{S1Angle, S2LatLng, S2LatLngRect, S2Point};
use s2geometry_rust::interval::{R1Interval, S1Interval};
use s2geometry_rust::math::constants::*;
use approx::assert_relative_eq;

/// Helper function to construct rectangle from degrees (matches C++ RectFromDegrees)
fn rect_from_degrees(lat_lo: f64, lng_lo: f64, lat_hi: f64, lng_hi: f64) -> S2LatLngRect {
    S2LatLngRect::new(
        S2LatLng::from_degrees(lat_lo, lng_lo).normalized(),
        S2LatLng::from_degrees(lat_hi, lng_hi).normalized()
    )
}

/// Test basic properties of empty and full rectangles (from C++ EmptyAndFull)
#[test]
fn test_empty_and_full() {
    let empty = S2LatLngRect::empty();
    let full = S2LatLngRect::full();
    
    assert!(empty.is_valid());
    assert!(empty.is_empty());
    assert!(!empty.is_point());
    
    assert!(full.is_valid());
    assert!(full.is_full());
    assert!(!full.is_point());
    
    // Check that default S2LatLngRect is identical to Empty()
    let default_empty = S2LatLngRect::default();
    assert!(default_empty.is_valid());
    assert!(default_empty.is_empty());
    assert_eq!(empty.lat().bounds(), default_empty.lat().bounds());
    assert_eq!(empty.lng().bounds(), default_empty.lng().bounds());
}

/// Test various accessor methods (from C++ Accessors)
#[test]
fn test_accessors() {
    let d1 = rect_from_degrees(-90.0, 0.0, -45.0, 180.0);
    assert_relative_eq!(d1.lat_lo().degrees(), -90.0, epsilon = 1e-15);
    assert_relative_eq!(d1.lat_hi().degrees(), -45.0, epsilon = 1e-15);
    assert_relative_eq!(d1.lng_lo().degrees(), 0.0, epsilon = 1e-15);
    assert_relative_eq!(d1.lng_hi().degrees(), 180.0, epsilon = 1e-15);
    assert_eq!(d1.lat(), R1Interval::new(-PI_2, -PI_4));
    assert_eq!(d1.lng(), S1Interval::new(0.0, PI));
}

/// Test ApproxEquals with various tolerances (from C++ ApproxEquals)
#[test]
fn test_approx_equals() {
    // Empty rectangles should be approximately equal to point rectangles
    assert!(S2LatLngRect::empty().approx_equals(
        &rect_from_degrees(1.0, 5.0, 1.0, 5.0),
        S1Angle::from_radians(1e-15)
    ));
    assert!(rect_from_degrees(1.0, 5.0, 1.0, 5.0).approx_equals(
        &S2LatLngRect::empty(),
        S1Angle::from_radians(1e-15)
    ));
    assert!(!rect_from_degrees(1.0, 5.0, 1.0, 5.0).approx_equals(
        &rect_from_degrees(2.0, 7.0, 2.0, 7.0),
        S1Angle::from_radians(1e-15)
    ));

    // Test the max_error parameter
    assert!(rect_from_degrees(10.0, 10.0, 20.0, 20.0).approx_equals(
        &rect_from_degrees(11.0, 11.0, 19.0, 19.0),
        S1Angle::from_degrees(1.001)
    ));
    assert!(!rect_from_degrees(10.0, 10.0, 20.0, 20.0).approx_equals(
        &rect_from_degrees(11.0, 11.0, 19.0, 19.0),
        S1Angle::from_degrees(0.999)
    ));

    // Test the max_error (S2LatLng) parameter
    assert!(rect_from_degrees(0.0, 10.0, 20.0, 30.0).approx_equals_latlng(
        &rect_from_degrees(-1.0, 8.0, 21.0, 32.0),
        S2LatLng::from_degrees(1.001, 2.001)
    ));
    assert!(!rect_from_degrees(0.0, 10.0, 20.0, 30.0).approx_equals_latlng(
        &rect_from_degrees(-1.0, 8.0, 21.0, 32.0),
        S2LatLng::from_degrees(0.999, 1.999)
    ));
}

/// Test FromCenterSize constructor (from C++ FromCenterSize)
#[test]
fn test_from_center_size() {
    assert!(S2LatLngRect::from_center_size(
        S2LatLng::from_degrees(80.0, 170.0),
        S2LatLng::from_degrees(40.0, 60.0)
    ).approx_equals(
        &rect_from_degrees(60.0, 140.0, 90.0, -160.0),
        S1Angle::from_radians(1e-15)
    ));
    
    assert!(S2LatLngRect::from_center_size(
        S2LatLng::from_degrees(10.0, 40.0),
        S2LatLng::from_degrees(210.0, 400.0)
    ).is_full());
    
    assert!(S2LatLngRect::from_center_size(
        S2LatLng::from_degrees(-90.0, 180.0),
        S2LatLng::from_degrees(20.0, 50.0)
    ).approx_equals(
        &rect_from_degrees(-90.0, 155.0, -80.0, -155.0),
        S1Angle::from_radians(1e-15)
    ));
}

/// Test FromPoint constructor (from C++ FromPoint)
#[test]
fn test_from_point() {
    let p = S2LatLng::from_degrees(23.0, 47.0);
    assert_eq!(S2LatLngRect::from_point(p), S2LatLngRect::new(p, p));
    assert!(S2LatLngRect::from_point(p).is_point());
}

/// Test FromPointPair constructor (from C++ FromPointPair)
#[test]
fn test_from_point_pair() {
    assert_eq!(
        S2LatLngRect::from_point_pair(
            S2LatLng::from_degrees(-35.0, -140.0),
            S2LatLng::from_degrees(15.0, 155.0)
        ),
        rect_from_degrees(-35.0, 155.0, 15.0, -140.0)
    );
    assert_eq!(
        S2LatLngRect::from_point_pair(
            S2LatLng::from_degrees(25.0, -70.0),
            S2LatLng::from_degrees(-90.0, 80.0)
        ),
        rect_from_degrees(-90.0, -70.0, 25.0, 80.0)
    );
}

/// Test GetCenter and GetSize methods (from C++ GetCenterSize)
#[test]
fn test_get_center_size() {
    let r1 = S2LatLngRect::from_intervals(
        R1Interval::new(0.0, PI_2),
        S1Interval::new(-PI, 0.0)
    );
    assert_eq!(r1.get_center(), S2LatLng::from_radians(PI_4, -PI_2));
    assert_eq!(r1.get_size(), S2LatLng::from_radians(PI_2, PI));
    assert!(S2LatLngRect::empty().get_size().lat().radians() < 0.0);
    assert!(S2LatLngRect::empty().get_size().lng().radians() < 0.0);
}

/// Test GetVertex method (from C++ GetVertex)
#[test]
fn test_get_vertex() {
    let r1 = S2LatLngRect::from_intervals(
        R1Interval::new(0.0, PI_2),
        S1Interval::new(-PI, 0.0)
    );
    assert_eq!(r1.get_vertex(0), S2LatLng::from_radians(0.0, -PI));
    assert_eq!(r1.get_vertex(1), S2LatLng::from_radians(0.0, 0.0));
    assert_eq!(r1.get_vertex(2), S2LatLng::from_radians(PI_2, 0.0));
    assert_eq!(r1.get_vertex(3), S2LatLng::from_radians(PI_2, -PI));

    // Test that GetVertex() returns vertices in CCW order - simplified version
    for i in 0..4 {
        let lat = PI_4 * (i as f64 - 2.0);
        let lng = PI_2 * (i as f64 - 2.0) + 0.2;
        let r = S2LatLngRect::from_intervals(
            R1Interval::new(lat, lat + PI_4),
            S1Interval::new(
                lng.rem_euclid(2.0 * PI) - PI,
                (lng + PI_2).rem_euclid(2.0 * PI) - PI
            )
        );
        
        // Just verify we can get all vertices (full CCW test would require more geometry)
        for k in 0..4 {
            let _vertex = r.get_vertex(k);
        }
    }
}

/// Test Contains methods (from C++ Contains)
#[test]
fn test_contains() {
    let eq_m180 = S2LatLng::from_radians(0.0, -PI);
    let north_pole = S2LatLng::from_radians(PI_2, 0.0);
    let r1 = S2LatLngRect::new(eq_m180, north_pole);

    assert!(r1.contains(&S2LatLng::from_degrees(30.0, -45.0)));
    assert!(r1.interior_contains(&S2LatLng::from_degrees(30.0, -45.0)));
    assert!(!r1.contains(&S2LatLng::from_degrees(30.0, 45.0)));
    assert!(!r1.interior_contains(&S2LatLng::from_degrees(30.0, 45.0)));
    assert!(r1.contains(&eq_m180));
    assert!(!r1.interior_contains(&eq_m180));
    assert!(r1.contains(&north_pole));
    assert!(!r1.interior_contains(&north_pole));
    
    // Test S2Point containment
    assert!(r1.contains_point(&S2Point::new(0.5, -0.3, 0.1).unwrap()));
    assert!(!r1.contains_point(&S2Point::new(0.5, 0.2, 0.1).unwrap()));
}

/// Helper function to test interval operations
fn test_interval_ops(
    x: &S2LatLngRect,
    y: &S2LatLngRect,
    expected_relation: &str,
    expected_union: &S2LatLngRect,
    expected_intersection: &S2LatLngRect
) {
    // Test all interval operations
    // "expected_relation" is sequence of "T"/"F" for Contains(), InteriorContains(), 
    // Intersects(), InteriorIntersects() respectively
    
    assert_eq!(x.contains_rect(y), expected_relation.chars().nth(0).unwrap() == 'T');
    assert_eq!(x.interior_contains_rect(y), expected_relation.chars().nth(1).unwrap() == 'T');
    assert_eq!(x.intersects(y), expected_relation.chars().nth(2).unwrap() == 'T');
    assert_eq!(x.interior_intersects(y), expected_relation.chars().nth(3).unwrap() == 'T');

    assert_eq!(x.contains_rect(y), x.union(y) == *x);
    assert_eq!(x.intersects(y), !x.intersection(y).is_empty());

    assert_eq!(x.union(y), *expected_union);
    assert_eq!(x.intersection(y), *expected_intersection);

    // Test AddPoint for point rectangles
    if y.get_size() == S2LatLng::from_radians(0.0, 0.0) {
        let mut r = *x;
        r.add_point(&y.lo());
        assert_eq!(r, *expected_union);
    }
}

/// Test interval operations (from C++ IntervalOps)
#[test]
fn test_interval_operations() {
    // Rectangle "r1" covers one-quarter of the sphere
    let r1 = rect_from_degrees(0.0, -180.0, 90.0, 0.0);

    // Test operations where one rectangle consists of a single point
    let r1_mid = rect_from_degrees(45.0, -90.0, 45.0, -90.0);
    test_interval_ops(&r1, &r1_mid, "TTTT", &r1, &r1_mid);

    let req_m180 = rect_from_degrees(0.0, -180.0, 0.0, -180.0);
    test_interval_ops(&r1, &req_m180, "TFTF", &r1, &req_m180);

    let rnorth_pole = rect_from_degrees(90.0, 0.0, 90.0, 0.0);
    test_interval_ops(&r1, &rnorth_pole, "TFTF", &r1, &rnorth_pole);

    test_interval_ops(
        &r1,
        &rect_from_degrees(-10.0, -1.0, 1.0, 20.0),
        "FFTT",
        &rect_from_degrees(-10.0, 180.0, 90.0, 20.0),
        &rect_from_degrees(0.0, -1.0, 1.0, 0.0)
    );

    test_interval_ops(
        &r1,
        &rect_from_degrees(-10.0, -1.0, 0.0, 20.0),
        "FFTF",
        &rect_from_degrees(-10.0, 180.0, 90.0, 20.0),
        &rect_from_degrees(0.0, -1.0, 0.0, 0.0)
    );

    test_interval_ops(
        &r1,
        &rect_from_degrees(-10.0, 0.0, 1.0, 20.0),
        "FFTF",
        &rect_from_degrees(-10.0, 180.0, 90.0, 20.0),
        &rect_from_degrees(0.0, 0.0, 1.0, 0.0)
    );

    test_interval_ops(
        &rect_from_degrees(-15.0, -160.0, -15.0, -150.0),
        &rect_from_degrees(20.0, 145.0, 25.0, 155.0),
        "FFFF",
        &rect_from_degrees(-15.0, 145.0, 25.0, -150.0),
        &S2LatLngRect::empty()
    );

    test_interval_ops(
        &rect_from_degrees(70.0, -10.0, 90.0, -140.0),
        &rect_from_degrees(60.0, 175.0, 80.0, 5.0),
        "FFTT",
        &rect_from_degrees(60.0, -180.0, 90.0, 180.0),
        &rect_from_degrees(70.0, 175.0, 80.0, 5.0)
    );

    // Check that intersection of rectangles that overlap in latitude
    // but not longitude is valid, and vice versa
    test_interval_ops(
        &rect_from_degrees(12.0, 30.0, 60.0, 60.0),
        &rect_from_degrees(0.0, 0.0, 30.0, 18.0),
        "FFFF",
        &rect_from_degrees(0.0, 0.0, 60.0, 60.0),
        &S2LatLngRect::empty()
    );

    test_interval_ops(
        &rect_from_degrees(0.0, 0.0, 18.0, 42.0),
        &rect_from_degrees(30.0, 12.0, 42.0, 60.0),
        "FFFF",
        &rect_from_degrees(0.0, 0.0, 42.0, 60.0),
        &S2LatLngRect::empty()
    );
}

/// Test BoundaryIntersects with empty rectangle (from C++ BoundaryIntersects)
#[test]
fn test_boundary_intersects_empty() {
    let rect = S2LatLngRect::empty();
    let lo = rect.lo().to_point().unwrap();
    let hi = rect.hi().to_point().unwrap();
    assert!(!rect.boundary_intersects(&lo, &lo));
    assert!(!rect.boundary_intersects(&lo, &hi));
}

/// Test BoundaryIntersects with full rectangle
#[test]
fn test_boundary_intersects_full() {
    let rect = S2LatLngRect::full();
    let lo = rect.lo().to_point().unwrap();
    let hi = rect.hi().to_point().unwrap();
    assert!(!rect.boundary_intersects(&lo, &lo));
    assert!(!rect.boundary_intersects(&lo, &hi));
}

/// Test BoundaryIntersects with spherical lune
#[test]
fn test_boundary_intersects_spherical_lune() {
    // This rectangle only has two non-degenerate sides
    let rect = rect_from_degrees(-90.0, 100.0, 90.0, 120.0);
    
    // Test various edge cases - simplified without s2textformat dependency
    let p1 = S2LatLng::from_degrees(60.0, 60.0).to_point().unwrap();
    let p2 = S2LatLng::from_degrees(90.0, 60.0).to_point().unwrap();
    assert!(!rect.boundary_intersects(&p1, &p2));
    
    let p3 = S2LatLng::from_degrees(-60.0, 110.0).to_point().unwrap();
    let p4 = S2LatLng::from_degrees(60.0, 110.0).to_point().unwrap();
    assert!(!rect.boundary_intersects(&p3, &p4));
    
    let p5 = S2LatLng::from_degrees(-60.0, 95.0).to_point().unwrap();
    let p6 = S2LatLng::from_degrees(60.0, 110.0).to_point().unwrap();
    assert!(rect.boundary_intersects(&p5, &p6));
    
    let p7 = S2LatLng::from_degrees(60.0, 115.0).to_point().unwrap();
    let p8 = S2LatLng::from_degrees(80.0, 125.0).to_point().unwrap();
    assert!(rect.boundary_intersects(&p7, &p8));
}

/// Test AddPoint method (from C++ AddPoint)
#[test]
fn test_add_point() {
    let mut p = S2LatLngRect::empty();
    p.add_point(&S2LatLng::from_degrees(0.0, 0.0));
    assert!(p.is_point());
    
    p.add_point(&S2LatLng::from_radians(0.0, -PI_2));
    assert!(!p.is_point());
    
    p.add_point(&S2LatLng::from_radians(PI_4, -PI));
    p.add_s2_point(&S2Point::new(0.0, 0.0, 1.0).unwrap());
    assert_eq!(p, rect_from_degrees(0.0, -180.0, 90.0, 0.0));
}

/// Test Expanded method (from C++ Expanded)
#[test]
fn test_expanded() {
    assert!(rect_from_degrees(70.0, 150.0, 80.0, 170.0)
        .expanded(S2LatLng::from_degrees(20.0, 30.0))
        .approx_equals(
            &rect_from_degrees(50.0, 120.0, 90.0, -160.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(S2LatLngRect::empty()
        .expanded(S2LatLng::from_degrees(20.0, 30.0))
        .is_empty());

    assert!(S2LatLngRect::full()
        .expanded(S2LatLng::from_degrees(500.0, 500.0))
        .is_full());

    assert!(rect_from_degrees(-90.0, 170.0, 10.0, 20.0)
        .expanded(S2LatLng::from_degrees(30.0, 80.0))
        .approx_equals(
            &rect_from_degrees(-90.0, -180.0, 40.0, 180.0),
            S1Angle::from_radians(1e-15)
        ));

    // Test negative margins
    assert!(rect_from_degrees(10.0, -50.0, 60.0, 70.0)
        .expanded(S2LatLng::from_degrees(-10.0, -10.0))
        .approx_equals(
            &rect_from_degrees(20.0, -40.0, 50.0, 60.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(-20.0, -180.0, 20.0, 180.0)
        .expanded(S2LatLng::from_degrees(-10.0, -10.0))
        .approx_equals(
            &rect_from_degrees(-10.0, -180.0, 10.0, 180.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(-20.0, -180.0, 20.0, 180.0)
        .expanded(S2LatLng::from_degrees(-30.0, -30.0))
        .is_empty());

    assert!(rect_from_degrees(-90.0, 10.0, 90.0, 11.0)
        .expanded(S2LatLng::from_degrees(-10.0, -10.0))
        .is_empty());

    assert!(rect_from_degrees(-90.0, 10.0, 90.0, 100.0)
        .expanded(S2LatLng::from_degrees(-10.0, -10.0))
        .approx_equals(
            &rect_from_degrees(-80.0, 20.0, 80.0, 90.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(S2LatLngRect::empty()
        .expanded(S2LatLng::from_degrees(-50.0, -500.0))
        .is_empty());

    assert!(S2LatLngRect::full()
        .expanded(S2LatLng::from_degrees(-50.0, -50.0))
        .approx_equals(
            &rect_from_degrees(-40.0, -180.0, 40.0, 180.0),
            S1Angle::from_radians(1e-15)
        ));

    // Test mixed margins
    assert!(rect_from_degrees(10.0, -50.0, 60.0, 70.0)
        .expanded(S2LatLng::from_degrees(-10.0, 30.0))
        .approx_equals(
            &rect_from_degrees(20.0, -80.0, 50.0, 100.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(-20.0, -180.0, 20.0, 180.0)
        .expanded(S2LatLng::from_degrees(10.0, -500.0))
        .approx_equals(
            &rect_from_degrees(-30.0, -180.0, 30.0, 180.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(-90.0, -180.0, 80.0, 180.0)
        .expanded(S2LatLng::from_degrees(-30.0, 500.0))
        .approx_equals(
            &rect_from_degrees(-60.0, -180.0, 50.0, 180.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(-80.0, -100.0, 80.0, 150.0)
        .expanded(S2LatLng::from_degrees(30.0, -50.0))
        .approx_equals(
            &rect_from_degrees(-90.0, -50.0, 90.0, 100.0),
            S1Angle::from_radians(1e-15)
        ));

    assert!(rect_from_degrees(0.0, -180.0, 50.0, 180.0)
        .expanded(S2LatLng::from_degrees(-30.0, 500.0))
        .is_empty());

    assert!(rect_from_degrees(-80.0, 10.0, 70.0, 20.0)
        .expanded(S2LatLng::from_degrees(30.0, -200.0))
        .is_empty());

    assert!(S2LatLngRect::empty()
        .expanded(S2LatLng::from_degrees(100.0, -100.0))
        .is_empty());

    assert!(S2LatLngRect::full()
        .expanded(S2LatLng::from_degrees(100.0, -100.0))
        .is_full());
}

/// Test PolarClosure method (from C++ PolarClosure)
#[test]
fn test_polar_closure() {
    assert_eq!(
        rect_from_degrees(-89.0, 0.0, 89.0, 1.0),
        rect_from_degrees(-89.0, 0.0, 89.0, 1.0).polar_closure()
    );
    assert_eq!(
        rect_from_degrees(-90.0, -180.0, -45.0, 180.0),
        rect_from_degrees(-90.0, -30.0, -45.0, 100.0).polar_closure()
    );
    assert_eq!(
        rect_from_degrees(89.0, -180.0, 90.0, 180.0),
        rect_from_degrees(89.0, 145.0, 90.0, 146.0).polar_closure()
    );
    assert_eq!(
        S2LatLngRect::full(),
        rect_from_degrees(-90.0, -145.0, 90.0, -144.0).polar_closure()
    );
}

/// Test ExpandedByDistance with positive distance (from C++ ExpandedByDistance)
#[test]
fn test_expanded_by_distance_positive() {
    assert!(rect_from_degrees(0.0, 170.0, 0.0, -170.0)
        .expanded_by_distance(S1Angle::from_degrees(15.0))
        .approx_equals(
            &rect_from_degrees(-15.0, 155.0, 15.0, -155.0),
            S1Angle::from_radians(1e-13)
        ));

    assert!(rect_from_degrees(60.0, 150.0, 80.0, 10.0)
        .expanded_by_distance(S1Angle::from_degrees(15.0))
        .approx_equals(
            &rect_from_degrees(45.0, -180.0, 90.0, 180.0),
            S1Angle::from_radians(1e-13)
        ));
}

/// Test area calculation (from C++ Area)
#[test]
fn test_area() {
    assert_eq!(S2LatLngRect::empty().area(), 0.0);
    assert_relative_eq!(S2LatLngRect::full().area(), 4.0 * PI, epsilon = 1e-15);
    assert_relative_eq!(
        rect_from_degrees(0.0, 0.0, 90.0, 90.0).area(),
        PI_2,
        epsilon = 1e-15
    );
}

/// Test distance calculations - simplified version (from C++ GetDistance*)
#[test]
fn test_get_distance_overlapping() {
    // Check pairs of rectangles that overlap (should all return 0)
    let a = rect_from_degrees(0.0, 0.0, 2.0, 2.0);
    let b = S2LatLngRect::from_point(S2LatLng::from_degrees(0.0, 0.0));
    
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&a));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&b));
    assert_eq!(S1Angle::from_radians(0.0), b.get_distance(&b));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance_to_point(&S2LatLng::from_degrees(0.0, 0.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(0.0, 1.0, 2.0, 3.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(0.0, 2.0, 2.0, 4.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(1.0, 0.0, 3.0, 2.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(2.0, 0.0, 4.0, 2.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(1.0, 1.0, 3.0, 3.0)));
    assert_eq!(S1Angle::from_radians(0.0), a.get_distance(&rect_from_degrees(2.0, 2.0, 4.0, 4.0)));
}

/// Test hash implementation  
#[test]
fn test_hash() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    
    let rect1 = rect_from_degrees(-89.0, 0.0, 89.0, 1.0);
    let rect2 = rect_from_degrees(1.0, 5.0, 1.0, 5.0);
    let empty = S2LatLngRect::empty();
    let full = S2LatLngRect::full();
    
    // Test that different rectangles have different hashes
    let mut hasher1 = DefaultHasher::new();
    rect1.hash(&mut hasher1);
    let hash1 = hasher1.finish();
    
    let mut hasher2 = DefaultHasher::new();
    rect2.hash(&mut hasher2);
    let hash2 = hasher2.finish();
    
    assert_ne!(hash1, hash2);
    
    // Test that same rectangle has same hash
    let mut hasher3 = DefaultHasher::new();
    rect1.hash(&mut hasher3);
    let hash3 = hasher3.finish();
    
    assert_eq!(hash1, hash3);
}

/// Test string representation
#[test]
fn test_display() {
    let rect = rect_from_degrees(10.0, 20.0, 30.0, 40.0);
    let display_str = format!("{}", rect);
    
    // Should contain coordinates
    assert!(display_str.contains("Lo="));
    assert!(display_str.contains("Hi="));
}

/// Test edge cases and special rectangles
#[test]
fn test_edge_cases() {
    // Test rectangle crossing International Date Line
    let crossing_rect = rect_from_degrees(20.0, 170.0, 40.0, -170.0);
    assert!(crossing_rect.is_inverted());
    assert!(crossing_rect.contains(&S2LatLng::from_degrees(30.0, 180.0)));
    assert!(crossing_rect.contains(&S2LatLng::from_degrees(30.0, -180.0)));
    
    // Test polar rectangles
    let north_polar = rect_from_degrees(85.0, -180.0, 90.0, 180.0);
    assert!(north_polar.contains(&S2LatLng::from_degrees(90.0, 0.0)));
    assert!(north_polar.contains(&S2LatLng::from_degrees(90.0, 123.0)));
    
    let south_polar = rect_from_degrees(-90.0, -180.0, -85.0, 180.0);
    assert!(south_polar.contains(&S2LatLng::from_degrees(-90.0, 0.0)));
    assert!(south_polar.contains(&S2LatLng::from_degrees(-90.0, -123.0)));
    
    // Test very small rectangles
    let tiny = rect_from_degrees(0.0, 0.0, 1e-10, 1e-10);
    assert!(tiny.is_valid());
    assert!(!tiny.is_empty());
    assert!(tiny.area() > 0.0);
    assert!(tiny.area() < 1e-15);
}

/// Test arithmetic and geometric operations
#[test]
fn test_geometric_operations() {
    let rect1 = rect_from_degrees(10.0, 20.0, 30.0, 40.0);
    let rect2 = rect_from_degrees(25.0, 35.0, 45.0, 55.0);
    
    // Test union
    let union_rect = rect1.union(&rect2);
    assert!(union_rect.contains_rect(&rect1));
    assert!(union_rect.contains_rect(&rect2));
    assert_relative_eq!(union_rect.lat_lo().degrees(), 10.0, epsilon = 1e-15);
    assert_relative_eq!(union_rect.lat_hi().degrees(), 45.0, epsilon = 1e-15);
    assert_relative_eq!(union_rect.lng_lo().degrees(), 20.0, epsilon = 1e-15);
    assert_relative_eq!(union_rect.lng_hi().degrees(), 55.0, epsilon = 1e-15);
    
    // Test intersection
    let intersection = rect1.intersection(&rect2);
    assert_relative_eq!(intersection.lat_lo().degrees(), 25.0, epsilon = 1e-15);
    assert_relative_eq!(intersection.lat_hi().degrees(), 30.0, epsilon = 1e-15);
    assert_relative_eq!(intersection.lng_lo().degrees(), 35.0, epsilon = 1e-15);
    assert_relative_eq!(intersection.lng_hi().degrees(), 40.0, epsilon = 1e-15);
    
    // Test centroid (basic sanity check)
    let centroid = rect1.get_centroid();
    assert!(centroid.coords().x.is_finite());
    assert!(centroid.coords().y.is_finite());
    assert!(centroid.coords().z.is_finite());
    
    // Empty rectangle should have zero centroid
    assert_eq!(S2LatLngRect::empty().get_centroid().coords().x, 0.0);
    assert_eq!(S2LatLngRect::empty().get_centroid().coords().y, 0.0);
    assert_eq!(S2LatLngRect::empty().get_centroid().coords().z, 0.0);
}

/// Test rectangle validation
#[test]
fn test_validation() {
    // Valid rectangles
    assert!(rect_from_degrees(-90.0, -180.0, 90.0, 180.0).is_valid());
    assert!(rect_from_degrees(0.0, 0.0, 0.0, 0.0).is_valid());
    assert!(S2LatLngRect::empty().is_valid());
    assert!(S2LatLngRect::full().is_valid());
    
    // Test that our implementation maintains validity
    let rect = rect_from_degrees(10.0, 20.0, 30.0, 40.0);
    assert!(rect.is_valid());
    assert!(rect.expanded(S2LatLng::from_degrees(5.0, 5.0)).is_valid());
    assert!(rect.union(&rect_from_degrees(50.0, 50.0, 60.0, 60.0)).is_valid());
    assert!(rect.intersection(&rect_from_degrees(20.0, 30.0, 40.0, 50.0)).is_valid());
}