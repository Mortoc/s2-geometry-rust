//! Port of r1interval_test.cc - Complete R1Interval tests
//!
//! This module ports all tests from C++ r1interval_test.cc to validate
//! our R1Interval implementation matches C++ behavior exactly.
//!
//! R1Interval represents closed intervals on the real line, supporting:
//! - Empty intervals and zero-length intervals (points)
//! - Linear interval operations (no wraparound like S1Interval)
//! - Interval operations (union, intersection, contains)
//! - Expansion and projection operations
//! - Edge cases and boundary conditions

use s2geometry_rust::R1Interval;
use std::f64::{EPSILON as DBL_EPSILON, INFINITY};

/// Test helper function that tests all interval operations
/// Matches the C++ TestIntervalOps function exactly
fn test_interval_ops(x: &R1Interval, y: &R1Interval, expected: &str) {
    // Test all of the interval operations on the given pair of intervals.
    // "expected" is a sequence of "T" and "F" characters corresponding to
    // the expected results of Contains(), InteriorContains(), Intersects(),
    // and InteriorIntersects() respectively.

    assert_eq!(expected.chars().nth(0).unwrap() == 'T', x.contains_interval(y));
    assert_eq!(expected.chars().nth(1).unwrap() == 'T', x.interior_contains_interval(y));
    assert_eq!(expected.chars().nth(2).unwrap() == 'T', x.intersects(y));
    assert_eq!(expected.chars().nth(3).unwrap() == 'T', x.interior_intersects(y));

    assert_eq!(x.contains_interval(y), x.union(y) == *x);
    assert_eq!(x.intersects(y), !x.intersection(y).is_empty());

    let mut z = *x;
    z.add_interval(y);
    assert_eq!(x.union(y), z);
}

#[test]
fn test_basic() {
    // Constructors and accessors.
    let unit = R1Interval::new(0.0, 1.0);
    let negunit = R1Interval::new(-1.0, 0.0);
    assert_eq!(0.0, unit.lo());
    assert_eq!(1.0, unit.hi());
    assert_eq!(-1.0, negunit[0]);
    assert_eq!(0.0, negunit[1]);
    
    let mut ten = R1Interval::new(0.0, 0.0);
    ten.set_hi(10.0);
    assert_eq!(R1Interval::new(0.0, 10.0), ten);
    ten[0] = -10.0;
    assert_eq!(R1Interval::new(-10.0, 10.0), ten);
    ten[1] = 0.0;
    assert_eq!([-10.0, 0.0], ten.bounds());
    *ten.bounds_mut() = [0.0, 10.0];
    assert_eq!(R1Interval::new(0.0, 10.0), ten);

    // is_empty()
    let half = R1Interval::new(0.5, 0.5);
    assert!(!unit.is_empty());
    assert!(!half.is_empty());
    let empty = R1Interval::empty();
    assert!(empty.is_empty());

    // == and !=
    assert_eq!(empty, empty);
    assert_eq!(unit, unit);
    assert_ne!(unit, empty);
    assert_ne!(R1Interval::new(1.0, 2.0), R1Interval::new(1.0, 3.0));

    // Check that the default R1Interval is identical to Empty().
    let default_empty = R1Interval::default();
    assert!(default_empty.is_empty());
    assert_eq!(empty.lo(), default_empty.lo());
    assert_eq!(empty.hi(), default_empty.hi());

    // GetCenter(), GetLength()
    assert_eq!(unit.get_center(), 0.5);
    assert_eq!(half.get_center(), 0.5);
    assert_eq!(negunit.get_length(), 1.0);
    assert_eq!(half.get_length(), 0.0);
    assert!(empty.get_length() < 0.0);

    // Contains(double), InteriorContains(double)
    assert!(unit.contains(0.5));
    assert!(unit.interior_contains(0.5));
    assert!(unit.contains(0.0));
    assert!(!unit.interior_contains(0.0));
    assert!(unit.contains(1.0));
    assert!(!unit.interior_contains(1.0));

    // Contains(R1Interval), InteriorContains(R1Interval)
    // Intersects(R1Interval), InteriorIntersects(R1Interval)
    test_interval_ops(&empty, &empty, "TTFF");
    test_interval_ops(&empty, &unit, "FFFF");
    test_interval_ops(&unit, &half, "TTTT");
    test_interval_ops(&unit, &unit, "TFTT");
    test_interval_ops(&unit, &empty, "TTFF");
    test_interval_ops(&unit, &negunit, "FFTF");
    test_interval_ops(&unit, &R1Interval::new(0.0, 0.5), "TFTT");
    test_interval_ops(&half, &R1Interval::new(0.0, 0.5), "FFTF");

    // AddPoint()
    let mut r = empty;
    r.add_point(5.0);
    assert_eq!(5.0, r.lo());
    assert_eq!(5.0, r.hi());
    r.add_point(-1.0);
    assert_eq!(-1.0, r.lo());
    assert_eq!(5.0, r.hi());
    r.add_point(0.0);
    assert_eq!(-1.0, r.lo());
    assert_eq!(5.0, r.hi());

    // Project()
    assert_eq!(0.3, R1Interval::new(0.1, 0.4).project(0.3));
    assert_eq!(0.1, R1Interval::new(0.1, 0.4).project(-7.0));
    assert_eq!(0.4, R1Interval::new(0.1, 0.4).project(0.6));

    // FromPointPair()
    assert_eq!(R1Interval::new(4.0, 4.0), R1Interval::from_point_pair(4.0, 4.0));
    assert_eq!(R1Interval::new(-2.0, -1.0), R1Interval::from_point_pair(-1.0, -2.0));
    assert_eq!(R1Interval::new(-5.0, 3.0), R1Interval::from_point_pair(-5.0, 3.0));

    // Expanded()
    assert_eq!(empty, empty.expanded(0.45));
    assert_eq!(R1Interval::new(-0.5, 1.5), unit.expanded(0.5));
    assert_eq!(R1Interval::new(0.5, 0.5), unit.expanded(-0.5));
    assert!(unit.expanded(-0.51).is_empty());
    assert!(unit.expanded(-0.51).expanded(0.51).is_empty());

    // Union(), Intersection()
    assert_eq!(R1Interval::new(99.0, 100.0), R1Interval::new(99.0, 100.0).union(&empty));
    assert_eq!(R1Interval::new(99.0, 100.0), empty.union(&R1Interval::new(99.0, 100.0)));
    assert!(R1Interval::new(5.0, 3.0).union(&R1Interval::new(0.0, -2.0)).is_empty());
    assert!(R1Interval::new(0.0, -2.0).union(&R1Interval::new(5.0, 3.0)).is_empty());
    assert_eq!(unit, unit.union(&unit));
    assert_eq!(R1Interval::new(-1.0, 1.0), unit.union(&negunit));
    assert_eq!(R1Interval::new(-1.0, 1.0), negunit.union(&unit));
    assert_eq!(unit, half.union(&unit));
    assert_eq!(half, unit.intersection(&half));
    assert_eq!(R1Interval::new(0.0, 0.0), unit.intersection(&negunit));
    assert!(negunit.intersection(&half).is_empty());
    assert!(unit.intersection(&empty).is_empty());
    assert!(empty.intersection(&unit).is_empty());
}

#[test]
fn test_approx_equals() {
    // Choose two values kLo and kHi such that it's okay to shift an endpoint by
    // kLo (i.e., the resulting interval is equivalent) but not by kHi.
    const K_LO: f64 = 4.0 * DBL_EPSILON;  // < max_error default
    const K_HI: f64 = 6.0 * DBL_EPSILON;  // > max_error default

    // Empty intervals.
    let empty = R1Interval::empty();
    assert!(empty.approx_equals(&empty, 1e-15));
    assert!(R1Interval::new(0.0, 0.0).approx_equals(&empty, 1e-15));
    assert!(empty.approx_equals(&R1Interval::new(0.0, 0.0), 1e-15));
    assert!(R1Interval::new(1.0, 1.0).approx_equals(&empty, 1e-15));
    assert!(empty.approx_equals(&R1Interval::new(1.0, 1.0), 1e-15));
    assert!(!empty.approx_equals(&R1Interval::new(0.0, 1.0), 1e-15));
    assert!(empty.approx_equals(&R1Interval::new(1.0, 1.0 + 2.0*K_LO), 1e-15));
    assert!(!empty.approx_equals(&R1Interval::new(1.0, 1.0 + 2.0*K_HI), 1e-15));

    // Singleton intervals.
    assert!(R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0, 1.0), 1e-15));
    assert!(R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0 - K_LO, 1.0 - K_LO), 1e-15));
    assert!(R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0 + K_LO, 1.0 + K_LO), 1e-15));
    assert!(!R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0 - K_HI, 1.0), 1e-15));
    assert!(!R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0, 1.0 + K_HI), 1e-15));
    assert!(R1Interval::new(1.0, 1.0).approx_equals(&R1Interval::new(1.0 - K_LO, 1.0 + K_LO), 1e-15));
    assert!(!R1Interval::new(0.0, 0.0).approx_equals(&R1Interval::new(1.0, 1.0), 1e-15));

    // Other intervals.
    assert!(R1Interval::new(1.0 - K_LO, 2.0 + K_LO).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
    assert!(R1Interval::new(1.0 + K_LO, 2.0 - K_LO).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
    assert!(!R1Interval::new(1.0 - K_HI, 2.0 + K_LO).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
    assert!(!R1Interval::new(1.0 + K_HI, 2.0 - K_LO).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
    assert!(!R1Interval::new(1.0 - K_LO, 2.0 + K_HI).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
    assert!(!R1Interval::new(1.0 + K_LO, 2.0 - K_HI).approx_equals(&R1Interval::new(1.0, 2.0), 1e-15));
}

#[test] 
fn test_constructors() {
    // Test various constructor methods
    let empty = R1Interval::empty();
    assert!(empty.is_empty());
    
    let point = R1Interval::from_point(5.0);
    assert_eq!(point.lo(), 5.0);
    assert_eq!(point.hi(), 5.0);
    assert!(!point.is_empty());
    assert_eq!(point.get_length(), 0.0);
    
    let pair1 = R1Interval::from_point_pair(1.0, 3.0);
    assert_eq!(pair1, R1Interval::new(1.0, 3.0));
    
    let pair2 = R1Interval::from_point_pair(3.0, 1.0);
    assert_eq!(pair2, R1Interval::new(1.0, 3.0));
}

#[test]
fn test_hausdorff_distance() {
    let empty = R1Interval::empty();
    let unit = R1Interval::new(0.0, 1.0);
    let half = R1Interval::new(0.5, 0.5);
    let neg = R1Interval::new(-1.0, 0.0);
    let far = R1Interval::new(10.0, 20.0);
    
    // Empty interval cases
    assert_eq!(empty.get_directed_hausdorff_distance(&empty), 0.0);
    assert_eq!(empty.get_directed_hausdorff_distance(&unit), 0.0);
    assert_eq!(unit.get_directed_hausdorff_distance(&empty), INFINITY);
    
    // Same intervals
    assert_eq!(unit.get_directed_hausdorff_distance(&unit), 0.0);
    assert_eq!(half.get_directed_hausdorff_distance(&half), 0.0);
    
    // Contained intervals  
    assert_eq!(unit.get_directed_hausdorff_distance(&half), 0.5);
    assert_eq!(half.get_directed_hausdorff_distance(&unit), 0.0);
    
    // Overlapping intervals
    assert_eq!(unit.get_directed_hausdorff_distance(&neg), 1.0);
    assert_eq!(neg.get_directed_hausdorff_distance(&unit), 1.0);
    
    // Disjoint intervals
    assert_eq!(unit.get_directed_hausdorff_distance(&far), 10.0);
    assert_eq!(far.get_directed_hausdorff_distance(&unit), 19.0);
}

#[test]
fn test_add_operations() {
    let mut interval = R1Interval::empty();
    
    // Adding to empty interval
    interval.add_point(5.0);
    assert_eq!(interval, R1Interval::from_point(5.0));
    
    // Expanding interval
    interval.add_point(3.0);
    assert_eq!(interval, R1Interval::new(3.0, 5.0));
    
    interval.add_point(7.0);
    assert_eq!(interval, R1Interval::new(3.0, 7.0));
    
    // Adding point inside interval
    interval.add_point(4.0);
    assert_eq!(interval, R1Interval::new(3.0, 7.0));
    
    // Test add_interval
    let mut base = R1Interval::new(0.0, 2.0);
    let other = R1Interval::new(1.0, 3.0);
    base.add_interval(&other);
    assert_eq!(base, R1Interval::new(0.0, 3.0));
    
    // Adding empty interval should not change anything
    base.add_interval(&R1Interval::empty());
    assert_eq!(base, R1Interval::new(0.0, 3.0));
    
    // Adding to empty interval
    let mut empty = R1Interval::empty();
    empty.add_interval(&other);
    assert_eq!(empty, other);
}

#[test]
fn test_projection() {
    let interval = R1Interval::new(2.0, 8.0);
    
    // Point inside interval
    assert_eq!(interval.project(5.0), 5.0);
    
    // Point below interval
    assert_eq!(interval.project(0.0), 2.0);
    
    // Point above interval
    assert_eq!(interval.project(10.0), 8.0);
    
    // Boundary points
    assert_eq!(interval.project(2.0), 2.0);
    assert_eq!(interval.project(8.0), 8.0);
    
    // Single point interval
    let point = R1Interval::from_point(5.0);
    assert_eq!(point.project(3.0), 5.0);
    assert_eq!(point.project(7.0), 5.0);
    assert_eq!(point.project(5.0), 5.0);
}

#[test]
fn test_expansion() {
    let interval = R1Interval::new(2.0, 5.0);
    
    // Positive expansion
    assert_eq!(interval.expanded(1.0), R1Interval::new(1.0, 6.0));
    assert_eq!(interval.expanded(0.5), R1Interval::new(1.5, 5.5));
    
    // No expansion
    assert_eq!(interval.expanded(0.0), interval);
    
    // Negative expansion (shrinking)
    assert_eq!(interval.expanded(-0.5), R1Interval::new(2.5, 4.5));
    assert_eq!(interval.expanded(-1.0), R1Interval::new(3.0, 4.0));
    
    // Shrinking to point
    assert_eq!(interval.expanded(-1.5), R1Interval::new(3.5, 3.5));
    
    // Over-shrinking (becomes empty)
    assert!(interval.expanded(-2.0).is_empty());
    
    // Empty interval expansion
    let empty = R1Interval::empty();
    assert_eq!(empty.expanded(1.0), empty);
    assert_eq!(empty.expanded(-1.0), empty);
}

#[test]
fn test_index_access() {
    let mut interval = R1Interval::new(1.0, 3.0);
    
    // Read access
    assert_eq!(interval[0], 1.0);
    assert_eq!(interval[1], 3.0);
    
    // Write access
    interval[0] = 0.5;
    interval[1] = 3.5;
    assert_eq!(interval, R1Interval::new(0.5, 3.5));
}

#[test]
fn test_display() {
    let interval = R1Interval::new(1.5, 2.5);
    let display_str = format!("{}", interval);
    assert_eq!(display_str, "[1.5, 2.5]");
    
    let empty = R1Interval::empty();
    let empty_str = format!("{}", empty);
    assert_eq!(empty_str, "[1, 0]");
}

#[test]
fn test_edge_cases() {
    // Test with very small intervals
    let tiny = R1Interval::new(1.0, 1.0 + DBL_EPSILON);
    assert!(!tiny.is_empty());
    assert!(tiny.get_length() > 0.0);
    
    // Test with very large numbers
    let huge = R1Interval::new(-1e100, 1e100);
    assert!(!huge.is_empty());
    assert_eq!(huge.get_length(), 2e100);
    
    // Test with infinity
    let inf_interval = R1Interval::new(-INFINITY, INFINITY);
    assert!(!inf_interval.is_empty());
    assert!(inf_interval.get_length().is_infinite());
    assert!(inf_interval.contains(0.0));
    assert!(inf_interval.contains(1e100));
    assert!(inf_interval.contains(-1e100));
    
    // Test intervals with NaN (should be handled gracefully)
    let nan_interval = R1Interval::new(f64::NAN, 1.0);
    // NaN comparisons are always false, so NaN > 1.0 is false, meaning not empty
    // This matches IEEE 754 behavior
    assert!(!nan_interval.is_empty());
    
    // However, contains() with NaN will behave strangely
    assert!(!nan_interval.contains(0.5)); // NaN >= 0.5 is false
    assert!(!nan_interval.contains(f64::NAN)); // NaN >= NaN is false
}

#[test]
fn test_set_operations_comprehensive() {
    let a = R1Interval::new(1.0, 3.0);
    let b = R1Interval::new(2.0, 4.0);
    let c = R1Interval::new(5.0, 7.0);
    let empty = R1Interval::empty();
    
    // Union tests
    assert_eq!(a.union(&b), R1Interval::new(1.0, 4.0));
    assert_eq!(a.union(&c), R1Interval::new(1.0, 7.0)); // Disjoint intervals
    assert_eq!(a.union(&empty), a);
    assert_eq!(empty.union(&a), a);
    
    // Intersection tests
    assert_eq!(a.intersection(&b), R1Interval::new(2.0, 3.0));
    assert!(a.intersection(&c).is_empty()); // Disjoint intervals
    assert!(a.intersection(&empty).is_empty());
    assert!(empty.intersection(&a).is_empty());
    
    // Self operations
    assert_eq!(a.union(&a), a);
    assert_eq!(a.intersection(&a), a);
}