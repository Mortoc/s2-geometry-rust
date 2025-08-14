//! Port of s1interval_test.cc - Complete S1Interval tests
//!
//! This module ports all tests from C++ s1interval_test.cc to validate
//! our S1Interval implementation matches C++ behavior exactly.
//!
//! S1Interval represents closed intervals on the unit circle, supporting:
//! - Empty and full intervals
//! - Wraparound handling (circular intervals)
//! - Interval operations (union, intersection, contains)
//! - Complement operations
//! - Edge cases and special values

use s2geometry_rust::S1Interval;
use s2geometry_rust::math::constants::*;
use approx::{assert_relative_eq, assert_abs_diff_eq};

/// Test base struct providing standard intervals for testing
/// 
/// These intervals cover various quadrants on the unit circle:
/// - quad1 == [0, π/2]
/// - quad2 == [π/2, π] 
/// - quad3 == [-π, -π/2]
/// - quad4 == [-π/2, 0]
struct S1IntervalTestBase {
    // Basic intervals
    empty: S1Interval,
    full: S1Interval,
    
    // Single-point intervals
    zero: S1Interval,
    pi2: S1Interval,    // π/2
    pi: S1Interval,     // π
    mipi: S1Interval,   // -π (normalized to π)
    mipi2: S1Interval,  // -π/2
    
    // Single quadrants
    quad1: S1Interval,  // [0, π/2]
    quad2: S1Interval,  // [π/2, π] (inverted: [π/2, -π])
    quad3: S1Interval,  // [π, -π/2] (inverted)
    quad4: S1Interval,  // [-π/2, 0]
    
    // Quadrant pairs
    quad12: S1Interval, // [0, π] (inverted: [0, -π])
    quad23: S1Interval, // [π/2, -π/2] (inverted)
    quad34: S1Interval, // [π, 0] (inverted: [-π, 0])
    quad41: S1Interval, // [-π/2, π/2]
    
    // Quadrant triples
    quad123: S1Interval, // [0, -π/2] (inverted)
    quad234: S1Interval, // [π/2, 0] (inverted)
    quad341: S1Interval, // [π, π/2] (inverted)
    quad412: S1Interval, // [-π/2, π] (inverted: [-π/2, -π])
    
    // Small intervals around quadrant midpoints
    mid12: S1Interval,   // Around π/2
    mid23: S1Interval,   // Around π
    mid34: S1Interval,   // Around -π/2
    mid41: S1Interval,   // Around 0
}

impl S1IntervalTestBase {
    fn new() -> Self {
        S1IntervalTestBase {
            empty: S1Interval::empty(),
            full: S1Interval::full(),
            
            // Single-point intervals
            zero: S1Interval::new(0.0, 0.0),
            pi2: S1Interval::new(PI_2, PI_2),
            pi: S1Interval::new(PI, PI),
            mipi: S1Interval::new(-PI, -PI),  // Normalized to [π, π]
            mipi2: S1Interval::new(-PI_2, -PI_2),
            
            // Single quadrants
            quad1: S1Interval::new(0.0, PI_2),
            quad2: S1Interval::new(PI_2, -PI),
            quad3: S1Interval::new(PI, -PI_2),
            quad4: S1Interval::new(-PI_2, 0.0),
            
            // Quadrant pairs  
            quad12: S1Interval::new(0.0, -PI),
            quad23: S1Interval::new(PI_2, -PI_2),
            quad34: S1Interval::new(-PI, 0.0),
            quad41: S1Interval::new(-PI_2, PI_2),
            
            // Quadrant triples
            quad123: S1Interval::new(0.0, -PI_2),
            quad234: S1Interval::new(PI_2, 0.0),
            quad341: S1Interval::new(PI, PI_2),
            quad412: S1Interval::new(-PI_2, -PI),
            
            // Small intervals around midpoints
            mid12: S1Interval::new(PI_2 - 0.01, PI_2 + 0.02),
            mid23: S1Interval::new(PI - 0.01, -PI + 0.02),
            mid34: S1Interval::new(-PI_2 - 0.01, -PI_2 + 0.02),
            mid41: S1Interval::new(-0.01, 0.02),
        }
    }
}

/// Test constructors and accessors - matches ConstructorsAndAccessors
#[test]
fn test_constructors_and_accessors() {
    let base = S1IntervalTestBase::new();
    
    // Spot-check the constructors and accessors
    assert_eq!(base.quad12.lo(), 0.0);
    assert_eq!(base.quad12.hi(), PI);
    assert_eq!(base.quad34.bounds()[0], PI);
    assert_eq!(base.quad34.bounds()[1], 0.0);
    assert_eq!(base.pi.lo(), PI);
    assert_eq!(base.pi.hi(), PI);
    
    // Check that [-π, -π] is normalized to [π, π]
    assert_eq!(base.mipi.lo(), PI);
    assert_eq!(base.mipi.hi(), PI);
    assert_eq!(base.quad23.lo(), PI_2);
    assert_eq!(base.quad23.hi(), -PI_2);
    
    // Check that the default S1Interval is identical to Empty()
    let default_empty = S1Interval::default();
    assert!(default_empty.is_valid());
    assert!(default_empty.is_empty());
    assert_eq!(base.empty.lo(), default_empty.lo());
    assert_eq!(base.empty.hi(), default_empty.hi());
}

/// Test simple predicates - matches SimplePredicates
#[test]
fn test_simple_predicates() {
    let base = S1IntervalTestBase::new();
    
    // is_valid(), is_empty(), is_full(), is_inverted()
    assert!(base.zero.is_valid() && !base.zero.is_empty() && !base.zero.is_full());
    assert!(base.empty.is_valid() && base.empty.is_empty() && !base.empty.is_full());
    assert!(base.empty.is_inverted());
    assert!(base.full.is_valid() && !base.full.is_empty() && base.full.is_full());
    assert!(!base.quad12.is_empty() && !base.quad12.is_full() && !base.quad12.is_inverted());
    assert!(!base.quad23.is_empty() && !base.quad23.is_full() && base.quad23.is_inverted());
    assert!(base.pi.is_valid() && !base.pi.is_empty() && !base.pi.is_inverted());
    assert!(base.mipi.is_valid() && !base.mipi.is_empty() && !base.mipi.is_inverted());
}

/// Test almost empty or full intervals - matches AlmostEmptyOrFull
#[test]
fn test_almost_empty_or_full() {
    // Test that rounding errors don't cause intervals that are almost empty or
    // full to be considered empty or full
    let k_almost_pi = PI - 2.0 * f64::EPSILON;
    assert!(!S1Interval::new(-k_almost_pi, PI).is_full());
    assert!(!S1Interval::new(-PI, k_almost_pi).is_full());
    assert!(!S1Interval::new(PI, -k_almost_pi).is_empty());
    assert!(!S1Interval::new(k_almost_pi, -PI).is_empty());
}

/// Test get_center method - matches GetCenter
#[test]
fn test_get_center() {
    let base = S1IntervalTestBase::new();
    
    assert_eq!(base.quad12.get_center(), PI_2);
    assert_abs_diff_eq!(S1Interval::new(3.1, 2.9).get_center(), 3.0 - PI, epsilon = 1e-15);
    assert_abs_diff_eq!(S1Interval::new(-2.9, -3.1).get_center(), PI - 3.0, epsilon = 1e-15);
    assert_abs_diff_eq!(S1Interval::new(2.1, -2.1).get_center(), PI, epsilon = 1e-15);
    assert_eq!(base.pi.get_center(), PI);
    assert_eq!(base.mipi.get_center(), PI);
    assert_eq!(base.quad23.get_center().abs(), PI);
    assert_abs_diff_eq!(base.quad123.get_center(), 0.75 * PI, epsilon = 1e-15);
}

/// Test get_length method - matches GetLength
#[test]
fn test_get_length() {
    let base = S1IntervalTestBase::new();
    
    assert_eq!(base.quad12.get_length(), PI);
    assert_eq!(base.pi.get_length(), 0.0);
    assert_eq!(base.mipi.get_length(), 0.0);
    assert_abs_diff_eq!(base.quad123.get_length(), 1.5 * PI, epsilon = 1e-15);
    assert_eq!(base.quad23.get_length().abs(), PI);
    assert_eq!(base.full.get_length(), 2.0 * PI);
    assert!(base.empty.get_length() < 0.0);
}

/// Test complement operation - matches Complement
#[test]
fn test_complement() {
    let base = S1IntervalTestBase::new();
    
    assert!(base.empty.complement().is_full());
    assert!(base.full.complement().is_empty());
    assert!(base.pi.complement().is_full());
    assert!(base.mipi.complement().is_full());
    assert!(base.zero.complement().is_full());
    assert!(base.quad12.complement().approx_equals(&base.quad34, 1e-15));
    assert!(base.quad34.complement().approx_equals(&base.quad12, 1e-15));
    assert!(base.quad123.complement().approx_equals(&base.quad4, 1e-15));
}

/// Test contains point operations - matches Contains
#[test]
fn test_contains() {
    let base = S1IntervalTestBase::new();
    
    // Contains(double), InteriorContains(double)
    assert!(!base.empty.contains_point(0.0) && !base.empty.contains_point(PI) &&
            !base.empty.contains_point(-PI));
    assert!(!base.empty.interior_contains_point(PI) && !base.empty.interior_contains_point(-PI));
    assert!(base.full.contains_point(0.0) && base.full.contains_point(PI) && base.full.contains_point(-PI));
    assert!(base.full.interior_contains_point(PI) && base.full.interior_contains_point(-PI));
    assert!(base.quad12.contains_point(0.0) && base.quad12.contains_point(PI) &&
            base.quad12.contains_point(-PI));
    assert!(base.quad12.interior_contains_point(PI_2) && !base.quad12.interior_contains_point(0.0));
    assert!(!base.quad12.interior_contains_point(PI) &&
            !base.quad12.interior_contains_point(-PI));
    assert!(base.quad23.contains_point(PI_2) && base.quad23.contains_point(-PI_2));
    assert!(base.quad23.contains_point(PI) && base.quad23.contains_point(-PI));
    assert!(!base.quad23.contains_point(0.0));
    assert!(!base.quad23.interior_contains_point(PI_2) &&
            !base.quad23.interior_contains_point(-PI_2));
    assert!(base.quad23.interior_contains_point(PI) && base.quad23.interior_contains_point(-PI));
    assert!(!base.quad23.interior_contains_point(0.0));
    assert!(base.pi.contains_point(PI) && base.pi.contains_point(-PI) && !base.pi.contains_point(0.0));
    assert!(!base.pi.interior_contains_point(PI) && !base.pi.interior_contains_point(-PI));
    assert!(base.mipi.contains_point(PI) && base.mipi.contains_point(-PI) && !base.mipi.contains_point(0.0));
    assert!(!base.mipi.interior_contains_point(PI) && !base.mipi.interior_contains_point(-PI));
    assert!(base.zero.contains_point(0.0) && !base.zero.interior_contains_point(0.0));
}

/// Test interval operations helper function - matches C++ TestIntervalOps
fn test_interval_operations(
    x: &S1Interval,
    y: &S1Interval, 
    expected_relation: &str,
    expected_union: &S1Interval,
    expected_intersection: &S1Interval,
) {
    // Test all interval operations on the given pair of intervals
    // expected_relation is a sequence of "T" and "F" characters corresponding
    // to the expected results of Contains(), InteriorContains(), Intersects(),
    // and InteriorIntersects() respectively

    let relation_chars: Vec<char> = expected_relation.chars().collect();
    assert_eq!(x.contains(y), relation_chars[0] == 'T');
    assert_eq!(x.interior_contains(y), relation_chars[1] == 'T');
    assert_eq!(x.intersects(y), relation_chars[2] == 'T');
    assert_eq!(x.interior_intersects(y), relation_chars[3] == 'T');

    let actual_union = x.union(y);
    let actual_intersection = x.intersection(y);
    
    assert_eq!(actual_union.bounds(), expected_union.bounds());
    assert_eq!(actual_intersection.bounds(), expected_intersection.bounds());

    assert_eq!(x.contains(y), x.union(y) == *x);
    assert_eq!(x.intersects(y), !x.intersection(y).is_empty());

    if y.lo() == y.hi() {
        let mut r = *x;
        r.add_point(y.lo());
        assert_eq!(r.bounds(), expected_union.bounds());
    }
}

/// Test interval operations - matches IntervalOps
#[test]
fn test_interval_ops() {
    let base = S1IntervalTestBase::new();
    
    test_interval_operations(&base.empty, &base.empty, "TTFF", &base.empty, &base.empty);
    test_interval_operations(&base.empty, &base.full, "FFFF", &base.full, &base.empty);
    test_interval_operations(&base.empty, &base.zero, "FFFF", &base.zero, &base.empty);
    test_interval_operations(&base.empty, &base.pi, "FFFF", &base.pi, &base.empty);
    test_interval_operations(&base.empty, &base.mipi, "FFFF", &base.mipi, &base.empty);

    test_interval_operations(&base.full, &base.empty, "TTFF", &base.full, &base.empty);
    test_interval_operations(&base.full, &base.full, "TTTT", &base.full, &base.full);
    test_interval_operations(&base.full, &base.zero, "TTTT", &base.full, &base.zero);
    test_interval_operations(&base.full, &base.pi, "TTTT", &base.full, &base.pi);
    test_interval_operations(&base.full, &base.mipi, "TTTT", &base.full, &base.mipi);
    test_interval_operations(&base.full, &base.quad12, "TTTT", &base.full, &base.quad12);
    test_interval_operations(&base.full, &base.quad23, "TTTT", &base.full, &base.quad23);

    test_interval_operations(&base.zero, &base.empty, "TTFF", &base.zero, &base.empty);
    test_interval_operations(&base.zero, &base.full, "FFTF", &base.full, &base.zero);
    test_interval_operations(&base.zero, &base.zero, "TFTF", &base.zero, &base.zero);
    test_interval_operations(&base.zero, &base.pi, "FFFF", &S1Interval::new(0.0, PI), &base.empty);
    test_interval_operations(&base.zero, &base.pi2, "FFFF", &base.quad1, &base.empty);
    test_interval_operations(&base.zero, &base.mipi, "FFFF", &base.quad12, &base.empty);
    test_interval_operations(&base.zero, &base.mipi2, "FFFF", &base.quad4, &base.empty);
    test_interval_operations(&base.zero, &base.quad12, "FFTF", &base.quad12, &base.zero);
    test_interval_operations(&base.zero, &base.quad23, "FFFF", &base.quad123, &base.empty);

    test_interval_operations(&base.pi2, &base.empty, "TTFF", &base.pi2, &base.empty);
    test_interval_operations(&base.pi2, &base.full, "FFTF", &base.full, &base.pi2);
    test_interval_operations(&base.pi2, &base.zero, "FFFF", &base.quad1, &base.empty);
    test_interval_operations(&base.pi2, &base.pi, "FFFF", &S1Interval::new(PI_2, PI), &base.empty);
    test_interval_operations(&base.pi2, &base.pi2, "TFTF", &base.pi2, &base.pi2);
    test_interval_operations(&base.pi2, &base.mipi, "FFFF", &base.quad2, &base.empty);
    test_interval_operations(&base.pi2, &base.mipi2, "FFFF", &base.quad23, &base.empty);
    test_interval_operations(&base.pi2, &base.quad12, "FFTF", &base.quad12, &base.pi2);
    test_interval_operations(&base.pi2, &base.quad23, "FFTF", &base.quad23, &base.pi2);

    test_interval_operations(&base.pi, &base.empty, "TTFF", &base.pi, &base.empty);
    test_interval_operations(&base.pi, &base.full, "FFTF", &base.full, &base.pi);
    test_interval_operations(&base.pi, &base.zero, "FFFF", &S1Interval::new(PI, 0.0), &base.empty);
    test_interval_operations(&base.pi, &base.pi, "TFTF", &base.pi, &base.pi);
    test_interval_operations(&base.pi, &base.pi2, "FFFF", &S1Interval::new(PI_2, PI), &base.empty);
    test_interval_operations(&base.pi, &base.mipi, "TFTF", &base.pi, &base.pi);
    test_interval_operations(&base.pi, &base.mipi2, "FFFF", &base.quad3, &base.empty);
    test_interval_operations(&base.pi, &base.quad12, "FFTF", &S1Interval::new(0.0, PI), &base.pi);
    test_interval_operations(&base.pi, &base.quad23, "FFTF", &base.quad23, &base.pi);

    test_interval_operations(&base.mipi, &base.empty, "TTFF", &base.mipi, &base.empty);
    test_interval_operations(&base.mipi, &base.full, "FFTF", &base.full, &base.mipi);
    test_interval_operations(&base.mipi, &base.zero, "FFFF", &base.quad34, &base.empty);
    test_interval_operations(&base.mipi, &base.pi, "TFTF", &base.mipi, &base.mipi);
    test_interval_operations(&base.mipi, &base.pi2, "FFFF", &base.quad2, &base.empty);
    test_interval_operations(&base.mipi, &base.mipi, "TFTF", &base.mipi, &base.mipi);
    test_interval_operations(&base.mipi, &base.mipi2, "FFFF", &S1Interval::new(-PI, -PI_2), &base.empty);
    test_interval_operations(&base.mipi, &base.quad12, "FFTF", &base.quad12, &base.mipi);
    test_interval_operations(&base.mipi, &base.quad23, "FFTF", &base.quad23, &base.mipi);

    test_interval_operations(&base.quad12, &base.empty, "TTFF", &base.quad12, &base.empty);
    test_interval_operations(&base.quad12, &base.full, "FFTT", &base.full, &base.quad12);
    test_interval_operations(&base.quad12, &base.zero, "TFTF", &base.quad12, &base.zero);
    test_interval_operations(&base.quad12, &base.pi, "TFTF", &base.quad12, &base.pi);
    test_interval_operations(&base.quad12, &base.mipi, "TFTF", &base.quad12, &base.mipi);
    test_interval_operations(&base.quad12, &base.quad12, "TFTT", &base.quad12, &base.quad12);
    test_interval_operations(&base.quad12, &base.quad23, "FFTT", &base.quad123, &base.quad2);
    test_interval_operations(&base.quad12, &base.quad34, "FFTF", &base.full, &base.quad12);

    test_interval_operations(&base.quad23, &base.empty, "TTFF", &base.quad23, &base.empty);
    test_interval_operations(&base.quad23, &base.full, "FFTT", &base.full, &base.quad23);
    test_interval_operations(&base.quad23, &base.zero, "FFFF", &base.quad234, &base.empty);
    test_interval_operations(&base.quad23, &base.pi, "TTTT", &base.quad23, &base.pi);
    test_interval_operations(&base.quad23, &base.mipi, "TTTT", &base.quad23, &base.mipi);
    test_interval_operations(&base.quad23, &base.quad12, "FFTT", &base.quad123, &base.quad2);
    test_interval_operations(&base.quad23, &base.quad23, "TFTT", &base.quad23, &base.quad23);
    test_interval_operations(&base.quad23, &base.quad34, "FFTT", &base.quad234, &S1Interval::new(-PI, -PI_2));

    test_interval_operations(&base.quad1, &base.quad23, "FFTF", &base.quad123, &S1Interval::new(PI_2, PI_2));
    
    // Extract individual quad2 and quad3 from the combined intervals for specific tests
    let quad2 = S1Interval::new(PI_2, PI);
    let quad3 = S1Interval::new(PI, -PI_2);
    
    test_interval_operations(&quad2, &quad3, "FFTF", &base.quad23, &base.mipi);
    test_interval_operations(&quad3, &quad2, "FFTF", &base.quad23, &base.pi);
    test_interval_operations(&quad2, &base.pi, "TFTF", &quad2, &base.pi);
    test_interval_operations(&quad2, &base.mipi, "TFTF", &quad2, &base.mipi);
    test_interval_operations(&quad3, &base.pi, "TFTF", &quad3, &base.pi);
    test_interval_operations(&quad3, &base.mipi, "TFTF", &quad3, &base.mipi);

    test_interval_operations(&base.quad12, &base.mid12, "TTTT", &base.quad12, &base.mid12);
    test_interval_operations(&base.mid12, &base.quad12, "FFTT", &base.quad12, &base.mid12);

    let quad12eps = S1Interval::new(base.quad12.lo(), base.mid23.hi());
    let quad2hi = S1Interval::new(base.mid23.lo(), base.quad12.hi());
    test_interval_operations(&base.quad12, &base.mid23, "FFTT", &quad12eps, &quad2hi);
    test_interval_operations(&base.mid23, &base.quad12, "FFTT", &quad12eps, &quad2hi);

    // Test union of disjoint intervals gives smallest containing interval
    let quad412eps = S1Interval::new(base.mid34.lo(), base.quad12.hi());
    test_interval_operations(&base.quad12, &base.mid34, "FFFF", &quad412eps, &base.empty);
    test_interval_operations(&base.mid34, &base.quad12, "FFFF", &quad412eps, &base.empty);

    let quadeps12 = S1Interval::new(base.mid41.lo(), base.quad12.hi());
    let quad1lo = S1Interval::new(base.quad12.lo(), base.mid41.hi());
    test_interval_operations(&base.quad12, &base.mid41, "FFTT", &quadeps12, &quad1lo);
    test_interval_operations(&base.mid41, &base.quad12, "FFTT", &quadeps12, &quad1lo);

    let quad2lo = S1Interval::new(base.quad23.lo(), base.mid12.hi());
    let quad3hi = S1Interval::new(base.mid34.lo(), base.quad23.hi());
    let quadeps23 = S1Interval::new(base.mid12.lo(), base.quad23.hi());
    let quad23eps = S1Interval::new(base.quad23.lo(), base.mid34.hi());
    let quadeps123 = S1Interval::new(base.mid41.lo(), base.quad23.hi());
    test_interval_operations(&base.quad23, &base.mid12, "FFTT", &quadeps23, &quad2lo);
    test_interval_operations(&base.mid12, &base.quad23, "FFTT", &quadeps23, &quad2lo);
    test_interval_operations(&base.quad23, &base.mid23, "TTTT", &base.quad23, &base.mid23);
    test_interval_operations(&base.mid23, &base.quad23, "FFTT", &base.quad23, &base.mid23);
    test_interval_operations(&base.quad23, &base.mid34, "FFTT", &quad23eps, &quad3hi);
    test_interval_operations(&base.mid34, &base.quad23, "FFTT", &quad23eps, &quad3hi);
    test_interval_operations(&base.quad23, &base.mid41, "FFFF", &quadeps123, &base.empty);
    test_interval_operations(&base.mid41, &base.quad23, "FFFF", &quadeps123, &base.empty);
}

/// Test add_point method - matches AddPoint
#[test]
fn test_add_point() {
    let base = S1IntervalTestBase::new();
    
    let mut r = base.empty;
    r.add_point(0.0);
    assert_eq!(r, base.zero);
    
    r = base.empty;
    r.add_point(PI);
    assert_eq!(r, base.pi);
    
    r = base.empty;
    r.add_point(-PI);
    assert_eq!(r, base.mipi);
    
    r = base.empty;
    r.add_point(PI);
    r.add_point(-PI);
    assert_eq!(r, base.pi);
    
    r = base.empty;
    r.add_point(-PI);
    r.add_point(PI);
    assert_eq!(r, base.mipi);
    
    r = base.empty;
    r.add_point(base.mid12.lo());
    r.add_point(base.mid12.hi());
    assert_eq!(r, base.mid12);
    
    r = base.empty;
    r.add_point(base.mid23.lo());
    r.add_point(base.mid23.hi());
    assert_eq!(r, base.mid23);
    
    r = base.quad1;
    r.add_point(-0.9 * PI);
    r.add_point(-PI_2);
    assert_eq!(r, base.quad123);
    
    r = base.full;
    r.add_point(0.0);
    assert!(r.is_full());
    
    r = base.full;
    r.add_point(PI);
    assert!(r.is_full());
    
    r = base.full;
    r.add_point(-PI);
    assert!(r.is_full());
}

/// Test project method - matches Project
#[test]
fn test_project() {
    let r = S1Interval::new(-PI, -PI);
    assert_eq!(PI, r.project(-PI));
    assert_eq!(PI, r.project(0.0));
    
    let r = S1Interval::new(0.0, PI);
    assert_eq!(0.1, r.project(0.1));
    assert_eq!(0.0, r.project(-PI_2 + 1e-15));
    assert_eq!(PI, r.project(-PI_2 - 1e-15));
    
    let r = S1Interval::new(PI - 0.1, -PI + 0.1);
    assert_eq!(PI, r.project(PI));
    assert_eq!(PI - 0.1, r.project(1e-15));
    assert_eq!(-PI + 0.1, r.project(-1e-15));
    
    assert_eq!(0.0, S1Interval::full().project(0.0));
    assert_eq!(PI, S1Interval::full().project(PI));
    assert_eq!(PI, S1Interval::full().project(-PI));
}

/// Test from_point_pair method - matches FromPointPair
#[test]
fn test_from_point_pair() {
    let base = S1IntervalTestBase::new();
    
    assert_eq!(S1Interval::from_point_pair(-PI, PI), base.pi);
    assert_eq!(S1Interval::from_point_pair(PI, -PI), base.pi);
    assert_eq!(S1Interval::from_point_pair(base.mid34.hi(), base.mid34.lo()), base.mid34);
    assert_eq!(S1Interval::from_point_pair(base.mid23.lo(), base.mid23.hi()), base.mid23);
}


/// Test expanded method - matches Expanded
#[test]
fn test_expanded() {
    let base = S1IntervalTestBase::new();
    
    assert_eq!(base.empty.expanded(1.0), base.empty);
    assert_eq!(base.full.expanded(1.0), base.full);
    assert_eq!(base.zero.expanded(1.0), S1Interval::new(-1.0, 1.0));
    assert_eq!(base.mipi.expanded(0.01), S1Interval::new(PI - 0.01, -PI + 0.01));
    assert_eq!(base.pi.expanded(27.0), base.full);
    assert_eq!(base.pi.expanded(PI_2), base.quad23);
    assert_eq!(base.pi2.expanded(PI_2), base.quad12);
    assert_eq!(base.mipi2.expanded(PI_2), base.quad34);

    assert_eq!(base.empty.expanded(-1.0), base.empty);
    assert_eq!(base.full.expanded(-1.0), base.full);
    assert_eq!(base.quad123.expanded(-27.0), base.empty);
    assert_eq!(base.quad234.expanded(-27.0), base.empty);
    assert_eq!(base.quad123.expanded(-PI_2), base.quad2);
    assert_eq!(base.quad341.expanded(-PI_2), base.quad4);
    assert_eq!(base.quad412.expanded(-PI_2), base.quad1);
}

/// Test approx_equals method - matches ApproxEquals
#[test]
fn test_approx_equals() {
    let base = S1IntervalTestBase::new();
    
    // Choose two values such that it's okay to shift an endpoint by
    // kLo but not by kHi
    let k_lo = 4.0 * f64::EPSILON;  // < max_error default
    let k_hi = 6.0 * f64::EPSILON;  // > max_error default

    // Empty intervals
    assert!(base.empty.approx_equals(&base.empty, 1e-15));
    assert!(base.zero.approx_equals(&base.empty, 1e-15) && base.empty.approx_equals(&base.zero, 1e-15));
    assert!(base.pi.approx_equals(&base.empty, 1e-15) && base.empty.approx_equals(&base.pi, 1e-15));
    assert!(base.mipi.approx_equals(&base.empty, 1e-15) && base.empty.approx_equals(&base.mipi, 1e-15));
    assert!(!base.empty.approx_equals(&base.full, 1e-15));
    assert!(base.empty.approx_equals(&S1Interval::new(1.0, 1.0 + 2.0 * k_lo), 1e-15));
    assert!(!base.empty.approx_equals(&S1Interval::new(1.0, 1.0 + 2.0 * k_hi), 1e-15));
    assert!(S1Interval::new(PI - k_lo, -PI + k_lo).approx_equals(&base.empty, 1e-15));

    // Full intervals
    assert!(base.full.approx_equals(&base.full, 1e-15));
    assert!(!base.full.approx_equals(&base.empty, 1e-15));
    assert!(!base.full.approx_equals(&base.zero, 1e-15));
    assert!(!base.full.approx_equals(&base.pi, 1e-15));
    assert!(base.full.approx_equals(&S1Interval::new(k_lo, -k_lo), 1e-15));
    assert!(!base.full.approx_equals(&S1Interval::new(2.0 * k_hi, 0.0), 1e-15));
    assert!(S1Interval::new(-PI + k_lo, PI - k_lo).approx_equals(&base.full, 1e-15));
    assert!(!S1Interval::new(-PI, PI - 2.0 * k_hi).approx_equals(&base.full, 1e-15));

    // Singleton intervals
    assert!(base.pi.approx_equals(&base.pi, 1e-15) && base.mipi.approx_equals(&base.pi, 1e-15));
    assert!(base.pi.approx_equals(&S1Interval::new(PI - k_lo, PI - k_lo), 1e-15));
    assert!(!base.pi.approx_equals(&S1Interval::new(PI - k_hi, PI - k_hi), 1e-15));
    assert!(base.pi.approx_equals(&S1Interval::new(PI - k_lo, -PI + k_lo), 1e-15));
    assert!(!base.pi.approx_equals(&S1Interval::new(PI - k_hi, -PI), 1e-15));
    assert!(!base.zero.approx_equals(&base.pi, 1e-15));
    assert!(base.pi.union(&base.mid12).union(&base.zero).approx_equals(&base.quad12, 1e-15));
    assert!(base.quad2.intersection(&base.quad3).approx_equals(&base.pi, 1e-15));
    assert!(base.quad3.intersection(&base.quad2).approx_equals(&base.pi, 1e-15));

    // Intervals with endpoints in opposite order (inverted intervals)
    assert!(!S1Interval::new(0.0, k_lo).approx_equals(&S1Interval::new(k_lo, 0.0), 1e-15));
    assert!(!S1Interval::new(PI - 0.5 * k_lo, -PI + 0.5 * k_lo)
            .approx_equals(&S1Interval::new(-PI + 0.5 * k_lo, PI - 0.5 * k_lo), 1e-15));

    // Other intervals
    assert!(S1Interval::new(1.0 - k_lo, 2.0 + k_lo).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(S1Interval::new(1.0 + k_lo, 2.0 - k_lo).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(S1Interval::new(2.0 - k_lo, 1.0 + k_lo).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
    assert!(S1Interval::new(2.0 + k_lo, 1.0 - k_lo).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
    assert!(!S1Interval::new(1.0 - k_hi, 2.0 + k_lo).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(!S1Interval::new(1.0 + k_hi, 2.0 - k_lo).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(!S1Interval::new(2.0 - k_hi, 1.0 + k_lo).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
    assert!(!S1Interval::new(2.0 + k_hi, 1.0 - k_lo).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
    assert!(!S1Interval::new(1.0 - k_lo, 2.0 + k_hi).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(!S1Interval::new(1.0 + k_lo, 2.0 - k_hi).approx_equals(&S1Interval::new(1.0, 2.0), 1e-15));
    assert!(!S1Interval::new(2.0 - k_lo, 1.0 + k_hi).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
    assert!(!S1Interval::new(2.0 + k_lo, 1.0 - k_hi).approx_equals(&S1Interval::new(2.0, 1.0), 1e-15));
}

/// Test equality operators - matches OperatorEquals
#[test]
fn test_operator_equals() {
    let base = S1IntervalTestBase::new();
    
    assert_eq!(base.empty, base.empty);
    assert_eq!(base.full, base.full);
    assert_ne!(base.full, base.empty);
}

/// Test get_directed_hausdorff_distance method - matches GetDirectedHausdorffDistance
#[test]
fn test_get_directed_hausdorff_distance() {
    let base = S1IntervalTestBase::new();
    
    assert_abs_diff_eq!(0.0, base.empty.get_directed_hausdorff_distance(&base.empty), epsilon = 1e-15);
    assert_abs_diff_eq!(0.0, base.empty.get_directed_hausdorff_distance(&base.mid12), epsilon = 1e-15);
    assert_abs_diff_eq!(PI, base.mid12.get_directed_hausdorff_distance(&base.empty), epsilon = 1e-15);

    assert_eq!(0.0, base.quad12.get_directed_hausdorff_distance(&base.quad123));
    
    let input = S1Interval::new(3.0, -3.0);  // interval whose complement center is 0
    assert_abs_diff_eq!(3.0, S1Interval::new(-0.1, 0.2).get_directed_hausdorff_distance(&input), epsilon = 1e-6);
    assert_abs_diff_eq!(3.0 - 0.1, S1Interval::new(0.1, 0.2).get_directed_hausdorff_distance(&input), epsilon = 1e-6);
    assert_abs_diff_eq!(3.0 - 0.1, S1Interval::new(-0.2, -0.1).get_directed_hausdorff_distance(&input), epsilon = 1e-6);
}

/// Additional comprehensive tests for edge cases and numerical precision
mod precision_tests {
    use super::*;

    #[test]
    fn test_boundary_precision() {
        // Test intervals very close to π/-π boundary (stay within valid range)
        let near_pi = S1Interval::new(PI - 1e-15, PI);
        assert!(near_pi.contains_point(PI));
        assert!(near_pi.contains_point(-PI));
        
        // Test very small intervals
        let tiny = S1Interval::new(0.0, 1e-15);
        assert_eq!(tiny.get_length(), 1e-15);
        assert!(tiny.contains_point(1e-16));
        
        // Test intervals spanning exactly π
        let half_circle = S1Interval::new(-PI_2, PI_2);
        assert_eq!(half_circle.get_length(), PI);
        assert!(!half_circle.is_inverted());
    }
    
    #[test]
    fn test_wraparound_operations() {
        // Test operations that cross the ±π boundary
        let cross_boundary = S1Interval::new(2.0, -2.0);
        assert!(cross_boundary.is_inverted());
        assert!(cross_boundary.contains_point(PI));
        assert!(cross_boundary.contains_point(-PI));
        assert!(!cross_boundary.contains_point(0.0));
        
        let complement = cross_boundary.complement();
        assert!(!complement.is_inverted());
        assert!(!complement.contains_point(PI));
        assert!(complement.contains_point(0.0));
    }
    
    #[test]
    fn test_interval_arithmetic() {
        // Test complex union/intersection patterns
        let i1 = S1Interval::new(-1.0, 1.0);
        let i2 = S1Interval::new(0.5, 2.0);
        
        let union = i1.union(&i2);
        assert_eq!(union, S1Interval::new(-1.0, 2.0));
        
        let intersection = i1.intersection(&i2);
        assert_eq!(intersection, S1Interval::new(0.5, 1.0));
        
        // Test with inverted intervals
        let i3 = S1Interval::new(2.0, -1.0); // Inverted
        let union_inv = i1.union(&i3);
        assert!(union_inv.is_full() || union_inv.get_length() > PI);
    }
    
    #[test]
    fn test_special_values() {
        // Test intervals at special angles
        let quarter_turn = S1Interval::new(0.0, PI_2);
        let half_turn = S1Interval::new(0.0, PI);
        let three_quarter_turn = S1Interval::new(0.0, -PI_2); // Inverted interval spanning 3π/2
        
        assert_eq!(quarter_turn.get_length(), PI_2);
        assert_eq!(half_turn.get_length(), PI);
        assert_eq!(three_quarter_turn.get_length(), 3.0 * PI_2);
        
        // Centers should be at expected positions
        assert_eq!(quarter_turn.get_center(), PI_4);
        assert_eq!(half_turn.get_center(), PI_2);
    }
}

/// Behavioral Driven Development (BDD) tests matching the style of other ports
mod bdd_tests {
    use super::*;
    
    #[test]
    fn given_empty_interval_when_check_predicates_then_correct_flags() {
        // Given
        let empty = S1Interval::empty();
        
        // When & Then
        assert!(empty.is_empty());
        assert!(!empty.is_full());
        assert!(empty.is_inverted());
        assert!(empty.is_valid());
    }
    
    #[test]
    fn given_full_interval_when_check_predicates_then_correct_flags() {
        // Given
        let full = S1Interval::full();
        
        // When & Then
        assert!(full.is_full());
        assert!(!full.is_empty());
        assert!(!full.is_inverted());
        assert!(full.is_valid());
    }
    
    #[test]
    fn given_normal_interval_when_contains_interior_point_then_true() {
        // Given
        let interval = S1Interval::new(0.0, PI_2);
        let interior_point = PI_4;
        
        // When
        let contains = interval.contains_point(interior_point);
        let interior_contains = interval.interior_contains_point(interior_point);
        
        // Then
        assert!(contains);
        assert!(interior_contains);
    }
    
    #[test]
    fn given_normal_interval_when_contains_boundary_point_then_contains_but_not_interior() {
        // Given
        let interval = S1Interval::new(0.0, PI_2);
        let boundary_point = 0.0;
        
        // When
        let contains = interval.contains_point(boundary_point);
        let interior_contains = interval.interior_contains_point(boundary_point);
        
        // Then
        assert!(contains);
        assert!(!interior_contains);
    }
    
    #[test]
    fn given_inverted_interval_when_check_wraparound_then_correct_containment() {
        // Given (interval that wraps around ±π)
        let inverted = S1Interval::new(PI_2, -PI_2);
        
        // When & Then
        assert!(inverted.is_inverted());
        assert!(inverted.contains_point(PI));        // Should contain π
        assert!(inverted.contains_point(-PI));       // Should contain -π 
        assert!(inverted.contains_point(0.9 * PI));  // Should contain near π
        assert!(!inverted.contains_point(0.0));      // Should not contain 0
    }
    
    #[test]
    fn given_two_overlapping_intervals_when_union_then_merged_interval() {
        // Given
        let i1 = S1Interval::new(0.0, PI_2);
        let i2 = S1Interval::new(PI_4, PI);
        
        // When  
        let union = i1.union(&i2);
        
        // Then
        assert_eq!(union, S1Interval::new(0.0, PI));
        assert!(union.contains_point(0.0));
        assert!(union.contains_point(PI_4));
        assert!(union.contains_point(PI_2));
        assert!(union.contains_point(PI));
    }
    
    #[test]
    fn given_two_overlapping_intervals_when_intersection_then_overlap_only() {
        // Given
        let i1 = S1Interval::new(0.0, PI_2);
        let i2 = S1Interval::new(PI_4, PI);
        
        // When
        let intersection = i1.intersection(&i2);
        
        // Then
        assert_eq!(intersection, S1Interval::new(PI_4, PI_2));
        assert!(!intersection.contains_point(0.0));
        assert!(intersection.contains_point(PI_4));
        assert!(intersection.contains_point(PI_2));
        assert!(!intersection.contains_point(PI));
    }
    
    #[test]
    fn given_interval_when_add_exterior_point_then_expanded_appropriately() {
        // Given
        let mut interval = S1Interval::new(PI_4, PI_2);
        let exterior_point = PI; // Outside the interval
        
        // When
        interval.add_point(exterior_point);
        
        // Then
        assert!(interval.contains_point(PI_4));  // Original bounds preserved
        assert!(interval.contains_point(PI_2));  // Original bounds preserved  
        assert!(interval.contains_point(PI));    // New point included
        // Should take shorter path: expand hi to π rather than lo to π
        assert_eq!(interval, S1Interval::new(PI_4, PI));
    }
    
    #[test]
    fn given_interval_when_project_interior_point_then_unchanged() {
        // Given
        let interval = S1Interval::new(0.0, PI_2);
        let interior_point = PI_4;
        
        // When
        let projected = interval.project(interior_point);
        
        // Then
        assert_eq!(projected, interior_point);
    }
    
    #[test]
    fn given_interval_when_project_exterior_point_then_nearest_endpoint() {
        // Given
        let interval = S1Interval::new(PI_4, PI_2);
        let exterior_point = 0.0; // Closer to PI_4 than PI_2
        
        // When
        let projected = interval.project(exterior_point);
        
        // Then
        assert_eq!(projected, PI_4); // Should project to nearest endpoint
    }
    
    #[test]
    fn given_intervals_when_approx_equals_with_small_differences_then_true() {
        // Given
        let i1 = S1Interval::new(0.0, PI_2);
        let i2 = S1Interval::new(1e-16, PI_2 + 1e-16); // Tiny differences
        
        // When
        let approx_equal = i1.approx_equals(&i2, 1e-15);
        
        // Then
        assert!(approx_equal);
    }
    
    #[test]
    fn given_intervals_when_approx_equals_with_large_differences_then_false() {
        // Given
        let i1 = S1Interval::new(0.0, PI_2);
        let i2 = S1Interval::new(0.1, PI_2 + 0.1); // Large differences
        
        // When
        let approx_equal = i1.approx_equals(&i2, 1e-15);
        
        // Then
        assert!(!approx_equal);
    }
}

/// Test compatibility with existing S2 types and operations
mod integration_tests {
    use super::*;
    
    #[test]
    fn test_interval_with_s1angle_values() {
        use s2geometry_rust::S1Angle;
        
        // Test that S1Interval works with S1Angle values
        let angle1 = S1Angle::from_degrees(45.0);
        let angle2 = S1Angle::from_degrees(90.0);
        
        let interval = S1Interval::new(angle1.radians(), angle2.radians());
        
        assert!(interval.contains_point(angle1.radians()));
        assert!(interval.contains_point(angle2.radians()));
        assert_relative_eq!(interval.get_length(), PI_4, epsilon = 1e-15);
    }
    
    #[test]
    fn test_display_formatting() {
        let interval = S1Interval::new(0.0, PI_2);
        let formatted = format!("{}", interval);
        
        assert!(formatted.contains("0"));
        assert!(formatted.contains(&format!("{}", PI_2)));
        assert!(formatted.starts_with('['));
        assert!(formatted.ends_with(']'));
    }
}