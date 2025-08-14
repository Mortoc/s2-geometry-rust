//! Port of s2latlng_test.cc - Critical S2LatLng tests
//!
//! This module ports the essential tests from C++ s2latlng_test.cc to validate
//! our S2LatLng implementation matches C++ behavior exactly.

use s2geometry_rust::{S1Angle, S2LatLng, S2Point, S2Error};
use s2geometry_rust::math::constants::*;
use s2geometry_rust::latlng::S2LatLngHash;
use approx::{assert_relative_eq, assert_abs_diff_eq};
use std::collections::HashMap;

/// Test basic construction and accessors (from C++ TestBasic)
#[test]
fn test_basic() {
    let ll_rad = S2LatLng::from_radians(PI_4, PI_2);
    assert_eq!(PI_4, ll_rad.lat().radians());
    assert_eq!(PI_2, ll_rad.lng().radians());
    assert!(ll_rad.is_valid());

    let ll_deg = S2LatLng::from_degrees(45.0, 90.0);
    assert_eq!(ll_rad, ll_deg);
    assert!(ll_deg.is_valid());

    assert!(!S2LatLng::from_degrees(-91.0, 0.0).is_valid());
    assert!(!S2LatLng::from_degrees(0.0, 181.0).is_valid());

    let bad = S2LatLng::from_degrees(120.0, 200.0);
    assert!(!bad.is_valid());
    let better = bad.normalized();
    assert!(better.is_valid());
    assert_eq!(S1Angle::from_degrees(90.0), better.lat());
    assert_abs_diff_eq!(
        S1Angle::from_degrees(-160.0).radians(), 
        better.lng().radians(), 
        epsilon = 1e-15
    );

    let bad2 = S2LatLng::from_degrees(-100.0, -360.0);
    assert!(!bad2.is_valid());
    let better2 = bad2.normalized();
    assert!(better2.is_valid());
    assert_eq!(S1Angle::from_degrees(-90.0), better2.lat());
    assert_abs_diff_eq!(0.0, better2.lng().radians(), epsilon = 1e-15);

    // Test arithmetic operations
    assert!((S2LatLng::from_degrees(10.0, 20.0) + S2LatLng::from_degrees(20.0, 30.0))
        .approx_equals(&S2LatLng::from_degrees(30.0, 50.0), S1Angle::from_radians(1e-15)));
    assert!((S2LatLng::from_degrees(10.0, 20.0) - S2LatLng::from_degrees(20.0, 30.0))
        .approx_equals(&S2LatLng::from_degrees(-10.0, -10.0), S1Angle::from_radians(1e-15)));
    assert!((0.5 * S2LatLng::from_degrees(10.0, 20.0))
        .approx_equals(&S2LatLng::from_degrees(5.0, 10.0), S1Angle::from_radians(1e-15)));

    // Check that Invalid() returns an invalid point
    let invalid = S2LatLng::invalid();
    assert!(!invalid.is_valid());

    // Check that the default constructor sets latitude and longitude to 0
    let default_ll = S2LatLng::default();
    assert!(default_ll.is_valid());
    assert_eq!(0.0, default_ll.lat().radians());
    assert_eq!(0.0, default_ll.lng().radians());
}

/// Test conversion between S2Point and S2LatLng (from C++ TestConversion)
#[test]
fn test_conversion() {
    // Test special cases: poles, "date line"
    assert_abs_diff_eq!(
        90.0,
        S2LatLng::from_point(S2LatLng::from_degrees(90.0, 65.0).to_point().unwrap()).lat().degrees(),
        epsilon = 1e-13
    );
    
    assert_abs_diff_eq!(
        -PI_2,
        S2LatLng::from_point(S2LatLng::from_radians(-PI_2, 1.0).to_point().unwrap()).lat().radians(),
        epsilon = 1e-15
    );
    
    assert_abs_diff_eq!(
        180.0,
        S2LatLng::from_point(S2LatLng::from_degrees(12.2, 180.0).to_point().unwrap()).lng().degrees().abs(),
        epsilon = 1e-13
    );
    
    assert_abs_diff_eq!(
        PI,
        S2LatLng::from_point(S2LatLng::from_radians(0.1, -PI).to_point().unwrap()).lng().radians().abs(),
        epsilon = 1e-15
    );

    // Test a bunch of random points (simplified version of C++ test)
    use rand::Rng;
    let mut rng = rand::thread_rng();
    
    for _ in 0..1000 {
        // Generate random unit vector
        let x: f64 = rng.gen_range(-1.0..1.0);
        let y: f64 = rng.gen_range(-1.0..1.0); 
        let z: f64 = rng.gen_range(-1.0..1.0);
        
        if let Ok(p) = S2Point::new(x, y, z) {
            let converted = S2LatLng::from_point(p).to_point().unwrap();
            
            // Check that conversion is approximately equal (allowing for numerical precision)
            let diff = (p.coords() - converted.coords()).length();
            assert!(diff < 1e-13, "Point conversion failed: original={:?}, converted={:?}, diff={}", 
                   p.coords(), converted.coords(), diff);
        }
    }
}

/// Helper function for exact bit-level comparison (matches C++ IsIdentical)
fn is_identical(x: f64, y: f64) -> bool {
    x.to_bits() == y.to_bits()
}

/// Test handling of negative zeros (from C++ NegativeZeros)
#[test]
fn test_negative_zeros() {
    // Create test points with specific coordinate patterns
    let p1 = S2Point::new(1.0, 0.0, -0.0).unwrap();
    assert!(is_identical(S2LatLng::latitude(p1).radians(), 0.0));
    
    let p2 = S2Point::new(1.0, -0.0, 0.0).unwrap();
    assert!(is_identical(S2LatLng::longitude(p2).radians(), 0.0));
    
    let p3 = S2Point::new(-1.0, -0.0, 0.0).unwrap();
    assert!(is_identical(S2LatLng::longitude(p3).radians(), PI));
    
    let p4 = S2Point::new(-0.0, 0.0, 1.0).unwrap();
    assert!(is_identical(S2LatLng::longitude(p4).radians(), 0.0));
    
    let p5 = S2Point::new(-0.0, -0.0, 1.0).unwrap();
    assert!(is_identical(S2LatLng::longitude(p5).radians(), 0.0));
}

/// Test infinity handling (from C++ InfIsInvalid)  
#[test]
fn test_inf_is_invalid() {
    assert!(!S2LatLng::from_degrees(f64::INFINITY, -122.0).is_valid());
    assert!(!S2LatLng::from_degrees(37.0, f64::INFINITY).is_valid());

    // Also check the results of normalized()
    assert!(!S2LatLng::from_degrees(f64::INFINITY, -122.0).normalized().is_valid());
    assert!(!S2LatLng::from_degrees(37.0, f64::INFINITY).normalized().is_valid());
}

/// Test NaN handling (from C++ NanIsInvalid)
#[test]
fn test_nan_is_invalid() {
    assert!(!S2LatLng::from_degrees(f64::NAN, -122.0).is_valid());
    assert!(!S2LatLng::from_degrees(37.0, f64::NAN).is_valid());

    // Also check the results of normalized()
    assert!(!S2LatLng::from_degrees(37.0, f64::NAN).normalized().is_valid());
    assert!(!S2LatLng::from_degrees(37.0, f64::NAN).normalized().is_valid());
}

/// Test distance calculation (from C++ TestDistance)
#[test]
fn test_distance() {
    assert_eq!(
        0.0,
        S2LatLng::from_degrees(90.0, 0.0)
            .get_distance(&S2LatLng::from_degrees(90.0, 0.0))
            .radians()
    );
    
    assert_abs_diff_eq!(
        77.0,
        S2LatLng::from_degrees(-37.0, 25.0)
            .get_distance(&S2LatLng::from_degrees(-66.0, -155.0))
            .degrees(),
        epsilon = 1e-13
    );
    
    assert_abs_diff_eq!(
        115.0,
        S2LatLng::from_degrees(0.0, 165.0)
            .get_distance(&S2LatLng::from_degrees(0.0, -80.0))
            .degrees(),
        epsilon = 1e-13
    );
    
    assert_abs_diff_eq!(
        180.0,
        S2LatLng::from_degrees(47.0, -127.0)
            .get_distance(&S2LatLng::from_degrees(-47.0, 53.0))
            .degrees(),
        epsilon = 2e-6
    );
}

/// Test string formatting (from C++ TestToString)
#[test]
fn test_to_string() {
    #[derive(Debug)]
    struct TestCase {
        lat: f64,
        lng: f64,
        expected_lat: f64,
        expected_lng: f64,
    }

    let test_cases = vec![
        TestCase { lat: 0.0, lng: 0.0, expected_lat: 0.0, expected_lng: 0.0 },
        TestCase { lat: 1.5, lng: 91.7, expected_lat: 1.5, expected_lng: 91.7 },
        TestCase { lat: 9.9, lng: -0.31, expected_lat: 9.9, expected_lng: -0.31 },
        TestCase { lat: 2.0_f64.sqrt(), lng: -5.0_f64.sqrt(), expected_lat: 1.414214, expected_lng: -2.236068 },
        TestCase { lat: 91.3, lng: 190.4, expected_lat: 90.0, expected_lng: -169.6 },
        TestCase { lat: -100.0, lng: -710.0, expected_lat: -90.0, expected_lng: 10.0 },
    ];

    for (i, test_case) in test_cases.iter().enumerate() {
        let p = S2LatLng::from_degrees(test_case.lat, test_case.lng);
        let output = p.to_string_in_degrees();

        // Parse the output string 
        let parts: Vec<&str> = output.split(',').collect();
        assert_eq!(parts.len(), 2, "Test case {}: Expected 2 parts, got {}", i, parts.len());

        let lat: f64 = parts[0].parse().expect(&format!("Test case {}: Failed to parse lat", i));
        let lng: f64 = parts[1].parse().expect(&format!("Test case {}: Failed to parse lng", i));

        if (test_case.expected_lat - lat).abs() > 1e-6 {
            panic!("Test case {}: lat mismatch. Expected {}, got {}", i, test_case.expected_lat, lat);
        }
        if (test_case.expected_lng - lng).abs() > 1e-6 {
            panic!("Test case {}: lng mismatch. Expected {}, got {}", i, test_case.expected_lng, lng);
        }
    }
}

/// Test hash functionality (from C++ TestHashCode)
#[test]
fn test_hash_code() {
    let mut map = HashMap::new();
    map.insert(S2LatLng::from_degrees(0.0, 10.0), 1);
    map.insert(S2LatLng::from_degrees(2.0, 12.0), 2);
    map.insert(S2LatLng::from_degrees(5.0, 15.0), 3);
    map.insert(S2LatLng::from_degrees(7.0, 17.0), 4);
    map.insert(S2LatLng::from_degrees(11.0, 19.0), 5);
    
    assert_eq!(map.len(), 5);
    assert_eq!(1, map[&S2LatLng::from_degrees(0.0, 10.0)]);
    assert_eq!(2, map[&S2LatLng::from_degrees(2.0, 12.0)]);
    assert_eq!(3, map[&S2LatLng::from_degrees(5.0, 15.0)]);
    assert_eq!(4, map[&S2LatLng::from_degrees(7.0, 17.0)]);
    assert_eq!(5, map[&S2LatLng::from_degrees(11.0, 19.0)]);
}

/// Test legacy hash function
#[test]
fn test_legacy_hash_function() {
    let hasher = S2LatLngHash::default();
    let ll = S2LatLng::from_degrees(47.6062, -122.3321);
    let hash_value = hasher.hash(&ll);
    
    // Just verify it produces a hash value (exact value may vary by implementation)
    assert_ne!(hash_value, 0);
    
    // Same input should produce same hash
    let hash_value2 = hasher.hash(&ll);
    assert_eq!(hash_value, hash_value2);
    
    // Different input should (very likely) produce different hash
    let ll2 = S2LatLng::from_degrees(47.6063, -122.3321);
    let hash_value3 = hasher.hash(&ll2);
    assert_ne!(hash_value, hash_value3);
}

/// Test E5/E6/E7 representations
#[test]
fn test_e5_e6_e7_representations() {
    // Test E5
    let ll_e5 = S2LatLng::from_e5(4500000, 9000000); // 45.0°, 90.0°
    assert_abs_diff_eq!(45.0, ll_e5.lat().degrees(), epsilon = 1e-10);
    assert_abs_diff_eq!(90.0, ll_e5.lng().degrees(), epsilon = 1e-10);
    
    // Test E6
    let ll_e6 = S2LatLng::from_e6(45000000, 90000000); // 45.0°, 90.0°
    assert_abs_diff_eq!(45.0, ll_e6.lat().degrees(), epsilon = 1e-11);
    assert_abs_diff_eq!(90.0, ll_e6.lng().degrees(), epsilon = 1e-11);
    
    // Test E7
    let ll_e7 = S2LatLng::from_e7(450000000, 900000000); // 45.0°, 90.0°
    assert_abs_diff_eq!(45.0, ll_e7.lat().degrees(), epsilon = 1e-12);
    assert_abs_diff_eq!(90.0, ll_e7.lng().degrees(), epsilon = 1e-12);
    
    // Test unsigned variants
    let ll_ue6 = S2LatLng::from_unsigned_e6(45000000_u32, 90000000_u32);
    assert_abs_diff_eq!(45.0, ll_ue6.lat().degrees(), epsilon = 1e-11);
    assert_abs_diff_eq!(90.0, ll_ue6.lng().degrees(), epsilon = 1e-11);
    
    let ll_ue7 = S2LatLng::from_unsigned_e7(450000000_u32, 900000000_u32);
    assert_abs_diff_eq!(45.0, ll_ue7.lat().degrees(), epsilon = 1e-12);
    assert_abs_diff_eq!(90.0, ll_ue7.lng().degrees(), epsilon = 1e-12);
}

/// Test coordinate access methods
#[test]
fn test_accessors() {
    let ll = S2LatLng::from_degrees(47.6062, -122.3321);
    
    // Test individual angle accessors
    let lat = ll.lat();
    let lng = ll.lng();
    
    assert_abs_diff_eq!(47.6062, lat.degrees(), epsilon = 1e-10);
    assert_abs_diff_eq!(-122.3321, lng.degrees(), epsilon = 1e-10);
    
    // Test coordinate access
    let coords = ll.coords();
    assert_abs_diff_eq!(lat.radians(), coords.x(), epsilon = 1e-15);
    assert_abs_diff_eq!(lng.radians(), coords.y(), epsilon = 1e-15);
}

/// Test static methods for extracting lat/lng from points
#[test]
fn test_static_lat_lng_methods() {
    let point = S2Point::new(1.0, 0.0, 0.0).unwrap();
    
    let lat = S2LatLng::latitude(point);
    let lng = S2LatLng::longitude(point);
    
    assert_abs_diff_eq!(0.0, lat.degrees(), epsilon = 1e-15);
    assert_abs_diff_eq!(0.0, lng.degrees(), epsilon = 1e-15);
    
    // Test point at north pole
    let north_pole = S2Point::new(0.0, 0.0, 1.0).unwrap();
    let lat_np = S2LatLng::latitude(north_pole);
    assert_abs_diff_eq!(90.0, lat_np.degrees(), epsilon = 1e-13);
}

/// Test approximate equality
#[test] 
fn test_approx_equals() {
    let ll1 = S2LatLng::from_degrees(45.0, 90.0);
    let ll2 = S2LatLng::from_degrees(45.0000001, 90.0000001);
    
    // Should not be exactly equal
    assert_ne!(ll1, ll2);
    
    // But should be approximately equal with reasonable tolerance
    assert!(ll1.approx_equals(&ll2, S1Angle::from_degrees(1e-5)));
    
    // Should not be approximately equal with tight tolerance
    assert!(!ll1.approx_equals(&ll2, S1Angle::from_degrees(1e-8)));
}

/// Test comparison operators
#[test]
fn test_comparison() {
    let ll1 = S2LatLng::from_degrees(45.0, 90.0);
    let ll2 = S2LatLng::from_degrees(45.0, 91.0); // Same lat, higher lng
    let ll3 = S2LatLng::from_degrees(46.0, 90.0); // Higher lat, same lng
    
    // Test lexicographic ordering (lat first, then lng)
    assert!(ll1 < ll2); // Same lat, ll1 has smaller lng
    assert!(ll1 < ll3); // ll1 has smaller lat
    assert!(ll2 > ll1);
    assert!(ll3 > ll1);
    
    // Test equality
    assert_eq!(ll1, ll1);
    assert_ne!(ll1, ll2);
}

/// Test edge cases and boundary conditions
#[test]
fn test_edge_cases() {
    // Test exact poles
    let north_pole = S2LatLng::from_degrees(90.0, 0.0);
    let south_pole = S2LatLng::from_degrees(-90.0, 0.0);
    
    assert!(north_pole.is_valid());
    assert!(south_pole.is_valid());
    
    // Test date line
    let date_line_east = S2LatLng::from_degrees(0.0, 180.0);
    let date_line_west = S2LatLng::from_degrees(0.0, -180.0);
    
    assert!(date_line_east.is_valid());
    assert!(date_line_west.is_valid());
    
    // Test equator and prime meridian
    let equator_prime = S2LatLng::from_degrees(0.0, 0.0);
    assert!(equator_prime.is_valid());
    
    // Test boundary values
    let max_valid_lat = S2LatLng::from_degrees(90.0, 0.0);
    let min_valid_lat = S2LatLng::from_degrees(-90.0, 0.0);
    let max_valid_lng = S2LatLng::from_degrees(0.0, 180.0);
    let min_valid_lng = S2LatLng::from_degrees(0.0, -180.0);
    
    assert!(max_valid_lat.is_valid());
    assert!(min_valid_lat.is_valid());
    assert!(max_valid_lng.is_valid());
    assert!(min_valid_lng.is_valid());
    
    // Test just outside boundaries
    let lat_too_high = S2LatLng::from_degrees(90.1, 0.0);
    let lat_too_low = S2LatLng::from_degrees(-90.1, 0.0);
    let lng_too_high = S2LatLng::from_degrees(0.0, 180.1);
    let lng_too_low = S2LatLng::from_degrees(0.0, -180.1);
    
    assert!(!lat_too_high.is_valid());
    assert!(!lat_too_low.is_valid());
    assert!(!lng_too_high.is_valid());
    assert!(!lng_too_low.is_valid());
}

/// Test normalization of complex cases
#[test]
fn test_complex_normalization() {
    // Test latitude over 90° (should reflect around pole and adjust longitude)
    let over_pole = S2LatLng::from_degrees(120.0, 30.0);
    let normalized = over_pole.normalized();
    assert!(normalized.is_valid());
    // 120° lat becomes 60° lat on other side, lng should be shifted by 180°
    assert_abs_diff_eq!(60.0, normalized.lat().degrees(), epsilon = 1e-13);
    assert_abs_diff_eq!(-150.0, normalized.lng().degrees(), epsilon = 1e-13);
    
    // Test negative latitude under -90°
    let under_pole = S2LatLng::from_degrees(-120.0, 30.0);
    let normalized2 = under_pole.normalized();
    assert!(normalized2.is_valid());
    // -120° lat becomes -60° lat on other side, lng should be shifted by 180°
    assert_abs_diff_eq!(-60.0, normalized2.lat().degrees(), epsilon = 1e-13);
    assert_abs_diff_eq!(-150.0, normalized2.lng().degrees(), epsilon = 1e-13);
    
    // Test longitude wrapping (multiple times around)
    let multi_wrap = S2LatLng::from_degrees(45.0, 720.0 + 30.0); // 2 full rotations + 30°
    let normalized3 = multi_wrap.normalized();
    assert!(normalized3.is_valid());
    assert_abs_diff_eq!(45.0, normalized3.lat().degrees(), epsilon = 1e-13);
    assert_abs_diff_eq!(30.0, normalized3.lng().degrees(), epsilon = 1e-13);
}

/// Test point conversion accuracy
#[test]
fn test_point_conversion_accuracy() {
    // Test various known coordinates
    let test_coords = vec![
        (0.0, 0.0),     // Equator, Prime Meridian
        (90.0, 0.0),    // North Pole
        (-90.0, 0.0),   // South Pole
        (45.0, 90.0),   // 45°N, 90°E
        (-45.0, -90.0), // 45°S, 90°W
        (0.0, 180.0),   // Equator, Date Line
        (60.0, 30.0),   // Arbitrary point
    ];
    
    for &(lat, lng) in &test_coords {
        let original = S2LatLng::from_degrees(lat, lng);
        
        // Convert to point and back
        let point = original.to_point().unwrap();
        let roundtrip = S2LatLng::from_point(point);
        
        // Should be very close (allowing for floating point precision)
        assert!(original.approx_equals(&roundtrip, S1Angle::from_degrees(1e-13)),
                "Round trip failed for ({}, {}): original={:?}, roundtrip={:?}", 
                lat, lng, original, roundtrip);
    }
}

/// Test arithmetic operations maintain precision  
#[test]
fn test_arithmetic_precision() {
    let a = S2LatLng::from_degrees(10.123456789, 20.987654321);
    let b = S2LatLng::from_degrees(5.111111111, 3.222222222);
    
    // Test addition
    let sum = a + b;
    assert_abs_diff_eq!(sum.lat().degrees(), 15.234567900, epsilon = 1e-9);
    assert_abs_diff_eq!(sum.lng().degrees(), 24.209876543, epsilon = 1e-9);
    
    // Test subtraction  
    let diff = a - b;
    assert_abs_diff_eq!(diff.lat().degrees(), 5.012345678, epsilon = 1e-9);
    assert_abs_diff_eq!(diff.lng().degrees(), 17.765432099, epsilon = 1e-9);
    
    // Test scalar multiplication
    let scaled = 2.5 * a;
    assert_abs_diff_eq!(scaled.lat().degrees(), 25.308641973, epsilon = 1e-9);
    assert_abs_diff_eq!(scaled.lng().degrees(), 52.469135803, epsilon = 1e-9);
    
    // Test commutative scalar multiplication
    let scaled2 = a * 2.5;
    assert_eq!(scaled, scaled2);
}

/// Test distance calculation accuracy with known values
#[test]
fn test_distance_accuracy() {
    // Test distance between New York and London (known distance ~5585 km ≈ 50.6°)
    let new_york = S2LatLng::from_degrees(40.7589, -73.9851);
    let london = S2LatLng::from_degrees(51.5074, -0.1278);
    
    let distance = new_york.get_distance(&london);
    
    // Should be approximately 50-51 degrees
    assert!(distance.degrees() > 48.0 && distance.degrees() < 53.0,
            "Distance between NYC and London should be ~50.6°, got {:.2}°", distance.degrees());
    
    // Test antipodal points (should be 180°)
    let point1 = S2LatLng::from_degrees(45.0, 0.0);
    let antipodal = S2LatLng::from_degrees(-45.0, 180.0);
    
    let antipodal_dist = point1.get_distance(&antipodal);
    assert_abs_diff_eq!(180.0, antipodal_dist.degrees(), epsilon = 1e-10);
    
    // Test nearby points (should be small distance)
    let seattle = S2LatLng::from_degrees(47.6062, -122.3321);
    let tacoma = S2LatLng::from_degrees(47.2529, -122.4443); // ~40 km south
    
    let nearby_dist = seattle.get_distance(&tacoma);
    assert!(nearby_dist.degrees() < 1.0, "Distance between nearby cities should be < 1°");
    assert!(nearby_dist.degrees() > 0.1, "Distance should be > 0.1°");
}

/// BDD-style tests for clear specification
mod bdd_tests {
    use super::*;

    /// BDD: Given valid lat/lng degrees, When I create S2LatLng, Then it should be valid
    #[test]
    fn given_valid_degrees_when_create_s2latlng_then_valid() {
        // Given
        let lat_deg = 47.6062;
        let lng_deg = -122.3321;
        
        // When  
        let latlng = S2LatLng::from_degrees(lat_deg, lng_deg);
        
        // Then
        assert!(latlng.is_valid());
        assert_abs_diff_eq!(lat_deg, latlng.lat().degrees(), epsilon = 1e-10);
        assert_abs_diff_eq!(lng_deg, latlng.lng().degrees(), epsilon = 1e-10);
    }
    
    /// BDD: Given invalid coordinates, When I normalize, Then result should be valid
    #[test]
    fn given_invalid_coordinates_when_normalize_then_valid() {
        // Given
        let invalid_latlng = S2LatLng::from_degrees(100.0, 270.0);
        assert!(!invalid_latlng.is_valid());
        
        // When
        let normalized = invalid_latlng.normalized();
        
        // Then
        assert!(normalized.is_valid());
    }
    
    /// BDD: Given two S2LatLng points, When I compute distance, Then result matches expected
    #[test]
    fn given_two_points_when_compute_distance_then_matches_expected() {
        // Given
        let origin = S2LatLng::from_degrees(0.0, 0.0);
        let point_90_degrees_away = S2LatLng::from_degrees(90.0, 0.0);
        
        // When
        let distance = origin.get_distance(&point_90_degrees_away);
        
        // Then
        assert_abs_diff_eq!(90.0, distance.degrees(), epsilon = 1e-13);
    }
    
    /// BDD: Given S2LatLng, When I convert to point and back, Then coordinates preserved  
    #[test]
    fn given_s2latlng_when_convert_to_point_and_back_then_coordinates_preserved() {
        // Given
        let original = S2LatLng::from_degrees(45.0, -122.0);
        
        // When
        let point = original.to_point().unwrap();
        let roundtrip = S2LatLng::from_point(point);
        
        // Then
        assert!(original.approx_equals(&roundtrip, S1Angle::from_degrees(1e-13)));
    }
    
    /// BDD: Given two S2LatLng values, When I add them, Then coordinates add componentwise
    #[test]
    fn given_two_s2latlng_when_add_then_coordinates_add_componentwise() {
        // Given
        let ll1 = S2LatLng::from_degrees(10.0, 20.0);
        let ll2 = S2LatLng::from_degrees(5.0, 15.0);
        
        // When
        let sum = ll1 + ll2;
        
        // Then
        assert_abs_diff_eq!(15.0, sum.lat().degrees(), epsilon = 1e-15);
        assert_abs_diff_eq!(35.0, sum.lng().degrees(), epsilon = 1e-15);
    }
    
    /// BDD: Given S2LatLng with high precision, When I format as string, Then precision preserved
    #[test]
    fn given_high_precision_s2latlng_when_format_string_then_precision_preserved() {
        // Given
        let precise_ll = S2LatLng::from_degrees(47.123456, -122.987654);
        
        // When  
        let formatted = precise_ll.to_string_in_degrees();
        
        // Then
        assert!(formatted.contains("47.123456"));
        assert!(formatted.contains("-122.987654"));
        assert!(formatted.contains(","));
    }
}

/// Integration tests with other S2 components
mod integration_tests {
    use super::*;

    /// Test integration with S1Angle
    #[test]
    fn test_s1angle_integration() {
        let angle_45 = S1Angle::from_degrees(45.0);
        let angle_90 = S1Angle::from_degrees(90.0);
        
        let latlng = S2LatLng::new(angle_45, angle_90);
        
        assert_eq!(angle_45, latlng.lat());
        assert_eq!(angle_90, latlng.lng());
        assert!(latlng.is_valid());
    }
    
    /// Test integration with S2Point
    #[test]
    fn test_s2point_integration() {
        // Create S2Point at known location
        let point = S2Point::new(1.0, 0.0, 0.0).unwrap(); // (1,0,0) -> (0°, 0°)
        
        // Convert to S2LatLng
        let latlng = S2LatLng::from_point(point);
        
        // Verify coordinates
        assert_abs_diff_eq!(0.0, latlng.lat().degrees(), epsilon = 1e-13);
        assert_abs_diff_eq!(0.0, latlng.lng().degrees(), epsilon = 1e-13);
        
        // Convert back to S2Point
        let point2 = latlng.to_point().unwrap();
        
        // Should be nearly identical
        let diff = (point.coords() - point2.coords()).length();
        assert!(diff < 1e-15);
    }
    
    /// Test with error handling
    #[test] 
    fn test_error_handling() {
        // Test conversion of invalid coordinates to point
        let invalid_ll = S2LatLng::from_degrees(f64::NAN, 0.0);
        let result = invalid_ll.to_point();
        
        assert!(result.is_err());
        if let Err(e) = result {
            assert!(format!("{}", e).contains("finite"));
        }
    }
}