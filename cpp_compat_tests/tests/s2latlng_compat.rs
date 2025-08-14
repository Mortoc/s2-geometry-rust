//! S2LatLng C++ Compatibility Tests
//!
//! Tests that verify functional equivalence between Rust and C++ S2LatLng implementations.

use s2geometry_rust::{S2Point, S2LatLng, S1Angle, math::DVec3};
use s2geometry_cpp_compat_tests::*;

const TOLERANCE: f64 = 1e-15;

#[test]
fn test_s2latlng_from_point_equivalence() {
    // Test conversion from S2Point to S2LatLng matches between implementations
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
        DVec3::new(0.5, 0.5, 0.7071).normalize(),
    ];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        // Test Rust implementation
        let rust_latlng = S2LatLng::from_point(&point);
        
        // Test C++ implementation
        let cpp_latlng = cpp_s2latlng_from_point(point.into());
        
        assert!(angles_equal_approx(rust_latlng.lat().radians(), cpp_latlng.lat_radians, TOLERANCE),
            "Latitude mismatch for point {}: Rust={}, C++={}, diff={}", 
            point_idx, rust_latlng.lat().radians(), cpp_latlng.lat_radians, 
            (rust_latlng.lat().radians() - cpp_latlng.lat_radians).abs());
            
        assert!(angles_equal_approx(rust_latlng.lng().radians(), cpp_latlng.lng_radians, TOLERANCE),
            "Longitude mismatch for point {}: Rust={}, C++={}, diff={}", 
            point_idx, rust_latlng.lng().radians(), cpp_latlng.lng_radians,
            (rust_latlng.lng().radians() - cpp_latlng.lng_radians).abs());
    }
}

#[test]
fn test_s2latlng_to_point_equivalence() {
    // Test conversion from S2LatLng to S2Point matches between implementations
    let test_latlngs = vec![
        // (lat_radians, lng_radians)
        (0.0, 0.0),                                    // Equator, Prime meridian
        (std::f64::consts::PI / 2.0, 0.0),            // North pole
        (-std::f64::consts::PI / 2.0, 0.0),           // South pole
        (0.0, std::f64::consts::PI),                  // Equator, antimeridian
        (std::f64::consts::PI / 4.0, std::f64::consts::PI / 4.0),  // 45°N, 45°E
        (-std::f64::consts::PI / 4.0, -std::f64::consts::PI / 4.0), // 45°S, 45°W
        (std::f64::consts::PI / 6.0, std::f64::consts::PI / 3.0),   // 30°N, 60°E
        (-std::f64::consts::PI / 3.0, -std::f64::consts::PI / 6.0), // 60°S, 30°W
        (0.5, 1.0),                                   // Arbitrary coordinates
        (-0.7, -2.5),                                 // More arbitrary coordinates
    ];
    
    for (latlng_idx, (lat_rad, lng_rad)) in test_latlngs.iter().enumerate() {
        let rust_latlng = S2LatLng::from_radians(*lat_rad, *lng_rad);
        
        // Test Rust implementation
        let rust_point = rust_latlng.to_point();
        
        // Test C++ implementation
        let cpp_point = cpp_s2latlng_to_point(rust_latlng.into());
        
        let rust_coords = rust_point.coords();
        assert!(points_equal_approx(
            S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
            cpp_point,
            TOLERANCE
        ), "Point mismatch for latlng {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
           latlng_idx,
           rust_coords.x, rust_coords.y, rust_coords.z,
           cpp_point.x, cpp_point.y, cpp_point.z);
    }
}

#[test]
fn test_s2latlng_distance_equivalence() {
    // Test distance calculations between S2LatLng points match
    let test_latlngs = vec![
        // (lat_radians, lng_radians)
        (0.0, 0.0),                                    // Equator, Prime meridian
        (std::f64::consts::PI / 2.0, 0.0),            // North pole
        (-std::f64::consts::PI / 2.0, 0.0),           // South pole
        (0.0, std::f64::consts::PI),                  // Equator, antimeridian
        (std::f64::consts::PI / 4.0, std::f64::consts::PI / 4.0),  // 45°N, 45°E
        (-std::f64::consts::PI / 4.0, -std::f64::consts::PI / 4.0), // 45°S, 45°W
    ];
    
    // Test all pairs of points
    for (i, (lat1, lng1)) in test_latlngs.iter().enumerate() {
        for (j, (lat2, lng2)) in test_latlngs.iter().enumerate() {
            if i >= j { continue; } // Avoid duplicate pairs and self-distance
            
            let rust_latlng1 = S2LatLng::from_radians(*lat1, *lng1);
            let rust_latlng2 = S2LatLng::from_radians(*lat2, *lng2);
            
            // Test Rust implementation
            let rust_distance = rust_latlng1.get_distance(&rust_latlng2).radians();
            
            // Test C++ implementation
            let cpp_distance = cpp_s2latlng_distance(rust_latlng1.into(), rust_latlng2.into());
            
            assert!(angles_equal_approx(rust_distance, cpp_distance, TOLERANCE),
                "Distance mismatch between points {} and {}: Rust={}, C++={}, diff={}", 
                i, j, rust_distance, cpp_distance, (rust_distance - cpp_distance).abs());
        }
    }
}

#[test]
fn test_s2latlng_round_trip_equivalence() {
    // Test that Point->LatLng->Point round trips produce equivalent results
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
        DVec3::new(-1.0, -1.0, -1.0).normalize(),
        DVec3::new(0.707, 0.707, 0.0),
    ];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let original_point = S2Point::from_normalized(*point_vec);
        
        // Rust round trip: Point -> LatLng -> Point
        let rust_latlng = S2LatLng::from_point(&original_point);
        let rust_round_trip = rust_latlng.to_point();
        
        // C++ round trip: Point -> LatLng -> Point
        let cpp_latlng = cpp_s2latlng_from_point(original_point.into());
        let cpp_round_trip = cpp_s2latlng_to_point(cpp_latlng);
        
        // Compare the round-trip results
        let rust_coords = rust_round_trip.coords();
        assert!(points_equal_approx(
            S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
            cpp_round_trip,
            TOLERANCE
        ), "Round trip mismatch for point {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
           point_idx,
           rust_coords.x, rust_coords.y, rust_coords.z,
           cpp_round_trip.x, cpp_round_trip.y, cpp_round_trip.z);
    }
}

#[test]
fn test_s2latlng_edge_cases() {
    // Test edge cases like poles, antimeridian crossings, etc.
    
    // North pole
    let north_pole = S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0));
    let rust_north_latlng = S2LatLng::from_point(&north_pole);
    let cpp_north_latlng = cpp_s2latlng_from_point(north_pole.into());
    
    assert!(angles_equal_approx(rust_north_latlng.lat().radians(), std::f64::consts::PI / 2.0, TOLERANCE),
        "North pole latitude incorrect in Rust: {}", rust_north_latlng.lat().radians());
    assert!(angles_equal_approx(cpp_north_latlng.lat_radians, std::f64::consts::PI / 2.0, TOLERANCE),
        "North pole latitude incorrect in C++: {}", cpp_north_latlng.lat_radians);
    
    // South pole
    let south_pole = S2Point::from_normalized(DVec3::new(0.0, 0.0, -1.0));
    let rust_south_latlng = S2LatLng::from_point(&south_pole);
    let cpp_south_latlng = cpp_s2latlng_from_point(south_pole.into());
    
    assert!(angles_equal_approx(rust_south_latlng.lat().radians(), -std::f64::consts::PI / 2.0, TOLERANCE),
        "South pole latitude incorrect in Rust: {}", rust_south_latlng.lat().radians());
    assert!(angles_equal_approx(cpp_south_latlng.lat_radians, -std::f64::consts::PI / 2.0, TOLERANCE),
        "South pole latitude incorrect in C++: {}", cpp_south_latlng.lat_radians);
    
    // Antimeridian crossing - test that longitude wrapping is consistent
    let antimeridian_point = S2Point::from_normalized(DVec3::new(-1.0, 0.0, 0.0));
    let rust_anti_latlng = S2LatLng::from_point(&antimeridian_point);
    let cpp_anti_latlng = cpp_s2latlng_from_point(antimeridian_point.into());
    
    // Longitude should be π (or close to it, considering normalization)
    let expected_lng = std::f64::consts::PI;
    assert!(angles_equal_approx(rust_anti_latlng.lng().radians().abs(), expected_lng, TOLERANCE),
        "Antimeridian longitude incorrect in Rust: {}", rust_anti_latlng.lng().radians());
    assert!(angles_equal_approx(cpp_anti_latlng.lng_radians.abs(), expected_lng, TOLERANCE),
        "Antimeridian longitude incorrect in C++: {}", cpp_anti_latlng.lng_radians);
}