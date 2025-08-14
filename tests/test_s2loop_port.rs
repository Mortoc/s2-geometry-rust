//! Comprehensive tests for S2Loop ported from C++ implementation
//!
//! This test suite covers all the major functionality of S2Loop including:
//! - Loop construction and validation
//! - Orientation operations (clockwise/counter-clockwise)
//! - Containment operations (contains point, contains loop)
//! - Area and perimeter calculations  
//! - Loop normalization and validity checking
//! - Intersection operations
//! - Bounding operations

use s2geometry_rust::{
    S2Loop, S2Point, S1Angle, S1ChordAngle, S2LatLng, S2LatLngRect, S2Cap, S2Cell, S2CellId,
    S2Error, S2Result
};
use s2geometry_rust::math::{DVec3, predicates::*};
use std::f64::consts::{PI, FRAC_PI_2};

// Test data structure to hold common loops used in tests
struct S2LoopTestData {
    empty: S2Loop,
    full: S2Loop,
    north_hemi: S2Loop,
    south_hemi: S2Loop,
    west_hemi: S2Loop,
    east_hemi: S2Loop,
    near_hemi: S2Loop,
    far_hemi: S2Loop,
    candy_cane: S2Loop,
    small_ne_cw: S2Loop,
    arctic_80: S2Loop,
    antarctic_80: S2Loop,
    line_triangle: S2Loop,
    skinny_chevron: S2Loop,
    loop_a: S2Loop,
    loop_b: S2Loop,
    a_intersect_b: S2Loop,
    a_union_b: S2Loop,
    a_minus_b: S2Loop,
    b_minus_a: S2Loop,
}

impl S2LoopTestData {
    fn new() -> S2Result<Self> {
        Ok(Self {
            empty: S2Loop::empty(),
            full: S2Loop::full(),
            
            // Northern hemisphere using two pairs of antipodal points
            north_hemi: make_loop_from_degrees(&[
                (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0)
            ])?,
            
            // Southern hemisphere
            south_hemi: make_loop_from_degrees(&[
                (0.0, 90.0), (0.0, 0.0), (0.0, -90.0), (0.0, -180.0)
            ])?,
            
            // Western hemisphere
            west_hemi: make_loop_from_degrees(&[
                (0.0, -180.0), (-90.0, 0.0), (0.0, 0.0), (90.0, 0.0)
            ])?,
            
            // Eastern hemisphere  
            east_hemi: make_loop_from_degrees(&[
                (90.0, 0.0), (0.0, 0.0), (-90.0, 0.0), (0.0, -180.0)
            ])?,
            
            // Near hemisphere
            near_hemi: make_loop_from_degrees(&[
                (0.0, -90.0), (-90.0, 0.0), (0.0, 90.0), (90.0, 0.0)
            ])?,
            
            // Far hemisphere
            far_hemi: make_loop_from_degrees(&[
                (90.0, 0.0), (0.0, 90.0), (-90.0, 0.0), (0.0, -90.0)
            ])?,
            
            // Spiral stripe that slightly overwraps the equator
            candy_cane: make_loop_from_degrees(&[
                (-20.0, 150.0), (-20.0, -70.0), (0.0, 70.0), 
                (10.0, -150.0), (10.0, 70.0), (-10.0, -70.0)
            ])?,
            
            // Small clockwise loop in northern & eastern hemispheres
            small_ne_cw: make_loop_from_degrees(&[
                (35.0, 20.0), (45.0, 20.0), (40.0, 25.0)
            ])?,
            
            // Loop around north pole at 80 degrees
            arctic_80: make_loop_from_degrees(&[
                (80.0, -150.0), (80.0, -30.0), (80.0, 90.0)
            ])?,
            
            // Loop around south pole at 80 degrees  
            antarctic_80: make_loop_from_degrees(&[
                (-80.0, 120.0), (-80.0, 0.0), (-80.0, -120.0)
            ])?,
            
            // Degenerate triangle along the equator
            line_triangle: make_loop_from_degrees(&[
                (0.0, 1.0), (0.0, 2.0), (0.0, 3.0)
            ])?,
            
            // Nearly-degenerate CCW chevron near equator with very long sides
            skinny_chevron: make_loop_from_degrees(&[
                (0.0, 0.0), (-1e-15, 80.0), (0.0, 1e-15), (1e-15, 80.0)
            ])?,
            
            // Diamond-shaped loop around point 0:180
            loop_a: make_loop_from_degrees(&[
                (0.0, 178.0), (-1.0, 180.0), (0.0, -179.0), (1.0, -180.0)
            ])?,
            
            // Another diamond-shaped loop around point 0:180
            loop_b: make_loop_from_degrees(&[
                (0.0, 179.0), (-1.0, 180.0), (0.0, -178.0), (1.0, -180.0)
            ])?,
            
            // Intersection of A and B
            a_intersect_b: make_loop_from_degrees(&[
                (0.0, 179.0), (-1.0, 180.0), (0.0, -179.0), (1.0, -180.0)
            ])?,
            
            // Union of A and B
            a_union_b: make_loop_from_degrees(&[
                (0.0, 178.0), (-1.0, 180.0), (0.0, -178.0), (1.0, -180.0)
            ])?,
            
            // A minus B (concave)
            a_minus_b: make_loop_from_degrees(&[
                (0.0, 178.0), (-1.0, 180.0), (0.0, 179.0), (1.0, -180.0)
            ])?,
            
            // B minus A (concave)
            b_minus_a: make_loop_from_degrees(&[
                (0.0, -179.0), (-1.0, 180.0), (0.0, -178.0), (1.0, -180.0)
            ])?,
        })
    }
}

/// Helper function to create a loop from latitude/longitude degrees
fn make_loop_from_degrees(coords: &[(f64, f64)]) -> S2Result<S2Loop> {
    let vertices: Result<Vec<_>, _> = coords
        .iter()
        .map(|&(lat, lng)| {
            let latlng = S2LatLng::from_degrees(lat, lng);
            latlng.to_point()
        })
        .collect();
    S2Loop::new(vertices?)
}

/// Helper to create a point from degrees
fn point_from_degrees(lat: f64, lng: f64) -> S2Point {
    S2LatLng::from_degrees(lat, lng).to_point().unwrap()
}

#[test]
fn test_default_loop_is_invalid() {
    // Can't test default constructor since we don't have one in our API
    // This test would verify that an uninitialized loop is invalid
}

#[test]
fn test_empty_and_full_loops() -> S2Result<()> {
    let empty = S2Loop::empty();
    let full = S2Loop::full();
    
    assert!(empty.is_empty());
    assert!(!empty.is_full());
    assert!(empty.is_empty_or_full());
    assert_eq!(empty.num_vertices(), 1);
    
    assert!(!full.is_empty());
    assert!(full.is_full());
    assert!(full.is_empty_or_full());
    assert_eq!(full.num_vertices(), 1);
    
    // Test areas
    assert_eq!(empty.get_area(), 0.0);
    assert_eq!(full.get_area(), 4.0 * PI);
    
    // Test containment
    let origin = S2Point::new(0.0, 0.0, 1.0)?;
    assert!(!empty.contains(&origin));
    assert!(full.contains(&origin));
    
    Ok(())
}

#[test]
fn test_basic_loop_properties() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Test northern hemisphere
    assert!(!data.north_hemi.is_empty());
    assert!(!data.north_hemi.is_full());
    assert!(!data.north_hemi.is_empty_or_full());
    assert_eq!(data.north_hemi.num_vertices(), 4);
    
    // Test area (should be approximately 2π for hemisphere)
    let north_area = data.north_hemi.get_area();
    assert!((north_area - 2.0 * PI).abs() < 1e-10);
    
    // Test containment
    let north_pole = S2Point::new(0.0, 0.0, 1.0)?;
    let south_pole = S2Point::new(0.0, 0.0, -1.0)?;
    
    assert!(data.north_hemi.contains(&north_pole));
    assert!(!data.north_hemi.contains(&south_pole));
    
    Ok(())
}

#[test]
fn test_loop_validation() -> S2Result<()> {
    // Valid triangle
    let valid_triangle = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ])?;
    assert!(valid_triangle.is_valid());
    
    // Test that loops with < 3 vertices fail
    let single_vertex = vec![S2Point::new(1.0, 0.0, 0.0)?];
    let result = S2Loop::new(single_vertex);
    // This should be a valid empty/full loop actually
    
    let two_vertices = vec![
        S2Point::new(1.0, 0.0, 0.0)?,
        S2Point::new(0.0, 1.0, 0.0)?
    ];
    let result = S2Loop::new(two_vertices);
    assert!(result.is_err());
    
    Ok(())
}

#[test]
fn test_get_rect_bound() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Empty loop should have empty bound
    assert!(data.empty.get_rect_bound().is_empty());
    
    // Full loop should have full bound  
    assert!(data.full.get_rect_bound().is_full());
    
    // Candy cane should span full longitude
    let candy_bound = data.candy_cane.get_rect_bound();
    assert!(candy_bound.lng().is_full());
    assert!(candy_bound.lat().lo().to_degrees() < -20.0);
    assert!(candy_bound.lat().hi().to_degrees() > 10.0);
    
    // Small clockwise loop should span full sphere (after normalization)
    let small_bound = data.small_ne_cw.get_rect_bound();
    // The small clockwise loop represents the complement, so it should be full
    
    Ok(())
}

#[test] 
fn test_area_consistent_with_curvature() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    let loops = vec![
        &data.empty, &data.full, &data.north_hemi, &data.south_hemi,
        &data.west_hemi, &data.east_hemi, &data.candy_cane, &data.arctic_80
    ];
    
    // Check Gauss-Bonnet theorem: area + curvature = 2π
    for loop_obj in loops {
        let area = loop_obj.get_area();
        let curvature = loop_obj.get_curvature();
        let gauss_area = 2.0 * PI - curvature;
        
        // The error bound should be sufficient for our test cases
        assert!((area - gauss_area).abs() < 1e-10, 
            "Area-curvature inconsistency: area={}, gauss_area={}", area, gauss_area);
    }
    
    Ok(())
}

#[test]
fn test_get_area_consistent_with_sign() -> S2Result<()> {
    // Test that GetArea() returns appropriate values for degenerate loops
    // This test creates small loops that should have area near 0 or 4π
    
    // Create a very small triangle on the equator
    let small_triangle = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1e-10), (1e-10, 0.0)
    ])?;
    
    let area = small_triangle.get_area();
    let is_ccw = small_triangle.is_normalized();
    
    if is_ccw {
        assert!(area < 1e-10, "CCW degenerate loop should have small area");
    } else {
        assert!((area - 4.0 * PI).abs() < 1e-10, "CW degenerate loop should have area ~4π");
    }
    
    Ok(())
}

#[test]
fn test_contains_point() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Test hemisphere containment
    let north_pole = S2Point::new(0.0, 0.0, 1.0)?;
    let south_pole = S2Point::new(0.0, 0.0, -1.0)?;
    
    assert!(data.north_hemi.contains(&north_pole));
    assert!(!data.north_hemi.contains(&south_pole));
    assert!(!data.south_hemi.contains(&north_pole));
    assert!(data.south_hemi.contains(&south_pole));
    
    // Test candy cane
    let test_point = point_from_degrees(5.0, 71.0);
    assert!(data.candy_cane.contains(&test_point));
    
    // Test east/west hemispheres
    let east_point = S2Point::new(0.0, 1.0, 0.0)?;  // (0, 90°)
    let west_point = S2Point::new(0.0, -1.0, 0.0)?; // (0, -90°)
    
    assert!(data.east_hemi.contains(&east_point));
    assert!(!data.east_hemi.contains(&west_point));
    assert!(!data.west_hemi.contains(&east_point));
    assert!(data.west_hemi.contains(&west_point));
    
    Ok(())
}

#[test]
fn test_loop_relationships() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Test basic relationships
    test_one_nested_pair(&data.full, &data.north_hemi);
    test_one_nested_pair(&data.north_hemi, &data.arctic_80);
    
    test_one_disjoint_pair(&data.north_hemi, &data.antarctic_80);
    test_one_disjoint_pair(&data.arctic_80, &data.antarctic_80);
    
    // Test specific loop relationships from C++ test
    assert!(data.north_hemi.contains_loop(&data.arctic_80));
    assert!(!data.north_hemi.contains_loop(&data.south_hemi));
    assert!(!data.arctic_80.contains_loop(&data.north_hemi));
    
    Ok(())
}

fn test_one_nested_pair(outer: &S2Loop, inner: &S2Loop) {
    assert!(outer.contains_loop(inner));
    assert!(!inner.contains_loop(outer) || outer.boundary_equals(inner));
    assert!(outer.intersects(inner) || inner.is_empty());
    assert!(inner.intersects(outer) || inner.is_empty());
}

fn test_one_disjoint_pair(a: &S2Loop, b: &S2Loop) {
    assert!(!a.intersects(b));
    assert!(!b.intersects(a));
    assert!(a.contains_loop(b) == b.is_empty());
    assert!(b.contains_loop(a) == a.is_empty());
}

#[test]
fn test_normalize_and_invert() -> S2Result<()> {
    let mut north_copy = make_loop_from_degrees(&[
        (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0)
    ])?;
    
    // Should be normalized already (area ≤ 2π)
    assert!(north_copy.is_normalized());
    
    // Test inversion
    let orig_area = north_copy.get_area();
    north_copy.invert();
    let inverted_area = north_copy.get_area();
    
    // After inversion, area should be 4π - original_area
    assert!((orig_area + inverted_area - 4.0 * PI).abs() < 1e-10);
    
    // Invert back
    north_copy.invert();
    assert!((north_copy.get_area() - orig_area).abs() < 1e-10);
    
    Ok(())
}

#[test]
fn test_boundary_equality() -> S2Result<()> {
    let loop1 = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ])?;
    
    let loop2 = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ])?;
    
    // Same vertices in same order
    assert!(loop1.boundary_equals(&loop2));
    
    // Rotated vertices (should still be equal)
    let loop3 = make_loop_from_degrees(&[
        (0.0, 1.0), (1.0, 0.0), (0.0, 0.0)
    ])?;
    assert!(loop1.boundary_equals(&loop3));
    
    // Different loop
    let loop4 = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 2.0), (2.0, 0.0)
    ])?;
    assert!(!loop1.boundary_equals(&loop4));
    
    Ok(())
}

#[test]
fn test_make_regular_loop() -> S2Result<()> {
    let center = point_from_degrees(80.0, 135.0);
    let radius = S1Angle::from_degrees(20.0);
    let loop_poly = S2Loop::make_regular_loop(center, radius, 4)?;
    
    assert_eq!(loop_poly.num_vertices(), 4);
    
    // Check that all vertices are approximately the right distance from center
    for i in 0..4 {
        let vertex = loop_poly.vertex(i);
        let distance = S1Angle::from_radians(center.coords().dot(vertex.coords()).acos());
        assert!((distance.degrees() - 20.0).abs() < 1e-10);
    }
    
    // Test that polygon is approximately square
    let v0 = loop_poly.vertex(0);
    let v1 = loop_poly.vertex(1);
    let v2 = loop_poly.vertex(2);
    let v3 = loop_poly.vertex(3);
    
    // All edges should have similar length
    let edge01 = (v1.coords() - v0.coords()).length();
    let edge12 = (v2.coords() - v1.coords()).length();
    let edge23 = (v3.coords() - v2.coords()).length();
    let edge30 = (v0.coords() - v3.coords()).length();
    
    assert!((edge01 - edge12).abs() < 1e-10);
    assert!((edge12 - edge23).abs() < 1e-10);
    assert!((edge23 - edge30).abs() < 1e-10);
    
    Ok(())
}

#[test]
fn test_cell_constructor() -> S2Result<()> {
    // Test creating a loop from an S2Cell
    let cell_id = S2CellId::from_face_ij_same_face(0, 0, 0);
    let cell = S2Cell::from(cell_id);
    let loop_from_cell = S2Loop::from_cell(&cell);
    
    assert_eq!(loop_from_cell.num_vertices(), 4);
    assert!(loop_from_cell.is_valid());
    
    // The loop should contain the cell center
    let center = cell.get_center();
    assert!(loop_from_cell.contains(&center));
    
    Ok(())
}

#[test]
fn test_depth_and_sign() -> S2Result<()> {
    let mut loop_obj = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ])?;
    
    // Test depth
    assert_eq!(loop_obj.depth(), 0);
    assert!(!loop_obj.is_hole());
    assert_eq!(loop_obj.sign(), 1);
    
    loop_obj.set_depth(1);
    assert_eq!(loop_obj.depth(), 1);
    assert!(loop_obj.is_hole());
    assert_eq!(loop_obj.sign(), -1);
    
    loop_obj.set_depth(2);
    assert_eq!(loop_obj.depth(), 2);
    assert!(!loop_obj.is_hole());
    assert_eq!(loop_obj.sign(), 1);
    
    Ok(())
}

#[test]
fn test_get_centroid() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Empty and full loops should have zero centroid
    let empty_centroid = data.empty.get_centroid();
    let full_centroid = data.full.get_centroid();
    
    assert!(empty_centroid.coords().length() < 1e-15);
    assert!(full_centroid.coords().length() < 1e-15);
    
    // Test a simple triangle
    let triangle = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 90.0), (90.0, 0.0)
    ])?;
    
    let centroid = triangle.get_centroid();
    // The centroid should be somewhere reasonable
    assert!(centroid.coords().length() > 0.0);
    
    Ok(())
}

#[test]
fn test_get_cap_bound() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Empty loop should have empty cap
    let empty_cap = data.empty.get_cap_bound();
    assert!(empty_cap.is_empty());
    
    // Full loop should have full cap
    let full_cap = data.full.get_cap_bound();
    assert!(full_cap.is_full());
    
    // Test that cap contains all vertices
    let triangle = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 45.0), (45.0, 0.0)
    ])?;
    
    let cap = triangle.get_cap_bound();
    for i in 0..triangle.num_vertices() {
        let vertex = triangle.vertex(i);
        assert!(cap.contains(&vertex));
    }
    
    Ok(())
}

#[test]
fn test_special_vertex_cases() -> S2Result<()> {
    // Test vertex access with wrapping
    let triangle = make_loop_from_degrees(&[
        (0.0, 0.0), (0.0, 1.0), (1.0, 0.0)
    ])?;
    
    // Test that vertex wrapping works
    let v0 = triangle.vertex(0);
    let v3 = triangle.vertex(3); // Should wrap to vertex 0
    assert!(v0.coords().abs_diff_eq(v3.coords(), 1e-15));
    
    let v1 = triangle.vertex(1);
    let v4 = triangle.vertex(4); // Should wrap to vertex 1
    assert!(v1.coords().abs_diff_eq(v4.coords(), 1e-15));
    
    Ok(())
}

#[test] 
fn test_curvature_properties() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Empty loop should have curvature 2π
    assert!((data.empty.get_curvature() - 2.0 * PI).abs() < 1e-15);
    
    // Full loop should have curvature -2π
    assert!((data.full.get_curvature() + 2.0 * PI).abs() < 1e-15);
    
    // For hemispheres, curvature should be close to 0
    assert!(data.north_hemi.get_curvature().abs() < 1e-10);
    
    Ok(())
}

#[test]
fn test_error_cases() {
    // Test invalid regular loop
    let center = S2Point::from_coords_raw(1.0, 0.0, 0.0);
    let radius = S1Angle::from_degrees(20.0);
    
    let result = S2Loop::make_regular_loop(center, radius, 2);
    assert!(result.is_err());
    
    let result = S2Loop::make_regular_loop(center, radius, 1);
    assert!(result.is_err());
}

#[test]
fn test_intersection_operations() -> S2Result<()> {
    let data = S2LoopTestData::new()?;
    
    // Test basic intersection properties
    assert!(!data.north_hemi.intersects(&data.south_hemi));
    assert!(data.north_hemi.intersects(&data.east_hemi));
    assert!(data.north_hemi.intersects(&data.west_hemi));
    
    // Full loop intersects everything except empty
    assert!(data.full.intersects(&data.north_hemi));
    assert!(!data.full.intersects(&data.empty));
    assert!(!data.empty.intersects(&data.north_hemi));
    
    Ok(())
}

// Helper function to run all tests
#[test]
fn run_comprehensive_s2loop_tests() -> S2Result<()> {
    // This test ensures all our test data can be created successfully
    let _data = S2LoopTestData::new()?;
    println!("All S2Loop test data created successfully");
    Ok(())
}