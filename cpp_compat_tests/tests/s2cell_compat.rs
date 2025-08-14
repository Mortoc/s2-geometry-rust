//! S2Cell C++ Compatibility Tests
//!
//! Tests that verify functional equivalence between Rust and C++ S2Cell implementations.

use s2geometry_rust::{S2Point, S2Cell, S2CellId, math::DVec3};
use s2geometry_cpp_compat_tests::*;

const TOLERANCE: f64 = 1e-15;

#[test]
fn test_s2cell_area_equivalence() {
    // Test cell area calculations match between implementations
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0), 
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
        DVec3::new(-1.0, 1.0, 0.0).normalize(),
        DVec3::new(0.5, -0.5, 0.7).normalize(),
    ];
    
    let test_levels = vec![0, 5, 10, 15, 20, 25, 30];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        for &level in &test_levels {
            let cell_id = S2CellId::from_point(&point).parent(level);
            let rust_cell = S2Cell::new(cell_id);
            
            // Test Rust implementation
            let rust_area = rust_cell.exact_area();
            
            // Test C++ implementation  
            let cpp_area = cpp_s2cell_area(cell_id.into());
            
            assert!(angles_equal_approx(rust_area, cpp_area, TOLERANCE),
                "Area mismatch for point {} level {}: Rust={}, C++={}, diff={}", 
                point_idx, level, rust_area, cpp_area, (rust_area - cpp_area).abs());
        }
    }
}

#[test]
fn test_s2cell_perimeter_equivalence() {
    // Test cell perimeter calculations match between implementations
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
    ];
    
    let test_levels = vec![0, 5, 10, 15, 20];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        for &level in &test_levels {
            let cell_id = S2CellId::from_point(&point).parent(level);
            let rust_cell = S2Cell::new(cell_id);
            
            // Calculate Rust perimeter by summing edge lengths
            let mut rust_perimeter = 0.0;
            for i in 0..4 {
                let v0 = rust_cell.vertex(i);
                let v1 = rust_cell.vertex((i + 1) % 4);
                rust_perimeter += v0.angle(&v1).radians();
            }
            
            // Test C++ implementation
            let cpp_perimeter = cpp_s2cell_perimeter(cell_id.into());
            
            assert!(angles_equal_approx(rust_perimeter, cpp_perimeter, TOLERANCE),
                "Perimeter mismatch for point {} level {}: Rust={}, C++={}, diff={}", 
                point_idx, level, rust_perimeter, cpp_perimeter, (rust_perimeter - cpp_perimeter).abs());
        }
    }
}

#[test]
fn test_s2cell_vertices_equivalence() {
    // Test that cell vertices match between implementations
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 0.0).normalize(),
    ];
    
    let test_levels = vec![0, 5, 10, 15, 20];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        for &level in &test_levels {
            let cell_id = S2CellId::from_point(&point).parent(level);
            let rust_cell = S2Cell::new(cell_id);
            
            // Get Rust vertices
            let rust_vertices: Vec<S2Point> = (0..4).map(|i| rust_cell.vertex(i)).collect();
            
            // Get C++ vertices
            let cpp_vertices = cpp_s2cell_vertices(cell_id.into());
            
            assert_eq!(rust_vertices.len(), cpp_vertices.len(),
                "Vertex count mismatch for point {} level {}: Rust={}, C++={}", 
                point_idx, level, rust_vertices.len(), cpp_vertices.len());
            
            // Compare each vertex
            for (vertex_idx, (rust_vertex, cpp_vertex)) in rust_vertices.iter().zip(cpp_vertices.iter()).enumerate() {
                let rust_coords = rust_vertex.coords();
                assert!(points_equal_approx(
                    S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
                    *cpp_vertex,
                    TOLERANCE
                ), "Vertex {} mismatch for point {} level {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
                   vertex_idx, point_idx, level, 
                   rust_coords.x, rust_coords.y, rust_coords.z,
                   cpp_vertex.x, cpp_vertex.y, cpp_vertex.z);
            }
        }
    }
}

#[test]
fn test_s2cellid_level_equivalence() {
    // Test that cell ID level calculations match
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
    ];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        for level in 0..=30 {
            let cell_id = S2CellId::from_point(&point).parent(level);
            
            // Test Rust implementation
            let rust_level = cell_id.level();
            
            // Test C++ implementation
            let cpp_level = cpp_s2cellid_level(cell_id.into());
            
            assert_eq!(rust_level, cpp_level as u32,
                "Level mismatch for point {} at level {}: Rust={}, C++={}", 
                point_idx, level, rust_level, cpp_level);
        }
    }
}

#[test]
fn test_s2cellid_parent_equivalence() {
    // Test that parent cell calculations match
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
    ];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        let leaf_cell_id = S2CellId::from_point(&point);
        
        for parent_level in 0..=25 {
            // Test Rust implementation
            let rust_parent = leaf_cell_id.parent(parent_level);
            
            // Test C++ implementation
            let cpp_parent = cpp_s2cellid_parent(leaf_cell_id.into(), parent_level as i32);
            
            assert_eq!(rust_parent.id(), cpp_parent.id,
                "Parent mismatch for point {} at level {}: Rust={}, C++={}", 
                point_idx, parent_level, rust_parent.id(), cpp_parent.id);
        }
    }
}

#[test]
fn test_s2cellid_to_point_equivalence() {
    // Test that cell ID to point conversion matches
    let test_points = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
    ];
    
    let test_levels = vec![0, 5, 10, 15, 20, 25, 30];
    
    for (point_idx, point_vec) in test_points.iter().enumerate() {
        let point = S2Point::from_normalized(*point_vec);
        
        for &level in &test_levels {
            let cell_id = S2CellId::from_point(&point).parent(level);
            
            // Test Rust implementation
            let rust_center = cell_id.to_point();
            
            // Test C++ implementation
            let cpp_center = cpp_s2cellid_to_point(cell_id.into());
            
            let rust_coords = rust_center.coords();
            assert!(points_equal_approx(
                S2PointCpp { x: rust_coords.x, y: rust_coords.y, z: rust_coords.z },
                cpp_center,
                TOLERANCE
            ), "Center point mismatch for point {} level {}: Rust=({}, {}, {}), C++=({}, {}, {})", 
               point_idx, level,
               rust_coords.x, rust_coords.y, rust_coords.z,
               cpp_center.x, cpp_center.y, cpp_center.z);
        }
    }
}