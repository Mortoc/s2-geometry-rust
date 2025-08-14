//! Basic S2Cell functionality tests
//! 
//! Simple tests to verify S2Cell implementation works correctly

use s2geometry_rust::{S2Cell, S2CellId, S2Point};

#[test]
fn test_s2cell_basic_creation() {
    // Test creating S2Cell from face
    for face in 0..6 {
        let cell = S2Cell::from_face(face);
        match cell {
            Ok(c) => {
                assert_eq!(c.face(), face);
                assert_eq!(c.level(), 0);
                assert!(!c.is_leaf());
                println!("Face {} cell created successfully", face);
            }
            Err(e) => {
                panic!("Failed to create face {} cell: {}", face, e);
            }
        }
    }
}

#[test]
fn test_s2cell_from_point() {
    // Test creating S2Cell from a point
    let point = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let cell_id = S2CellId::from_point(&point);
    
    match S2Cell::new(cell_id) {
        Ok(cell) => {
            println!("Cell from point: face={}, level={}", cell.face(), cell.level());
            assert!(cell.level() >= 0);
            assert!(cell.face() >= 0 && cell.face() < 6);
        }
        Err(e) => {
            panic!("Failed to create cell from point: {}", e);
        }
    }
}

#[test]
fn test_s2cell_vertices() {
    // Test vertex computation
    let cell = S2Cell::from_face(0).unwrap();
    
    for k in 0..4 {
        let vertex = cell.get_vertex(k);
        let vertex_raw = cell.get_vertex_raw(k);
        
        // Vertices should be valid points
        assert!(vertex.coords().length() > 0.0);
        assert!(vertex_raw.coords().length() > 0.0);
        
        println!("Vertex {}: ({}, {}, {})", k, vertex.x(), vertex.y(), vertex.z());
    }
}

#[test]
fn test_s2cell_edges() {
    // Test edge computation
    let cell = S2Cell::from_face(0).unwrap();
    
    for k in 0..4 {
        let edge = cell.get_edge(k);
        let edge_raw = cell.get_edge_raw(k);
        
        // Edges should be valid vectors
        assert!(edge.coords().length() > 0.0);
        assert!(edge_raw.coords().length() > 0.0);
        
        println!("Edge {}: ({}, {}, {})", k, edge.x(), edge.y(), edge.z());
    }
}

#[test]
fn test_s2cell_area() {
    // Test area calculations
    let cell = S2Cell::from_face(0).unwrap();
    
    let avg_area = cell.get_average_area();
    let approx_area = cell.approx_area();
    let exact_area = cell.exact_area();
    
    // All areas should be positive
    assert!(avg_area > 0.0);
    assert!(approx_area > 0.0);
    assert!(exact_area > 0.0);
    
    println!("Areas - avg: {}, approx: {}, exact: {}", avg_area, approx_area, exact_area);
}

#[test]
fn test_s2cell_containment() {
    // Test basic containment
    let cell = S2Cell::from_face(0).unwrap();
    let center = cell.get_center();
    
    // Cell should contain its center
    assert!(cell.contains(&center));
    println!("Cell contains its center: OK");
}