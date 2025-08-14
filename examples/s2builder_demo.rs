//! S2Builder Demonstration
//!
//! This example demonstrates the core functionality of S2Builder including:
//! - Building polygons from edge loops
//! - Building polylines from edge chains  
//! - Using different snap functions for robustness
//! - Handling multiple output layers simultaneously
//! - Error handling and validation

use s2geometry_rust::builder::*;
use s2geometry_rust::{S2Point, S2Polyline, S2Loop, S2Builder, EdgeType, S2LatLng};
use s2geometry_rust::math::DVec3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("S2Builder Demonstration");
    println!("=====================");
    
    // Demonstrate basic builder usage
    basic_triangle_example()?;
    
    // Demonstrate snap functions
    snap_function_examples()?;
    
    // Demonstrate polyline construction
    polyline_construction_example()?;
    
    // Demonstrate multiple layers
    multiple_layers_example()?;
    
    println!("\nAll demonstrations completed successfully!");
    Ok(())
}

/// Demonstrate building a simple triangle polygon
fn basic_triangle_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n1. Basic Triangle Construction");
    println!("------------------------------");
    
    // Create builder with default options
    let mut builder = S2Builder::default();
    
    // Define triangle vertices
    let vertices = vec![
        S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),   // Point on +X axis
        S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),   // Point on +Y axis  
        S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),   // Point on +Z axis
    ];
    
    println!("Triangle vertices:");
    for (i, vertex) in vertices.iter().enumerate() {
        println!("  Vertex {}: {:?}", i, vertex.coords());
    }
    
    // Add triangle edges as undirected loop (for polygon construction)
    builder.add_loop(&vertices, EdgeType::Undirected)?;
    
    // Set up polygon layer
    let mut loops = Vec::new();
    let polygon_options = PolygonLayerOptions::default()
        .with_edge_type(EdgeType::Undirected)
        .with_validate(true);
    let layer = S2PolygonLayer::new(&mut loops, polygon_options);
    builder.add_layer(Box::new(layer))?;
    
    println!("Added {} edges to builder", builder.num_input_edges());
    println!("Added {} layers to builder", builder.num_layers());
    
    // Note: Full build() would require complete implementation
    // This demonstrates the setup and API usage
    
    Ok(())
}

/// Demonstrate different snap function behaviors
fn snap_function_examples() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n2. Snap Function Examples");
    println!("-------------------------");
    
    let test_point = S2Point::from_normalized(DVec3::new(1.0, 0.1, 0.1));
    println!("Original point: {:?}", test_point.coords());
    
    // Identity snap function (no snapping)
    let identity_snap = IdentitySnapFunction::new(1e-15);
    let identity_result = identity_snap.snap_point(test_point);
    println!("Identity snap result: {:?}", identity_result.coords());
    println!("  Snap radius: {:.2e}", identity_snap.snap_radius());
    
    // S2CellId snap function (snap to cell centers)
    let cell_snap = S2CellIdSnapFunction::new(10); // Level 10 cells
    let cell_result = cell_snap.snap_point(test_point);
    println!("S2CellId snap result: {:?}", cell_result.coords());
    println!("  Snap radius: {:.2e}", cell_snap.snap_radius());
    println!("  Min vertex separation: {:.2e}", cell_snap.min_vertex_separation());
    
    // Integer lat/lng snap function
    let latlng_snap = IntLatLngSnapFunction::degrees(); // 1 degree precision
    let latlng_result = latlng_snap.snap_point(test_point);
    println!("Lat/Lng snap result: {:?}", latlng_result.coords());
    println!("  Snap radius: {:.2e}", latlng_snap.snap_radius());
    
    // Demonstrate coordinate changes
    let original_latlng = S2LatLng::from_point(test_point);
    let snapped_latlng = S2LatLng::from_point(latlng_result);
    println!("Original lat/lng: {:.6}°, {:.6}°", 
             original_latlng.lat().degrees(), original_latlng.lng().degrees());
    println!("Snapped lat/lng: {:.6}°, {:.6}°", 
             snapped_latlng.lat().degrees(), snapped_latlng.lng().degrees());
    
    Ok(())
}

/// Demonstrate polyline construction from edge chains
fn polyline_construction_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n3. Polyline Construction");
    println!("------------------------");
    
    // Create builder with S2CellId snapping for robustness
    let snap_function = S2CellIdSnapFunction::new(15); // Fine-grained snapping
    let options = Options::default()
        .with_snap_function(snap_function)
        .with_validate_input(true)
        .with_validate_output(true);
    let mut builder = S2Builder::new(options);
    
    // Create a path along the equator
    let mut path_vertices = Vec::new();
    for i in 0..=4 {
        let lng = (i as f64) * 30.0; // 0°, 30°, 60°, 90°, 120° longitude
        let latlng = S2LatLng::from_degrees(0.0, lng);
        path_vertices.push(latlng.to_point()?);
    }
    
    println!("Polyline path vertices:");
    for (i, vertex) in path_vertices.iter().enumerate() {
        let latlng = S2LatLng::from_point(*vertex);
        println!("  Vertex {}: {:.1}°, {:.1}°", i, 
                latlng.lat().degrees(), latlng.lng().degrees());
    }
    
    // Add as directed polyline
    builder.add_polyline(&path_vertices, EdgeType::Directed)?;
    
    // Set up polyline layer
    let mut polyline_result = None;
    let polyline_options = PolylineLayerOptions::default()
        .with_edge_type(EdgeType::Directed)
        .with_validate(true);
    let layer = S2PolylineLayer::new(&mut polyline_result, polyline_options);
    builder.add_layer(Box::new(layer))?;
    
    println!("Added {} edges for polyline", builder.num_input_edges());
    
    Ok(())
}

/// Demonstrate building multiple output types simultaneously
fn multiple_layers_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n4. Multiple Layers Example");
    println!("--------------------------");
    
    // Create builder with integer lat/lng snapping
    let snap_function = IntLatLngSnapFunction::microdegrees(); // Microdegree precision
    let options = Options::default()
        .with_snap_function(snap_function)
        .with_split_crossing_edges(true)
        .with_validate_input(true);
    let mut builder = S2Builder::new(options);
    
    // Add a square polygon
    let square_vertices = vec![
        S2LatLng::from_degrees(0.0, 0.0).to_point()?,   // Bottom-left
        S2LatLng::from_degrees(0.0, 1.0).to_point()?,   // Bottom-right  
        S2LatLng::from_degrees(1.0, 1.0).to_point()?,   // Top-right
        S2LatLng::from_degrees(1.0, 0.0).to_point()?,   // Top-left
    ];
    
    // Add as polygon (undirected edges)
    builder.add_loop(&square_vertices, EdgeType::Undirected)?;
    
    // Add a diagonal polyline across the square
    let diagonal_vertices = vec![
        S2LatLng::from_degrees(0.0, 0.0).to_point()?,   // Bottom-left
        S2LatLng::from_degrees(1.0, 1.0).to_point()?,   // Top-right
    ];
    builder.add_polyline(&diagonal_vertices, EdgeType::Directed)?;
    
    // Set up multiple output layers
    
    // Layer 1: Polygon output
    let mut polygon_loops = Vec::new();
    let polygon_layer = S2PolygonLayer::with_defaults(&mut polygon_loops);
    builder.add_layer(Box::new(polygon_layer))?;
    
    // Layer 2: All polylines
    let mut all_polylines = Vec::new();
    let polyline_layer = S2PolylineVectorLayer::with_defaults(&mut all_polylines);
    builder.add_layer(Box::new(polyline_layer))?;
    
    println!("Input summary:");
    println!("  Total edges: {}", builder.num_input_edges());
    println!("  Total layers: {}", builder.num_layers());
    println!("  Polygon edges: 4 (square boundary)");
    println!("  Polyline edges: 1 (diagonal)");
    
    println!("Expected outputs:");
    println!("  Polygon layer: 1 square loop");
    println!("  Polyline layer: Multiple polylines from all edges");
    
    Ok(())
}

/// Demonstrate error handling and validation
fn error_handling_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n5. Error Handling Examples");
    println!("--------------------------");
    
    let mut builder = S2Builder::default();
    
    // Try to add invalid edge (antipodal points)
    let point_a = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    let point_b = S2Point::from_normalized(DVec3::new(-1.0, 0.0, 0.0)); // Antipodal
    
    match builder.add_edge(point_a, point_b) {
        Ok(_) => println!("ERROR: Should have failed for antipodal points"),
        Err(e) => println!("Correctly rejected antipodal edge: {}", e),
    }
    
    // Try to add vertex after building
    builder.built = true; // Simulate build() call
    let vertex = S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0));
    
    match builder.add_vertex(vertex) {
        Ok(_) => println!("ERROR: Should have failed after build"),
        Err(e) => println!("Correctly rejected modification after build: {}", e),
    }
    
    Ok(())
}

/// Demonstrate advanced configuration options
fn advanced_configuration_example() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n6. Advanced Configuration");
    println!("-------------------------");
    
    // Create highly customized builder options
    let base_snap = IdentitySnapFunction::new(1e-12);
    let min_edge_snap = MinEdgeLengthSnapFunction::new(Box::new(base_snap), 1e-7);
    
    let options = Options::new()
        .with_snap_function(min_edge_snap)
        .with_split_crossing_edges(false)
        .with_intersection_tolerance(1e-10)
        .with_simplify_edge_chains(true)
        .with_idempotent(true)
        .with_validate_input(false)  // Skip for performance
        .with_validate_output(true); // But validate results
    
    println!("Configuration summary:");
    println!("  Snap function: MinEdgeLength(Identity, 1e-7)");
    println!("  Split crossing edges: {}", options.split_crossing_edges());
    println!("  Intersection tolerance: {:.1e}", options.intersection_tolerance());
    println!("  Simplify edge chains: {}", options.simplify_edge_chains());
    println!("  Idempotent mode: {}", options.idempotent());
    println!("  Validate input: {}", options.validate_input());
    println!("  Validate output: {}", options.validate_output());
    
    let builder = S2Builder::new(options);
    println!("Builder created with custom configuration");
    
    Ok(())
}