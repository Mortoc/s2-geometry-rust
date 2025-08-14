//! Demo of S2ShapeIndex functionality
//!
//! This example demonstrates how to use the S2ShapeIndex system to create
//! and query spatial indexes for various types of geometric shapes.

use s2geometry_rust::{
    S2Point, MutableS2ShapeIndex, S2PointShape, S2MultiPointShape, S2PolylineShape,
    S2LoopShape, S2Polyline, S2Loop, S2ShapeIndex, S2ShapeIndexIterator, 
    MutableS2ShapeIndexTrait, math::DVec3
};
use std::sync::Arc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("S2 Shape Index Demo");
    println!("==================");

    // Create a mutable shape index
    let mut index = MutableS2ShapeIndex::new();
    
    // 1. Add some point shapes
    println!("\n1. Adding point shapes...");
    
    let point1 = S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0));
    let point_shape1 = Arc::new(S2PointShape::new(point1));
    let point_id1 = index.add_shape(point_shape1);
    println!("Added single point shape with ID: {}", point_id1);
    
    let points = vec![
        S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
        S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),
        S2Point::from_normalized(DVec3::new(1.0, 1.0, 0.0).normalize()),
    ];
    let multi_point_shape = Arc::new(S2MultiPointShape::new(points));
    let multi_point_id = index.add_shape(multi_point_shape);
    println!("Added multi-point shape with ID: {}", multi_point_id);
    
    // 2. Add a polyline shape
    println!("\n2. Adding polyline shape...");
    
    let polyline_vertices = vec![
        S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),
        S2Point::from_normalized(DVec3::new(0.8, 0.6, 0.0).normalize()),
        S2Point::from_normalized(DVec3::new(0.6, 0.8, 0.0).normalize()),
        S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),
    ];
    let polyline = Arc::new(S2Polyline::new(polyline_vertices)?);
    let polyline_shape = Arc::new(S2PolylineShape::new(polyline));
    let polyline_id = index.add_shape(polyline_shape);
    println!("Added polyline shape with ID: {}", polyline_id);
    
    // 3. Add a polygon shape (loop)
    println!("\n3. Adding polygon shape...");
    
    let loop_vertices = vec![
        S2Point::from_normalized(DVec3::new(0.0, 0.0, 1.0)),  // North pole
        S2Point::from_normalized(DVec3::new(1.0, 0.0, 0.0)),  // Point on equator
        S2Point::from_normalized(DVec3::new(0.0, 1.0, 0.0)),  // Another point on equator
    ];
    let loop_shape = Arc::new(S2LoopShape::new(Arc::new(S2Loop::new(loop_vertices)?)));
    let loop_id = index.add_shape(loop_shape);
    println!("Added polygon (loop) shape with ID: {}", loop_id);
    
    // 4. Query the index
    println!("\n4. Querying the index...");
    
    println!("Total number of shape IDs: {}", index.num_shape_ids());
    println!("Index memory usage: {} bytes", index.space_used());
    
    // List all shapes in the index
    for shape_id in 0..index.num_shape_ids() {
        if let Some(shape) = index.shape(shape_id) {
            println!("Shape ID {}: dimension={:?}, edges={}", 
                     shape_id, shape.dimension(), shape.num_edges());
        }
    }
    
    // 5. Create an iterator and traverse the index
    println!("\n5. Traversing the index...");
    
    let mut iter = index.iter();
    let mut cell_count = 0;
    
    while !iter.done() {
        let cell_id = iter.cell_id();
        let cell = iter.cell();
        
        println!("Cell {:?}: {} clipped shapes", cell_id, cell.num_clipped());
        
        // List shapes in this cell
        for i in 0..cell.num_clipped() {
            if let Some(clipped) = cell.clipped(i) {
                println!("  - Shape ID {} with {} edges", 
                         clipped.shape_id(), clipped.num_edges());
            }
        }
        
        cell_count += 1;
        iter.next();
    }
    
    println!("Total index cells: {}", cell_count);
    
    // 6. Remove a shape
    println!("\n6. Removing a shape...");
    
    let removed_shape = index.remove_shape(point_id1);
    if removed_shape.is_some() {
        println!("Successfully removed shape with ID: {}", point_id1);
    }
    
    println!("Remaining shapes: {}", index.len());
    
    println!("\nDemo completed successfully!");
    
    Ok(())
}