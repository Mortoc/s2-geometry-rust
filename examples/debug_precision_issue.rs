//! Debug script to investigate S2Cap intersection precision issue
//!
//! This script reproduces the exact case mentioned:
//! - Cap center: (0, -1, 0) - face 3 (negative Y axis)  
//! - Testing intersection with face 0 cell
//! - Vertex: (0.5774, -0.5774, -0.5774) has distance vs radius difference of -1.55e-15
//! - Rust says contained (true), C++ says not contained (false)

use s2geometry_rust::{S2Point, S2Cap, S2Cell, S2CellId, S1Angle, math::DVec3};

fn main() {
    println!("=== S2Cap Precision Issue Investigation ===\n");
    
    // 1. Create the cap with center at (0, -1, 0) - face 3
    let cap_center = S2Point::from_normalized(DVec3::new(0.0, -1.0, 0.0));
    println!("Cap center: ({:.10}, {:.10}, {:.10})", 
             cap_center.coords().x, cap_center.coords().y, cap_center.coords().z);
    
    // Let's use a small radius to create a precision boundary case  
    // The issue occurs at very precise boundary conditions
    let cap_radius = S1Angle::from_degrees(54.735610317245346); // Precise angle for the boundary case
    let cap = S2Cap::from_center_angle(&cap_center, cap_radius);
    
    println!("Cap radius: {:.10} degrees", cap_radius.degrees());
    println!("Cap radius squared (chord): {:.15}", cap.radius().length2());
    println!();
    
    // 2. Create face 0 cell and examine its vertices
    let face0_cell = S2Cell::from_face(0).unwrap();
    println!("Face 0 cell vertices:");
    for i in 0..4 {
        let vertex = face0_cell.get_vertex(i);
        let coords = vertex.coords();
        println!("  Vertex {}: ({:.10}, {:.10}, {:.10})", 
                 i, coords.x, coords.y, coords.z);
        
        // Test containment in cap
        let contained = cap.contains(&vertex);
        println!("    Contained in cap: {}", contained);
        
        // Calculate distance from cap center to vertex
        let distance_squared = cap.center().coords().distance_squared(coords);
        let radius_squared = cap.radius().length2();
        let diff = distance_squared - radius_squared;
        
        println!("    Distance² - Radius²: {:.15e}", diff);
        
        // Check if this matches the problematic vertex
        if (coords.x - 0.5774).abs() < 0.001 && 
           (coords.y - (-0.5774)).abs() < 0.001 && 
           (coords.z - (-0.5774)).abs() < 0.001 {
            println!("    *** This matches the problematic vertex! ***");
        }
        println!();
    }
    
    // 3. Let's also test the specific vertex mentioned
    let problem_vertex_raw = DVec3::new(0.5774, -0.5774, -0.5774).normalize();
    let problem_vertex = S2Point::from_normalized(problem_vertex_raw);
    println!("Testing specific problem vertex: ({:.15}, {:.15}, {:.15})", 
             problem_vertex.coords().x, problem_vertex.coords().y, problem_vertex.coords().z);
    
    let contained_specific = cap.contains(&problem_vertex);
    println!("Contained in cap: {}", contained_specific);
    
    let distance_squared_specific = cap.center().coords().distance_squared(problem_vertex.coords());
    let diff_specific = distance_squared_specific - cap.radius().length2();
    println!("Distance² - Radius²: {:.15e}", diff_specific);
    println!();
    
    // 4. Let's examine the coordinate transformation functions
    println!("=== Coordinate Transformation Analysis ===");
    
    // Check face determination
    let face = s2geometry_rust::math::coords::get_face(problem_vertex.coords());
    println!("Face for problem vertex: {}", face);
    
    // Get UV coordinates for face 0
    let (face_calc, u, v) = s2geometry_rust::math::coords::xyz_to_face_uv(problem_vertex.coords());
    println!("Face from xyz_to_face_uv: {}", face_calc);
    println!("UV coordinates: ({:.15}, {:.15})", u, v);
    
    // Convert back to XYZ
    let back_to_xyz = s2geometry_rust::math::coords::face_uv_to_xyz(face_calc, u, v);
    println!("Back to XYZ: ({:.15}, {:.15}, {:.15})", back_to_xyz.x, back_to_xyz.y, back_to_xyz.z);
    
    // Check roundtrip accuracy
    let roundtrip_error = (back_to_xyz - problem_vertex.coords()).length();
    println!("Roundtrip error: {:.15e}", roundtrip_error);
    println!();
    
    // 5. Let's try to find the actual face 0 cell vertex coordinates precisely
    println!("=== Face 0 Cell Vertex Generation Analysis ===");
    
    // A face 0 cell should have UV bounds [-1, 1] x [-1, 1]
    let uv_corners = [
        (-1.0, -1.0), // bottom-left
        (1.0, -1.0),  // bottom-right  
        (1.0, 1.0),   // top-right
        (-1.0, 1.0),  // top-left
    ];
    
    for (i, (u, v)) in uv_corners.iter().enumerate() {
        let xyz = s2geometry_rust::math::coords::face_uv_to_xyz(0, *u, *v);
        let normalized = xyz.normalize();
        println!("UV ({:.1}, {:.1}) -> XYZ ({:.15}, {:.15}, {:.15})", 
                 u, v, normalized.x, normalized.y, normalized.z);
        
        // Check if this matches our problem vertex
        let matches = (normalized.x - 0.5774).abs() < 0.001 && 
                     (normalized.y - (-0.5774)).abs() < 0.001 && 
                     (normalized.z - (-0.5774)).abs() < 0.001;
        if matches {
            println!("    *** This is likely the problematic vertex! ***");
        }
    }
}