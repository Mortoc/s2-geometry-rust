/// Example showing WASM-compatible S2 geometry operations
/// 
/// Compile for WASM with:
/// cargo build --target wasm32-unknown-unknown --features wasm --no-default-features --example wasm_example

use s2geometry_rust::{S1Angle, S2Point, math::DVec3};

/// Demonstrate basic S2 operations that work in WASM
fn main() {
    // Create points on the sphere
    let north_pole = S2Point::new(0.0, 0.0, 1.0).unwrap();
    let equator_point = S2Point::new(1.0, 0.0, 0.0).unwrap();
    
    // Calculate angle between points
    let angle = S1Angle::from_points(north_pole, equator_point);
    
    // This would print π/2 radians (90 degrees)
    println!("Angle between north pole and equator: {:.6} radians", angle.radians());
    println!("Angle in degrees: {:.2}°", angle.degrees());
    
    // Demonstrate coordinate transformations
    let coords = north_pole.coords();
    println!("North pole coordinates: ({:.2}, {:.2}, {:.2})", coords.x, coords.y, coords.z);
    
    // Demonstrate exact arithmetic predicates
    use s2geometry_rust::math::predicates::exact_sign;
    let a = DVec3::new(1.0, 0.0, 0.0);
    let b = DVec3::new(0.0, 1.0, 0.0);
    let c = DVec3::new(0.0, 0.0, 1.0);
    
    let orientation = exact_sign(a, b, c);
    println!("Orientation test result: {}", orientation); // Should be 1 (positive)
}