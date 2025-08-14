//! Integration tests for glam f64 math types
//! 
//! This demonstrates that our math foundation using Bevy's glam library
//! provides the f64 precision and operations we need for S2 geometry.

use s2geometry_rust::math::{DVec3, DMat3, DQuat};
use approx::assert_relative_eq;

#[test]
fn test_glam_dvec3_operations() {
    // Test basic DVec3 operations that S2 geometry will require
    let v1 = DVec3::new(1.0, 0.0, 0.0);
    let v2 = DVec3::new(0.0, 1.0, 0.0);
    
    // Unit vector properties
    assert_relative_eq!(v1.length(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(v2.length(), 1.0, epsilon = 1e-15);
    
    // Cross product (essential for S2 orientation calculations)
    let cross = v1.cross(v2);
    let expected_cross = DVec3::new(0.0, 0.0, 1.0);
    assert_relative_eq!(cross.x, expected_cross.x, epsilon = 1e-15);
    assert_relative_eq!(cross.y, expected_cross.y, epsilon = 1e-15);
    assert_relative_eq!(cross.z, expected_cross.z, epsilon = 1e-15);
    
    // Dot product (essential for S2 distance calculations)
    let dot = v1.dot(v2);
    assert_relative_eq!(dot, 0.0, epsilon = 1e-15);
}

#[test]
fn test_glam_dmat3_rotations() {
    // Test 3x3 matrix rotations that S2 coordinate transforms will use
    let axis = DVec3::new(0.0, 0.0, 1.0); // Z-axis
    let angle = std::f64::consts::FRAC_PI_2; // 90 degrees
    
    // Create rotation matrix around Z-axis
    let rotation = DMat3::from_axis_angle(axis, angle);
    
    // Test rotation of X-axis vector to Y-axis
    let x_axis = DVec3::new(1.0, 0.0, 0.0);
    let rotated = rotation * x_axis;
    
    let expected_y = DVec3::new(0.0, 1.0, 0.0);
    assert_relative_eq!(rotated.x, expected_y.x, epsilon = 1e-15);
    assert_relative_eq!(rotated.y, expected_y.y, epsilon = 1e-15);
    assert_relative_eq!(rotated.z, expected_y.z, epsilon = 1e-15);
}

#[test]
fn test_glam_dquat_unit_quaternions() {
    // Test quaternion operations for S2 rotations
    let q1 = DQuat::IDENTITY;
    let q2 = DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2);
    
    // Quaternions should be normalized
    assert_relative_eq!(q1.length(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(q2.length(), 1.0, epsilon = 1e-15);
    
    // Test quaternion multiplication for composition
    let composed = q1 * q2;
    assert_relative_eq!(composed.length(), 1.0, epsilon = 1e-15);
}

#[test]
fn test_spherical_geometry_precision() {
    // Test that glam f64 provides adequate precision for spherical calculations
    let lat = std::f64::consts::FRAC_PI_4; // 45 degrees
    let lng = std::f64::consts::FRAC_PI_3; // 60 degrees
    
    // Convert spherical coordinates to Cartesian (this will be core S2 functionality)
    let x = lat.cos() * lng.cos();
    let y = lat.cos() * lng.sin();
    let z = lat.sin();
    
    let point = DVec3::new(x, y, z);
    let length = point.length();
    
    // Should be very close to 1.0 for a point on the unit sphere
    assert_relative_eq!(length, 1.0, epsilon = 1e-15);
    
    // Test that we maintain precision in normalization
    let normalized = point.normalize();
    assert_relative_eq!(normalized.length(), 1.0, epsilon = 1e-16);
}

#[test]
fn test_high_precision_geometric_predicates() {
    // Test the precision needed for robust geometric predicates
    let a = DVec3::new(1.0, 0.0, 0.0);
    let b = DVec3::new(0.0, 1.0, 0.0);
    let c = DVec3::new(0.0, 0.0, 1.0);
    
    // Triple scalar product for orientation testing
    let orientation = a.cross(b).dot(c);
    
    // Should be exactly 1.0 for this right-handed coordinate system
    assert_relative_eq!(orientation, 1.0, epsilon = 1e-15);
    
    // Test with nearly collinear points (this will stress exact arithmetic needs)
    let epsilon = 1e-14;
    let nearly_a = DVec3::new(1.0, epsilon, 0.0).normalize();
    let nearly_b = DVec3::new(1.0, -epsilon, 0.0).normalize();
    let nearly_c = DVec3::new(1.0, 0.0, epsilon).normalize();
    
    let small_orientation = nearly_a.cross(nearly_b).dot(nearly_c);
    
    // This should be very small but non-zero - demonstrates where we need exact arithmetic
    assert!(small_orientation.abs() > 0.0);
    assert!(small_orientation.abs() < 1e-12);
}