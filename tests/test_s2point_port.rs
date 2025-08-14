//! Port of s2point_test.cc - Critical S2Point tests
//!
//! This module ports the essential tests from C++ s2point_test.cc to validate
//! our S2Point implementation matches C++ behavior exactly.

use s2geometry_rust::{S2Point, S2Error};
use s2geometry_rust::math::{DVec3, constants::*};
use approx::{assert_relative_eq, assert_abs_diff_eq};

/// Test basic S2Point operations match C++ behavior
#[test]
fn test_point_basic_operations() {
    // Test from C++ SubtractionWorks
    let mut a = S2Point::from_coords_raw(1.0, 2.0, 3.0);
    let b = S2Point::from_coords_raw(1.0, 1.0, 1.0);
    a -= b;
    
    let expected = S2Point::from_coords_raw(0.0, 1.0, 2.0);
    assert_relative_eq!(a.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(a.y(), expected.y(), epsilon = 1e-10); 
    assert_relative_eq!(a.z(), expected.z(), epsilon = 1e-10);
}

/// Test component-wise division operation
#[test]
fn test_element_wise_division() {
    // Test from C++ ElementWiseDivisionWorks
    let a = S2Point::from_coords_raw(4.0, 8.0, 16.0);
    let b = S2Point::from_coords_raw(2.0, 2.0, 2.0);
    let result = a.div_components(&b);
    
    let expected = S2Point::from_coords_raw(2.0, 4.0, 8.0);
    assert_relative_eq!(result.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(result.y(), expected.y(), epsilon = 1e-10);
    assert_relative_eq!(result.z(), expected.z(), epsilon = 1e-10);
}

/// Test component-wise square root operation
#[test]
fn test_sqrt_works() {
    // Test from C++ SqrtWorks
    let a = S2Point::from_coords_raw(4.0, 9.0, 16.0);
    let result = a.sqrt();
    
    let expected = S2Point::from_coords_raw(2.0, 3.0, 4.0);
    assert_relative_eq!(result.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(result.y(), expected.y(), epsilon = 1e-10);
    assert_relative_eq!(result.z(), expected.z(), epsilon = 1e-10);
}

/// Test component-wise floor operation
#[test] 
fn test_floor_works() {
    // Test from C++ FloorWorks
    let a = S2Point::from_coords_raw(1.4, 1.5, 1.6);
    let result = a.floor();
    
    let expected = S2Point::from_coords_raw(1.0, 1.0, 1.0);
    assert_relative_eq!(result.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(result.y(), expected.y(), epsilon = 1e-10);
    assert_relative_eq!(result.z(), expected.z(), epsilon = 1e-10);
}

/// Test component-wise ceil operation
#[test]
fn test_ceil_works() {
    // Test from C++ CeilWorks
    let a = S2Point::from_coords_raw(1.4, 1.5, 1.6);
    let result = a.ceil();
    
    let expected = S2Point::from_coords_raw(2.0, 2.0, 2.0);
    assert_relative_eq!(result.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(result.y(), expected.y(), epsilon = 1e-10);
    assert_relative_eq!(result.z(), expected.z(), epsilon = 1e-10);
}

/// Test component-wise round operation
#[test]
fn test_fround_works() {
    // Test from C++ FRoundWorks  
    let a = S2Point::from_coords_raw(1.4, 1.5, 1.6);
    let result = a.fround();
    
    let expected = S2Point::from_coords_raw(1.0, 2.0, 2.0);
    assert_relative_eq!(result.x(), expected.x(), epsilon = 1e-10);
    assert_relative_eq!(result.y(), expected.y(), epsilon = 1e-10);
    assert_relative_eq!(result.z(), expected.z(), epsilon = 1e-10);
}

/// Test that S2Point maintains compatibility with Vector3_d size
#[test]  
fn test_is_a_vector() {
    // Test from C++ IsAVector - ensure S2Point has same memory layout as DVec3
    assert_eq!(std::mem::size_of::<S2Point>(), std::mem::size_of::<DVec3>());
    assert_eq!(std::mem::align_of::<S2Point>(), std::mem::align_of::<DVec3>());
}

/// Test S2Point normalization behavior
#[test]
fn test_point_normalization() {
    // Test normalization preserves direction and ensures unit length
    let coords = DVec3::new(3.0, 4.0, 0.0);
    let point = S2Point::from_vec3(coords).unwrap();
    
    // Should be normalized to unit length
    assert_relative_eq!(point.coords().length(), 1.0, epsilon = 1e-15);
    
    // Should preserve ratios
    let expected_x = 3.0 / 5.0; // 3/sqrt(9+16) = 3/5
    let expected_y = 4.0 / 5.0; // 4/sqrt(9+16) = 4/5
    assert_relative_eq!(point.x(), expected_x, epsilon = 1e-15);
    assert_relative_eq!(point.y(), expected_y, epsilon = 1e-15);
    assert_relative_eq!(point.z(), 0.0, epsilon = 1e-15);
}

/// Test angle calculation between points
#[test]
fn test_point_angles() {
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let p2 = S2Point::new(0.0, 1.0, 0.0).unwrap(); 
    let p3 = S2Point::new(0.0, 0.0, 1.0).unwrap();
    
    // Test orthogonal points have π/2 angle
    assert_relative_eq!(p1.angle(&p2), PI_2, epsilon = 1e-15);
    assert_relative_eq!(p1.angle(&p3), PI_2, epsilon = 1e-15);
    assert_relative_eq!(p2.angle(&p3), PI_2, epsilon = 1e-15);
    
    // Test point with itself has zero angle
    assert_relative_eq!(p1.angle(&p1), 0.0, epsilon = 1e-15);
    
    // Test antipodal points have π angle
    let p1_neg = -p1;
    assert_relative_eq!(p1.angle(&p1_neg), PI, epsilon = 1e-15);
}

/// Test spherical interpolation (slerp)
#[test]
fn test_point_interpolation() {
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let p2 = S2Point::new(0.0, 1.0, 0.0).unwrap();
    
    // Test endpoints
    let interp0 = p1.interpolate(&p2, 0.0);
    assert_relative_eq!(interp0.x(), p1.x(), epsilon = 1e-15);
    assert_relative_eq!(interp0.y(), p1.y(), epsilon = 1e-15);
    assert_relative_eq!(interp0.z(), p1.z(), epsilon = 1e-15);
    
    let interp1 = p1.interpolate(&p2, 1.0);
    assert_relative_eq!(interp1.x(), p2.x(), epsilon = 1e-15);
    assert_relative_eq!(interp1.y(), p2.y(), epsilon = 1e-15);
    assert_relative_eq!(interp1.z(), p2.z(), epsilon = 1e-15);
    
    // Test midpoint interpolation
    let mid = p1.interpolate(&p2, 0.5);
    assert_relative_eq!(mid.coords().length(), 1.0, epsilon = 1e-15);
    
    // Midpoint should be equidistant from both endpoints
    let dist1 = mid.angle(&p1);
    let dist2 = mid.angle(&p2);
    assert_relative_eq!(dist1, dist2, epsilon = 1e-15);
}

/// Test distance calculations  
#[test]
fn test_point_distances() {
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let p2 = S2Point::new(0.0, 1.0, 0.0).unwrap();
    
    // Test distance consistency with angle
    let angle = p1.angle(&p2);
    let chord_length = p1.distance(&p2);
    
    // For unit sphere: chord_length = 2 * sin(angle/2)
    let expected_chord = 2.0 * (angle / 2.0).sin();
    assert_relative_eq!(chord_length, expected_chord, epsilon = 1e-15);
    
    // Test distance squared
    let dist_sq = p1.distance_squared(&p2);
    assert_relative_eq!(dist_sq, chord_length * chord_length, epsilon = 1e-15);
}

/// Test cross product operations
#[test] 
fn test_cross_product() {
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let p2 = S2Point::new(0.0, 1.0, 0.0).unwrap();
    
    let cross = p1.cross(&p2);
    
    // Cross product of x and y unit vectors should be z direction
    assert_relative_eq!(cross.x, 0.0, epsilon = 1e-15);
    assert_relative_eq!(cross.y, 0.0, epsilon = 1e-15);
    assert_relative_eq!(cross.z, 1.0, epsilon = 1e-15);
    
    // Cross product magnitude should equal sin of angle
    let angle = p1.angle(&p2);
    let expected_magnitude = angle.sin();
    assert_relative_eq!(cross.length(), expected_magnitude, epsilon = 1e-15);
}

/// Test arithmetic operations
#[test]
fn test_arithmetic_operations() {
    let p1 = S2Point::from_coords_raw(1.0, 2.0, 3.0);
    let p2 = S2Point::from_coords_raw(0.5, 1.0, 1.5);
    
    // Test addition (raw vector addition)
    let sum = p1 + p2;
    assert_relative_eq!(sum.x(), 1.5, epsilon = 1e-15);
    assert_relative_eq!(sum.y(), 3.0, epsilon = 1e-15);
    assert_relative_eq!(sum.z(), 4.5, epsilon = 1e-15);
    
    // Test subtraction
    let diff = p1 - p2;  
    assert_relative_eq!(diff.x(), 0.5, epsilon = 1e-15);
    assert_relative_eq!(diff.y(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(diff.z(), 1.5, epsilon = 1e-15);
    
    // Test scalar multiplication
    let scaled = p1 * 2.0;
    assert_relative_eq!(scaled.x(), 2.0, epsilon = 1e-15);
    assert_relative_eq!(scaled.y(), 4.0, epsilon = 1e-15);
    assert_relative_eq!(scaled.z(), 6.0, epsilon = 1e-15);
    
    // Test scalar division
    let divided = p1 / 2.0;
    assert_relative_eq!(divided.x(), 0.5, epsilon = 1e-15);
    assert_relative_eq!(divided.y(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(divided.z(), 1.5, epsilon = 1e-15);
    
    // Test negation
    let negated = -p1;
    assert_relative_eq!(negated.x(), -p1.x(), epsilon = 1e-15);
    assert_relative_eq!(negated.y(), -p1.y(), epsilon = 1e-15); 
    assert_relative_eq!(negated.z(), -p1.z(), epsilon = 1e-15);
}

/// Test error conditions
#[test]
fn test_error_conditions() {
    // Test zero vector error
    let result = S2Point::from_vec3(DVec3::new(0.0, 0.0, 0.0));
    assert!(result.is_err());
    
    match result {
        Err(S2Error::InvalidPoint { reason: _ }) => {}, // Expected
        _ => panic!("Expected InvalidPoint error"),
    }
    
    // Test very small vector (should normalize successfully)
    let tiny = DVec3::new(1e-100, 0.0, 0.0);
    let point = S2Point::from_vec3(tiny).unwrap();
    assert_relative_eq!(point.x(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(point.coords().length(), 1.0, epsilon = 1e-15);
}

/// Test hash function properties (ported from C++ HashSpreads test concept)  
#[test]
fn test_hash_consistency() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    
    let point = S2Point::new(1.0, 2.0, 3.0).unwrap();
    
    // Test hash consistency
    let mut hasher1 = DefaultHasher::new();
    let mut hasher2 = DefaultHasher::new();
    
    point.hash(&mut hasher1);  
    point.hash(&mut hasher2);
    
    assert_eq!(hasher1.finish(), hasher2.finish());
    
    // Test that different points have different hashes (with high probability)
    let point2 = S2Point::new(3.0, 2.0, 1.0).unwrap();
    let mut hasher3 = DefaultHasher::new();
    point2.hash(&mut hasher3);
    
    // Should be different (though not guaranteed, very likely)
    assert_ne!(hasher1.finish(), hasher3.finish());
}

/// Test numerical precision edge cases
#[test] 
fn test_numerical_precision() {
    // Test points very close to each other
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let tiny_offset = f64::EPSILON;
    let p2 = S2Point::new(1.0 + tiny_offset, 0.0, 0.0).unwrap();
    
    // Should still be properly normalized
    assert_relative_eq!(p1.coords().length(), 1.0, epsilon = 1e-15);
    assert_relative_eq!(p2.coords().length(), 1.0, epsilon = 1e-15);
    
    // Test that very small angles are handled correctly
    let angle = p1.angle(&p2);
    assert!(angle >= 0.0);
    assert!(angle < 1e-10); // Should be very small but positive
    
    // Test dot product precision
    let dot = p1.dot(&p2);
    assert!(dot <= 1.0 + f64::EPSILON); // Should not exceed 1 due to normalization
    assert!(dot >= 1.0 - 1e-10); // Should be very close to 1
}

/// Behavioral Driven Development (BDD) test structure matching other ports
mod bdd_tests {
    use super::*;
    
    /// BDD: Given two orthogonal unit vectors, When I compute their cross product, Then I get a unit vector perpendicular to both
    #[test]
    fn given_orthogonal_vectors_when_cross_product_then_perpendicular_unit_vector() {
        // Given
        let x_axis = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let y_axis = S2Point::new(0.0, 1.0, 0.0).unwrap();
        
        // When  
        let cross = x_axis.cross(&y_axis);
        
        // Then
        assert_relative_eq!(cross.length(), 1.0, epsilon = 1e-15);
        assert_relative_eq!(cross.dot(x_axis.coords()), 0.0, epsilon = 1e-15);
        assert_relative_eq!(cross.dot(y_axis.coords()), 0.0, epsilon = 1e-15);
        assert_relative_eq!(cross.z, 1.0, epsilon = 1e-15);
    }
    
    /// BDD: Given identical points, When I compute angle between them, Then angle should be zero
    #[test]
    fn given_identical_points_when_compute_angle_then_zero() {
        // Given
        let point = S2Point::new(0.6, 0.8, 0.0).unwrap();
        
        // When
        let angle = point.angle(&point);
        
        // Then  
        assert_abs_diff_eq!(angle, 0.0, epsilon = 1e-15);
    }
    
    /// BDD: Given antipodal points, When I compute angle between them, Then angle should be π
    #[test]
    fn given_antipodal_points_when_compute_angle_then_pi() {
        // Given
        let point = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let antipodal = -point;
        
        // When
        let angle = point.angle(&antipodal);
        
        // Then
        assert_relative_eq!(angle, PI, epsilon = 1e-15);
    }
    
    /// BDD: Given a point and normalization, When I check length, Then it should be exactly 1
    #[test]
    fn given_point_when_normalized_then_unit_length() {
        // Given
        let unnormalized = DVec3::new(3.0, 4.0, 12.0); // Length = 13
        
        // When  
        let point = S2Point::from_vec3(unnormalized).unwrap();
        
        // Then
        assert_relative_eq!(point.coords().length(), 1.0, epsilon = 1e-15);
        
        // And direction should be preserved
        let expected_dir = unnormalized.normalize();
        assert_relative_eq!(point.coords().x, expected_dir.x, epsilon = 1e-15);
        assert_relative_eq!(point.coords().y, expected_dir.y, epsilon = 1e-15);  
        assert_relative_eq!(point.coords().z, expected_dir.z, epsilon = 1e-15);
    }
}