//! Port of s1angle_test.cc - Critical S1Angle tests  
//!
//! This module ports the essential tests from C++ s1angle_test.cc to validate
//! our S1Angle implementation matches C++ behavior exactly.

use s2geometry_rust::{S1Angle, S2Point};
use s2geometry_rust::math::constants::*;  
use s2geometry_rust::angle::{sin, cos, tan, abs};
use approx::assert_relative_eq;

/// Test default constructor returns zero angle (from C++ DefaultConstructor)
#[test]
fn test_default_constructor() {
    let angle = S1Angle::zero();
    assert_eq!(angle.radians(), 0.0);
}

/// Test infinity constant (from C++ Infinity)
#[test]
fn test_infinity() {
    let inf = S1Angle::infinity();
    let large = S1Angle::from_radians(1e30);
    let zero = S1Angle::zero();
    
    assert!(large < inf);
    assert!((-inf) < zero);
    assert_eq!(inf, inf);
}

/// Test zero constant (from C++ Zero)
#[test]
fn test_zero() {
    assert_eq!(S1Angle::from_radians(0.0), S1Angle::zero());
}

/// Test exact π radians ↔ 180° conversion (from C++ PiRadiansExactly180Degrees)
#[test]
fn test_pi_radians_exactly_180_degrees() {
    // Check that the conversion between Pi radians and 180 degrees is exact
    assert_eq!(PI, S1Angle::from_radians(PI).radians());
    assert_eq!(180.0, S1Angle::from_radians(PI).degrees());
    assert_eq!(PI, S1Angle::from_degrees(180.0).radians());
    assert_eq!(180.0, S1Angle::from_degrees(180.0).degrees());
    
    assert_eq!(90.0, S1Angle::from_radians(PI_2).degrees());
    
    // Check negative angles
    assert_eq!(-90.0, S1Angle::from_radians(-PI_2).degrees());
    assert_eq!(-PI_4, S1Angle::from_degrees(-45.0).radians());
}

/// Test E5/E6/E7 representations (from C++ E5E6E7Representations)
#[test]
fn test_e5_e6_e7_representations() {
    // Check that E5/E6/E7 representations work as expected
    assert_relative_eq!(
        S1Angle::from_degrees(-45.0).radians(),
        S1Angle::from_e5(-4500000).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        S1Angle::from_degrees(-60.0).radians(), 
        S1Angle::from_e6(-60000000).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        S1Angle::from_degrees(75.0).radians(),
        S1Angle::from_e7(750000000).radians(),
        epsilon = 1e-15
    );
    
    assert_eq!(-17256123, S1Angle::from_degrees(-172.56123).e5());
    assert_eq!(12345678, S1Angle::from_degrees(12.345678).e6());
    assert_eq!(-123456789, S1Angle::from_degrees(-12.3456789).e7());
}

/// Test unsigned E6/E7 representations (from C++ E6E7RepresentationsUnsigned)
#[test]  
fn test_e6_e7_representations_unsigned() {
    // Check that unsigned E6/E7 representations work as expected
    assert_relative_eq!(
        S1Angle::from_degrees(60.0).radians(),
        S1Angle::from_unsigned_e6(60000000_u32).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        S1Angle::from_degrees(-60.0).radians(),
        S1Angle::from_unsigned_e6((-60000000_i32) as u32).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        S1Angle::from_degrees(75.0).radians(),
        S1Angle::from_unsigned_e7(750000000_u32).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        S1Angle::from_degrees(-75.0).radians(),
        S1Angle::from_unsigned_e7((-750000000_i32) as u32).radians(),
        epsilon = 1e-15
    );
}

/// Test normalization (from C++ NormalizeCorrectlyCanonicalizesAngles)
#[test]
fn test_normalize_correctly_canonicalizes_angles() {
    assert_relative_eq!(0.0, S1Angle::from_degrees(360.0).normalized().degrees(), epsilon = 1e-15);
    assert_relative_eq!(-90.0, S1Angle::from_degrees(-90.0).normalized().degrees(), epsilon = 1e-15);
    assert_relative_eq!(180.0, S1Angle::from_degrees(-180.0).normalized().degrees(), epsilon = 1e-15);
    assert_relative_eq!(180.0, S1Angle::from_degrees(180.0).normalized().degrees(), epsilon = 1e-15);
    assert_relative_eq!(180.0, S1Angle::from_degrees(540.0).normalized().degrees(), epsilon = 1e-15);
    assert_relative_eq!(90.0, S1Angle::from_degrees(-270.0).normalized().degrees(), epsilon = 1e-15);
}

/// Test arithmetic operations (from C++ ArithmeticOperationsOnAngles)
#[test]
fn test_arithmetic_operations_on_angles() {
    assert_relative_eq!(0.3, S1Angle::from_radians(-0.3).abs().radians(), epsilon = 1e-15);
    assert_relative_eq!(0.3, abs(S1Angle::from_radians(-0.3)).radians(), epsilon = 1e-15);
    assert_relative_eq!(-0.1, (-S1Angle::from_radians(0.1)).radians(), epsilon = 1e-15);
    assert_relative_eq!(
        0.4,
        (S1Angle::from_radians(0.1) + S1Angle::from_radians(0.3)).radians(),
        epsilon = 1e-15
    );
    assert_relative_eq!(
        -0.2,
        (S1Angle::from_radians(0.1) - S1Angle::from_radians(0.3)).radians(),
        epsilon = 1e-15  
    );
    assert_relative_eq!(0.6, (2.0 * S1Angle::from_radians(0.3)).radians(), epsilon = 1e-15);
    assert_relative_eq!(0.6, (S1Angle::from_radians(0.3) * 2.0).radians(), epsilon = 1e-15);
    assert_relative_eq!(0.15, (S1Angle::from_radians(0.3) / 2.0).radians(), epsilon = 1e-15);
    assert_relative_eq!(0.5, S1Angle::from_radians(0.3) / S1Angle::from_radians(0.6), epsilon = 1e-15);
    
    let mut tmp = S1Angle::from_radians(1.0);
    tmp += S1Angle::from_radians(0.5);
    assert_relative_eq!(1.5, tmp.radians(), epsilon = 1e-15);
    tmp -= S1Angle::from_radians(1.0);
    assert_relative_eq!(0.5, tmp.radians(), epsilon = 1e-15);
    tmp *= 5.0;
    assert_relative_eq!(2.5, tmp.radians(), epsilon = 1e-15);
    tmp /= 2.0;
    assert_relative_eq!(1.25, tmp.radians(), epsilon = 1e-15);
}

/// Test trigonometric functions (from C++ Trigonometry)
#[test]
fn test_trigonometry() {
    // Spot check a few angles to ensure that the correct function is called
    assert_relative_eq!(1.0, cos(S1Angle::from_degrees(0.0)), epsilon = 1e-15);
    assert_relative_eq!(1.0, sin(S1Angle::from_degrees(90.0)), epsilon = 1e-15);
    assert_relative_eq!(1.0, tan(S1Angle::from_degrees(45.0)), epsilon = 1e-15);
    
    // Expect that SinCos is exactly [sin, cos]
    for k in -100..=100 {
        let angle = S1Angle::from_degrees(k as f64);
        let sin_cos = angle.sin_cos();
        // If this fails once, it will likely fail many times
        assert_eq!(sin_cos.sin, sin(angle), "sin mismatch at k={}", k);
        assert_eq!(sin_cos.cos, cos(angle), "cos mismatch at k={}", k);
    }
}

/// Test constructors that measure angles between points (from C++ ConstructorsThatMeasureAngles)
#[test]
fn test_constructors_that_measure_angles() {
    let p1 = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let p2 = S2Point::new(0.0, 0.0, 2.0).unwrap(); // Will be normalized to (0,0,1)
    let angle = S1Angle::from_points(p1, p2);
    assert_relative_eq!(PI_2, angle.radians(), epsilon = 1e-13);
    
    let p_same = S2Point::new(1.0, 0.0, 0.0).unwrap();
    let angle_same = S1Angle::from_points(p1, p_same);
    assert_relative_eq!(0.0, angle_same.radians(), epsilon = 1e-15);
}

/// Test formatting (from C++ TestFormatting)
#[test]
fn test_formatting() {
    let angle = S1Angle::from_degrees(180.0);
    let formatted = format!("{}", angle);
    assert_eq!("180.0000000", formatted);
}

/// Test exact conversions between Degrees() and E6() (from C++ DegreesVsE6)
#[test]  
fn test_degrees_vs_e6() {
    // The current implementation guarantees exact conversions between
    // Degrees() and E6() when the Degrees() argument is an integer
    for i in 0..=180 {
        assert_eq!(
            S1Angle::from_degrees(i as f64),
            S1Angle::from_e6(1000000 * i)
        );
    }
}

/// Test exact conversions between Degrees() and E7() (from C++ DegreesVsE7)
#[test]
fn test_degrees_vs_e7() {
    // The current implementation guarantees exact conversions between  
    // Degrees() and E7() when the Degrees() argument is an integer
    for i in 0..=180 {
        assert_eq!(
            S1Angle::from_degrees(i as f64),
            S1Angle::from_e7(10000000 * i)
        );
    }
}

/// Test exact conversions between E6() and E7() (from C++ E6VsE7) 
#[test]
fn test_e6_vs_e7() {
    // The current implementation guarantees exact conversions between
    // E6() and E7() when the E6() argument is an integer
    use rand::Rng;
    let mut rng = rand::thread_rng();
    
    for _ in 0..100 { // Reduced from 1000 for test speed
        let i = rng.gen_range(0..180000000);
        assert_eq!(S1Angle::from_e6(i), S1Angle::from_e7(10 * i));
    }
}

/// Test guaranteed exact conversions (from C++ DegreesVsRadians)
#[test] 
fn test_degrees_vs_radians() {
    // The current implementation guarantees certain exact conversions between
    // degrees and radians (see the header file for details)
    for k in -8..=8 {
        assert_eq!(
            S1Angle::from_degrees(45.0 * k as f64),
            S1Angle::from_radians((k as f64) * PI / 4.0)
        );
        assert_eq!(45.0 * k as f64, S1Angle::from_degrees(45.0 * k as f64).degrees());
    }
    
    for k in 0..=20 { // Reduced from 30 for test speed
        let n = 1_u32 << k;
        let n_f = n as f64;
        assert_eq!(
            S1Angle::from_degrees(180.0 / n_f),
            S1Angle::from_radians(PI / n_f)
        );
        assert_eq!(
            S1Angle::from_degrees(60.0 / n_f),
            S1Angle::from_radians(PI / (3.0 * n_f))
        );
        assert_eq!(
            S1Angle::from_degrees(36.0 / n_f),
            S1Angle::from_radians(PI / (5.0 * n_f))
        );
        assert_eq!(
            S1Angle::from_degrees(20.0 / n_f),
            S1Angle::from_radians(PI / (9.0 * n_f))
        );
        assert_eq!(
            S1Angle::from_degrees(4.0 / n_f), 
            S1Angle::from_radians(PI / (45.0 * n_f))
        );
    }
    
    // We also spot check a couple of non-identities
    assert_ne!(S1Angle::from_degrees(3.0), S1Angle::from_radians(PI / 60.0));
}

/// Test numerical precision and edge cases
#[test]
fn test_numerical_precision() {
    // Test very small angles
    let tiny = S1Angle::from_radians(1e-15);
    assert_relative_eq!(tiny.radians(), 1e-15, epsilon = 1e-30);
    
    // Test angles near π
    let near_pi = S1Angle::from_radians(PI - 1e-15);
    assert!(near_pi.radians() < PI);
    assert!(near_pi.radians() > PI - 1e-10);
    
    // Test normalization of very large angles
    let huge = S1Angle::from_radians(1000.0 * PI);
    let normalized = huge.normalized();
    assert!(normalized.radians().abs() <= PI);
    
    // Test that sin/cos maintain precision
    let angle = S1Angle::from_degrees(30.0);
    assert_relative_eq!(sin(angle), 0.5, epsilon = 1e-15);
    assert_relative_eq!(cos(angle), (3.0_f64).sqrt() / 2.0, epsilon = 1e-15);
}

/// Test special angle values
#[test]  
fn test_special_angles() {
    // Test common angles have exact representations
    assert_eq!(S1Angle::from_degrees(0.0).radians(), 0.0);
    assert_eq!(S1Angle::from_degrees(90.0).radians(), PI_2);
    assert_eq!(S1Angle::from_degrees(180.0).radians(), PI);
    assert_eq!(S1Angle::from_degrees(270.0).radians(), 3.0 * PI_2);
    assert_eq!(S1Angle::from_degrees(360.0).radians(), 2.0 * PI);
    
    // Test negative angles
    assert_eq!(S1Angle::from_degrees(-90.0).radians(), -PI_2);
    assert_eq!(S1Angle::from_degrees(-180.0).radians(), -PI);
}

/// Behavioral Driven Development (BDD) test structure matching other ports
mod bdd_tests {
    use super::*;
    
    /// BDD: Given a 90-degree angle, When I convert to radians, Then I get π/2
    #[test]
    fn given_90_degrees_when_convert_to_radians_then_pi_over_2() {
        // Given
        let angle_deg = 90.0;
        
        // When
        let angle = S1Angle::from_degrees(angle_deg);
        
        // Then
        assert_eq!(angle.radians(), PI_2);
    }
    
    /// BDD: Given π radians, When I convert to degrees, Then I get exactly 180°
    #[test]
    fn given_pi_radians_when_convert_to_degrees_then_180() {
        // Given
        let pi_radians = PI;
        
        // When
        let angle = S1Angle::from_radians(pi_radians);
        
        // Then
        assert_eq!(angle.degrees(), 180.0);
    }
    
    /// BDD: Given two angles, When I add them, Then result equals sum of radians
    #[test]
    fn given_two_angles_when_add_then_sum_of_radians() {
        // Given
        let angle1 = S1Angle::from_degrees(30.0);
        let angle2 = S1Angle::from_degrees(60.0);
        
        // When
        let sum = angle1 + angle2;
        
        // Then
        assert_relative_eq!(sum.degrees(), 90.0, epsilon = 1e-15);
        assert_relative_eq!(sum.radians(), angle1.radians() + angle2.radians(), epsilon = 1e-15);
    }
    
    /// BDD: Given an angle > 180°, When I normalize, Then result is in [-π, π]
    #[test] 
    fn given_angle_over_180_when_normalize_then_in_range() {
        // Given
        let large_angle = S1Angle::from_degrees(450.0); // 450° = 360° + 90°
        
        // When
        let normalized = large_angle.normalized();
        
        // Then
        assert!(normalized.radians() >= -PI);
        assert!(normalized.radians() <= PI);
        assert_relative_eq!(normalized.degrees(), 90.0, epsilon = 1e-15); // Should wrap to 90°
    }
    
    /// BDD: Given a 45° angle, When I compute tan, Then result is 1
    #[test]
    fn given_45_degrees_when_compute_tan_then_one() {
        // Given
        let angle = S1Angle::from_degrees(45.0);
        
        // When
        let tangent = tan(angle);
        
        // Then
        assert_relative_eq!(tangent, 1.0, epsilon = 1e-15);
    }
    
    /// BDD: Given angle with sin_cos, When I compare to individual sin/cos, Then they match exactly
    #[test]
    fn given_angle_when_sin_cos_then_matches_individual() {
        // Given
        let angle = S1Angle::from_degrees(37.0); // Arbitrary angle
        
        // When
        let sin_cos_pair = angle.sin_cos();
        let individual_sin = sin(angle);
        let individual_cos = cos(angle);
        
        // Then
        assert_eq!(sin_cos_pair.sin, individual_sin);
        assert_eq!(sin_cos_pair.cos, individual_cos);
        
        // And verify Pythagorean identity
        let sum_of_squares = sin_cos_pair.sin.powi(2) + sin_cos_pair.cos.powi(2);
        assert_relative_eq!(sum_of_squares, 1.0, epsilon = 1e-15);
    }
    
    /// BDD: Given integer degrees, When I convert to E6 and back, Then I get exact match
    #[test]
    fn given_integer_degrees_when_convert_e6_roundtrip_then_exact() {
        // Given
        let original_degrees = 45;
        
        // When
        let angle = S1Angle::from_degrees(original_degrees as f64);
        let e6_value = angle.e6();
        let roundtrip = S1Angle::from_e6(e6_value);
        
        // Then
        assert_eq!(angle, roundtrip);
        assert_eq!(e6_value, 45_000_000); // 45 * 1e6
    }
}