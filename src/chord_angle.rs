//! S1ChordAngle - Efficient chord angle representation
//!
//! S1ChordAngle represents the angle subtended by a chord (i.e., the straight
//! line segment connecting two points on the sphere). Its representation
//! makes it very efficient for computing and comparing distances, but unlike
//! S1Angle it is only capable of representing angles between 0 and Pi radians.

use crate::angle::S1Angle;
use crate::point::S2Point;
use std::fmt;

/// An angle measurement using squared chord length for efficient comparison.
///
/// S1ChordAngle provides an efficient representation for angles that allows
/// fast comparison operations without trigonometric function calls.
/// The internal representation is the squared chord length, which ranges
/// from 0 to 4.
///
/// # Mathematical Background
///
/// For an angle θ, the chord length is 2*sin(θ/2), so the squared chord
/// length is 4*sin²(θ/2). This provides:
/// - θ = 0° → length² = 0
/// - θ = 90° → length² = 2  
/// - θ = 180° → length² = 4
///
/// # Example
/// ```rust,ignore
/// use s2geometry_rust::S1ChordAngle;
///
/// let angle = S1ChordAngle::from_degrees(90.0);
/// assert_eq!(angle.length2(), 2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct S1ChordAngle {
    length2: f64,
}

impl S1ChordAngle {
    /// Maximum relative error when summing two S1ChordAngles together.
    /// From C++: kRelativeSumError = 2.02 * DBL_EPSILON
    pub const RELATIVE_SUM_ERROR: f64 = 2.02 * f64::EPSILON;

    /// Zero chord angle constant
    #[inline]
    pub const fn zero() -> Self {
        Self { length2: 0.0 }
    }

    /// Right angle (90 degrees) constant
    #[inline]
    pub const fn right() -> Self {
        Self { length2: 2.0 }
    }

    /// Straight angle (180 degrees) constant - maximum finite chord angle
    #[inline]
    pub const fn straight() -> Self {
        Self { length2: 4.0 }
    }

    /// Infinity constant - larger than any finite chord angle
    #[inline]
    pub const fn infinity() -> Self {
        Self { length2: f64::INFINITY }
    }

    /// Negative constant - smaller than zero (special value)
    #[inline]
    pub const fn negative() -> Self {
        Self { length2: -1.0 }
    }

    /// Create S1ChordAngle from squared chord length
    /// The argument is automatically clamped to [0, 4] to handle roundoff errors
    #[inline]
    pub const fn from_length2(length2: f64) -> Self {
        Self { 
            length2: if length2 > 4.0 { 4.0 } else { length2 }
        }
    }

    /// Create S1ChordAngle from squared length (alias for from_length2)
    #[inline]
    pub const fn from_length_squared(length_squared: f64) -> Self {
        Self::from_length2(length_squared)
    }

    /// Create S1ChordAngle from two points on the sphere
    pub fn from_points(x: S2Point, y: S2Point) -> Self {
        let diff = x.coords() - y.coords();
        Self::from_length2(diff.length_squared().min(4.0))
    }

    /// Create S1ChordAngle between two points on the sphere (alias for from_points)
    #[inline]
    pub fn between_points(x: &S2Point, y: &S2Point) -> Self {
        Self::from_points(*x, *y)
    }

    /// Create S1ChordAngle from S1Angle
    pub fn from_angle(angle: S1Angle) -> Self {
        if angle.radians() < 0.0 {
            Self::negative()
        } else if angle.radians() == f64::INFINITY {
            Self::infinity()
        } else {
            // The chord length is 2 * sin(angle / 2)
            let half_angle = (0.5 * angle.radians().min(std::f64::consts::PI)).sin();
            let length = 2.0 * half_angle;
            Self { length2: length * length }
        }
    }

    /// Create S1ChordAngle from radians
    #[inline]
    pub fn from_radians(radians: f64) -> Self {
        Self::from_angle(S1Angle::from_radians(radians))
    }

    /// Create S1ChordAngle from degrees
    #[inline]
    pub fn from_degrees(degrees: f64) -> Self {
        Self::from_angle(S1Angle::from_degrees(degrees))
    }

    /// Create S1ChordAngle from E5 representation
    #[inline]
    pub fn from_e5(e5: i32) -> Self {
        Self::from_angle(S1Angle::from_e5(e5))
    }

    /// Create S1ChordAngle from E6 representation
    #[inline]
    pub fn from_e6(e6: i32) -> Self {
        Self::from_angle(S1Angle::from_e6(e6))
    }

    /// Create S1ChordAngle from E7 representation
    #[inline]
    pub fn from_e7(e7: i32) -> Self {
        Self::from_angle(S1Angle::from_e7(e7))
    }

    /// Fast upper bound construction from S1Angle
    /// Returns an S1ChordAngle that is guaranteed to be >= the input angle
    /// Accurate to within 1% for distances up to ~3100km on Earth
    #[inline]
    pub fn fast_upper_bound_from(angle: S1Angle) -> Self {
        let radians = angle.radians();
        Self::from_length2(radians * radians)
    }

    /// Convert to S1Angle
    pub fn to_angle(self) -> S1Angle {
        if self.is_negative() {
            S1Angle::from_radians(-1.0)
        } else if self.is_infinity() {
            S1Angle::infinity()
        } else {
            // Angle = 2 * asin(sqrt(length2) / 2)
            S1Angle::from_radians(2.0 * (0.5 * self.length2.sqrt()).asin())
        }
    }

    /// Get angle in radians (convenience method)
    #[inline]
    pub fn radians(self) -> f64 {
        self.to_angle().radians()
    }

    /// Get angle in degrees (convenience method)  
    #[inline]
    pub fn degrees(self) -> f64 {
        self.to_angle().degrees()
    }

    /// Get E5 representation (convenience method)
    #[inline]
    pub fn e5(self) -> i32 {
        self.to_angle().e5()
    }

    /// Get E6 representation (convenience method)
    #[inline]
    pub fn e6(self) -> i32 {
        self.to_angle().e6()
    }

    /// Get E7 representation (convenience method)
    #[inline]
    pub fn e7(self) -> i32 {
        self.to_angle().e7()
    }

    /// Get the squared chord length
    #[inline]
    pub fn length2(self) -> f64 {
        self.length2
    }

    /// Test if this is the zero angle
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.length2 == 0.0
    }

    /// Test if this is negative (special value)
    #[inline]  
    pub const fn is_negative(self) -> bool {
        self.length2 < 0.0
    }

    /// Test if this is infinity
    #[inline]
    pub const fn is_infinity(self) -> bool {
        self.length2 == f64::INFINITY
    }

    /// Test if this is a special value (negative or infinity)
    #[inline]
    pub const fn is_special(self) -> bool {
        self.is_negative() || self.is_infinity()
    }

    /// Test if this is a valid S1ChordAngle
    #[inline]
    pub const fn is_valid(self) -> bool {
        (self.length2 >= 0.0 && self.length2 <= 4.0) || self.is_special()
    }

    /// Get the next larger representable S1ChordAngle
    pub fn successor(self) -> Self {
        if self.length2 >= 4.0 {
            Self::infinity()
        } else if self.length2 < 0.0 {
            Self::zero()
        } else {
            Self { length2: next_after(self.length2, 10.0) }
        }
    }

    /// Get the next smaller representable S1ChordAngle  
    pub fn predecessor(self) -> Self {
        if self.length2 <= 0.0 {
            Self::negative()
        } else if self.length2 > 4.0 {
            Self::straight()
        } else {
            Self { length2: next_after(self.length2, -10.0) }
        }
    }

    /// Add error bound to this S1ChordAngle
    pub fn plus_error(self, error: f64) -> Self {
        if self.is_special() {
            self
        } else {
            Self::from_length2((self.length2 + error).clamp(0.0, 4.0))
        }
    }

    /// Get maximum error for S2Point constructor
    pub fn get_s2_point_constructor_max_error(self) -> f64 {
        // From C++: relative error of 2.5 * DBL_EPSILON when computing squared distance,
        // plus relative error of 2 * DBL_EPSILON and absolute error of 16 * DBL_EPSILON²
        4.5 * f64::EPSILON * self.length2 + 16.0 * f64::EPSILON * f64::EPSILON
    }

    /// Get maximum error for S1Angle constructor
    pub fn get_s1_angle_constructor_max_error(self) -> f64 {
        // From C++: sin() and multiply each have relative error of 0.5 * DBL_EPSILON
        1.5 * f64::EPSILON * self.length2
    }

    /// Return the smaller of two chord angles
    pub fn min(self, other: Self) -> Self {
        if self.length2 <= other.length2 {
            self
        } else {
            other
        }
    }

    /// Return the larger of two chord angles  
    pub fn max(self, other: Self) -> Self {
        if self.length2 >= other.length2 {
            self
        } else {
            other
        }
    }

    /// Get sine and cosine pair efficiently
    pub fn sin_cos(self) -> crate::angle::SinCosPair {
        crate::angle::SinCosPair {
            sin: sin(self),
            cos: cos(self),
        }
    }
}

/// Implement Default to return zero angle (matches C++ default constructor)
impl Default for S1ChordAngle {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

/// Comparison operations based on squared chord length
impl PartialOrd for S1ChordAngle {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.length2.partial_cmp(&other.length2)
    }
}

/// Arithmetic operations (addition and subtraction only)
impl std::ops::Add for S1ChordAngle {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        debug_assert!(!self.is_special() && !other.is_special());
        
        let a2 = self.length2;
        let b2 = other.length2;
        
        // Optimization for zero
        if b2 == 0.0 {
            return self;
        }
        
        // Clamp to straight angle
        if a2 + b2 >= 4.0 {
            return Self::straight();
        }
        
        // Use the formula from C++ implementation
        let x = a2 * (1.0 - 0.25 * b2);
        let y = b2 * (1.0 - 0.25 * a2);
        Self::from_length2((x + y + 2.0 * (x * y).sqrt()).min(4.0))
    }
}

impl std::ops::Sub for S1ChordAngle {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        debug_assert!(!self.is_special() && !other.is_special());
        
        let a2 = self.length2;
        let b2 = other.length2;
        
        if b2 == 0.0 {
            return self;
        }
        if a2 <= b2 {
            return Self::zero();
        }
        
        let x = a2 * (1.0 - 0.25 * b2);
        let y = b2 * (1.0 - 0.25 * a2);
        
        // Use two square roots to avoid cancellation error
        let c = (x.sqrt() - y.sqrt()).max(0.0);
        Self { length2: c * c }
    }
}

impl std::ops::AddAssign for S1ChordAngle {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::ops::SubAssign for S1ChordAngle {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

/// Trigonometric functions - more accurate than converting to S1Angle first
pub fn sin(angle: S1ChordAngle) -> f64 {
    debug_assert!(!angle.is_special());
    sin2(angle).sqrt()
}

pub fn cos(angle: S1ChordAngle) -> f64 {
    debug_assert!(!angle.is_special());
    1.0 - 0.5 * angle.length2()
}

pub fn tan(angle: S1ChordAngle) -> f64 {
    sin(angle) / cos(angle)
}

/// Compute sin²(angle) efficiently
pub fn sin2(angle: S1ChordAngle) -> f64 {
    debug_assert!(!angle.is_special());
    // sin²(2A) = 4 * sin²(A) * cos²(A) = 4 * sin²(A) * (1 - sin²(A))
    // where sin²(A) = length2/4, so:
    // sin²(2A) = length2 * (1 - length2/4)
    angle.length2() * (1.0 - 0.25 * angle.length2())
}

/// Helper function for next representable floating-point value
fn next_after(x: f64, direction: f64) -> f64 {
    // Simple implementation - in production, use proper nextafter function
    if x == direction {
        x
    } else if direction > x {
        if x == f64::MAX {
            f64::INFINITY
        } else {
            let bits = x.to_bits();
            f64::from_bits(bits + 1)
        }
    } else {
        if x == f64::MIN {
            -f64::INFINITY  
        } else {
            let bits = x.to_bits();
            f64::from_bits(bits - 1)
        }
    }
}

/// Conversion from S1Angle
impl From<S1Angle> for S1ChordAngle {
    #[inline]
    fn from(angle: S1Angle) -> Self {
        Self::from_angle(angle)
    }
}

/// Display implementation
impl fmt::Display for S1ChordAngle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_angle())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_constants() {
        assert_eq!(S1ChordAngle::zero().length2(), 0.0);
        assert_eq!(S1ChordAngle::right().length2(), 2.0);
        assert_eq!(S1ChordAngle::straight().length2(), 4.0);
        assert!(S1ChordAngle::infinity().is_infinity());
        assert!(S1ChordAngle::negative().is_negative());
    }

    #[test]
    fn test_from_length2() {
        assert_eq!(S1ChordAngle::from_length2(0.0).degrees(), 0.0);
        assert_relative_eq!(S1ChordAngle::from_length2(1.0).degrees(), 60.0, epsilon = 1e-10);
        assert_relative_eq!(S1ChordAngle::from_length2(2.0).degrees(), 90.0, epsilon = 1e-10);
        assert_eq!(S1ChordAngle::from_length2(4.0).degrees(), 180.0);
        assert_eq!(S1ChordAngle::from_length2(5.0).degrees(), 180.0); // Clamped
    }

    #[test]
    fn test_predicates() {
        assert!(S1ChordAngle::zero().is_zero());
        assert!(!S1ChordAngle::zero().is_negative());
        assert!(!S1ChordAngle::zero().is_special());
        assert!(!S1ChordAngle::straight().is_special());
        assert!(S1ChordAngle::negative().is_negative());
        assert!(S1ChordAngle::negative().is_special());
        assert!(S1ChordAngle::infinity().is_infinity());
        assert!(S1ChordAngle::infinity().is_special());
    }

    #[test]
    fn test_arithmetic() {
        let zero = S1ChordAngle::zero();
        let degree30 = S1ChordAngle::from_degrees(30.0);
        let degree60 = S1ChordAngle::from_degrees(60.0);
        let degree90 = S1ChordAngle::from_degrees(90.0);
        
        assert_eq!((zero + zero).degrees(), 0.0);
        assert_eq!((zero - zero).degrees(), 0.0);
        assert_eq!((degree60 - degree60).degrees(), 0.0);
        assert_relative_eq!((degree30 + degree60).degrees(), 90.0, epsilon = 1e-10);
        assert_relative_eq!((degree90 - degree30).degrees(), 60.0, epsilon = 1e-10);
    }

    #[test]
    fn test_trigonometry() {
        let angle90 = S1ChordAngle::from_length2(2.0);
        let angle180 = S1ChordAngle::from_length2(4.0);
        
        assert_eq!(sin(angle90), 1.0);
        assert_eq!(cos(angle90), 0.0);
        assert_eq!(tan(angle90), f64::INFINITY);
        assert_eq!(sin(angle180), 0.0);
        assert_eq!(cos(angle180), -1.0);
        assert_eq!(tan(angle180), 0.0);
    }
}