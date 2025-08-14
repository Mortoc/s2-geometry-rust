//! Interval implementations matching Google's C++ S2 library
//!
//! This module contains both S1Interval (circular intervals) and R1Interval (linear intervals).
//!
//! ## S1Interval
//! An S1Interval represents a closed interval on a unit circle (1-dimensional sphere).
//! It can represent the empty interval, full interval, and zero-length intervals 
//! (containing a single point).
//!
//! Points are represented by the angle they make with the positive x-axis in
//! the range [-π, π]. An interval is represented by its lower and upper bounds
//! (both inclusive). The lower bound may be greater than the upper bound, in
//! which case the interval is "inverted" (passes through the point (-1, 0)).
//!
//! Special intervals:
//! - Empty interval: [π, -π]  
//! - Full interval: [-π, π]
//! - The point (-1, 0) has two representations: π and -π
//!
//! ## R1Interval  
//! An R1Interval represents a closed, bounded interval on the real line.
//! It is capable of representing the empty interval (containing no points)
//! and zero-length intervals (containing a single point). Unlike S1Interval,
//! R1Interval is linear with no wraparound behavior.

use crate::math::constants::*;
use std::fmt;
use std::f64::EPSILON as DBL_EPSILON;

/// IEEE 754 remainder function (like C++ remainder)
fn remainder(x: f64, y: f64) -> f64 {
    // Use round-to-even (banker's rounding) for IEEE remainder
    let quotient = x / y;
    let n = if quotient.fract().abs() == 0.5 {
        // Exactly halfway - round to even
        let trunc = quotient.trunc();
        if trunc as i64 % 2 == 0 {
            trunc
        } else {
            trunc + quotient.signum()
        }
    } else {
        quotient.round()
    };
    x - n * y
}

/// S1Interval represents a closed interval on the unit circle
///
/// This matches the C++ S1Interval implementation exactly, including
/// all edge cases and numerical precision handling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct S1Interval {
    bounds: [f64; 2], // [lo, hi]
}

impl S1Interval {
    /// Create a new S1Interval with given bounds
    ///
    /// Both endpoints must be in the range [-π, π]. The value -π is converted
    /// internally to π except for the Full() and Empty() intervals.
    pub fn new(lo: f64, hi: f64) -> S1Interval {
        debug_assert!(lo.abs() <= PI, "S1Interval lo bound {} exceeds PI", lo);
        debug_assert!(hi.abs() <= PI, "S1Interval hi bound {} exceeds PI", hi);
        
        let mut bounds = [lo, hi];
        
        // Normalize -π to π except for special cases
        if lo == -PI && hi != PI {
            bounds[0] = PI;
        }
        if hi == -PI && lo != PI {
            bounds[1] = PI;
        }
        
        let result = S1Interval { bounds };
        debug_assert!(result.is_valid());
        result
    }
    
    /// Create interval with pre-validated bounds (internal use)
    fn new_unchecked(lo: f64, hi: f64) -> S1Interval {
        let result = S1Interval { bounds: [lo, hi] };
        debug_assert!(result.is_valid());
        result
    }
    
    /// Create the empty interval
    pub fn empty() -> S1Interval {
        S1Interval { bounds: [PI, -PI] }
    }
    
    /// Create the full interval  
    pub fn full() -> S1Interval {
        S1Interval { bounds: [-PI, PI] }
    }
    
    /// Create an interval containing a single point
    pub fn from_point(p: f64) -> S1Interval {
        let normalized_p = if p == -PI { PI } else { p };
        S1Interval::new_unchecked(normalized_p, normalized_p)
    }
    
    /// Create the minimal interval containing the two given points
    pub fn from_point_pair(p1: f64, p2: f64) -> S1Interval {
        debug_assert!(p1.abs() <= PI + DBL_EPSILON);
        debug_assert!(p2.abs() <= PI + DBL_EPSILON);
        
        let norm_p1 = if p1 == -PI { PI } else { p1 };
        let norm_p2 = if p2 == -PI { PI } else { p2 };
        
        if positive_distance(norm_p1, norm_p2) <= PI {
            S1Interval::new_unchecked(norm_p1, norm_p2)
        } else {
            S1Interval::new_unchecked(norm_p2, norm_p1)
        }
    }
    
    /// Get the lower bound
    #[inline]
    pub fn lo(&self) -> f64 {
        self.bounds[0]
    }
    
    /// Get the upper bound
    #[inline] 
    pub fn hi(&self) -> f64 {
        self.bounds[1]
    }
    
    /// Get bounds as array [lo, hi]
    #[inline]
    pub fn bounds(&self) -> [f64; 2] {
        self.bounds
    }
    
    /// Check if interval is valid
    pub fn is_valid(&self) -> bool {
        (self.lo().abs() <= PI + DBL_EPSILON) &&
        (self.hi().abs() <= PI + DBL_EPSILON) &&
        !(self.lo() == -PI && self.hi() != PI) &&
        !(self.hi() == -PI && self.lo() != PI)
    }
    
    /// Check if interval is full (contains all points)
    #[inline]
    pub fn is_full(&self) -> bool {
        self.lo() == -PI && self.hi() == PI
    }
    
    /// Check if interval is empty (contains no points)
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.lo() == PI && self.hi() == -PI
    }
    
    /// Check if interval is inverted (lo > hi)
    #[inline]
    pub fn is_inverted(&self) -> bool {
        self.lo() > self.hi()
    }
    
    /// Get the center point of the interval
    pub fn get_center(&self) -> f64 {
        let center = 0.5 * (self.lo() + self.hi());
        if !self.is_inverted() {
            center
        } else {
            // Return center in range (-π, π]
            if center <= 0.0 {
                center + PI
            } else {
                center - PI
            }
        }
    }
    
    /// Get the length of the interval
    pub fn get_length(&self) -> f64 {
        let mut length = self.hi() - self.lo();
        if length >= 0.0 {
            length
        } else {
            length += 2.0 * PI;
            // Empty intervals have negative length
            if length > 0.0 { length } else { -1.0 }
        }
    }
    
    /// Get the complement of the interval
    pub fn complement(&self) -> S1Interval {
        if self.lo() == self.hi() {
            // Singleton interval -> full interval
            S1Interval::full()
        } else {
            // Swap bounds (handles empty and full correctly)
            S1Interval::new_unchecked(self.hi(), self.lo())
        }
    }
    
    /// Get the center of the complement
    pub fn get_complement_center(&self) -> f64 {
        if self.lo() != self.hi() {
            self.complement().get_center()
        } else {
            // Singleton case - return antipodal point
            if self.hi() <= 0.0 {
                self.hi() + PI
            } else {
                self.hi() - PI
            }
        }
    }
    
    /// Check if the interval contains a point (internal, assumes normalized point)
    fn fast_contains(&self, p: f64) -> bool {
        if self.is_inverted() {
            (p >= self.lo() || p <= self.hi()) && !self.is_empty()
        } else {
            p >= self.lo() && p <= self.hi()
        }
    }
    
    /// Check if the interval contains a point
    pub fn contains_point(&self, p: f64) -> bool {
        debug_assert!(p.abs() <= PI + DBL_EPSILON);
        let normalized_p = if p == -PI { PI } else { p };
        self.fast_contains(normalized_p)
    }
    
    /// Check if the interior of the interval contains a point
    pub fn interior_contains_point(&self, p: f64) -> bool {
        debug_assert!(p.abs() <= PI + DBL_EPSILON);
        let normalized_p = if p == -PI { PI } else { p };
        
        if self.is_inverted() {
            normalized_p > self.lo() || normalized_p < self.hi()
        } else {
            (normalized_p > self.lo() && normalized_p < self.hi()) || self.is_full()
        }
    }
    
    /// Check if this interval contains another interval
    pub fn contains(&self, other: &S1Interval) -> bool {
        if self.is_inverted() {
            if other.is_inverted() {
                other.lo() >= self.lo() && other.hi() <= self.hi()
            } else {
                (other.lo() >= self.lo() || other.hi() <= self.hi()) && !self.is_empty()
            }
        } else {
            if other.is_inverted() {
                self.is_full() || other.is_empty()
            } else {
                other.lo() >= self.lo() && other.hi() <= self.hi()
            }
        }
    }
    
    /// Check if the interior of this interval contains another interval
    pub fn interior_contains(&self, other: &S1Interval) -> bool {
        if self.is_inverted() {
            if !other.is_inverted() {
                other.lo() > self.lo() || other.hi() < self.hi()
            } else {
                (other.lo() > self.lo() && other.hi() < self.hi()) || other.is_empty()
            }
        } else {
            if other.is_inverted() {
                self.is_full() || other.is_empty()
            } else {
                (other.lo() > self.lo() && other.hi() < self.hi()) || self.is_full()
            }
        }
    }
    
    /// Check if this interval intersects with another
    pub fn intersects(&self, other: &S1Interval) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }
        
        if self.is_inverted() {
            // Every non-empty inverted interval contains π
            other.is_inverted() || other.lo() <= self.hi() || other.hi() >= self.lo()
        } else {
            if other.is_inverted() {
                other.lo() <= self.hi() || other.hi() >= self.lo()
            } else {
                other.lo() <= self.hi() && other.hi() >= self.lo()
            }
        }
    }
    
    /// Check if the interior of this interval intersects with another
    pub fn interior_intersects(&self, other: &S1Interval) -> bool {
        if self.is_empty() || other.is_empty() || self.lo() == self.hi() {
            return false;
        }
        
        if self.is_inverted() {
            other.is_inverted() || other.lo() < self.hi() || other.hi() > self.lo()
        } else {
            if other.is_inverted() {
                other.lo() < self.hi() || other.hi() > self.lo()
            } else {
                (other.lo() < self.hi() && other.hi() > self.lo()) || self.is_full()
            }
        }
    }
    
    /// Get the union of this interval with another
    pub fn union(&self, other: &S1Interval) -> S1Interval {
        if other.is_empty() {
            return *self;
        }
        
        if self.fast_contains(other.lo()) {
            if self.fast_contains(other.hi()) {
                // Either self contains other, or union is full
                if self.contains(other) {
                    *self
                } else {
                    S1Interval::full()
                }
            } else {
                S1Interval::new_unchecked(self.lo(), other.hi())
            }
        } else if self.fast_contains(other.hi()) {
            S1Interval::new_unchecked(other.lo(), self.hi())
        } else {
            // Self contains neither endpoint of other
            if self.is_empty() || other.fast_contains(self.lo()) {
                *other
            } else {
                // Check which pair of endpoints are closer
                let dlo = positive_distance(other.hi(), self.lo());
                let dhi = positive_distance(self.hi(), other.lo());
                if dlo < dhi {
                    S1Interval::new_unchecked(other.lo(), self.hi())
                } else {
                    S1Interval::new_unchecked(self.lo(), other.hi())
                }
            }
        }
    }
    
    /// Get the intersection of this interval with another
    pub fn intersection(&self, other: &S1Interval) -> S1Interval {
        if other.is_empty() {
            return S1Interval::empty();
        }
        
        if self.fast_contains(other.lo()) {
            if self.fast_contains(other.hi()) {
                // Either self contains other, or intersection consists of two disjoint parts
                if other.get_length() < self.get_length() {
                    *other
                } else {
                    *self
                }
            } else {
                S1Interval::new_unchecked(other.lo(), self.hi())
            }
        } else if self.fast_contains(other.hi()) {
            S1Interval::new_unchecked(self.lo(), other.hi())
        } else {
            // Self contains neither endpoint of other
            if other.fast_contains(self.lo()) {
                *self
            } else {
                debug_assert!(!self.intersects(other));
                S1Interval::empty()
            }
        }
    }
    
    /// Add a point to the interval, expanding if necessary
    pub fn add_point(&mut self, p: f64) {
        debug_assert!(p.abs() <= PI + DBL_EPSILON);
        let normalized_p = if p == -PI { PI } else { p };
        
        if self.fast_contains(normalized_p) {
            return;
        }
        
        if self.is_empty() {
            self.bounds[0] = normalized_p;
            self.bounds[1] = normalized_p;
        } else {
            // Compute distance from p to each endpoint
            let dlo = positive_distance(normalized_p, self.lo());
            let dhi = positive_distance(self.hi(), normalized_p);
            if dlo < dhi {
                self.bounds[0] = normalized_p;
            } else {
                self.bounds[1] = normalized_p;
            }
        }
    }
    
    /// Project a point to the closest point in the interval
    pub fn project(&self, p: f64) -> f64 {
        debug_assert!(!self.is_empty());
        debug_assert!(p.abs() <= PI + DBL_EPSILON);
        
        let normalized_p = if p == -PI { PI } else { p };
        if self.fast_contains(normalized_p) {
            return normalized_p;
        }
        
        // Compute distance from p to each endpoint
        let dlo = positive_distance(normalized_p, self.lo());
        let dhi = positive_distance(self.hi(), normalized_p);
        if dlo < dhi { self.lo() } else { self.hi() }
    }
    
    /// Expand the interval by the given margin on each side
    pub fn expanded(&self, margin: f64) -> S1Interval {
        if margin >= 0.0 {
            if self.is_empty() {
                return *self;
            }
            
            // Check if interval becomes full after expansion
            if self.get_length() + 2.0 * margin + 2.0 * DBL_EPSILON >= 2.0 * PI {
                return S1Interval::full();
            }
        } else {
            if self.is_full() {
                return *self;
            }
            
            // Check if interval becomes empty after expansion
            if self.get_length() + 2.0 * margin - 2.0 * DBL_EPSILON <= 0.0 {
                return S1Interval::empty();
            }
        }
        
        // Use remainder function like C++ to get results in [-π, π]
        let new_lo = remainder(self.lo() - margin, 2.0 * PI);
        let new_hi = remainder(self.hi() + margin, 2.0 * PI);
        
        // Create result with internal constructor to avoid automatic normalization
        let mut result = S1Interval {
            bounds: [new_lo, new_hi]
        };
        
        // Apply the same normalization logic as the C++ constructor
        if result.lo() <= -PI {
            result.set_lo(PI);
        }
        if result.hi() == -PI && result.lo() != PI {
            result.set_hi(PI);
        }
        
        result
    }
    
    /// Set the lower bound (internal use)
    fn set_lo(&mut self, p: f64) {
        self.bounds[0] = p;
        debug_assert!(self.is_valid());
    }
    
    /// Set the upper bound (internal use)  
    fn set_hi(&mut self, p: f64) {
        self.bounds[1] = p;
        debug_assert!(self.is_valid());
    }
    
    /// Get the directed Hausdorff distance to another interval
    pub fn get_directed_hausdorff_distance(&self, other: &S1Interval) -> f64 {
        if other.contains(self) {
            return 0.0; // Includes case where self is empty
        }
        if other.is_empty() {
            return PI; // Maximum possible distance on S1
        }
        
        let other_complement_center = other.get_complement_center();
        if self.contains_point(other_complement_center) {
            positive_distance(other.hi(), other_complement_center)
        } else {
            // Hausdorff distance is realized by either two hi() endpoints or two lo() endpoints
            let hi_hi = if S1Interval::new_unchecked(other.hi(), other_complement_center)
                .contains_point(self.hi()) {
                positive_distance(other.hi(), self.hi())
            } else {
                0.0
            };
            
            let lo_lo = if S1Interval::new_unchecked(other_complement_center, other.lo())
                .contains_point(self.lo()) {
                positive_distance(self.lo(), other.lo())
            } else {
                0.0
            };
            
            debug_assert!(hi_hi > 0.0 || lo_lo > 0.0);
            hi_hi.max(lo_lo)
        }
    }
    
    /// Check if this interval approximately equals another
    pub fn approx_equals(&self, other: &S1Interval, max_error: f64) -> bool {
        // Full and empty intervals require special cases
        if self.is_empty() {
            return other.get_length() <= 2.0 * max_error;
        }
        if other.is_empty() {
            return self.get_length() <= 2.0 * max_error;
        }
        if self.is_full() {
            return other.get_length() >= 2.0 * (PI - max_error);
        }
        if other.is_full() {
            return self.get_length() >= 2.0 * (PI - max_error);
        }
        
        // Check endpoints and length differences using remainder like C++
        let lo_diff = remainder(other.lo() - self.lo(), 2.0 * PI).abs();
        let hi_diff = remainder(other.hi() - self.hi(), 2.0 * PI).abs();
        let length_diff = (self.get_length() - other.get_length()).abs();
        
        lo_diff <= max_error &&
        hi_diff <= max_error &&
        length_diff <= 2.0 * max_error
    }
}

impl Default for S1Interval {
    fn default() -> Self {
        S1Interval::empty()
    }
}

impl fmt::Display for S1Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.lo(), self.hi())
    }
}

/// Compute positive distance from a to b in range [0, 2π)
fn positive_distance(a: f64, b: f64) -> f64 {
    let d = b - a;
    if d >= 0.0 {
        d
    } else {
        // Ensure proper handling of boundary cases
        (b + PI) - (a - PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_constructors() {
        let empty = S1Interval::empty();
        assert!(empty.is_empty());
        assert!(empty.is_valid());
        
        let full = S1Interval::full();
        assert!(full.is_full());
        assert!(full.is_valid());
        
        let point = S1Interval::from_point(0.5);
        assert_eq!(point.lo(), 0.5);
        assert_eq!(point.hi(), 0.5);
        assert!(!point.is_empty());
        assert!(!point.is_full());
    }
    
    #[test]
    fn test_predicates() {
        let empty = S1Interval::empty();
        let full = S1Interval::full();
        let zero = S1Interval::from_point(0.0);
        
        assert!(empty.is_empty() && !empty.is_full() && empty.is_inverted());
        assert!(full.is_full() && !full.is_empty() && !full.is_inverted());
        assert!(!zero.is_empty() && !zero.is_full() && !zero.is_inverted());
    }
    
    #[test]
    fn test_contains() {
        let interval = S1Interval::new(0.0, PI_2);
        assert!(interval.contains_point(0.0));
        assert!(interval.contains_point(PI_4));
        assert!(interval.contains_point(PI_2));
        assert!(!interval.contains_point(-0.1));
        assert!(!interval.contains_point(PI_2 + 0.1));
    }
    
    #[test]
    fn test_inverted_interval() {
        let inverted = S1Interval::new(PI_2, -PI_2);
        assert!(inverted.is_inverted());
        assert!(inverted.contains_point(PI));
        assert!(inverted.contains_point(-PI));
        assert!(inverted.contains_point(0.9 * PI));
        assert!(!inverted.contains_point(0.0));
    }
}

/// R1Interval represents a closed, bounded interval on the real line
///
/// This matches the C++ R1Interval implementation exactly, including
/// all edge cases and numerical precision handling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct R1Interval {
    bounds: [f64; 2], // [lo, hi]
}

impl R1Interval {
    /// Create a new R1Interval with given bounds
    /// If lo > hi, the interval is empty.
    pub fn new(lo: f64, hi: f64) -> R1Interval {
        R1Interval { bounds: [lo, hi] }
    }

    /// Create the empty interval
    /// Any interval where lo > hi is considered empty
    pub fn empty() -> R1Interval {
        R1Interval { bounds: [1.0, 0.0] }
    }

    /// Create an interval containing a single point
    pub fn from_point(p: f64) -> R1Interval {
        R1Interval { bounds: [p, p] }
    }

    /// Create the minimal interval containing the two given points
    /// This is equivalent to starting with an empty interval and calling add_point() twice
    pub fn from_point_pair(p1: f64, p2: f64) -> R1Interval {
        if p1 <= p2 {
            R1Interval { bounds: [p1, p2] }
        } else {
            R1Interval { bounds: [p2, p1] }
        }
    }

    /// Get the lower bound
    #[inline]
    pub fn lo(&self) -> f64 {
        self.bounds[0]
    }

    /// Get the upper bound
    #[inline]
    pub fn hi(&self) -> f64 {
        self.bounds[1]
    }

    /// Set the lower bound
    pub fn set_lo(&mut self, p: f64) {
        self.bounds[0] = p;
    }

    /// Set the upper bound
    pub fn set_hi(&mut self, p: f64) {
        self.bounds[1] = p;
    }

    /// Access bounds by index (0 = lo, 1 = hi)
    pub fn get(&self, i: usize) -> f64 {
        debug_assert!(i < 2);
        self.bounds[i]
    }

    /// Mutably access bounds by index (0 = lo, 1 = hi)
    pub fn get_mut(&mut self, i: usize) -> &mut f64 {
        debug_assert!(i < 2);
        &mut self.bounds[i]
    }

    /// Get bounds as array [lo, hi]
    #[inline]
    pub fn bounds(&self) -> [f64; 2] {
        self.bounds
    }

    /// Get mutable reference to bounds array
    #[inline]
    pub fn bounds_mut(&mut self) -> &mut [f64; 2] {
        &mut self.bounds
    }

    /// Check if interval is empty (contains no points)
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.lo() > self.hi()
    }

    /// Get the center of the interval
    /// For empty intervals, the result is arbitrary
    pub fn get_center(&self) -> f64 {
        0.5 * (self.lo() + self.hi())
    }

    /// Get the length of the interval
    /// The length of an empty interval is negative
    pub fn get_length(&self) -> f64 {
        self.hi() - self.lo()
    }

    /// Check if the interval contains a point
    pub fn contains(&self, p: f64) -> bool {
        p >= self.lo() && p <= self.hi()
    }

    /// Check if the interior of the interval contains a point
    pub fn interior_contains(&self, p: f64) -> bool {
        p > self.lo() && p < self.hi()
    }

    /// Check if this interval contains another interval
    pub fn contains_interval(&self, other: &R1Interval) -> bool {
        if other.is_empty() {
            return true;
        }
        other.lo() >= self.lo() && other.hi() <= self.hi()
    }

    /// Check if the interior of this interval contains another interval
    pub fn interior_contains_interval(&self, other: &R1Interval) -> bool {
        if other.is_empty() {
            return true;
        }
        other.lo() > self.lo() && other.hi() < self.hi()
    }

    /// Check if this interval intersects another interval
    pub fn intersects(&self, other: &R1Interval) -> bool {
        if self.lo() <= other.lo() {
            other.lo() <= self.hi() && other.lo() <= other.hi()
        } else {
            self.lo() <= other.hi() && self.lo() <= self.hi()
        }
    }

    /// Check if the interior of this interval intersects another interval
    pub fn interior_intersects(&self, other: &R1Interval) -> bool {
        other.lo() < self.hi() && self.lo() < other.hi() && self.lo() < self.hi() && other.lo() <= other.hi()
    }

    /// Get the directed Hausdorff distance to another interval
    pub fn get_directed_hausdorff_distance(&self, other: &R1Interval) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        if other.is_empty() {
            return f64::INFINITY;
        }
        f64::max(0.0, f64::max(self.hi() - other.hi(), other.lo() - self.lo()))
    }

    /// Expand the interval to contain the given point
    pub fn add_point(&mut self, p: f64) {
        if self.is_empty() {
            self.bounds[0] = p;
            self.bounds[1] = p;
        } else if p < self.lo() {
            self.bounds[0] = p;
        } else if p > self.hi() {
            self.bounds[1] = p;
        }
    }

    /// Expand the interval to contain another interval
    pub fn add_interval(&mut self, other: &R1Interval) {
        if other.is_empty() {
            return;
        }
        if self.is_empty() {
            *self = *other;
            return;
        }
        if other.lo() < self.lo() {
            self.bounds[0] = other.lo();
        }
        if other.hi() > self.hi() {
            self.bounds[1] = other.hi();
        }
    }

    /// Project a point to the closest point in the interval
    /// The interval must be non-empty
    pub fn project(&self, p: f64) -> f64 {
        debug_assert!(!self.is_empty());
        p.clamp(self.lo(), self.hi())
    }

    /// Return an interval expanded by the given margin on each side
    /// If margin is negative, shrink the interval instead
    pub fn expanded(&self, margin: f64) -> R1Interval {
        if self.is_empty() {
            return *self;
        }
        R1Interval::new(self.lo() - margin, self.hi() + margin)
    }

    /// Return the union of this interval with another
    pub fn union(&self, other: &R1Interval) -> R1Interval {
        if self.is_empty() {
            return *other;
        }
        if other.is_empty() {
            return *self;
        }
        R1Interval::new(
            f64::min(self.lo(), other.lo()),
            f64::max(self.hi(), other.hi())
        )
    }

    /// Return the intersection of this interval with another
    pub fn intersection(&self, other: &R1Interval) -> R1Interval {
        R1Interval::new(
            f64::max(self.lo(), other.lo()),
            f64::min(self.hi(), other.hi())
        )
    }

    /// Check if this interval approximately equals another
    pub fn approx_equals(&self, other: &R1Interval, max_error: f64) -> bool {
        if self.is_empty() {
            return other.get_length() <= 2.0 * max_error;
        }
        if other.is_empty() {
            return self.get_length() <= 2.0 * max_error;
        }
        (self.lo() - other.lo()).abs() <= max_error && 
        (self.hi() - other.hi()).abs() <= max_error
    }
}

impl Default for R1Interval {
    fn default() -> Self {
        R1Interval::empty()
    }
}

impl fmt::Display for R1Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.lo(), self.hi())
    }
}

// Implement index access like C++
impl std::ops::Index<usize> for R1Interval {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.bounds[index]
    }
}

impl std::ops::IndexMut<usize> for R1Interval {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.bounds[index]
    }
}