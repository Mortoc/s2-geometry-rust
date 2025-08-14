//! R2 geometry types - Points and rectangles in 2D Euclidean space
//!
//! This module provides 2D geometric primitives used throughout S2 geometry:
//! - R2Point: A point in 2D Euclidean space  
//! - R2Rect: An axis-aligned rectangle in 2D space
//!
//! These types are used for projections, bounding boxes, and 2D geometric operations
//! within the S2 geometry system.

use crate::interval::R1Interval;
use std::fmt;
use std::ops::{Add, Sub, Mul};

/// A point in 2D Euclidean space
///
/// R2Point represents a point in the Euclidean plane with x and y coordinates.
/// This matches the C++ Vector2_d functionality used in S2Geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct R2Point {
    x: f64,
    y: f64,
}

impl R2Point {
    /// Create a new R2Point with the given coordinates
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        R2Point { x, y }
    }

    /// Get the x coordinate
    #[inline]
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Get the y coordinate
    #[inline]
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Set the x coordinate
    #[inline]
    pub fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    /// Set the y coordinate
    #[inline]
    pub fn set_y(&mut self, y: f64) {
        self.y = y;
    }

    /// Get coordinates as array [x, y]
    #[inline]
    pub fn coords(&self) -> [f64; 2] {
        [self.x, self.y]
    }

    /// Compute the dot product with another point
    #[inline]
    pub fn dot_prod(&self, other: &R2Point) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// Compute the cross product with another point (returns scalar)
    #[inline]
    pub fn cross_prod(&self, other: &R2Point) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Get a vector orthogonal to this one (rotated 90 degrees counterclockwise)
    #[inline]
    pub fn ortho(&self) -> R2Point {
        R2Point::new(-self.y, self.x)
    }

    /// Compute the squared distance to another point
    #[inline]
    pub fn distance_squared(&self, other: &R2Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Compute the distance to another point
    #[inline]
    pub fn distance(&self, other: &R2Point) -> f64 {
        self.distance_squared(other).sqrt()
    }

    /// Get the squared norm (length squared) of this point as a vector
    #[inline]
    pub fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    /// Get the norm (length) of this point as a vector
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Check if approximately equal to another point within max_error
    /// This matches the C++ aequal function used in S2LatLng::ApproxEquals
    #[inline]
    pub fn aequal(&self, other: R2Point, max_error: f64) -> bool {
        (self.x - other.x).abs() <= max_error && (self.y - other.y).abs() <= max_error
    }
}

impl Default for R2Point {
    fn default() -> Self {
        R2Point::new(0.0, 0.0)
    }
}

impl Add for R2Point {
    type Output = R2Point;

    fn add(self, other: R2Point) -> R2Point {
        R2Point::new(self.x + other.x, self.y + other.y)
    }
}

impl Sub for R2Point {
    type Output = R2Point;

    fn sub(self, other: R2Point) -> R2Point {
        R2Point::new(self.x - other.x, self.y - other.y)
    }
}

impl Mul<f64> for R2Point {
    type Output = R2Point;

    fn mul(self, scalar: f64) -> R2Point {
        R2Point::new(self.x * scalar, self.y * scalar)
    }
}

impl fmt::Display for R2Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

/// An axis-aligned rectangle in 2D Euclidean space
///
/// R2Rect represents a closed, axis-aligned rectangle in the (x,y) plane.
/// It's implemented using two R1Intervals for the x and y dimensions.
/// This matches the C++ R2Rect implementation exactly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct R2Rect {
    bounds: [R1Interval; 2], // [x_interval, y_interval]
}

impl R2Rect {
    /// Construct a rectangle from the given lower-left and upper-right points
    pub fn new(lo: R2Point, hi: R2Point) -> Self {
        R2Rect {
            bounds: [
                R1Interval::new(lo.x(), hi.x()),
                R1Interval::new(lo.y(), hi.y()),
            ],
        }
    }

    /// Construct a rectangle from the given intervals in x and y
    /// The two intervals must either be both empty or both non-empty
    pub fn from_intervals(x: R1Interval, y: R1Interval) -> Self {
        let result = R2Rect { bounds: [x, y] };
        debug_assert!(result.is_valid());
        result
    }

    /// Create an empty rectangle
    pub fn empty() -> Self {
        R2Rect {
            bounds: [R1Interval::empty(), R1Interval::empty()],
        }
    }

    /// Construct a rectangle from a center point and size in each dimension
    /// Both components of size should be non-negative
    pub fn from_center_size(center: R2Point, size: R2Point) -> Self {
        let half_x = size.x() * 0.5;
        let half_y = size.y() * 0.5;
        R2Rect::new(
            R2Point::new(center.x() - half_x, center.y() - half_y),
            R2Point::new(center.x() + half_x, center.y() + half_y),
        )
    }

    /// Construct a rectangle containing a single point
    pub fn from_point(p: R2Point) -> Self {
        R2Rect::new(p, p)
    }

    /// Construct the minimal bounding rectangle containing the two given points
    pub fn from_point_pair(p1: R2Point, p2: R2Point) -> Self {
        R2Rect {
            bounds: [
                R1Interval::from_point_pair(p1.x(), p2.x()),
                R1Interval::from_point_pair(p1.y(), p2.y()),
            ],
        }
    }

    /// Get the x interval
    #[inline]
    pub fn x(&self) -> R1Interval {
        self.bounds[0]
    }

    /// Get the y interval
    #[inline]
    pub fn y(&self) -> R1Interval {
        self.bounds[1]
    }

    /// Get the lower-left corner
    #[inline]
    pub fn lo(&self) -> R2Point {
        R2Point::new(self.x().lo(), self.y().lo())
    }

    /// Get the upper-right corner
    #[inline]
    pub fn hi(&self) -> R2Point {
        R2Point::new(self.x().hi(), self.y().hi())
    }

    /// Access bounds by index (0 = x, 1 = y)
    #[inline]
    pub fn get_bound(&self, i: usize) -> R1Interval {
        debug_assert!(i < 2);
        self.bounds[i]
    }

    /// Mutably access bounds by index (0 = x, 1 = y)
    #[inline]
    pub fn get_bound_mut(&mut self, i: usize) -> &mut R1Interval {
        debug_assert!(i < 2);
        &mut self.bounds[i]
    }

    /// Check if the rectangle is valid
    /// A rectangle is valid if both intervals are either empty or non-empty
    pub fn is_valid(&self) -> bool {
        self.x().is_empty() == self.y().is_empty()
    }

    /// Check if the rectangle is empty (contains no points)
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.x().is_empty()
    }

    /// Get the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order
    /// Vertex 0 is in the lower-left corner
    pub fn get_vertex(&self, k: i32) -> R2Point {
        // Twiddle bits to return points in CCW order (lower left, lower right, upper right, upper left)
        let k = k & 3; // Reduce modulo 4
        let j = (k >> 1) & 1;
        self.get_vertex_ij(j ^ (k & 1), j)
    }

    /// Get the vertex at position (i, j) where i=0/1 is left/right and j=0/1 is bottom/top
    pub fn get_vertex_ij(&self, i: i32, j: i32) -> R2Point {
        debug_assert!(i >= 0 && i <= 1);
        debug_assert!(j >= 0 && j <= 1);
        R2Point::new(
            self.bounds[0].get(i as usize),
            self.bounds[1].get(j as usize),
        )
    }

    /// Get the center of the rectangle
    #[inline]
    pub fn get_center(&self) -> R2Point {
        R2Point::new(self.x().get_center(), self.y().get_center())
    }

    /// Get the size (width, height) of the rectangle
    #[inline]
    pub fn get_size(&self) -> R2Point {
        R2Point::new(self.x().get_length(), self.y().get_length())
    }

    /// Check if the rectangle contains the given point
    #[inline]
    pub fn contains(&self, p: R2Point) -> bool {
        self.x().contains(p.x()) && self.y().contains(p.y())
    }

    /// Check if the interior of the rectangle contains the given point
    #[inline]
    pub fn interior_contains(&self, p: R2Point) -> bool {
        self.x().interior_contains(p.x()) && self.y().interior_contains(p.y())
    }

    /// Check if this rectangle contains another rectangle
    pub fn contains_rect(&self, other: &R2Rect) -> bool {
        self.x().contains_interval(&other.x()) && self.y().contains_interval(&other.y())
    }

    /// Check if the interior of this rectangle contains another rectangle
    pub fn interior_contains_rect(&self, other: &R2Rect) -> bool {
        self.x().interior_contains_interval(&other.x()) && self.y().interior_contains_interval(&other.y())
    }

    /// Check if this rectangle intersects another rectangle
    pub fn intersects(&self, other: &R2Rect) -> bool {
        self.x().intersects(&other.x()) && self.y().intersects(&other.y())
    }

    /// Check if the interior of this rectangle intersects another rectangle
    pub fn interior_intersects(&self, other: &R2Rect) -> bool {
        self.x().interior_intersects(&other.x()) && self.y().interior_intersects(&other.y())
    }

    /// Expand the rectangle to include the given point
    pub fn add_point(&mut self, p: R2Point) {
        self.bounds[0].add_point(p.x());
        self.bounds[1].add_point(p.y());
    }

    /// Expand the rectangle to include another rectangle
    pub fn add_rect(&mut self, other: &R2Rect) {
        self.bounds[0].add_interval(&other.x());
        self.bounds[1].add_interval(&other.y());
    }

    /// Project a point to the closest point in the rectangle
    /// The rectangle must be non-empty
    pub fn project(&self, p: R2Point) -> R2Point {
        debug_assert!(!self.is_empty());
        R2Point::new(
            self.x().project(p.x()),
            self.y().project(p.y()),
        )
    }

    /// Return a rectangle expanded by the given margin on each side
    pub fn expanded(&self, margin: R2Point) -> R2Rect {
        let xx = self.x().expanded(margin.x());
        let yy = self.y().expanded(margin.y());
        if xx.is_empty() || yy.is_empty() {
            return R2Rect::empty();
        }
        R2Rect { bounds: [xx, yy] }
    }

    /// Return a rectangle expanded by the given margin on all sides
    #[inline]
    pub fn expanded_by_margin(&self, margin: f64) -> R2Rect {
        self.expanded(R2Point::new(margin, margin))
    }

    /// Return the union of this rectangle with another
    pub fn union(&self, other: &R2Rect) -> R2Rect {
        R2Rect {
            bounds: [
                self.x().union(&other.x()),
                self.y().union(&other.y()),
            ],
        }
    }

    /// Return the intersection of this rectangle with another
    pub fn intersection(&self, other: &R2Rect) -> R2Rect {
        let xx = self.x().intersection(&other.x());
        let yy = self.y().intersection(&other.y());
        if xx.is_empty() || yy.is_empty() {
            return R2Rect::empty();
        }
        R2Rect { bounds: [xx, yy] }
    }

    /// Check if this rectangle approximately equals another
    pub fn approx_equals(&self, other: &R2Rect, max_error: f64) -> bool {
        self.x().approx_equals(&other.x(), max_error) && 
        self.y().approx_equals(&other.y(), max_error)
    }
}

impl Default for R2Rect {
    fn default() -> Self {
        R2Rect::empty()
    }
}

// Implement indexing like C++
impl std::ops::Index<usize> for R2Rect {
    type Output = R1Interval;

    fn index(&self, index: usize) -> &Self::Output {
        &self.bounds[index]
    }
}

impl std::ops::IndexMut<usize> for R2Rect {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.bounds[index]
    }
}

impl fmt::Display for R2Rect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "EMPTY")
        } else {
            write!(f, "[{}, {}] x [{}, {}]", 
                   self.x().lo(), self.x().hi(),
                   self.y().lo(), self.y().hi())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r2point_basic() {
        let p = R2Point::new(1.0, 2.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
        
        let p2 = R2Point::new(3.0, 4.0);
        let sum = p + p2;
        assert_eq!(sum.x(), 4.0);
        assert_eq!(sum.y(), 6.0);
    }

    #[test]
    fn test_r2rect_basic() {
        let rect = R2Rect::new(R2Point::new(0.0, 1.0), R2Point::new(2.0, 3.0));
        assert!(!rect.is_empty());
        assert!(rect.is_valid());
        
        assert_eq!(rect.lo(), R2Point::new(0.0, 1.0));
        assert_eq!(rect.hi(), R2Point::new(2.0, 3.0));
        
        assert!(rect.contains(R2Point::new(1.0, 2.0)));
        assert!(!rect.contains(R2Point::new(3.0, 2.0)));
    }

    #[test]
    fn test_r2rect_empty() {
        let empty = R2Rect::empty();
        assert!(empty.is_empty());
        assert!(empty.is_valid());
    }
}