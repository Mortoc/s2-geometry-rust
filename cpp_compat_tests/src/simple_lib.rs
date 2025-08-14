//! Simple C++ Compatibility Test Bridge
//!
//! This is a simplified version focusing on basic S2 operations without complex vector returns

use std::fmt;

#[cxx::bridge]
mod ffi {
    // Simple shared types
    struct S2PointCpp {
        x: f64,
        y: f64,
        z: f64,
    }

    struct S2LatLngCpp {
        lat_radians: f64,
        lng_radians: f64,
    }

    struct S2CellIdCpp {
        id: u64,
    }

    unsafe extern "C++" {
        include!("simple_bridge.h");

        // Simple S2Point functions
        fn simple_s2point_normalize(point: S2PointCpp) -> S2PointCpp;
        fn simple_s2point_angle(a: S2PointCpp, b: S2PointCpp) -> f64;
        fn simple_s2point_cross_prod(a: S2PointCpp, b: S2PointCpp) -> S2PointCpp;

        // Simple S2LatLng functions
        fn simple_s2latlng_from_point(point: S2PointCpp) -> S2LatLngCpp;
        fn simple_s2latlng_to_point(latlng: S2LatLngCpp) -> S2PointCpp;
        fn simple_s2latlng_distance(a: S2LatLngCpp, b: S2LatLngCpp) -> f64;

        // Simple S2CellId functions
        fn simple_s2cellid_from_point(point: S2PointCpp) -> S2CellIdCpp;
        fn simple_s2cellid_level(cell_id: S2CellIdCpp) -> i32;
        fn simple_s2cellid_parent(cell_id: S2CellIdCpp, level: i32) -> S2CellIdCpp;
        fn simple_s2cellid_to_point(cell_id: S2CellIdCpp) -> S2PointCpp;

        // Simple S2Cap/S2Cell functions  
        fn simple_s2cap_contains_point(center: S2PointCpp, radius_radians: f64, point: S2PointCpp) -> bool;
        fn simple_s2cap_area(center: S2PointCpp, radius_radians: f64) -> f64;
        fn simple_s2cell_area(cell_id: S2CellIdCpp) -> f64;
    }
}

pub use ffi::*;

// Helper functions to convert between Rust and C++ types
impl From<s2geometry_rust::S2Point> for S2PointCpp {
    fn from(point: s2geometry_rust::S2Point) -> Self {
        let coords = point.coords();
        S2PointCpp {
            x: coords.x,
            y: coords.y,
            z: coords.z,
        }
    }
}

impl From<S2PointCpp> for s2geometry_rust::S2Point {
    fn from(point: S2PointCpp) -> Self {
        s2geometry_rust::S2Point::from_normalized(
            s2geometry_rust::math::DVec3::new(point.x, point.y, point.z)
        )
    }
}

impl From<s2geometry_rust::S2LatLng> for S2LatLngCpp {
    fn from(latlng: s2geometry_rust::S2LatLng) -> Self {
        S2LatLngCpp {
            lat_radians: latlng.lat().radians(),
            lng_radians: latlng.lng().radians(),
        }
    }
}

impl From<S2LatLngCpp> for s2geometry_rust::S2LatLng {
    fn from(latlng: S2LatLngCpp) -> Self {
        s2geometry_rust::S2LatLng::from_radians(latlng.lat_radians, latlng.lng_radians)
    }
}

impl From<s2geometry_rust::S2CellId> for S2CellIdCpp {
    fn from(cell_id: s2geometry_rust::S2CellId) -> Self {
        S2CellIdCpp {
            id: cell_id.id(),
        }
    }
}

impl From<S2CellIdCpp> for s2geometry_rust::S2CellId {
    fn from(cell_id: S2CellIdCpp) -> Self {
        s2geometry_rust::S2CellId::new(cell_id.id)
    }
}

impl fmt::Debug for S2PointCpp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "S2PointCpp({}, {}, {})", self.x, self.y, self.z)
    }
}

impl fmt::Debug for S2LatLngCpp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "S2LatLngCpp({}, {})", self.lat_radians, self.lng_radians)
    }
}

impl fmt::Debug for S2CellIdCpp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "S2CellIdCpp({})", self.id)
    }
}

// Comparison helpers with tolerance
pub fn points_equal_approx(a: S2PointCpp, b: S2PointCpp, tolerance: f64) -> bool {
    (a.x - b.x).abs() < tolerance &&
    (a.y - b.y).abs() < tolerance &&
    (a.z - b.z).abs() < tolerance
}

pub fn angles_equal_approx(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() < tolerance
}