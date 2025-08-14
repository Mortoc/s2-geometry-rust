#include "s2geometry-cpp-compat-tests/src/simple_lib.rs.h"
#include "simple_bridge.h"
#include "s2/s2point.h"
#include "s2/s2latlng.h"
#include "s2/s2cell_id.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"

// Helper conversion function implementations
S2Point simple_to_s2point(const S2PointCpp& p) {
    return S2Point(p.x, p.y, p.z);
}

S2PointCpp simple_from_s2point(const S2Point& p) {
    return {p.x(), p.y(), p.z()};
}

S2LatLng simple_to_s2latlng(const S2LatLngCpp& ll) {
    return S2LatLng::FromRadians(ll.lat_radians, ll.lng_radians);
}

S2LatLngCpp simple_from_s2latlng(const S2LatLng& ll) {
    return {ll.lat().radians(), ll.lng().radians()};
}

S2CellId simple_to_s2cellid(const S2CellIdCpp& c) {
    return S2CellId(c.id);
}

S2CellIdCpp simple_from_s2cellid(const S2CellId& c) {
    return {c.id()};
}

// S2Point functions
S2PointCpp simple_s2point_normalize(S2PointCpp point) {
    S2Point p = simple_to_s2point(point);
    return simple_from_s2point(p.Normalize());
}

double simple_s2point_angle(S2PointCpp a, S2PointCpp b) {
    return simple_to_s2point(a).Angle(simple_to_s2point(b)).radians();
}

S2PointCpp simple_s2point_cross_prod(S2PointCpp a, S2PointCpp b) {
    return simple_from_s2point(simple_to_s2point(a).CrossProd(simple_to_s2point(b)));
}

// S2LatLng functions
S2LatLngCpp simple_s2latlng_from_point(S2PointCpp point) {
    return simple_from_s2latlng(S2LatLng(simple_to_s2point(point)));
}

S2PointCpp simple_s2latlng_to_point(S2LatLngCpp latlng) {
    return simple_from_s2point(simple_to_s2latlng(latlng).ToPoint());
}

double simple_s2latlng_distance(S2LatLngCpp a, S2LatLngCpp b) {
    return simple_to_s2latlng(a).GetDistance(simple_to_s2latlng(b)).radians();
}

// S2CellId functions
S2CellIdCpp simple_s2cellid_from_point(S2PointCpp point) {
    return simple_from_s2cellid(S2CellId(simple_to_s2point(point)));
}

int simple_s2cellid_level(S2CellIdCpp cell_id) {
    return simple_to_s2cellid(cell_id).level();
}

S2CellIdCpp simple_s2cellid_parent(S2CellIdCpp cell_id, int level) {
    return simple_from_s2cellid(simple_to_s2cellid(cell_id).parent(level));
}

S2PointCpp simple_s2cellid_to_point(S2CellIdCpp cell_id) {
    return simple_from_s2point(simple_to_s2cellid(cell_id).ToPoint());
}

// S2Cap functions
bool simple_s2cap_contains_point(S2PointCpp center, double radius_radians, S2PointCpp point) {
    S2Cap cap(simple_to_s2point(center), S1Angle::Radians(radius_radians));
    return cap.Contains(simple_to_s2point(point));
}

double simple_s2cap_area(S2PointCpp center, double radius_radians) {
    S2Cap cap(simple_to_s2point(center), S1Angle::Radians(radius_radians));
    return cap.GetArea();
}

// S2Cell functions
double simple_s2cell_area(S2CellIdCpp cell_id) {
    S2Cell cell(simple_to_s2cellid(cell_id));
    return cell.ExactArea();
}