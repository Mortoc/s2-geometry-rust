#include "cpp_bridge.h"
#include "s2/s2point.h"
#include "s2/s2latlng.h"
#include "s2/s2cell_id.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2region_coverer.h"

// Helper conversion function implementations
S2Point to_s2point(const S2PointCpp& p) {
    return S2Point(p.x, p.y, p.z);
}

S2PointCpp from_s2point(const S2Point& p) {
    return {p.x(), p.y(), p.z()};
}

S2LatLng to_s2latlng(const S2LatLngCpp& ll) {
    return S2LatLng::FromRadians(ll.lat_radians, ll.lng_radians);
}

S2LatLngCpp from_s2latlng(const S2LatLng& ll) {
    return {ll.lat().radians(), ll.lng().radians()};
}

S2CellId to_s2cellid(const S2CellIdCpp& c) {
    return S2CellId(c.id);
}

S2CellIdCpp from_s2cellid(const S2CellId& c) {
    return {c.id()};
}

// S2Point functions
S2PointCpp cpp_s2point_normalize(S2PointCpp point) {
    S2Point p = to_s2point(point);
    return from_s2point(p.Normalize());
}

double cpp_s2point_angle(S2PointCpp a, S2PointCpp b) {
    return to_s2point(a).Angle(to_s2point(b)).radians();
}

S2PointCpp cpp_s2point_cross_prod(S2PointCpp a, S2PointCpp b) {
    return from_s2point(to_s2point(a).CrossProd(to_s2point(b)));
}

// S2LatLng functions
S2LatLngCpp cpp_s2latlng_from_point(S2PointCpp point) {
    return from_s2latlng(S2LatLng(to_s2point(point)));
}

S2PointCpp cpp_s2latlng_to_point(S2LatLngCpp latlng) {
    return from_s2point(to_s2latlng(latlng).ToPoint());
}

double cpp_s2latlng_distance(S2LatLngCpp a, S2LatLngCpp b) {
    return to_s2latlng(a).GetDistance(to_s2latlng(b)).radians();
}

// S2CellId functions
S2CellIdCpp cpp_s2cellid_from_point(S2PointCpp point) {
    return from_s2cellid(S2CellId(to_s2point(point)));
}

int cpp_s2cellid_level(S2CellIdCpp cell_id) {
    return to_s2cellid(cell_id).level();
}

S2CellIdCpp cpp_s2cellid_parent(S2CellIdCpp cell_id, int level) {
    return from_s2cellid(to_s2cellid(cell_id).parent(level));
}

S2PointCpp cpp_s2cellid_to_point(S2CellIdCpp cell_id) {
    return from_s2point(to_s2cellid(cell_id).ToPoint());
}

// S2Cap functions
void cpp_s2cap_from_point_angle(S2PointCpp center, double radius_radians, std::vector<S2CellIdCpp>& result) {
    S2Cap cap(to_s2point(center), S1Angle::Radians(radius_radians));
    
    // Use S2RegionCoverer to get cell covering
    S2RegionCoverer coverer;
    coverer.mutable_options()->set_max_level(10);
    coverer.mutable_options()->set_max_cells(100);
    
    std::vector<S2CellId> covering;
    coverer.GetCovering(cap, &covering);
    
    result.clear();
    result.reserve(covering.size());
    for (const auto& cell_id : covering) {
        result.push_back(from_s2cellid(cell_id));
    }
}

bool cpp_s2cap_contains_point(S2PointCpp center, double radius_radians, S2PointCpp point) {
    S2Cap cap(to_s2point(center), S1Angle::Radians(radius_radians));
    return cap.Contains(to_s2point(point));
}

double cpp_s2cap_area(S2PointCpp center, double radius_radians) {
    S2Cap cap(to_s2point(center), S1Angle::Radians(radius_radians));
    return cap.GetArea();
}

// S2Cell functions
double cpp_s2cell_area(S2CellIdCpp cell_id) {
    S2Cell cell(to_s2cellid(cell_id));
    return cell.ExactArea();
}

double cpp_s2cell_perimeter(S2CellIdCpp cell_id) {
    S2Cell cell(to_s2cellid(cell_id));
    double perimeter = 0.0;
    for (int i = 0; i < 4; ++i) {
        S2Point v0 = cell.GetVertex(i);
        S2Point v1 = cell.GetVertex((i + 1) % 4);
        perimeter += v0.Angle(v1).radians();
    }
    return perimeter;
}

void cpp_s2cell_vertices(S2CellIdCpp cell_id, std::vector<S2PointCpp>& result) {
    S2Cell cell(to_s2cellid(cell_id));
    result.clear();
    result.reserve(4);
    for (int i = 0; i < 4; ++i) {
        result.push_back(from_s2point(cell.GetVertex(i)));
    }
}