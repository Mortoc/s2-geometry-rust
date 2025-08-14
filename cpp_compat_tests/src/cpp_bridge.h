#pragma once

#include "s2/s2point.h"
#include "s2/s2latlng.h"
#include "s2/s2cell_id.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2region_coverer.h"
#include <vector>

// Forward declarations - actual definitions come from cxx bridge
struct S2PointCpp;
struct S2LatLngCpp;  
struct S2CellIdCpp;

// Helper conversion function declarations
S2Point to_s2point(const S2PointCpp& p);
S2PointCpp from_s2point(const S2Point& p);
S2LatLng to_s2latlng(const S2LatLngCpp& ll);
S2LatLngCpp from_s2latlng(const S2LatLng& ll);
S2CellId to_s2cellid(const S2CellIdCpp& c);
S2CellIdCpp from_s2cellid(const S2CellId& c);

// S2Point function declarations
S2PointCpp cpp_s2point_normalize(S2PointCpp point);
double cpp_s2point_angle(S2PointCpp a, S2PointCpp b);
S2PointCpp cpp_s2point_cross_prod(S2PointCpp a, S2PointCpp b);

// S2LatLng function declarations
S2LatLngCpp cpp_s2latlng_from_point(S2PointCpp point);
S2PointCpp cpp_s2latlng_to_point(S2LatLngCpp latlng);
double cpp_s2latlng_distance(S2LatLngCpp a, S2LatLngCpp b);

// S2CellId function declarations
S2CellIdCpp cpp_s2cellid_from_point(S2PointCpp point);
int cpp_s2cellid_level(S2CellIdCpp cell_id);
S2CellIdCpp cpp_s2cellid_parent(S2CellIdCpp cell_id, int level);
S2PointCpp cpp_s2cellid_to_point(S2CellIdCpp cell_id);

// S2Cap function declarations
void cpp_s2cap_from_point_angle(S2PointCpp center, double radius_radians, std::vector<S2CellIdCpp>& result);
bool cpp_s2cap_contains_point(S2PointCpp center, double radius_radians, S2PointCpp point);
double cpp_s2cap_area(S2PointCpp center, double radius_radians);

// S2Cell function declarations
double cpp_s2cell_area(S2CellIdCpp cell_id);
double cpp_s2cell_perimeter(S2CellIdCpp cell_id);
void cpp_s2cell_vertices(S2CellIdCpp cell_id, std::vector<S2PointCpp>& result);