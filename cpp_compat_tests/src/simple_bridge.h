#pragma once

#include "s2/s2point.h"
#include "s2/s2latlng.h"
#include "s2/s2cell_id.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"

// Forward declarations - actual definitions come from cxx bridge
struct S2PointCpp;
struct S2LatLngCpp;  
struct S2CellIdCpp;

// Helper conversion function declarations
S2Point simple_to_s2point(const S2PointCpp& p);
S2PointCpp simple_from_s2point(const S2Point& p);
S2LatLng simple_to_s2latlng(const S2LatLngCpp& ll);
S2LatLngCpp simple_from_s2latlng(const S2LatLng& ll);
S2CellId simple_to_s2cellid(const S2CellIdCpp& c);
S2CellIdCpp simple_from_s2cellid(const S2CellId& c);

// Simple S2Point function declarations
S2PointCpp simple_s2point_normalize(S2PointCpp point);
double simple_s2point_angle(S2PointCpp a, S2PointCpp b);
S2PointCpp simple_s2point_cross_prod(S2PointCpp a, S2PointCpp b);

// Simple S2LatLng function declarations
S2LatLngCpp simple_s2latlng_from_point(S2PointCpp point);
S2PointCpp simple_s2latlng_to_point(S2LatLngCpp latlng);
double simple_s2latlng_distance(S2LatLngCpp a, S2LatLngCpp b);

// Simple S2CellId function declarations
S2CellIdCpp simple_s2cellid_from_point(S2PointCpp point);
int simple_s2cellid_level(S2CellIdCpp cell_id);
S2CellIdCpp simple_s2cellid_parent(S2CellIdCpp cell_id, int level);
S2PointCpp simple_s2cellid_to_point(S2CellIdCpp cell_id);

// Simple S2Cap/S2Cell function declarations
bool simple_s2cap_contains_point(S2PointCpp center, double radius_radians, S2PointCpp point);
double simple_s2cap_area(S2PointCpp center, double radius_radians);
double simple_s2cell_area(S2CellIdCpp cell_id);