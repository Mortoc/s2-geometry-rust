#include "s2geometry-cpp/src/s2/s2cap.h"
#include "s2geometry-cpp/src/s2/s2cell.h"
#include "s2geometry-cpp/src/s2/s2cell_id.h"
#include "s2geometry-cpp/src/s2/s1angle.h"
#include "s2geometry-cpp/src/s2/s2point.h"
#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    std::cout << std::fixed << std::setprecision(15);
    
    // Test the specific case from the failing Rust test:
    // Cap center: (0, 1, 0) - exactly on +Y axis (face 2)
    // Cap radius: atan(sqrt(2)) + 1e-15 (face_radius + EPS)
    // Cell: S2Cell::FromFace(0) - the entire +X face (face 0)
    
    S2Point cap_center(0.0, 1.0, 0.0);  // +Y axis (face 2)
    
    // Calculate face_radius = atan(sqrt(2)) + EPS
    double face_radius = std::atan(std::sqrt(2.0)) + 1e-15;
    
    std::cout << "=== S2Cap::MayIntersect Test Case ===" << std::endl;
    std::cout << "Cap center: (" << cap_center.x() << ", " << cap_center.y() << ", " << cap_center.z() << ")" << std::endl;
    std::cout << "Face radius (atan(sqrt(2))): " << std::atan(std::sqrt(2.0)) << std::endl;
    std::cout << "EPS: " << 1e-15 << std::endl;
    std::cout << "Cap radius: " << face_radius << " radians" << std::endl;
    std::cout << "Cap radius: " << (face_radius * 180.0 / M_PI) << " degrees" << std::endl;
    
    // Create the cap
    S2Cap cap(cap_center, S1Angle::Radians(face_radius));
    
    // Create S2Cell for face 0 (+X face)
    S2Cell cell = S2Cell::FromFace(0);
    
    std::cout << "\nCell (face 0) details:" << std::endl;
    std::cout << "Cell ID: " << cell.id().ToString() << std::endl;
    std::cout << "Cell level: " << cell.level() << std::endl;
    std::cout << "Cell vertices:" << std::endl;
    for (int i = 0; i < 4; i++) {
        S2Point vertex = cell.GetVertex(i);
        std::cout << "  Vertex " << i << ": (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ")" << std::endl;
    }
    
    // Test MayIntersect
    bool may_intersect = cap.MayIntersect(cell);
    
    std::cout << "\n=== RESULT ===" << std::endl;
    std::cout << "cap.MayIntersect(cell): " << (may_intersect ? "true" : "false") << std::endl;
    
    // Additional debugging: test distance from cap center to cell vertices
    std::cout << "\n=== DEBUG INFO ===" << std::endl;
    std::cout << "Distance from cap center to cell vertices:" << std::endl;
    for (int i = 0; i < 4; i++) {
        S2Point vertex = cell.GetVertex(i);
        double distance = cap_center.Angle(vertex);
        std::cout << "  Vertex " << i << " distance: " << distance << " radians (" 
                  << (distance * 180.0 / M_PI) << " degrees)" << std::endl;
        std::cout << "  Vertex " << i << " within cap: " << (distance <= face_radius ? "true" : "false") << std::endl;
    }
    
    // Test distance to cell edges
    std::cout << "\nDistance from cap center to cell edges:" << std::endl;
    for (int i = 0; i < 4; i++) {
        S2Point v0 = cell.GetVertex(i);
        S2Point v1 = cell.GetVertex((i + 1) % 4);
        
        // Simple approach: test midpoint of edge
        S2Point edge_midpoint = (v0 + v1).Normalize();
        double distance_to_edge = cap_center.Angle(edge_midpoint);
        std::cout << "  Edge " << i << " midpoint distance: " << distance_to_edge << " radians (" 
                  << (distance_to_edge * 180.0 / M_PI) << " degrees)" << std::endl;
        std::cout << "  Edge " << i << " midpoint within cap: " << (distance_to_edge <= face_radius ? "true" : "false") << std::endl;
    }
    
    // Also test the other cap methods for comparison
    std::cout << "\n=== OTHER CAP METHODS ===" << std::endl;
    std::cout << "cap.Contains(cell): " << (cap.Contains(cell) ? "true" : "false") << std::endl;
    
    return 0;
}