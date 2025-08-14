# S2 Geometry Rust - Missing C++ Test Coverage

**Status**: 5/121 C++ tests ported (4% coverage) - 116 tests remaining

This document tracks all missing C++ test files that need to be ported to complete the S2 Geometry Rust implementation.

## ✅ Already Ported (5 tests)

- [x] `s1angle_test.cc` → `test_s1angle_port.rs` (24 tests)
- [x] `s2coords_test.cc` → `test_s2coords_port.rs` (23 tests) 
- [x] `s2edge_crossings_test.cc` → `test_s2edge_crossings_port.rs` (9 tests)
- [x] `s2point_test.cc` → `test_s2point_port.rs` (20 tests)
- [x] `s2predicates_test.cc` → `test_s2predicates_port.rs` (7 tests)

---

## 🔥 TIER 1: FOUNDATIONAL TESTS (Critical Dependencies)

### Priority A: Core Data Types (MUST IMPLEMENT FIRST)

- [ ] **`s2cell_id_test.cc`** (820 lines) - **CRITICAL**
  - Hierarchical cell addressing system
  - Blocks: All spatial indexing, cell operations, geographic queries
  
- [ ] **`s2cell_test.cc`** (790 lines) - **CRITICAL**
  - Individual S2 cells on sphere
  - Blocks: All region operations, covering algorithms, spatial analysis

- [ ] **`s1chord_angle_test.cc`** (279 lines) - **HIGH**
  - Efficient angle representation used throughout distance calculations
  - Blocks: Distance calculations, geometric measurements

### Priority B: Core Measurements & Utilities

- [ ] **`s2measures_test.cc`** (158 lines) - **HIGH**
  - Area, length, angle calculations - core geometric operations

- [ ] **`s2metrics_test.cc`** (141 lines) - **HIGH**
  - Cell size and projection metrics for accurate spatial analysis

- [ ] **`s1interval_test.cc`** (478 lines) - **HIGH**
  - 1D interval operations - foundation for 2D rectangles

- [ ] **`s2earth_test.cc`** (166 lines) - **MEDIUM**
  - Earth-specific constants and calculations for geo applications

- [ ] **`r1interval_test.cc`** (189 lines) - **MEDIUM**
  - Real number intervals - supporting utility

- [ ] **`r2rect_test.cc`** (239 lines) - **MEDIUM**
  - 2D rectangles in plane for planar projections

---

## 🏗️ TIER 2: CORE GEOMETRY TESTS (Essential Structures)

### Priority A: Core Geometric Structures

- [ ] **`s2loop_test.cc`** (1,359 lines) - **CRITICAL**
  - Closed polygonal loops - foundation for all polygon operations

- [ ] **`s2cap_test.cc`** (407 lines) - **HIGH**
  - Spherical caps (circular regions) - basic region type

- [ ] **`s2cell_union_test.cc`** (879 lines) - **HIGH**
  - Collections of cells for efficient region representation

- [ ] **`s2latlng_test.cc`** (247 lines) - **HIGH**
  - Latitude/longitude coordinates - standard geographic interface

- [ ] **`s2latlng_rect_test.cc`** (1,063 lines) - **HIGH**
  - Geographic bounding rectangles for geographic queries

### Priority B: Advanced Shapes

- [ ] **`s2polygon_test.cc`** (3,276 lines) - **HIGH**
  - Complex polygons with holes - most complex geometric primitive

- [ ] **`s2polyline_test.cc`** (729 lines) - **MEDIUM**
  - Connected line segments for linear features

- [ ] **`s2region_test.cc`** (349 lines) - **MEDIUM**
  - Abstract region interface for polymorphic regions

### Priority C: Spatial Indexing Foundation

- [ ] **`mutable_s2shape_index_test.cc`** (825 lines) - **HIGH**
  - Dynamic spatial index - core indexing structure

- [ ] **`s2shape_index_test.cc`** (48 lines) - **HIGH**
  - Read-only spatial index for query performance

- [ ] **`s2region_coverer_test.cc`** (560 lines) - **HIGH**
  - Cell covering algorithms - approximation system

---

## 🧩 TIER 3: ADVANCED FUNCTIONALITY (Complex Operations)

### Priority A: Boolean Operations

- [ ] **`s2boolean_operation_test.cc`** (2,348 lines) - **MEDIUM**
  - Union, intersection, difference operations

- [ ] **`s2builder_test.cc`** (1,819 lines) - **MEDIUM**
  - Robust polygon construction handling degenerate cases

- [ ] **`s2buffer_operation_test.cc`** (585 lines) - **MEDIUM**
  - Buffer operations for geometric morphology

### Priority B: Query Operations

- [ ] **`s2closest_edge_query_test.cc`** (878 lines) - **MEDIUM**
  - Nearest edge queries for spatial queries

- [ ] **`s2closest_point_query_test.cc`** (339 lines) - **MEDIUM**
  - Nearest point queries

- [ ] **`s2contains_point_query_test.cc`** (204 lines) - **MEDIUM**
  - Point-in-shape queries

- [ ] **`s2closest_edge_query_base_test.cc`** (127 lines) - **LOW**
  - Base functionality for edge queries

- [ ] **`s2closest_point_query_base_test.cc`** (82 lines) - **LOW**
  - Base functionality for point queries

- [ ] **`s2closest_cell_query_test.cc`** (291 lines) - **LOW**
  - Nearest cell queries

- [ ] **`s2closest_cell_query_base_test.cc`** (45 lines) - **LOW**
  - Base functionality for cell queries

### Priority C: Specialized Operations

- [ ] **`s2edge_distances_test.cc`** (753 lines) - **LOW**
  - Distance calculations between edges

- [ ] **`s2edge_tessellator_test.cc`** (577 lines) - **LOW**
  - Edge subdivision algorithms

- [ ] **`s2fractal_test.cc`** (160 lines) - **LOW**
  - Fractal curve generation for testing

---

## 🔧 TIER 4: UTILITIES & SPECIALIZED FEATURES

### Builder Utilities

- [ ] `s2builderutil_closed_set_normalizer_test.cc` (186 lines)
- [ ] `s2builderutil_find_polygon_degeneracies_test.cc` (201 lines)
- [ ] `s2builderutil_get_snapped_winding_delta_test.cc` (89 lines)
- [ ] `s2builderutil_lax_polygon_layer_test.cc` (391 lines)
- [ ] `s2builderutil_lax_polyline_layer_test.cc` (209 lines)
- [ ] `s2builderutil_s2point_vector_layer_test.cc` (118 lines)
- [ ] `s2builderutil_s2polygon_layer_test.cc` (376 lines)
- [ ] `s2builderutil_s2polyline_layer_test.cc` (312 lines)
- [ ] `s2builderutil_s2polyline_vector_layer_test.cc` (157 lines)
- [ ] `s2builderutil_snap_functions_test.cc` (198 lines)
- [ ] `s2builderutil_testing_test.cc` (65 lines)

### Shape Utilities

- [ ] `s2shapeutil_build_polygon_boundaries_test.cc` (278 lines)
- [ ] `s2shapeutil_coding_test.cc` (89 lines)
- [ ] `s2shapeutil_contains_brute_force_test.cc` (67 lines)
- [ ] `s2shapeutil_conversion_test.cc` (213 lines)
- [ ] `s2shapeutil_count_edges_test.cc` (44 lines)
- [ ] `s2shapeutil_count_vertices_test.cc` (44 lines)
- [ ] `s2shapeutil_edge_iterator_test.cc` (86 lines)
- [ ] `s2shapeutil_edge_wrap_test.cc` (49 lines)
- [ ] `s2shapeutil_get_reference_point_test.cc` (128 lines)
- [ ] `s2shapeutil_shape_edge_id_test.cc` (84 lines)
- [ ] `s2shapeutil_visit_crossing_edge_pairs_test.cc` (197 lines)

### Shape Types

- [ ] `s2edge_vector_shape_test.cc` (169 lines)
- [ ] `s2lax_loop_shape_test.cc` (198 lines)
- [ ] `s2lax_polygon_shape_test.cc` (412 lines)
- [ ] `s2lax_polyline_shape_test.cc` (263 lines)
- [ ] `s2point_vector_shape_test.cc` (131 lines)
- [ ] `s2wrapped_shape_test.cc` (66 lines)

### Indexing & Data Structures

- [ ] `s2cell_index_test.cc` (203 lines)
- [ ] `s2cell_iterator_join_test.cc` (219 lines)
- [ ] `s2cell_iterator_testing_test.cc` (49 lines)
- [ ] `s2cell_range_iterator_test.cc` (145 lines)
- [ ] `s2density_tree_test.cc` (165 lines)
- [ ] `s2padded_cell_test.cc` (124 lines)

### Measurements & Analysis

- [ ] `s2centroids_test.cc` (204 lines)
- [ ] `s2loop_measures_test.cc` (422 lines)
- [ ] `s2polyline_measures_test.cc` (131 lines)
- [ ] `s2shape_measures_test.cc` (188 lines)
- [ ] `s2shape_index_measures_test.cc` (94 lines)

### Geographic & Projections

- [ ] `s2latlng_rect_bounder_test.cc` (135 lines)
- [ ] `s2point_compression_test.cc` (143 lines)
- [ ] `s2projections_test.cc` (295 lines)
- [ ] `s2r2rect_test.cc` (89 lines)

### Advanced Queries

- [ ] `s2chain_interpolation_query_test.cc` (188 lines)
- [ ] `s2contains_vertex_query_test.cc` (157 lines)
- [ ] `s2convex_hull_query_test.cc` (203 lines)
- [ ] `s2crossing_edge_query_test.cc` (236 lines)
- [ ] `s2furthest_edge_query_test.cc` (370 lines)
- [ ] `s2hausdorff_distance_query_test.cc` (134 lines)
- [ ] `s2max_distance_targets_test.cc` (89 lines)
- [ ] `s2min_distance_targets_test.cc` (267 lines)
- [ ] `s2shape_nesting_query_test.cc` (178 lines)
- [ ] `s2validation_query_test.cc` (82 lines)

### Regions & Indexing

- [ ] `s2point_index_test.cc` (267 lines)
- [ ] `s2point_region_test.cc` (87 lines)
- [ ] `s2region_sharder_test.cc` (189 lines)
- [ ] `s2region_term_indexer_test.cc` (267 lines)
- [ ] `s2region_union_test.cc` (67 lines)
- [ ] `s2shape_index_buffered_region_test.cc` (123 lines)
- [ ] `s2shape_index_region_test.cc` (89 lines)

### Specialized Algorithms

- [ ] `s2edge_clipping_test.cc` (467 lines)
- [ ] `s2edge_crosser_test.cc` (179 lines)
- [ ] `s2polyline_alignment_test.cc` (278 lines)
- [ ] `s2polyline_simplifier_test.cc` (188 lines)
- [ ] `s2wedge_relations_test.cc` (389 lines)
- [ ] `s2winding_operation_test.cc` (134 lines)

### Encoding & Compression

- [ ] `encoded_s2cell_id_vector_test.cc` (378 lines)
- [ ] `encoded_s2point_vector_test.cc` (668 lines)
- [ ] `encoded_s2shape_index_test.cc` (211 lines)
- [ ] `encoded_string_vector_test.cc` (178 lines)
- [ ] `encoded_uint_vector_test.cc` (189 lines)

### Testing & Development Tools

- [ ] `s2error_test.cc` (45 lines)
- [ ] `s2memory_tracker_test.cc` (67 lines)
- [ ] `s2pointutil_test.cc` (123 lines)
- [ ] `s2random_test.cc` (78 lines)
- [ ] `s2text_format_test.cc` (267 lines)
- [ ] `gmock_matchers_test.cc` (89 lines)

### Internal Components

- [ ] `internal/s2disjoint_set_test.cc` (67 lines)
- [ ] `internal/s2index_cell_data_test.cc` (44 lines)

### Utility Data Structures

- [ ] `id_set_lexicon_test.cc` (145 lines)
- [ ] `sequence_lexicon_test.cc` (167 lines)
- [ ] `value_lexicon_test.cc` (134 lines)

---

## 📊 Summary Statistics

- **Total C++ tests**: 121 files
- **Already ported**: 5 files (4%)
- **Remaining**: 116 files (96%)
- **Estimated total lines**: ~47,000 lines of test code
- **Critical path blocking**: 4 foundational tests (s2cell_id, s2cell, s1chord_angle, s2measures)

## 🎯 Implementation Strategy

### Phase 1: Unblock Development (TIER 1)
Focus on the 9 foundational tests that enable all higher-level functionality.

### Phase 2: Core Geometry (TIER 2) 
Implement essential geometric structures (loops, caps, regions).

### Phase 3: Advanced Features (TIER 3+)
Add complex operations, queries, and specialized utilities.

### Phase 4: Complete Coverage
Fill in remaining utilities, encodings, and specialized features.

---

*Last updated: 2025-01-11*
*Progress tracking: See main TODO list in project*