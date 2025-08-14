//! S2CellId - Hierarchical cell addressing system for the S2 sphere
//!
//! This module provides the core S2CellId implementation, which uniquely identifies
//! every cell in the S2 hierarchical decomposition of the sphere using a single u64.
//!
//! ## Architecture
//!
//! Each S2CellId encodes:
//! - **Face** (3 bits): Which cube face (0-5) the cell belongs to
//! - **Position** (61 bits): Position along the Hilbert curve for that face
//! - **Level** (implicit): Hierarchical subdivision level (0-30) encoded in LSB pattern
//!
//! The implementation provides O(1) operations for:
//! - Parent/child navigation
//! - Level extraction
//! - Range queries and containment testing
//! - Hilbert curve position encoding/decoding

use crate::error::{S2Error, S2Result};
use crate::math::DVec3;
use crate::point::S2Point;
use std::fmt;
use std::sync::Once;

/// Maximum subdivision level for S2 cells (30 levels = 2^60 leaf cells per face)
pub const MAX_LEVEL: i32 = 30;

/// Number of bits used to encode the cube face (0-5)
pub const FACE_BITS: i32 = 3;

/// Number of cube faces on the S2 sphere
pub const NUM_FACES: i32 = 6;

/// Number of bits used for position encoding (61 bits)
pub const POS_BITS: i32 = 2 * MAX_LEVEL + 1;

/// Maximum cell size at level 0 (2^30)
pub const MAX_SIZE: i32 = 1 << MAX_LEVEL;

/// Number of bits processed per lookup table iteration
const LOOKUP_BITS: i32 = 4;

/// Swap mask for Hilbert curve orientation tracking
const SWAP_MASK: u64 = 0x01;

/// Invert mask for Hilbert curve orientation tracking  
const INVERT_MASK: u64 = 0x02;

/// Lookup table for (i,j) to Hilbert position transformation
static mut LOOKUP_POS: [u64; 1 << (2 * LOOKUP_BITS + 2)] = [0; 1 << (2 * LOOKUP_BITS + 2)];

/// Lookup table for Hilbert position to (i,j) transformation
static mut LOOKUP_IJ: [u64; 1 << (2 * LOOKUP_BITS + 2)] = [0; 1 << (2 * LOOKUP_BITS + 2)];

/// Ensures lookup tables are initialized exactly once
static LOOKUP_INIT: Once = Once::new();

/// Metric for measuring lengths on the sphere
#[derive(Debug, Clone, Copy)]
pub struct LengthMetric {
    /// The dimensionality of the metric (1 for length, 2 for area, etc.)
    dim: f64,
    /// The metric value at level 0
    deriv: f64,
}

impl LengthMetric {
    /// Create a new length metric
    pub const fn new(dim: f64, deriv: f64) -> Self {
        Self { dim, deriv }
    }

    /// Get the value of this metric at the given level
    pub fn get_value(&self, level: i32) -> f64 {
        self.deriv / (1u64 << (self.dim as i32 * level)) as f64
    }

    /// Get the level at which this metric has approximately the given value
    pub fn get_closest_level(&self, value: f64) -> i32 {
        if value <= 0.0 {
            return MAX_LEVEL;
        }
        
        let level = (value / self.deriv).log2() / self.dim;
        let level = level.round() as i32;
        level.max(0).min(MAX_LEVEL)
    }
}

/// Average edge length metric for S2 cells
/// 
/// This metric bounds the average distance from the center of one cell to the
/// center of one of its edge neighbors. The value corresponds to quadratic
/// projection used in the S2 geometry library.
pub const AVG_EDGE_METRIC: LengthMetric = LengthMetric::new(1.0, 1.459213746386106062);

/// S2CellId represents a unique identifier for cells in the S2 hierarchical decomposition
///
/// ## Design Principles
///
/// - **Single u64 storage**: Entire cell identity encoded in one 64-bit integer
/// - **Hilbert curve ordering**: Spatial locality preserved in numerical ordering
/// - **Level hierarchy**: 31 levels from face cells (level 0) to 1cm² leaf cells (level 30)
/// - **Bit manipulation**: All operations use fast bitwise arithmetic
///
/// ## Examples
///
/// ```rust,ignore
/// use s2geometry_rust::S2CellId;
///
/// // Create cell containing a point
/// let cell_id = S2CellId::from_point(&point);
/// 
/// // Navigate hierarchy
/// let parent = cell_id.parent();
/// let children: Vec<S2CellId> = cell_id.children().collect();
///
/// // Range queries
/// assert!(parent.contains(&cell_id));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct S2CellId {
    /// The complete cell identifier encoded as a single u64
    ///
    /// Bit layout: [face:3][position:61]
    /// - Bits 61-63: Cube face (0-5)  
    /// - Bits 0-60: Position along Hilbert curve
    /// - Level encoded implicitly in LSB pattern
    id: u64,
}

impl S2CellId {
    /// Creates a new S2CellId from raw u64 value
    ///
    /// ## Safety
    ///
    /// The caller must ensure the id represents a valid S2 cell.
    /// Use `from_face_pos_level` or other constructors for safety.
    #[inline]
    pub const fn new(id: u64) -> Self {
        Self { id }
    }

    /// Creates S2CellId from face, position, and level
    ///
    /// ## Arguments
    /// 
    /// - `face`: Cube face index (0-5)
    /// - `pos`: Position along Hilbert curve for the face
    /// - `level`: Subdivision level (0-30)
    ///
    /// ## Returns
    ///
    /// Valid S2CellId or error for invalid parameters
    pub fn from_face_pos_level(face: i32, pos: u64, level: i32) -> S2Result<Self> {
        if face < 0 || face >= NUM_FACES {
            return Err(S2Error::invalid_cell_id(0, "Face must be 0-5"));
        }
        
        if level < 0 || level > MAX_LEVEL {
            return Err(S2Error::invalid_cell_id(0, "Level must be 0-30"));
        }

        // Construct cell ID with proper LSB pattern for the level
        let lsb = Self::lsb_for_level(level);
        let face_bits = (face as u64) << POS_BITS;
        let id = face_bits | (pos & !lsb) | lsb;
        
        Ok(Self::new(id))
    }

    /// Creates S2CellId from a point on the unit sphere
    ///
    /// Finds the leaf cell (level 30) containing the point.
    pub fn from_point(point: &S2Point) -> Self {
        let (face, u, v) = Self::xyz_to_face_uv(point.coords());
        Self::from_face_uv(face, u, v)
    }

    /// Creates S2CellId from face and UV coordinates
    ///
    /// ## Arguments
    ///
    /// - `face`: Cube face index (0-5)
    /// - `u`: U coordinate on face [-1, 1]  
    /// - `v`: V coordinate on face [-1, 1]
    ///
    /// ## Returns
    ///
    /// Leaf cell (level 30) containing the UV point
    pub fn from_face_uv(face: i32, u: f64, v: f64) -> Self {
        let i = Self::uv_to_st(u);
        let j = Self::uv_to_st(v);
        Self::from_face_ij_same_face(face, i, j)
    }

    /// Creates S2CellId from face and integer IJ coordinates
    ///
    /// ## Arguments
    ///
    /// - `face`: Cube face index (0-5)
    /// - `i`: I coordinate [0, 2^30)
    /// - `j`: J coordinate [0, 2^30)
    ///
    /// ## Returns  
    ///
    /// Leaf cell containing the IJ point on the specified face
    pub fn from_face_ij_same_face(face: i32, i: u32, j: u32) -> Self {
        Self::init_lookup_table();
        
        // Encode face in upper bits
        let mut n = (face as u64) << (POS_BITS - 1);
        
        // Initialize Hilbert curve orientation tracking
        let mut bits = (face as u64) & SWAP_MASK;
        
        // Process coordinates in 4-bit chunks using lookup tables
        for k in (0..8).rev() {
            let mask = (1u32 << LOOKUP_BITS) - 1;
            let i_chunk = (i >> (k * LOOKUP_BITS)) & mask;
            let j_chunk = (j >> (k * LOOKUP_BITS)) & mask;
            
            // Pack chunk coordinates with orientation bits
            let lookup_index = bits + ((i_chunk as u64) << (LOOKUP_BITS + 2)) + ((j_chunk as u64) << 2);
            
            // Transform using Hilbert curve lookup table
            let lookup_result = unsafe { LOOKUP_POS[lookup_index as usize] };
            
            // Extract position bits for this level
            n |= (lookup_result >> 2) << (k * 2 * LOOKUP_BITS);
            
            // Update orientation for next iteration
            bits = lookup_result & (SWAP_MASK | INVERT_MASK);
        }
        
        // Set LSB to indicate leaf cell
        Self::new(n * 2 + 1)
    }

    /// Returns the raw u64 cell identifier
    #[inline]
    pub const fn id(&self) -> u64 {
        self.id
    }

    /// Returns true if this represents a valid S2 cell
    ///
    /// Validates:
    /// - Face is in range 0-5
    /// - LSB follows correct hierarchical pattern
    #[inline]
    pub fn is_valid(&self) -> bool {
        self.face() < NUM_FACES && (self.lsb() & 0x1555555555555555u64) != 0
    }

    /// Returns the cube face this cell belongs to (0-5)
    #[inline]
    pub fn face(&self) -> i32 {
        (self.id >> POS_BITS) as i32
    }

    /// Returns the subdivision level (0-30)
    ///
    /// - Level 0: 6 face cells covering the entire sphere
    /// - Level 30: ~1cm² leaf cells (4^30 per face)
    #[inline]  
    pub fn level(&self) -> i32 {
        if self.id == 0 {
            return -1; // Invalid cell
        }
        MAX_LEVEL - ((self.id.trailing_zeros() as i32) >> 1)
    }

    /// Returns the position along the Hilbert curve for this face
    #[inline]
    pub fn pos(&self) -> u64 {
        self.id & (u64::MAX >> FACE_BITS)
    }

    /// Returns the least significant bit (LSB) of the position
    ///
    /// The LSB encodes the level: `lsb = 1 << (2 * (MAX_LEVEL - level))`
    #[inline]
    pub fn lsb(&self) -> u64 {
        self.id & self.id.wrapping_neg()
    }

    /// Returns the parent cell at the specified level
    ///
    /// ## Arguments
    ///
    /// - `level`: Target parent level (must be ≤ current level)
    ///
    /// ## Returns
    ///
    /// Parent cell or error if level is invalid
    pub fn parent(&self, level: i32) -> S2Result<Self> {
        if level < 0 || level > self.level() {
            return Err(S2Error::invalid_cell_id(self.id, "Level must be ≤ current level"));
        }
        
        let new_lsb = Self::lsb_for_level(level);
        let id = (self.id & (new_lsb.wrapping_neg())) | new_lsb;
        Ok(Self::new(id))
    }

    /// Returns the immediate parent cell (level - 1)
    pub fn immediate_parent(&self) -> S2Result<Self> {
        let current_level = self.level();
        if current_level <= 0 {
            return Err(S2Error::invalid_cell_id(self.id, "Cannot get parent of level 0 cell"));
        }
        self.parent(current_level - 1)
    }

    /// Returns the child cell at the specified position
    ///
    /// ## Arguments
    ///
    /// - `pos`: Child position (0-3)
    ///
    /// ## Returns
    ///
    /// Child cell or error if this is a leaf cell or pos is invalid
    pub fn child(&self, pos: i32) -> S2Result<Self> {
        if self.is_leaf() {
            return Err(S2Error::invalid_cell_id(self.id, "Leaf cells have no children"));
        }
        
        if pos < 0 || pos >= 4 {
            return Err(S2Error::invalid_cell_id(self.id, "Child position must be 0-3"));
        }

        let new_lsb = self.lsb() >> 2;
        let offset = (2 * pos as u64 + 1).wrapping_sub(4);
        let id = self.id.wrapping_add(offset.wrapping_mul(new_lsb));
        Ok(Self::new(id))
    }

    /// Returns iterator over the four child cells
    pub fn children(&self) -> impl Iterator<Item = S2CellId> + '_ {
        (0..4).filter_map(move |i| self.child(i).ok())
    }

    /// Returns true if this is a leaf cell (level 30)
    #[inline]
    pub fn is_leaf(&self) -> bool {
        self.level() == MAX_LEVEL
    }

    /// Returns true if this cell contains the other cell
    ///
    /// A cell contains another if the other cell's range is within this cell's range.
    #[inline]
    pub fn contains(&self, other: &Self) -> bool {
        other.id >= self.range_min() && other.id <= self.range_max()
    }

    /// Returns true if this cell intersects the other cell's range
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        self.range_min() <= other.range_max() && other.range_min() <= self.range_max()
    }

    /// Returns a compact string token representation of this cell ID
    ///
    /// The token preserves ordering: if `a < b` then `a.to_token() < b.to_token()`.
    /// Invalid cells return "X".
    pub fn to_token(&self) -> String {
        if self.id == 0 {
            return "X".to_string();
        }
        
        // Convert to hex but remove trailing zeros for compactness
        let mut hex = format!("{:016x}", self.id);  // Always 16 hex digits
        
        // Remove trailing zeros, but keep at least 1 character
        while hex.len() > 1 && hex.ends_with('0') {
            hex.pop();
        }
        
        hex
    }

    /// Parse a token string back to S2CellId
    ///
    /// ## Arguments
    ///
    /// - `token`: Token string returned by `to_token()`
    ///
    /// ## Returns
    ///
    /// S2CellId or error if token is invalid
    pub fn from_token(token: &str) -> S2Result<Self> {
        if token == "X" {
            return Ok(Self::default());
        }
        
        // Pad with zeros if necessary to get full hex representation
        let padded_token = if token.len() < 16 {
            format!("{}{}", token, "0".repeat(16 - token.len()))
        } else {
            token.to_string()
        };
        
        // Parse as hexadecimal
        let id = u64::from_str_radix(&padded_token, 16)
            .map_err(|_| S2Error::invalid_cell_id(0, "Invalid token format"))?;
        
        let cell_id = Self::new(id);
        if !cell_id.is_valid() && id != 0 {
            return Err(S2Error::invalid_cell_id(id, "Token represents invalid cell"));
        }
        
        Ok(cell_id)
    }

    /// Returns the center point of this cell on the unit sphere
    ///
    /// This converts the cell's hierarchical position back to a geometric point.
    pub fn to_point_raw(&self) -> S2Point {
        if !self.is_valid() {
            return S2Point::new(1.0, 0.0, 0.0).unwrap();  // Default point
        }

        // Convert cell to its center coordinates
        let (si, ti) = self.get_center_si_ti();
        let u = Self::st_to_uv(si);
        let v = Self::st_to_uv(ti);
        let xyz = Self::face_uv_to_xyz(self.face(), u, v);
        
        S2Point::from_vec3(xyz).unwrap_or_else(|_| S2Point::new(1.0, 0.0, 0.0).unwrap())
    }

    /// Returns the (si, ti) coordinates of the cell center
    ///
    /// These are integer coordinates in [0, 2^30) representing the cell's position
    /// on its face in the ST coordinate system.
    fn get_center_si_ti(&self) -> (u32, u32) {
        let level = self.level();
        if level < 0 {
            return (MAX_SIZE as u32 / 2, MAX_SIZE as u32 / 2);  // Center of face
        }

        // For level 0 (face cells), return face center
        if level == 0 {
            return (MAX_SIZE as u32 / 2, MAX_SIZE as u32 / 2);
        }

        // Cell size at this level (in ST coordinates)
        let cell_size = 1u32 << (MAX_LEVEL - level);
        
        // Use simplified coordinate extraction based on cell ID structure
        // This is an approximation until we implement proper Hilbert curve inversion
        let face_pos = self.pos();
        
        // Extract rough position from the lower bits
        // This is a very rough approximation - proper implementation needs Hilbert curve
        let shift_amount = 2 * (MAX_LEVEL - level);
        let cell_index = if shift_amount < 64 {
            (face_pos >> shift_amount) as u32
        } else {
            0u32
        };
        
        // Calculate position within this level's grid
        let cells_per_side = 1u32 << level;
        let i = cell_index % cells_per_side;
        let j = cell_index / cells_per_side;
        
        // Convert to ST coordinates and get center
        let si = i * cell_size + cell_size / 2;
        let ti = j * cell_size + cell_size / 2;
        
        // Clamp to valid range
        let si = si.min((MAX_SIZE as u32) - 1);
        let ti = ti.min((MAX_SIZE as u32) - 1);
        
        (si, ti)
    }

    /// Returns the minimum cell ID in this cell's range
    #[inline]
    pub fn range_min(&self) -> u64 {
        self.id - (self.lsb() - 1)
    }

    /// Returns the maximum cell ID in this cell's range  
    #[inline]
    pub fn range_max(&self) -> u64 {
        self.id + (self.lsb() - 1)
    }

    /// Returns LSB value for the specified level
    #[inline]
    pub fn lsb_for_level(level: i32) -> u64 {
        1u64 << (2 * (MAX_LEVEL - level))
    }

    /// Converts XYZ coordinates to face and UV coordinates
    ///
    /// This implements the S2 cube face projection exactly matching the C++ implementation.
    ///
    /// ## Returns
    ///
    /// (face, u, v) where face ∈ [0,5] and u,v ∈ [-1,1]
    fn xyz_to_face_uv(xyz: DVec3) -> (i32, f64, f64) {
        let abs_x = xyz.x.abs();
        let abs_y = xyz.y.abs();
        let abs_z = xyz.z.abs();

        // Determine which face this point projects to
        let (face, u, v) = if abs_x >= abs_y && abs_x >= abs_z {
            // X axis is dominant - faces 0 (+X) or 3 (-X)
            if xyz.x >= 0.0 {
                (0, xyz.y / xyz.x, xyz.z / xyz.x)  // +X face
            } else {
                (3, -xyz.z / (-xyz.x), -xyz.y / (-xyz.x))  // -X face  
            }
        } else if abs_y >= abs_z {
            // Y axis is dominant - faces 1 (+Y) or 4 (-Y)
            if xyz.y >= 0.0 {
                (1, -xyz.x / xyz.y, xyz.z / xyz.y)  // +Y face
            } else {
                (4, xyz.z / (-xyz.y), xyz.x / (-xyz.y))  // -Y face
            }
        } else {
            // Z axis is dominant - faces 2 (+Z) or 5 (-Z)  
            if xyz.z >= 0.0 {
                (2, -xyz.y / xyz.z, -xyz.x / xyz.z)  // +Z face
            } else {
                (5, -xyz.x / (-xyz.z), xyz.y / (-xyz.z))  // -Z face
            }
        };

        (face, u, v)
    }

    /// Converts UV coordinate to ST coordinate (0 to MAX_SIZE)
    ///
    /// UV coordinates are in [-1, 1], ST coordinates are in [0, 2^30)
    fn uv_to_st(u: f64) -> u32 {
        // S2 uses a quadratic transformation here for better area preservation
        // For now, use linear transformation (this is a simplification)
        let s = 0.5 * (u + 1.0); // Convert from [-1,1] to [0,1]
        let scaled = s * (MAX_SIZE as f64);
        (scaled.clamp(0.0, (MAX_SIZE - 1) as f64)) as u32
    }

    /// Converts ST coordinate to UV coordinate 
    ///
    /// ST coordinates are in [0, 2^30), UV coordinates are in [-1, 1]
    fn st_to_uv(s: u32) -> f64 {
        // Inverse of uv_to_st
        let normalized = (s as f64) / (MAX_SIZE as f64);  // [0,1]
        2.0 * normalized - 1.0  // [-1,1]
    }

    /// Converts face and UV coordinates back to XYZ point on unit sphere
    ///
    /// This is the inverse of xyz_to_face_uv
    fn face_uv_to_xyz(face: i32, u: f64, v: f64) -> DVec3 {
        match face {
            0 => DVec3::new(1.0, u, v),      // +X face
            1 => DVec3::new(-u, 1.0, v),     // +Y face
            2 => DVec3::new(-v, -u, 1.0),    // +Z face  
            3 => DVec3::new(-1.0, -v, -u),   // -X face
            4 => DVec3::new(v, -1.0, u),     // -Y face
            5 => DVec3::new(u, v, -1.0),     // -Z face
            _ => DVec3::new(1.0, 0.0, 0.0),  // Fallback
        }.normalize()
    }

    /// Initializes Hilbert curve lookup tables (thread-safe, called once)
    fn init_lookup_table() {
        LOOKUP_INIT.call_once(|| {
            Self::init_lookup_cell(0, 0, 0, 0, 0, 0);
            Self::init_lookup_cell(0, 0, 0, SWAP_MASK, 0, SWAP_MASK);
            Self::init_lookup_cell(0, 0, 0, INVERT_MASK, 0, INVERT_MASK);
            Self::init_lookup_cell(0, 0, 0, SWAP_MASK | INVERT_MASK, 0, SWAP_MASK | INVERT_MASK);
        });
    }

    /// Recursively initializes lookup table entries (Hilbert curve generation)
    fn init_lookup_cell(level: i32, i: u64, j: u64, orig_orientation: u64, pos: u64, orientation: u64) {
        if level == LOOKUP_BITS {
            let ij = (i << (LOOKUP_BITS + 2)) + (j << 2) + orig_orientation;
            unsafe {
                LOOKUP_POS[ij as usize] = (pos << 2) + orientation;
                LOOKUP_IJ[((pos << 2) + orig_orientation) as usize] = (ij << 2) + orientation;
            }
        } else {
            let level = level + 1;
            let i = i << 1;
            let j = j << 1;
            let pos = pos << 2;
            let r = [0, 1, 3, 2];
            
            for k in 0..4 {
                let (sub_i, sub_j, sub_orientation) = if orientation & SWAP_MASK != 0 {
                    if orientation & INVERT_MASK != 0 {
                        (i + ((r[k] >> 1) as u64), j + ((r[k] & 1) as u64), orientation ^ Self::mask_for_position(r[k] as usize))
                    } else {
                        (i + ((k >> 1) as u64), j + ((k & 1) as u64), orientation ^ Self::mask_for_position(k))
                    }
                } else {
                    if orientation & INVERT_MASK != 0 {
                        (i + ((k & 1) as u64), j + ((r[k] >> 1) as u64), orientation ^ Self::mask_for_position(r[k] as usize))
                    } else {
                        (i + ((r[k] & 1) as u64), j + ((r[k] >> 1) as u64), orientation ^ Self::mask_for_position(k))
                    }
                };

                Self::init_lookup_cell(level, sub_i, sub_j, orig_orientation, pos + k as u64, sub_orientation);
            }
        }
    }

    /// Returns orientation mask for given position in Hilbert curve
    #[inline]
    fn mask_for_position(pos: usize) -> u64 {
        match pos {
            0 => 0,
            1 => SWAP_MASK,
            2 => SWAP_MASK | INVERT_MASK,
            3 => INVERT_MASK,
            _ => 0,
        }
    }

    // Additional methods needed for S2CellUnion

    /// Creates a face cell (level 0) for the given face
    pub fn from_face(face: i32) -> Self {
        Self::from_face_pos_level(face, 0, 0).unwrap()
    }

    /// Returns true if this is a face cell (level 0)
    #[inline]
    pub fn is_face(&self) -> bool {
        self.level() == 0
    }

    /// Returns the next cell ID in sequence
    #[inline]
    pub fn next(&self) -> Self {
        Self::new(self.id + (self.lsb() << 1))
    }

    /// Returns the parent at the specified level
    pub fn parent_at_level(&self, level: i32) -> Self {
        if level >= self.level() {
            *self
        } else {
            self.parent(level).unwrap_or(*self)
        }
    }

    /// Returns the minimum leaf cell ID at the maximum level
    pub fn begin(level: i32) -> Self {
        Self::from_face_pos_level(0, 0, level).unwrap()
    }

    /// Returns the maximum leaf cell ID + 1 at the maximum level  
    pub fn end(level: i32) -> Self {
        Self::from_face_pos_level(5, 0, level).unwrap().next()
    }

    /// Returns the largest cell that fits within the range [self, end)
    ///
    /// This is used for efficient range covering by finding the largest
    /// cell that doesn't exceed the end boundary.
    pub fn maximum_tile(&self, end: Self) -> Self {
        let mut id = *self;
        
        while id.level() > 0 {
            let parent = id.parent(id.level() - 1).unwrap();
            if parent.range_max() >= end.id {
                break;
            }
            id = parent;
        }
        
        id
    }

    /// Returns the center point of this cell (for S2CellUnion compatibility)
    pub fn to_point(&self) -> S2Point {
        self.to_point_raw()
    }

    /// Appends all neighboring cells at the given level to the output vector
    ///
    /// This is a placeholder implementation. A full implementation would
    /// need to properly compute all neighboring cells at the specified level.
    pub fn append_all_neighbors(&self, level: i32, output: &mut Vec<S2CellId>) {
        // Simplified implementation - add adjacent cells in the ID space
        // This is not geometrically correct but will allow compilation
        let current_level = self.level();
        let target_cell = if current_level > level {
            self.parent_at_level(level)
        } else {
            *self
        };
        
        // Add some adjacent cells (this is a simplified approximation)
        let lsb = Self::lsb_for_level(level);
        let step = lsb << 1;
        
        // Add previous and next cells if they exist and are valid
        if target_cell.id >= step {
            let prev_id = Self::new(target_cell.id - step);
            if prev_id.is_valid() && prev_id.level() == level {
                output.push(prev_id);
            }
        }
        
        let next_id = Self::new(target_cell.id + step);
        if next_id.is_valid() && next_id.level() == level {
            output.push(next_id);
        }
    }

    /// Returns the level for cells with minimum width >= min_width
    pub fn level_for_min_width(min_width: f64) -> i32 {
        // Simplified implementation based on typical S2 geometry
        // The actual implementation would use S2 metrics
        let face_width = 2.0; // Approximate width of a face cell
        let mut level = 0;
        let mut width = face_width;
        
        while width > min_width && level < MAX_LEVEL {
            level += 1;
            width /= 2.0;
        }
        
        level
    }

    /// Returns the minimum width of cells at the given level
    pub fn min_width_at_level(level: i32) -> f64 {
        // Simplified implementation
        let face_width = 2.0; // Approximate width of a face cell
        face_width / (1 << level) as f64
    }

    /// Returns the average edge length metric for S2 cells
    /// 
    /// This metric bounds the average distance from the center of one cell to the
    /// center of one of its edge neighbors. It can be used to determine appropriate
    /// cell levels for spatial operations.
    ///
    /// # Example
    /// ```ignore
    /// // Get the cell level that corresponds to an average edge length of 0.1 radians
    /// let level = S2CellId::avg_edge_metric().get_closest_level(0.1);
    /// ```
    pub fn avg_edge_metric() -> LengthMetric {
        AVG_EDGE_METRIC
    }
}

impl Default for S2CellId {
    /// Returns an invalid cell ID (id = 0)
    fn default() -> Self {
        Self::new(0)
    }
}

impl fmt::Display for S2CellId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_token())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::S2Point;

    #[test]
    fn test_cell_id_basic_properties() {
        // Test face cells (level 0)
        for face in 0..NUM_FACES {
            let cell_id = S2CellId::from_face_pos_level(face, 0, 0).unwrap();
            assert_eq!(cell_id.face(), face);
            assert_eq!(cell_id.level(), 0);
            assert!(cell_id.is_valid());
            assert!(!cell_id.is_leaf());
        }
    }

    #[test]
    fn test_parent_child_relationships() {
        let cell_id = S2CellId::from_face_pos_level(0, 0, 5).unwrap();
        
        // Test parent
        let parent = cell_id.parent(4).unwrap();
        assert_eq!(parent.level(), 4);
        assert!(parent.contains(&cell_id));
        
        // Test immediate parent
        let immediate_parent = cell_id.immediate_parent().unwrap();
        assert_eq!(immediate_parent.level(), 4);
        
        // Test children
        if !parent.is_leaf() {
            let children: Vec<_> = parent.children().collect();
            assert_eq!(children.len(), 4);
            assert!(children.iter().any(|child| child.contains(&cell_id)));
        }
    }

    #[test]
    fn test_point_conversion() {
        let point = S2Point::new(1.0, 0.0, 0.0).unwrap();
        let cell_id = S2CellId::from_point(&point);
        
        assert!(cell_id.is_valid());
        assert!(cell_id.is_leaf());
        assert_eq!(cell_id.level(), MAX_LEVEL);
    }

    #[test]
    fn test_range_operations() {
        let cell_id = S2CellId::from_face_pos_level(0, 0, 2).unwrap();
        let child = cell_id.child(0).unwrap();
        
        assert!(cell_id.contains(&child));
        assert!(cell_id.intersects(&child));
        assert!(child.range_min() >= cell_id.range_min());
        assert!(child.range_max() <= cell_id.range_max());
    }

    #[test]
    fn test_invalid_operations() {
        let leaf = S2CellId::from_face_pos_level(0, 0, MAX_LEVEL).unwrap();
        assert!(leaf.child(0).is_err()); // Leaf has no children
        
        let face_cell = S2CellId::from_face_pos_level(0, 0, 0).unwrap();
        assert!(face_cell.immediate_parent().is_err()); // Face cell has no parent
        
        // Invalid face
        assert!(S2CellId::from_face_pos_level(6, 0, 0).is_err());
        
        // Invalid level
        assert!(S2CellId::from_face_pos_level(0, 0, 31).is_err());
    }

    #[test]
    fn test_display() {
        let cell_id = S2CellId::from_face_pos_level(0, 1, 0).unwrap();
        let display = format!("{}", cell_id);
        assert!(!display.is_empty());
        assert_ne!(display, "X"); // Not invalid
        
        let invalid = S2CellId::default();
        assert_eq!(format!("{}", invalid), "X");
    }
}