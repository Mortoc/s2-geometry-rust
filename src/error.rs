//! Error types for S2 geometry operations
//!
//! This module defines the comprehensive error handling for the S2 geometry library,
//! covering invalid input, numerical precision issues, and geometric impossibilities.

use thiserror::Error;

/// Result type alias for S2 geometry operations
pub type S2Result<T> = Result<T, S2Error>;

/// Comprehensive error types for S2 geometry operations
#[derive(Error, Debug, Clone, PartialEq)]
pub enum S2Error {
    /// Invalid point coordinates (not normalized)
    #[error("Invalid point: {reason}")]
    InvalidPoint { 
        /// Description of why the point is invalid
        reason: String 
    },

    /// Invalid latitude value (must be in [-π/2, π/2])
    #[error("Invalid latitude: {value} (must be in [-π/2, π/2])")]
    InvalidLatitude { 
        /// The invalid latitude value
        value: f64 
    },

    /// Invalid longitude value (must be in [-π, π])
    #[error("Invalid longitude: {value} (must be in [-π, π])")]
    InvalidLongitude { 
        /// The invalid longitude value
        value: f64 
    },

    /// Invalid angle value
    #[error("Invalid angle: {reason}")]
    InvalidAngle { 
        /// Description of why the angle is invalid
        reason: String 
    },

    /// Invalid cell ID
    #[error("Invalid S2CellId: {cell_id:#018x} - {reason}")]
    InvalidCellId { 
        /// The invalid cell ID value
        cell_id: u64, 
        /// Description of why the cell ID is invalid
        reason: String 
    },

    /// Invalid cell level
    #[error("Invalid cell level: {level} (must be in [0, {max_level}])")]
    InvalidCellLevel { 
        /// The invalid level value
        level: i32, 
        /// The maximum allowed level
        max_level: i32 
    },

    /// Invalid face number
    #[error("Invalid face: {face} (must be in [0, 5])")]
    InvalidFace { 
        /// The invalid face number
        face: i32 
    },

    /// Invalid loop (self-intersecting, duplicate vertices, etc.)
    #[error("Invalid loop: {reason}")]
    InvalidLoop { 
        /// Description of why the loop is invalid
        reason: String 
    },

    /// Invalid polygon (invalid loops, incorrect nesting, etc.)
    #[error("Invalid polygon: {reason}")]
    InvalidPolygon { 
        /// Description of why the polygon is invalid
        reason: String 
    },

    /// Invalid polyline
    #[error("Invalid polyline: {reason}")]
    InvalidPolyline { 
        /// Description of why the polyline is invalid
        reason: String 
    },

    /// Insufficient numerical precision for computation
    #[error("Computation failed: insufficient precision for {operation}")]
    InsufficientPrecision { 
        /// The operation that failed due to precision limits
        operation: String 
    },

    /// Geometric degeneracy (e.g., zero-area polygon, coincident points)
    #[error("Geometric degeneracy: {reason}")]
    GeometricDegeneracy { 
        /// Description of the geometric degeneracy
        reason: String 
    },

    /// Index construction failure
    #[error("Index construction failed: {reason}")]
    IndexError { 
        /// Description of the index error
        reason: String 
    },

    /// Query operation failure
    #[error("Query failed: {reason}")]
    QueryError { 
        /// Description of the query error
        reason: String 
    },

    /// Boolean operation failure
    #[error("Boolean operation failed: {operation} - {reason}")]
    BooleanOperationError { 
        /// The boolean operation that failed
        operation: String, 
        /// Description of why the operation failed
        reason: String 
    },

    /// Builder operation failure
    #[error("Builder operation failed: {reason}")]
    BuilderError { 
        /// Description of the builder error
        reason: String 
    },

    /// Internal library error (should not occur in correct usage)
    #[error("Internal error: {reason} - please report this bug")]
    InternalError { 
        /// Description of the internal error
        reason: String 
    },
}

impl S2Error {
    /// Create an invalid point error
    ///
    /// # Example
    /// ```rust,ignore
    /// let error = S2Error::invalid_point("coordinates must be normalized");
    /// ```
    pub fn invalid_point(reason: impl Into<String>) -> Self {
        S2Error::InvalidPoint {
            reason: reason.into(),
        }
    }

    /// Create an invalid latitude error
    pub fn invalid_latitude(value: f64) -> Self {
        S2Error::InvalidLatitude { value }
    }

    /// Create an invalid longitude error
    pub fn invalid_longitude(value: f64) -> Self {
        S2Error::InvalidLongitude { value }
    }

    /// Create an invalid cell ID error
    pub fn invalid_cell_id(cell_id: u64, reason: impl Into<String>) -> Self {
        S2Error::InvalidCellId {
            cell_id,
            reason: reason.into(),
        }
    }

    /// Create an invalid cell level error
    pub fn invalid_cell_level(level: i32, max_level: i32) -> Self {
        S2Error::InvalidCellLevel { level, max_level }
    }

    /// Create an invalid face error
    pub fn invalid_face(face: i32) -> Self {
        S2Error::InvalidFace { face }
    }

    /// Create an insufficient precision error
    pub fn insufficient_precision(operation: impl Into<String>) -> Self {
        S2Error::InsufficientPrecision {
            operation: operation.into(),
        }
    }

    /// Create a geometric degeneracy error
    pub fn geometric_degeneracy(reason: impl Into<String>) -> Self {
        S2Error::GeometricDegeneracy {
            reason: reason.into(),
        }
    }

    /// Create an internal error (for impossible conditions)
    pub fn internal_error(reason: impl Into<String>) -> Self {
        S2Error::InternalError {
            reason: reason.into(),
        }
    }

    /// Create an invalid loop error
    pub fn invalid_loop(reason: impl Into<String>) -> Self {
        S2Error::InvalidLoop {
            reason: reason.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_creation() {
        let error = S2Error::invalid_point("not normalized");
        assert!(matches!(error, S2Error::InvalidPoint { .. }));
        assert!(error.to_string().contains("not normalized"));
    }

    #[test]
    fn test_error_display() {
        let error = S2Error::invalid_latitude(1.8);
        let message = error.to_string();
        assert!(message.contains("Invalid latitude"));
        assert!(message.contains("1.8"));
        assert!(message.contains("[-π/2, π/2]"));
    }

    #[test]
    fn test_s2result_type() {
        let success: S2Result<i32> = Ok(42);
        let failure: S2Result<i32> = Err(S2Error::invalid_face(7));
        
        assert!(success.is_ok());
        assert!(failure.is_err());
    }
}