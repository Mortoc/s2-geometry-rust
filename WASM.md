# WebAssembly (WASM) Support

S2 Geometry Rust provides full support for WebAssembly compilation, enabling spherical geometry operations in web browsers and other WASM environments.

## Quick Start

### Compilation for WASM

```bash
# Compile library for WASM
cargo build --target wasm32-unknown-unknown --features wasm --no-default-features

# Compile example for WASM
cargo build --target wasm32-unknown-unknown --features wasm --no-default-features --example wasm_example
```

### Features Configuration

The library uses conditional compilation to ensure WASM compatibility:

- **`wasm` feature**: Enables WASM-specific dependencies (like `getrandom/js`)
- **`--no-default-features`**: Excludes non-WASM-compatible features
- **Excluded features for WASM**:
  - `proptest-support`: Property-based testing (uses incompatible rand version)
  - `bevy-support`: Bevy math integration (ahash/getrandom conflict)
  - `bench-support`: Benchmarking with criterion (uses rayon)

## Available Functionality in WASM

All core S2 geometry operations work in WASM:

### ✅ Fully Supported
- **S2Point**: Points on the unit sphere
- **S1Angle**: Angle measurements and conversions
- **Coordinate transformations**: UV ↔ XYZ, face projections
- **Geometric predicates**: Exact orientation tests, robust sign computation
- **Mathematical utilities**: Exact arithmetic with BigRational fallback
- **Error handling**: Complete S2Error enum and S2Result types

### ❌ Not Available in WASM
- Property-based testing (proptest)
- Bevy math integration (bevy_math)
- Benchmarking (criterion)

## Example Usage

```rust
use s2geometry_rust::{S1Angle, S2Point, math::DVec3};
use s2geometry_rust::math::predicates::exact_sign;

// Create points on the sphere
let north_pole = S2Point::new(0.0, 0.0, 1.0).unwrap();
let equator_point = S2Point::new(1.0, 0.0, 0.0).unwrap();

// Calculate angle between points (π/2 radians = 90°)
let angle = S1Angle::from_points(north_pole, equator_point);
println!("Angle: {:.2}°", angle.degrees());

// Exact geometric predicates
let a = DVec3::new(1.0, 0.0, 0.0);
let b = DVec3::new(0.0, 1.0, 0.0);
let c = DVec3::new(0.0, 0.0, 1.0);
let orientation = exact_sign(a, b, c); // Returns 1 (positive orientation)
```

## Integration with Web Applications

### With wasm-bindgen

```rust
use wasm_bindgen::prelude::*;
use s2geometry_rust::{S1Angle, S2Point};

#[wasm_bindgen]
pub fn calculate_spherical_distance(
    lat1: f64, lng1: f64,
    lat2: f64, lng2: f64
) -> f64 {
    // Convert to S2Points and calculate angle
    // (Implementation would include coordinate conversion)
    let angle = S1Angle::from_radians(0.0); // Placeholder
    angle.degrees()
}
```

### Performance Notes

- **Exact arithmetic**: Uses BigRational for deterministic results
- **Memory usage**: No heap allocations for core point/angle operations
- **Size optimization**: Use `wee_alloc` allocator for smaller WASM binaries

```rust
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;
```

## Browser Compatibility

The WASM build is compatible with:
- ✅ Modern browsers (Chrome 57+, Firefox 52+, Safari 11+)
- ✅ Node.js with WASM support
- ✅ Deno runtime
- ✅ Web Workers and Service Workers

## Size Optimization

For production use, compile with optimizations:

```bash
cargo build --target wasm32-unknown-unknown --features wasm --no-default-features --release
```

Additional size reduction techniques:
- Use `wasm-opt` for further optimization
- Enable LTO in `Cargo.toml` (already configured)
- Consider `wasm-snip` to remove unused code

## Limitations

1. **Random number generation**: Limited to browser APIs via `getrandom/js`
2. **Threading**: Single-threaded execution model
3. **File I/O**: No direct file system access
4. **Debugging**: Use browser developer tools

## Testing

Run tests for WASM compatibility:

```bash
# Test compilation without WASM-incompatible features
cargo check --target wasm32-unknown-unknown --features wasm --no-default-features

# Run tests with WASM feature enabled (for compatibility verification)
cargo test --features wasm --no-default-features
```

## Future Enhancements

- WASM-specific optimizations
- Integration with web mapping libraries
- Streaming geometry processing
- WebGL visualization examples