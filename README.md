# S2 Geometry Rust

An experimental Rust port of the Google S2 geometry library, created with Claude Code.

## Status

This project provides a Rust implementation of core S2 geometry functionality including:

- S2Point: 3D points on the unit sphere
- S1Angle & S1ChordAngle: Spherical angle representations  
- S2LatLng: Latitude/longitude coordinates
- S2Cap: Spherical caps with precision-aware intersection testing
- S2Cell & S2CellId: Hierarchical cell decomposition system

## Test Coverage

Currently achieving 93% test coverage for S2Cap implementation (13/14 tests passing) with systematic C++/Rust compatibility validation.

## Building

```bash
cargo build
cargo test
```

## C++ Compatibility Tests (Optional)

This project includes optional C++ compatibility tests that verify identical behavior between the Rust and original C++ S2 implementations. These tests are automatically skipped if the C++ submodule is not available.

**Default behavior**: Rust-only tests run, C++ tests are skipped.

**To enable C++ compatibility tests**:
```bash
git submodule update --init --recursive
cd s2geometry-cpp && mkdir build && cd build && cmake .. && make
cd ../../cpp_compat_tests && cargo test
```

## Continuous Integration

This project uses GitHub Actions for automated testing:

- **Rust Tests**: Tests the core Rust library across stable/beta/nightly toolchains without C++ dependencies
- **C++ Verification**: Builds the C++ submodule and runs compatibility tests to verify identical behavior

Both workflows automatically handle the conditional testing system - the C++ verification workflow initializes the submodule, while the Rust-only workflow demonstrates graceful degradation when C++ is unavailable.

## Architecture

The implementation follows a three-tier mathematical precision approach (fast/stable/exact paths) while maintaining memory safety and C++ compatibility.

🤖 Generated with [Claude Code](https://claude.ai/code)