# S2 Geometry C++ Compatibility Test Suite

This directory contains a **conditionally enabled** test suite that validates functional equivalence between the Rust S2 implementation and the original Google C++ S2 library when available as a git submodule.

## Conditional Testing System

The C++ compatibility tests are **automatically enabled/disabled** based on submodule availability:

- **✅ Submodule Available**: Full C++ compatibility tests run, comparing Rust and C++ implementations
- **⚠️ Submodule Missing**: Tests are gracefully skipped with informative messages
- **🔧 Build Failed**: Clear instructions provided for building the C++ library

This design allows the main Rust library to be used independently without requiring C++ dependencies.

## Quick Start

### Option 1: Use Rust-Only (Default)
No setup required! The library works independently and C++ tests are automatically skipped.

### Option 2: Enable C++ Compatibility Tests

```bash
# 1. Initialize the submodule
git submodule update --init --recursive

# 2. Build the C++ library  
cd s2geometry-cpp
mkdir -p build && cd build
cmake .. && make

# 3. Run tests
cd ../../cpp_compat_tests
cargo test
```

## Test Architecture

The compatibility test suite:
- **Runtime Detection**: Automatically detects if C++ submodule is available
- **Conditional Compilation**: FFI bridge only compiled when C++ library is present
- **Graceful Fallback**: Clear messaging when tests are skipped
- **Status Reporting**: Shows exactly which tests are running/skipped

## Test Files

- `tests/cpp_compat.rs` - S2Point operations (normalize, angle, cross product)
- `tests/s2latlng_compat.rs` - Latitude/longitude conversions and distances
- `tests/s2cell_compat.rs` - Cell area, perimeter, vertices, and hierarchy
- `tests/s2cap_compat.rs` - Spherical cap operations and cell covering
- `src/lib.rs` - FFI bridge definitions and type conversions
- `src/cpp_bridge.h` - C++ header declarations for bridge functions
- `src/cpp_bridge.cpp` - C++ implementation calling Google S2 library

## Setup Requirements

### Dependencies

The C++ compatibility tests require:
1. **Google S2 C++ Library** - Built and available in `../s2geometry-cpp/build/`
2. **Abseil C++ Library** - Required by Google S2
3. **C++ Build Tools** - cmake, g++, make

### Building Dependencies

#### Option 1: Using System Package Manager

On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install cmake build-essential libabsl-dev
```

On Fedora/CentOS:
```bash
sudo dnf install cmake gcc-c++ make abseil-cpp-devel
```

#### Option 2: Building from Source

1. Build Abseil:
```bash
git clone https://github.com/abseil/abseil-cpp.git
cd abseil-cpp
mkdir build && cd build
cmake .. -DCMAKE_CXX_STANDARD=17 -DBUILD_SHARED_LIBS=ON
make -j$(nproc)
sudo make install
```

2. Build Google S2:
```bash
cd ../s2geometry-cpp
mkdir build && cd build  
cmake .. -DCMAKE_CXX_STANDARD=17
make -j$(nproc)
```

#### Option 3: NixOS Development Shell

Update `../shell.nix` to include Abseil:
```nix
{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    # Existing dependencies...
    rustc cargo rustfmt clippy rust-analyzer
    cmake gcc
    
    # Add Abseil for C++ S2 compatibility tests
    abseil-cpp
    
    # Documentation tools
    mdbook
  ];
  
  shellHook = ''
    echo "S2 Geometry Rust Development Environment"
    echo "Ready for geometric computations!"
  '';
}
```

Then build S2:
```bash
nix-shell --run "cd s2geometry-cpp && mkdir -p build && cd build && cmake .. && make -j$(nproc)"
```

## Running Tests

Once dependencies are built:

```bash
# Run all compatibility tests
cargo test

# Run specific test suites
cargo test s2cap_compat
cargo test s2cell_compat  
cargo test s2latlng_compat
cargo test cpp_compat
```

## Test Design

### Equivalence Testing Methodology

Each test compares Rust and C++ implementations by:
1. Creating identical input data
2. Calling both Rust and C++ functions
3. Comparing results within numerical tolerance (1e-15)
4. Testing edge cases and boundary conditions

### FFI Bridge Architecture

The FFI bridge uses the `cxx` crate for safe Rust-C++ interop:
- Shared struct definitions for data exchange
- Conversion traits between Rust and C++ types  
- Error-safe function call interfaces
- Memory-safe data transfer

### Test Coverage

- **S2Point**: Normalization, angle calculations, cross products
- **S2LatLng**: Point conversions, distance calculations, round trips
- **S2CellId**: Level calculations, parent hierarchy, point conversions
- **S2Cell**: Area calculations, perimeter, vertex extraction
- **S2Cap**: Point containment, area calculations, cell covering

## Current Status

✅ Test framework implemented  
✅ All test files created  
✅ FFI bridge designed  
⚠️  Requires external dependency setup (Abseil)  
⚠️  C++ library build needed  

The compatibility test suite is ready to run once the required C++ dependencies are installed and built.

## Troubleshooting

**Build Error: "absl/base/attributes.h: No such file or directory"**
- Install Abseil C++ library (see Setup Requirements above)

**Link Error: "cannot find -ls2"**  
- Build the Google S2 C++ library in `../s2geometry-cpp/build/`

**Test Failures with Large Numerical Differences**
- Check compiler flags and optimization levels match between Rust and C++
- Verify both implementations use the same mathematical constants

## Integration with Main Test Suite

These compatibility tests complement the main Rust test suite by:
- Validating mathematical correctness against the reference implementation
- Catching algorithmic differences that unit tests might miss
- Ensuring consistent behavior across different platforms
- Providing confidence in the Rust port's accuracy