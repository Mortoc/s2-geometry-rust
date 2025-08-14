# S2 Geometry Rust Development Commands

# Default recipe - show available commands
default:
    @just --list

# Run all tests
test:
    cargo test

# Run tests in release mode for performance
test-release:
    cargo test --release

# Run a specific test
test-one TEST:
    cargo test {{TEST}}

# Run tests with output
test-verbose:
    cargo test -- --nocapture

# Check code formatting
fmt:
    cargo fmt --all

# Check code with clippy
lint:
    cargo clippy --all-targets --all-features

# Fix common clippy issues
fix:
    cargo clippy --all-targets --all-features --fix --allow-dirty

# Build the project
build:
    cargo build

# Build in release mode
build-release:
    cargo build --release

# Generate documentation
docs:
    cargo doc --open

# Generate test coverage report
coverage:
    cargo tarpaulin --out html --output-dir target/coverage

# Run benchmarks
bench:
    cargo bench

# Clean build artifacts
clean:
    cargo clean

# Port a specific C++ test file to Rust
port-test CPP_TEST:
    echo "Porting {{CPP_TEST}} to Rust..."
    # This would contain the logic to analyze and port a C++ test

# Run all checks (fmt, lint, test)
check: fmt lint test

# Setup development environment
setup:
    cargo install cargo-tarpaulin cargo-watch cargo-edit cargo-expand