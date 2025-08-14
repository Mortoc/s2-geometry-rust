{
  description = "S2 Geometry Rust Port Development Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, rust-overlay, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        
        rustVersion = pkgs.rust-bin.stable.latest.default.override {
          extensions = [ "rust-src" "rust-analyzer" ];
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Rust toolchain
            rustVersion
            cargo-watch
            cargo-edit
            cargo-expand
            cargo-tarpaulin  # Code coverage
            
            # C++ toolchain for reference tests
            cmake
            gcc
            gdb
            abseil-cpp
            
            # Development tools
            git
            just  # Command runner
            
            # Libraries needed for S2 geometry
            openssl
            pkg-config
            
            # Testing and benchmarking
            valgrind
            hyperfine
          ];

          shellHook = ''
            echo "🦀 S2 Geometry Rust Development Environment"
            echo "Rust version: $(rustc --version)"
            echo "Cargo version: $(cargo --version)"
            echo ""
            echo "Available commands:"
            echo "  cargo test          - Run all tests"
            echo "  cargo test --release - Run tests in release mode"
            echo "  cargo bench         - Run benchmarks"
            echo "  cargo tarpaulin     - Generate test coverage"
            echo "  just --list         - Show available just commands"
            echo ""
          '';

          # Environment variables
          RUST_BACKTRACE = "1";
          RUST_LOG = "debug";
        };
      });
}