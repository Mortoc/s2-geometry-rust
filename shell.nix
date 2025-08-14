{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    # Rust toolchain (will use rust-toolchain.toml if present)
    rustc
    cargo
    rustfmt
    clippy
    
    # Development tools
    rust-analyzer
    
    # For C++ reference comparison (optional)
    cmake
    gcc
    abseil-cpp
    openssl
    pkg-config
    
    # Documentation tools
    mdbook
  ];
  
  # Environment variables for optimal performance
  RUSTFLAGS = "-C target-cpu=native";
  
  shellHook = ''
    echo "S2 Geometry Rust Development Environment"
    echo "Rust version: $(rustc --version)"
    echo "Cargo version: $(cargo --version)"
    echo ""
    echo "Ready for geometric computations!"
  '';
}