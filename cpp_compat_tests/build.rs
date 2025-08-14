use std::path::PathBuf;

fn main() {
    // Build the simple C++ bridge
    cxx_build::bridge("src/simple_lib.rs")
        .file("src/simple_bridge.cpp")
        .std("c++17")
        .include("src") // Include the local src directory for simple_bridge.h
        .include("../s2geometry-cpp/src")
        .flag_if_supported("-w") // Suppress warnings from C++ code
        .compile("s2-simple-bridge");

    // Link against the S2 C++ library
    let s2_lib_path = PathBuf::from("../s2geometry-cpp/build");
    println!("cargo:rustc-link-search=native={}", s2_lib_path.display());
    println!("cargo:rustc-link-lib=dylib=s2");
    
    // System libraries that S2 depends on
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=m");

    // Tell cargo to invalidate the built crate whenever the C++ code changes
    println!("cargo:rerun-if-changed=src/simple_bridge.cpp");
    println!("cargo:rerun-if-changed=src/simple_bridge.h");
}