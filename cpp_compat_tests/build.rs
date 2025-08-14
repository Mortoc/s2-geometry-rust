use std::path::{Path, PathBuf};

fn main() {
    // Check if the s2geometry-cpp submodule exists and is built
    let submodule_path = Path::new("../s2geometry-cpp");
    let build_path = submodule_path.join("build");
    let lib_path = build_path.join("libs2.so");
    
    let cpp_available = submodule_path.exists() 
        && submodule_path.join("src/s2").exists()
        && lib_path.exists();
    
    if cpp_available {
        println!("cargo:rustc-cfg=feature=\"cpp-compat\"");
        println!("cargo:warning=Building with C++ compatibility tests enabled");
        
        // Build the simple C++ bridge
        cxx_build::bridge("src/lib.rs")
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
    } else {
        println!("cargo:warning=s2geometry-cpp submodule not available - C++ compatibility tests disabled");
        println!("cargo:warning=To enable: git submodule update --init && cd s2geometry-cpp && mkdir build && cd build && cmake .. && make");
    }
    
    // Always rerun if submodule status changes
    println!("cargo:rerun-if-changed=../s2geometry-cpp");
    println!("cargo:rerun-if-changed=../s2geometry-cpp/build/libs2.so");
}