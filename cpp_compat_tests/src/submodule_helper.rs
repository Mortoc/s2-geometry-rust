//! Helper functions to detect if the s2geometry-cpp submodule is available

use std::path::Path;

/// Check if the s2geometry-cpp submodule is available and built
pub fn is_cpp_submodule_available() -> bool {
    let submodule_path = Path::new("../s2geometry-cpp");
    let source_path = submodule_path.join("src/s2");
    let build_path = submodule_path.join("build");
    let lib_path = build_path.join("libs2.so");
    
    // Check if submodule exists and contains source files
    if !submodule_path.exists() || !source_path.exists() {
        return false;
    }
    
    // Check if the library has been built
    if !build_path.exists() || !lib_path.exists() {
        eprintln!("Warning: s2geometry-cpp submodule exists but library not built.");
        eprintln!("Run: cd s2geometry-cpp && mkdir build && cd build && cmake .. && make");
        return false;
    }
    
    true
}

/// Print helpful messages about the C++ submodule status
pub fn print_submodule_status() {
    if is_cpp_submodule_available() {
        println!("✓ s2geometry-cpp submodule available - C++ compatibility tests enabled");
    } else {
        println!("⚠ s2geometry-cpp submodule not available - C++ compatibility tests skipped");
        println!("  To enable C++ tests:");
        println!("  1. git submodule update --init --recursive");
        println!("  2. cd s2geometry-cpp && mkdir -p build && cd build");
        println!("  3. cmake .. && make");
    }
}

/// Macro to conditionally run tests only when C++ submodule is available
#[macro_export]
macro_rules! cpp_test {
    ($test_name:ident, $test_body:block) => {
        #[test]
        fn $test_name() {
            if !$crate::submodule_helper::is_cpp_submodule_available() {
                println!("Skipping {} - C++ submodule not available", stringify!($test_name));
                return;
            }
            $test_body
        }
    };
}

/// Conditional compilation helper for modules that depend on C++ FFI
#[macro_export]
macro_rules! cpp_module {
    ($($item:item)*) => {
        $(
            #[cfg(feature = "cpp-compat")]
            $item
        )*
    };
}