use std::env;
use std::path::PathBuf;

/// `uno-lightgbm-sys`'s own build.rs sets up the link search path + link
/// library for the vendored `lib_lightgbm.so` (that part propagates fine via
/// `cargo:rustc-link-search`), but a raw `-rpath` linker arg from a
/// non-binary dependency does NOT propagate to a downstream binary crate --
/// only a build script belonging to the crate that actually produces the
/// binary can emit `rustc-link-arg` for it. So the final rpath has to be set
/// here, in the `uno` bin crate's own build script.
fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor_dir = manifest_dir.join("../uno-lightgbm-sys/vendor");
    let vendor_dir = std::fs::canonicalize(&vendor_dir).unwrap_or(vendor_dir);
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", vendor_dir.display());
}
