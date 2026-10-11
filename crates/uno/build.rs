use std::env;
use std::path::PathBuf;

/// `uno-ms-lightgbm-sys`'s own build.rs sets up the link search path + link
/// library for the vendored `lib_lightgbm.so` (that part propagates fine via
/// `cargo:rustc-link-search`), but a raw `-rpath` linker arg from a
/// non-binary dependency does NOT propagate to a downstream binary crate --
/// only a build script belonging to the crate that actually produces the
/// binary can emit `rustc-link-arg` for it. So the final rpath has to be set
/// here, in the `uno` bin crate's own build script.
///
/// The directory comes from the sys crate's `links = "lightgbm"` metadata, so
/// it is right both in this workspace and when built from crates.io (where the
/// two crates are unpacked in separate registry directories).
fn main() {
    let vendor_dir = env::var_os("DEP_LIGHTGBM_VENDOR_DIR")
        .map(PathBuf::from)
        .expect("uno-ms-lightgbm-sys did not report its vendor directory");
    let vendor_dir = std::fs::canonicalize(&vendor_dir).unwrap_or(vendor_dir);
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", vendor_dir.display());
}
