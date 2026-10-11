use std::env;
use std::path::PathBuf;

/// Links against the vendored LightGBM shared library (extracted from the
/// official `lightgbm` PyPI wheel, see vendor/README.md) and bakes an absolute
/// rpath into the built binary so it finds `lib_lightgbm.so` at runtime without
/// requiring `LD_LIBRARY_PATH`. No cmake / C++ build needed — this crate never
/// compiles LightGBM itself, only calls its stable C API against a prebuilt .so.
fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor_dir = manifest_dir.join("vendor");
    println!("cargo:rustc-link-search=native={}", vendor_dir.display());
    // file is `lib_lightgbm.so` -> link name `_lightgbm` (linker prepends "lib").
    println!("cargo:rustc-link-lib=dylib=_lightgbm");
    // NOTE: a raw `-rpath` linker arg can't be emitted from here -- Cargo only
    // applies `rustc-link-arg` to binaries/cdylibs/tests OWNED BY THIS PACKAGE,
    // and this package only has a lib target (this crate's own `cargo test`
    // picks up the rpath fine; a downstream binary like `uno` would not). The
    // actual `uno` binary's build.rs sets its own rpath to this vendor dir.
    // This crate's own test binaries need the rpath too (cargo only puts
    // native search paths inside target/ on the loader path).
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", vendor_dir.display());
    // Exposed to direct dependents as DEP_LIGHTGBM_VENDOR_DIR.
    println!("cargo:vendor_dir={}", vendor_dir.display());
    println!("cargo:rerun-if-changed=vendor/lib_lightgbm.so");
}
