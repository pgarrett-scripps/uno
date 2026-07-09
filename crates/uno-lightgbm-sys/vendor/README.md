`lib_lightgbm.so` extracted from the official `lightgbm` PyPI wheel (v4.6.0,
manylinux x86_64): `pip download lightgbm==4.6.0 --no-deps` then unzip the wheel
and take `lightgbm/lib/lib_lightgbm.so`. MIT-licensed (Microsoft/LightGBM).

Vendored rather than built from source because this environment has no cmake,
and the C API this crate binds against (`LGBM_Dataset*` / `LGBM_Booster*`) has
been stable for years — linking a prebuilt library is simpler and faster than
vendoring + building the C++ sources. If you need a different platform/arch,
replace this file with the matching `lib_lightgbm.so` from the same wheel family.
