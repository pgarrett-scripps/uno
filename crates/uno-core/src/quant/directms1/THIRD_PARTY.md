# Attribution

The DirectMS1Quant algorithm and reference source are by Mark Ivanov and
contributors: https://github.com/markmipt/ms1searchpy,
commit 9df4a941a901aac657b1c9ab5229cb49671e3c7d.

Copyright 2021 Mark Ivanov. Licensed under Apache-2.0; see
`licenses/Apache-2.0.txt`. The Python file under `tools/quant-reference/reference/`
is unmodified and retains its upstream license/provenance. This Rust crate is
an adapted reimplementation, not an endorsement by or release from those authors.

Citation: Ivanov et al., DirectMS1Quant: Ultrafast Quantitative Proteomics with
MS/MS-Free Mass Spectrometry. https://doi.org/10.1021/acs.analchem.2c02255.

`src/order.rs` adapts NumPy's scalar introsort ordering, originally by
Charles R. Harris and maintained by NumPy developers (NumPy 1.26.4,
`numpy/core/src/npysort/quicksort.cpp`). Its BSD-3-Clause notice is preserved in
`licenses/NUMPY-BSD.txt`. Other Rust dependencies retain their respective licenses.

The Rust port uses independent Rust numerical implementations and a CPython-
compatible MT19937 stream. Statistical behavior and explicitly documented
unsupported cases are described in README.md; outputs remain experimental.
