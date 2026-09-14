# UNO publication freeze — 13 September 2026

The `paper-2026-09-13` tag consolidates the eight manuscript/scoring commits
formerly on `scaffold-migration`, the validated KOTH 0.5.0 long-matrix search,
and the complete current manuscript and reproduction scripts. The paper's
primary and matrix studies retain their original pinned versions and frozen
source/binary archives; this tag identifies the consolidated submission tree,
not a claim that every historical study ran the same executable.

The dormant IM gate and the new IM trainer/model-loader extensions were removed
from this tree and preserved for `codex/im-timstof`. The existing LFQ worktree
retains exploratory scripts, seed-pooling results, and its original files.
Both dirty checkouts were backed up with checksums before integration in
`experiments/consolidation-20260913`.

Validation passed: 73 Rust workspace tests, eight matrix benchmark tests, and
`just verify` for manuscript formatting, extractors, prose, provenance,
freshness and PDF overflow. The freshly rebuilt executable reproduces all
protein, LFQ and abundance tables byte for byte in four seed-zero controls:
an Astral 100 ng individual search, its joint matrix search, its joint search
with cell covariates, and the IonStar 20-run joint search. No study thresholds
or manuscript statistics changed during consolidation.

`benchmark/paper_freeze_20260913.json` records source checksums, the tested
binary digest, and the four parity receipts. Full logs and the parity driver
are under `benchmark/results/paper_freeze_20260913`. The test environment needs
`LD_LIBRARY_PATH=crates/uno-lightgbm-sys/vendor` for the standalone LightGBM
library test; the matrix Python tests require the packages in
`benchmark/matrix_study_requirements.txt`, including PyArrow.

The separate IM experiment is development work and does not update this tag
or the manuscript. No UNO GitHub release or package publication is implied.
