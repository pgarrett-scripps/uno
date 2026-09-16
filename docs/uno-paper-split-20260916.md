# Manuscript and study separation — 2026-09-16

The identification manuscript, its study drivers, figures, tables, analysis,
and local benchmark evidence live in `../uno-paper/`. LFQ and cohort workflow
studies live in `../uno-lfq-paper/`, including their drivers and local evidence.
Both repositories use paper-scaffold for shared compilation and formatting.

Uno retains reusable software, runtime models, and software tests. Test-only
configuration fixtures are under `tests/fixtures/`; optional large regression
inputs are under ignored `tests/data/`. The `paper` compatibility link points
to the identification manuscript and does not own its files. Migration maps
and original source snapshots are recorded in each paper’s `provenance/`.
