# KOTH 0.6.0 LFQ update

Use release v0.6.0 (7cc80fb84625), whose only algorithmic change is narrower LFQ
extraction defaults: RT half-window 0.005 of run span and IM half-window 0.015
absolute 1/K0. Detection, alignment, grouping and output schemas are unchanged.
Both manuscript datasets are Orbitrap data, so only the RT change is active.

Reuse the 32 archived 0.5.0 feature/hill tables and all 96 separate searches.
Rebuild all eight LFQ bundles and run both joint arms across all three seeds
(48 searches). Keep the original UNO executable, database, RT model, cohorts,
grouping settings and score configuration. Audit observation-coordinate parity,
all reused output hashes, and all accession counts at 1%, 3% and 5% empirical
entrapment FDP. Repeat the four serialized warm-cache timing controls.

Shorten the main-text LFQ paragraph and supporting prose while retaining the
Astral and IonStar comparisons, seed ranges, quantitative limitations and
reproducibility settings. Keep complete evidence in generated SI tables and
the versioned result archive. Do not mix the separate IM experiment into this
rerun or modify the KOTH manuscript. Preserve paper-2026-09-13 as an older freeze.
