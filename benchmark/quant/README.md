# MS1 quant prototype (Steps 1-4)

Dev/validation tools for uno's multi-sample quant, validated on the IonStar
benchmark (PXD003881: human constant, E. coli spiked A=3.3% .. E=90%).

- `make_consensus.py` — STAND-IN for koth_align: builds a koth-schema
  `intensity_matrix.tsv` from per-run koth `features.parquet` by mass+RT+charge
  clustering. Detected-only (NO MBR) because the raw mzML/hills were deleted.
  Replace with real `koth_align` output (same schema) once raw data is available.
- `de_proto.py` — DirectMS1Quant-style protein DE test (co-direction k +
  decoy-calibrated FDR). Prototype for the eventual Rust `--out-quant`.

## Validation status (A vs E, 3 reps each)
- Step 1 (`--matrix` search): WORKS. min_runs>=3 -> 1011 proteins @1%
  (860 human + 151 ecoli), decoy FDR 0.9%.
- Step 3 (DE): APPROACH VALIDATED IN PRINCIPLE — confident ecoli proteins get
  mean co-direction k=3.77 vs human 1.72 (frac k>=3: 0.67 vs 0.19).
- KNOWN GAPS before production:
  1. `--out-peptides` emits the best-MASS feature (no RT gate) -> mass-coincident
     human contamination flattens many ecoli peptides. Needs the RT-gated feature
     uno actually used (thread `pep_best` out of score.rs).
  2. No MBR (stand-in consensus) -> on/off DE peptides (present in E, absent in A)
     are missing-value-dominated. Real koth_align MBR fixes this.
  3. Normalize to the constant species / median-of-non-changing, not equal total
     (IonStar: equal-total normalization makes human appear to drift).
