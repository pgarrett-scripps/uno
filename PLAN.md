# uno — MS1-only search engine (Rust)

A Rust MS1-only protein identification engine, in the spirit of DirectMS1 /
ms1searchpy. `uno` does **not** read raw data or find features — it consumes a
**single-run biosaur2 feature file** (`*.features.tsv`), the same input
ms1searchpy uses, assigns peptide/protein identities, and controls FDR.
Single-run input is required: it carries the per-feature quality columns that
separate target from decoy, which cross-run alignment strips out.

It is, in effect, a Rust port of the [blitzms1](../blitzms1) (Python) search
backend, with digestion from [sage-core](../sage) and RT prediction distilled
from DeepLC (per [sherlock](../sherlock)).

## Decisions (locked)

- **RT model:** `gbdt` crate (pure Rust), distilled from DeepLC labels, Klammer
  composition features, sherlock's monotone piecewise-linear calibrator.
- **Scope (v1):** ID + target/decoy FDR + quant rollup. No differential stats yet.
- **Scorer:** GBDT discriminant (ms1searchpy/LightGBM-style), **not** a QDA port.
- **Digestion:** depend on `sage-core` (path dep), use enzyme/peptide/mass/fasta only.

## Data flow

```
FASTA ─► sage-core digest ─► mass-sorted peptide DB (+reversed decoys)
                                     │  binary search ±ppm, ±isotope-err
single-run features ──► Feature ─────┘─► Peptide-Feature Matches (PFMs)
 (massCalib, charge, rtApex, im,                 │
  combined_score, cosine/isotope_score, ...)     ▼
                              per-PFM features ─► distilled RT model
                                                 ─► bootstrap calibration
                                                 ─► GBDT discriminant (CV'd)
                                                 ─► target/decoy FDR (PSM/peptide)
                                                 ─► protein inference + protein FDR
                                                 ─► roll up koth intensities
                                                 ─► peptides.tsv / proteins.tsv / quant matrices
```

## Crates

- `uno-core` — library: input parsing, DB build, matching, scoring, FDR, rollup.
- `uno` — CLI over uno-core (TOML config, like koth/sage).
- `uno-rt-train` — offline binary: DeepLC label TSV ─► train `gbdt` RT model ─► serialize.

Dependencies: `sage-core` (path), `gbdt` (RT model + scorer), `rayon`, `serde`,
`csv`, `toml`, `clap`. Reuse sherlock's calibration + Klammer feature code.

## Input contract (the seam — lock first)

Read a single-run biosaur2 `*.features.tsv` (columns located by name):

```
Feature { mass: f64, charge: u8, rt: f64, im: Option<f64>,
          combined_score: f64, cosine_score/isotope_score/ppm_error/... : Option<f64>,
          intensities: Vec<f64> }   // single-run intensity (1 element)
```

RT is already in koth's reference-run space → single calibration, not per-file.

## Phases

1. **Input contract + skeleton.** Parse the three koth TSVs into `Feature`. Lock schema.
2. **Peptide DB.** sage-core digest (enzyme, missed cleavages, len bounds, fixed/var
   mods, reversed decoys) → mass-sorted `Vec<Peptide>`. TOML config mirrors sage's enzyme block.
3. **Matching.** Binary-search PFMs within ppm; enumerate isotope-error offsets
   (±1, ±2 Da/charge). Emit PFM records + raw features.
4. **RT distillation (parallel track).** sherlock `deeplc_distill.py` on a large diverse
   digest → label TSV → train `gbdt` in `uno-rt-train` → serialize. **Measure student-vs-DeepLC residual.**
5. **Scoring.** Per-PFM features: mass_ppm, rt_diff (post-calibration), im/CCS diff,
   n_isotopes, charge, koth combined_score, log-intensity, n_contributing_runs.
   GBDT discriminant (target=1/decoy=0). Bootstrap: loose first pass → calibrate RT →
   rescore. **3-fold CV by peptide to prevent decoy overfit (see below).** Target/decoy FDR.
6. **Protein inference + protein FDR.** Greedy peptide→protein assignment; protein-level
   aggregate score; protein-level target/decoy FDR.
   (Quant rollup was removed — uno is identification-only. Single-run input,
   protein output.)

## Critical notes

- **GBDT scorer must be cross-validated.** A flexible model trained to separate targets
  from decoys on the same data it scores will overfit and inflate IDs. Use 3-fold CV split
  **by peptide sequence** (train on 2 folds, score the held-out fold), as Percolator /
  mokapot / ms1searchpy do. Reuse sherlock's group-k-fold pattern.
- **MS1-only is information-poor** (mass + RT). The `im` column (timsTOF CCS) is a major
  extra discriminator — predict CCS and feed it in where available. Without IM/FAIMS,
  expect more modest depth than ms1searchpy's FAIMS headline numbers.
- **Distilled RT ≤ DeepLC, and composition-GBDT loses positional info.** Acceptable because
  RT is a soft score, not a hard filter — but validate the residual (phase 4). Fallback:
  burn CNN (sherlock `rt_net`) if residual too high.
- **Don't re-implement quant.** koth owns alignment + XIC + LFQ. uno = identity + rollup.
- **Calibration bootstrap is fragile** with few confident IDs — loose-then-tighten first pass.
- **sage-core is beta/MS2-oriented.** Depending pulls in unused fragment-index code; fine to
  ignore. If it bites, vendor enzyme/peptide/mass/fasta.

## Findings from phase-5 validation (human+ecoli benchmark, no ion mobility)

Built and verified phases 1–4 end-to-end. Phase 5 surfaced the real science of MS1-only:

- **RT distillation works.** gbdt student R² 0.93 vs DeepLC; **0.94 vs observed RT** on
  5000 Sage-confident peptides. koth rtApex is a valid RT axis (0.75 vs Sage rt, many exact).
- **Calibration is the hard bootstrap.** The tightest-mass match per feature is only ~5-6%
  correct at this DB density — mass barely constrains identity for low-mass charge-2 features.
  Binned-median / windowed-mode / Theil-Sen all fail at >70% contamination. **RANSAC works**
  (calibrate.rs::fit_ransac): isolates the ~1900-point true-RT diagonal, median residual ~1.1 min.
- **Composition gbdt cannot separate a peptide from its reverse** (median Δpred ~1.8 min < RT
  noise). So mass-preserving (reversed) decoys are indistinguishable twins. Fix applied:
  **mass-shifted decoys** (database.rs `decoy_mass_shift`, default 11 Da, like koth's LFQ TDC) —
  decoys match only by random coincidence, RT regains full discriminating power.
- **PSM-level FDR has an inherent floor.** Net target-excess over decoys ≈ 0.5–1k features on a
  ~24k random-coincidence floor. Random peptides matching in BOTH mass (tight) and RT (close)
  are indistinguishable from real IDs using mass+RT alone. This is the known MS1-only limit;
  DirectMS1 overcomes it with **protein-level aggregation** (+ FAIMS/ion mobility, absent here).

### Implication for the design
- Build **phase 6 protein-level aggregation** (binomial / multi-peptide score) — this is where
  the ~1k real-peptide signal concentrates and decoys (singletons, random proteins) fall away.
  Evaluate FDR at protein level via the human+ecoli two-proteome entrapment.
- The user's real **rapidms1 dataset has ion mobility** — predicting CCS and adding an IM-match
  term is a large orthogonal gain the IM-free benchmark cannot show. Strongly worth using.
- Consider tighter precursor tolerance only with leaner DB; on its own it cut signal faster than noise.

### Current dev invocation
```
uno --features .../run.features.tsv \
    --fasta .../human_ecoli.fasta --rt-model models/rt_*.urtc --config <toml> --out-proteins out.tsv
```
RT model: `uno-rt-train --labels rt_deeplc.tsv --out model.gbdt` (distill from DeepLC labels).
`uno-rt-train --predict-seqs <file>` validates a model against known peptides.

## Reference reading before phase 5

`blitzms1/src/blitzms1/searcher.py` — copy the exact PFM/protein feature lists even though
we swap QDA → GBDT.
