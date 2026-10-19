# Changelog

## Unreleased — effective-configuration report

- Every run now prints its complete effective configuration to stderr, between
  `# uno effective configuration` and `# end of effective configuration`. No
  output file changes, so archived results stay byte-identical.
- The report covers what the config format cannot: the tuned defaults no config
  field exposes (`[tuned]`: 300 boosting rounds, depth 7, learning rate 0.05,
  sequoia `min_child_weight` 80, the 14 tightest of 40 cumulative protein
  inclusion bands, deviation cut max 8), the per-stage seeds as resolved, the
  Bruker `.d` feature-quality gate when it is applied, and any `UNO_*`
  environment override in effect.
- Each line carries its provenance: `(tuned)`, `(derived)`, `(env)`,
  `(overridden)` for a clamped field, `(inactive)` for a field the config format
  accepts and no code path reads. Setting an inactive or clamped field from a
  config file now also emits a warning naming the value the engine uses instead.
- `example_config.toml` marks its own inactive and clamped fields, and points at
  the stderr report as the record to keep with a result.

## Unreleased — KOTH 0.7.0 compatibility

- Require the sibling KOTH dependency at version 0.7.0 and record the released source revision. Pin additive isotope-evidence defaults in both raw profiles and update the standalone feature configuration.
- Accept long-bundle schemas 2/3 alongside legacy schema 1, validate group confidence and exclusive-signal provenance, and retain suppressed-cell zeros. New group admission replaces the removed seed floor; original anchor-quality checks remain.
- Add the current LFQ configuration using `max_group_qvalue = 0.05`. The approved publication rerun is 23 extractions, two alignments and 75 searches; current manuscript results remain the archived 0.6.0 evidence until rerun.

## Unreleased — experimental RT-first IM sequence CNN

- Added charge-aware IM training and fine-tuning using the shared RT sequence CNN, with Rust inspection, prediction, and independent evaluation.
- The bundled IM CNN supplies a default soft classifier feature (`--no-im` disables it; `--im-model` overrides it) after RT/mass shortlisting; there is no IM cutoff. Incompatible legacy feature schemas are rejected.
- IM search precomputes the target/decoy library at observed IM charges and caches raw predictions across runs and calibration attempts, keyed by model, sequence, and charge. `UNO_IM_CACHE_DIR` sets storage; `UNO_IM_CACHE=0` disables caching. Shortlist calibration remains run-specific.
- IM caches now use Zstd-compressed Parquet chunks in reusable, per-user OS temporary storage. `UNO_IM_CACHE_DIR` can select persistent storage. Atomic publication keeps incomplete writes invisible; old TSV caches remain untouched.
- The [27-search trial](https://github.com/pgarrett-scripps/uno-paper/blob/main/docs/im-cnn-results-20260914.md) improved six of nine CNN/off pairs, lost two, and tied one. Median depth increased in all three acquisitions, but the predefined consistency rule was not met; the default promotion retains this limitation.

## Unreleased — top_k throttle fix + ms1searchpy-comparable peptide metric

Head-to-head on ms1searchpy's exact single-run input (koth HeLa `.features.tsv`,
same target FASTA, matched digest: mc=0, len 7–30, no variable mods, 10 ppm,
shuffled decoys) showed uno at protein parity (149 vs 155) but reporting far
fewer "peptides" — a metric artifact plus a candidate-depth throttle.

- **`top_k` default 25 → 50.** The per-feature candidate cap was silently
  dropping real peptides before RT could rescue them. Human+yeast entrapment on
  the single-run HeLa data: 25→105 IDs, 50→208, 100→275, all at ~1% entrapment
  FDP. Raising it ≈2× the honest IDs at negligible cost. This also overturns the
  previously-recorded "~218 single-run ceiling" — that ceiling was the throttle,
  not the data.
- **Matched-peptides metric** (`--out-matched-peptides`, always logged): distinct
  target peptides mapping to a confident (q≤1%) protein with ANY feature inside
  3σ mass + 3σ RT — ms1searchpy's "Matched peptides for top-scored proteins".
  Lets the two tools be compared on the same (lenient, non-peptide-FDR) number
  instead of uno's stricter 1%-peptide-FDR count. On the matched input uno now
  reports 268 proteins / 2631 matched peptides vs ms1searchpy 155 / 1700.
- RT model (pure-Rust CNN vs feeding ms1searchpy's own DeepLC predictions) makes
  no material difference to the result; the throttle was the lever.

## Unreleased — DirectMS1 single-run + mass-anchor rewrite (major)

uno was ~20x behind ms1searchpy on single runs and weak on aligned matrices (62
proteins on koth human+E.coli). Root cause: decoys were mass-twins (interior
permuted in place → same mass), inflating the random-match rate 4x and breaking
both the bootstrap and FDR. Reworked the scoring around ms1searchpy's algorithm,
now the DEFAULT:

- **Per-protein-pool shuffled decoys** (`decoy_shuffle`, default on): different
  mass per decoy → honest random-match null. Deterministic (sorted dealing).
- **Mass-anchor identifier** (`mass_anchor`, default on): iterative size-adjusted
  binomial bootstrap (mass stage at >=3 isotopes), then a **min-combined-deviation
  any-feature** pass — each peptide scored by the feature minimizing
  `mass_z² + (rt_diff/rt_sig)²`, a deviation-cutoff scan, and a target/decoy FDR.
- **Two-stage RT with histogram-peak (mode) calibration**: recalibrate RT on the
  deviation-confident peptides using a per-bin MODE fit (`Calibration::fit_robust`).
  The confident proteins' peptide sets are real-RT-diagonal + random-RT background;
  the mode tracks the diagonal where a median PWL is pulled off it (~5 min → toward
  DeepLC's true ~1-2 min), recovering each protein's lower-confidence peptides.
- **PAVA isotonic monotonicity** in `fit_robust`: the monotone constraint on the mode
  knots was a running-max upper envelope, which biases the curve UP at every noisy
  dip. Replaced with pool-adjacent-violators isotonic regression (least-squares
  monotone fit). Single-run 806 → 847 honest, aligned 469 → 604.
- **Greedy parsimony + target/decoy pair deletion** in the mass-anchor protein
  scoring (ms1searchpy `final_iteration`): each peptide is assigned to one protein
  and removed from competitors before rescoring, and the higher-scoring of each
  target/decoy accession pair is kept. Counting every shared/homologous peptide for
  every protein inflated spurious proteins and made the FDR optimistic. Tightens the
  entrapment-validated true FDR (3.6% → 3.0% at q≤0.01).
- A CV PFM rescorer (ms1searchpy `-ml 1`) is implemented (`ml_rescore`, OFF by
  default), using the fast closed-form CV-LDA (the `gbdt` crate is minutes/run and
  separates no better). On the current koth-derived input the mass+RT-gated candidate
  set is ~50/50 target/decoy and decoys (shuffled) share the targets' composition, so
  no classifier beats the min-deviation score at honest FDR — the lever is a purer
  gate (tighter RT), not the classifier. Adds amino-acid composition + `id_count`
  (the features ms1searchpy's ML uses) for when the gate is purified.
- `top_k` default 5 → 25.

**FDR honesty note:** the entrapment must use an RT library covering the entrapment
(yeast) peptides; otherwise they lack RT, are dropped from the deviation stage, and
"0 yeast" is guaranteed rather than earned. With full coverage, uno's target/decoy
q≤0.01 corresponds to ~3% true FDR; at a true 1% FDR uno identifies ~600 single-run
proteins vs ms1searchpy's 1510 (which is genuinely honest, 0 yeast). The remaining
gap is upstream feature richness + decoy coincidence, not the protein scoring.

Result (koth human+E.coli, 1% FDR, **entrapment-validated honest**, deterministic):
single-run **66 → 847** (vs ms1searchpy 1510); aligned **62 → 604**. The mass stage
reaches ms1searchpy parity (226 vs 204). Remaining ~1.8x to SOTA is the RT-recovery
multiplier. (An on-data MLP-readout fine-tune was tried in the mass-anchor stage and
did NOT help — the robust calibration already captures the run-level RT shift and a
129-d readout refit overfits the confident subset; left out.)

## Unreleased

### Added
- ~~**On-data RT fine-tuning** (`scoring.ondata_rt_finetune`, on by default). Freezes
  the distilled MLP's penultimate embedding and refits only its linear readout on
  the run's confident hits (ridge, closed-form, deterministic). Keeps the pretrained
  representation — unlike the from-scratch composition ridge, which degrades depth —
  and is FDR-honest under two-proteome entrapment (0 entrapment false positives).
  +51% protein depth on the koth human+E. coli benchmark (39 → 59). Falls back to
  the frozen RT path automatically when the RT model exposes no embedding (gbdt).~~
  **REMOVED in `a91c6a9` (discriminant-path purge).** The implementation no longer
  exists; `scoring.ondata_rt*` are inert config keys pending removal from the
  struct. See `crates/uno-core/src/score.rs` ("was tried here … not a validated win
  anywhere"). The +51% figure was measured on the earlier development configuration
  and is retired to the manuscript's Supporting Information.

### Performance
- Protein inference reworked: parallel threshold-optimization grid (rayon),
  memoized binomial tail, and best-first greedy parsimony with lazy invalidation.
  Protein phase ~40 s → ~1 s single-threaded (~40× algorithmic); full search
  43 s → 3.4 s (~13×). **Output byte-identical** to the prior implementation.
- Release profile: thin LTO + single codegen unit.

### Docs
- Paper ported to Typst (`paper/paper.typ`, arkheion template); adds RT-fine-tune
  and performance results.

## 0.1.0 (unreleased) — initial MS1-only search engine

A complete DirectMS1-style MS1-only protein identification pipeline in Rust,
consuming a single-run biosaur2 feature file (the same input ms1searchpy uses).

### Pipeline (crates/uno-core)
- **input** — reads a single-run biosaur2 `*.features.tsv` into a `FeatureMatrix`.
- **database** — in-silico digest via `sage-core` (enzyme, missed cleavages, mods, decoys);
  f64 mass recompute; mass-sorted binary-search candidate lookup. Mass-shifted decoys
  (`decoy_mass_shift`, default 11 Da) so decoys are not mass twins.
- **matching** — ppm mass matching with isotope-error offsets → Peptide-Feature Matches.
- **rt** — Klammer + positional-hydrophobicity features; gbdt RT model distilled from DeepLC.
- **calibrate** — monotone PWL calibration; **RANSAC** fit for heavily-contaminated bootstrap.
- **score** — top-K reduction, semi-supervised (mokapot-style) CV gbdt discriminant split by
  feature, RANSAC-seeded bootstrap RT calibration, target/decoy q-values.
- **protein** — binomial (Poisson-tail) protein scoring over searchable peptides; protein FDR.
- **quant** — rollup of koth per-sample intensities to peptide and protein matrices.

### Tools
- **uno** (CLI) — runs the full pipeline; TOML config; TSV outputs.
- **uno-rt-train** — distills a DeepLC label TSV into a gbdt RT model; `--predict-seqs` for validation.

### Validation highlights
- RT model: R² 0.93 vs DeepLC; **Pearson 0.94 vs observed RT** on 5,000 Sage-confirmed peptides.
- RANSAC calibration: median residual ~1.1 min (from contaminated best-mass seeds).
- Protein FDR controlled (0 decoys at the top); the most abundant proteins identified.

### Ion mobility (CCS) — added
- **im** module: gbdt 1/K0 model distilled from **im2deep** (charge-aware features).
  338k-label library generated by im2deep; distilled student R² **0.992** vs im2deep.
- `uno-rt-train --im-labels` trains it; `uno --im-model` enables IM scoring.
- Scoring: per-PFM 1/K0 prediction, **per-charge** RANSAC calibration (im2deep→instrument),
  `im_diff` discriminant feature, IM-purified positives, and an orthogonal RT∩IM
  protein "seen" filter. `matching.min_charge` skips singly-charged contaminants.
- On the rapidms1 timsTOF data IM measurably cut the random hit rate (e.g. 4.2%→1.9%,
  up to ~6×) — the orthogonality works.

### Known limitation
PSM- and protein-level depth is constrained by MS1's inherent random mass+RT coincidence
floor on this no-ion-mobility, two-proteome benchmark — the hardest case. (Was documented in
`PLAN.md`, deleted 2026-07-31 as superseded; the manuscript's Limitations section carries
this now.) The largest expected gains are ion mobility / CCS (timsTOF), FAIMS, and a leaner
search space; the user's `rapidms1` timsTOF dataset (with IM) is the natural next target.
