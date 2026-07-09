# Porting ms1searchpy's PFM ML rescorer into uno

## Why (the gap)
On run10 entrapment, at matched TRUE FDR: uno (binomial/deviation, honest) = 660 @1% /
896 @2%; ms1searchpy (with `-ml 1`) = **1827 @1% / 1881 @2%** — entrapment-honest (1 yeast
in top 1600). uno's binomial *beats* ms1searchpy's binomial stage (660 vs 220); the entire
~2.3-2.8x final gap is ms1searchpy's LightGBM PFM rescorer, which uno lacks an effective
version of. uno's prior ML (LDA, CV-GBDT) either inflated FDR or did nothing (gave 38
confident proteins, discarded). See memory uno-ml-protein-scoring.

## What ms1searchpy actually does (source: ~/Repos/ms1searchpy/ms1searchpy/main.py)
A semi-supervised, Percolator/mokapot-style PFM classifier:

1. **Label scheme (THE honesty mechanism), main.py:1876-1878**
   - `decoy2` = the TRUE target/decoy flag (peptide maps only to decoy proteins). Reserved
     ONLY for FDR estimation; NEVER a training label.
   - `decoy` (the ML LABEL) is REDEFINED = "peptide is NOT from a top-25%-FDR confident
     target protein". So training positives = peptides of confident target proteins;
     negatives = dim target peptides.
   - **Training set EXCLUDES true decoys** (`train_df[~train_df['decoy2']]`, line 1954).
   - After predicting, the TRUE decoy flag is restored (`df1['decoy']=df1['decoy2']`,
     line 1967) and FDR is estimated on those held-out decoys (lines 1965-1990).
   => Decoys are a pristine null, never seen in training. This is exactly the leakage uno's
   LDA had (decoys were training negatives AND set the FDR -> optimistic).

2. **Peptide-disjoint 3-fold CV, lines 1917-1954**
   - Fold key `G` is assigned per PEPTIDE (`seq_gmap`), and train excludes any peptide id
     present in the test fold. uno currently folds BY FEATURE (`fold_of(feature_idx)`) —
     different and a likely leak.

3. **Model**: LightGBM binary classifier, ~100-500 boost rounds, early stopping on AUC.

4. **Rich feature set (~40+), get_features_pfms = all columns minus banned. Built 1600-1730:**
   - Per-PFM: `mass_diff`, `rt_diff`, `ch`(charge), `im`, `Isotopes`(n), `Scans`(n),
     `mzraw`, `plen`, intensity.
   - **Composition: 20 AA counts `c_<aa>` + `c_DP`/`c_KP`/`c_RP`** (proline dipeptides).
   - **Per-peptide-group consistency (the big lever uno lacks):**
     - `id_count` = # features matching this peptide; `seq_count`, `charge_count` (distinct
       charges), `im_count` (distinct IM) — multi-evidence support.
     - within-group normalized deviations: `rt_diff_abs_pdiff` (= |rt_diff| − group median),
       `rt_diff_abs_pnorm` (= |rt_diff| / group sum); same for `mass_diff`.
   - `mass diff q median` (RT-binned mass calibration residual).

5. **PFM score -> peptide -> protein, lines 1965-1990 + final_iteration (line 140)**
   - Bin `preds` into 50 quantiles, walk bins, find the qpred cutoff where the (decoy2)
     decoy ratio is acceptable; keep peptides below it, re-bin to 10, map seq->qpred.
   - `final_iteration` scores proteins by a binomial/hypergeometric over their kept
     peptides + target/decoy protein FDR with pair deletion.

## What uno already has (reuse)
- gbdt crate (GBM) — `crates/uno-core/src/score.rs` `make_config`/`featurize` (~14 feats).
- Greedy-parsimony protein scorer + target/decoy pair deletion + binomial:
  `crates/uno-core/src/protein.rs` `score_proteins_greedy` — this is ms1searchpy's
  `final_iteration` equivalent; KEEP it, just feed it better per-peptide calls.
- Confident-protein set from the mass-anchor/deviation bootstrap (the `conf` HashSet in
  score.rs ~line 586) — exactly the "top-25% targets" positive labels.
- Full per-peptide PFM lists (`pep_pfms`/`pep_prots`) for group aggregates.

## What to BUILD (the port, in order)
1. **Group-aware featurizer** (extend `featurize`): add 20 AA counts + DP/KP/RP; and
   per-peptide-group aggregates (id_count, distinct-charge/IM counts, within-group
   |rt_diff|/|mass| pdiff & pnorm). Needs a pass grouping PFMs by `peptide_idx` first.
2. **Semi-supervised labels**: positive = PFM whose peptide ∈ confident proteins (`conf`);
   negative = target PFM not in conf. EXCLUDE decoy PFMs from the training matrix.
3. **Peptide-disjoint CV**: fold by `peptide_idx` (swap `fold_of(feature_idx)` ->
   `fold_of(peptide_idx)`), and exclude test-fold peptides from train (already implied by
   peptide-disjoint folds; add an id-exclusion guard for shared peptides).
4. **Predict all PFMs (incl decoys)**, take per-peptide best (or sum) pred.
5. **FDR on held-out decoys**: q-threshold the per-peptide preds using the TRUE decoy
   flag, then feed kept peptides to `score_proteins_greedy` (unchanged).
6. **Gate**: keep ML result only if entrapment-honest in tests; default behind a config
   flag until validated.

## Validation (NON-NEGOTIABLE — uno's ML inflated before)
Use the yeast-covered entrapment harness every iteration:
`crates/uno-rt-train/scripts/validate_rt_in_search.sh` pattern (rt_cnn_chrono_ch32 +
/tmp/entrap_db.fasta). Success = real-protein depth rises toward ~1800 AND TRUE entrap
FDR stays ≈ nominal (yeast-scaled). If yeast leak in -> the leakage isn't fixed; do NOT ship.

## IMPLEMENTATION RESULT (2026-06-28) — port done, but blocked on INPUT FEATURES
Ported the recipe into `score.rs` (mass-anchor ML block): semi-supervised labels
(conf-target=1 / dim-target=0 / decoy=None held out for FDR), peptide-disjoint CV
(`fold_of(peptide_idx)`), nonlinear gbdt (`cv_train_predict`), + composition/DP-KP-RP/
group-consistency features. It is FDR-honest and now actually separates conf from the
rest (was ~0 with the old target-vs-decoy LDA), BUT it does NOT beat uno's deviation
(370-464 confident proteins vs 672) so it stays discarded.

ROOT CAUSE (diagnosed via UNO_DUMP_ML feature dump): on uno's input, **dim-real peptides
and decoys are statistically identical on every feature** — means: rt_abs 12.19 vs 12.66,
mass 0.61 vs 0.63, ln(n_scans) 2.51 vs 2.49, ln(n_iso) 1.195 vs 1.193, intensity 6.78 vs
6.77. The only separation is conf-vs-rest = mass+RT tightness, which deviation already uses
optimally. So DECOY_mean(0.200) ≈ dim_mean(0.210); no ML signal beyond mass+RT exists.

WHY: uno's `--intensity-matrix` "shape" format lacks the discriminative biosaur2 columns
ms1searchpy's ML relies on. ms1searchpy input (run10_ent_features.tsv) has: `isoerror`,
`isoerror2` (isotope-pattern fit error — the key real-vs-random signal), `rtStart`/`rtEnd`
(peak width/shape), `mono_hills_scan_lists`/`mono_hills_intensity_list` (per-scan isotope
profiles), separate `intensityApex`/`intensitySum`, `scanApex`, FAIMS. uno's matrix has
only massCalib, mz, charge, rtApex, im, combined_score(=1 const), n_contributing_runs(=1),
n_isotopes, n_scans, intensity. Missing exactly the features that separate real from decoy.
ms1searchpy also genuinely matches ~2x more peptides (1357 vs 670 proteins with 5+ peptides
— NOT a parsimony artifact). The 2x gap is real and feature-driven.

NEXT STEP (the real fix): extend uno's input + FeatureMatrix to carry the rich biosaur2
per-feature descriptors (start with `isoerror`/`isoerror2` + peak width), thread them into
`featurize`, THEN the semi-supervised ML (already built) should separate real from decoy.
This is an input/feature-extraction change, not an ML-algorithm change. The ML code is in
place and waiting for the features.

## Risk / effort
- Medium. The gbdt crate ≈ LightGBM for this; the real work is the group-aware features
  (#1) and the label/CV/FDR plumbing (#2,#3,#5). ~1-2 focused sessions.
- Main risk: reproducing the honesty. The decoy2-holdout + semi-supervised-labels +
  peptide-disjoint-CV combination is what keeps ms1searchpy honest; uno failed before by
  training target-vs-decoy with feature-folds. Implement all three together or it will
  inflate again.
- gbdt crate may underperform LightGBM (no early-stopping/regularization parity). If depth
  stalls, consider binding LightGBM via FFI or porting its inference only.

## ENSEMBLE PORT SHIPPED (2026-07-07)

The piece described in `uno-ml-protein-scoring` memory as "the CORRECTED PLAN" — ms1searchpy's
`final_iteration` 10-step cumulative qpred-bin threshold sweep + per-protein score AVERAGING —
was validated only in a Python side-harness (`decisive.py`, never checked into this repo) and
then never actually ported to Rust. It sat as a single-cutoff peptide filter in `score.rs`
(bins 0..last_bin, ONE pass into `score_proteins_feature_compete`), which is exactly the
"collapsed" ablation the Python test showed loses most of the lever (1073 -> 632 when
qpred bins were collapsed to one cutoff).

Ported properly now (`score.rs`, the `cfg.ml_rescore` block): rank peptides by the existing
CV-gbdt qpred, sweep 10 CUMULATIVE inclusion thresholds (top 10%, 20%, ..., 100%), run
`score_proteins_feature_compete` at each, accumulate each protein's binomial score across all
10 passes (missing pass = 0, same accumulator pattern as the existing `UNO_ENSEMBLE` deviation
path), average, pair-compete target/decoy, (d+1)/t FDR. Replace the deviation result only if
`ml_n > best_n` (a safe, automatic no-regression guard — no env var needed anymore).

VALIDATED on PXD015669 (QEHF1_09744_JB), human+yeast entrapment, run through koth_align's
native biosaur2 feature file (isoerror/isoerror2/rt_width present, converted to uno's
intensity_matrix.tsv format — see /tmp/cmp/biosaur_intensity_matrix.tsv build recipe: map
biosaur2's massCalib/mz/charge/rtApex/im/nIsotopes/nScans/isoerror/isoerror2 directly,
combined_score=1.0, seed_run=0, n_contributing_runs=1, rt_width=rtEnd-rtStart,
intensity=intensityApex):

| | real proteins @1% nominal FDR | yeast (false positives) |
|---|---|---|
| deviation only (baseline) | 223 | 4 |
| **+ qpred ensemble (this port)** | **249-252** (~1-2 protein run-to-run gbdt noise) | **1** |

+12-13% depth AND better honesty (fewer yeast leaks), not a tradeoff — the ensemble
structure itself is the lever, confirmed again here: `ML diag` still shows
`dim_mean ~= DECOY_mean` (0.070 vs 0.068) — uno's pure-Rust gbdt classifier still barely
separates a dim-real peptide from a decoy per-PFM, exactly as before. All of the gain is
the 10-step sweep-and-average structure, not classifier quality.

On a koth-ALIGNED matrix (this same run, but through `koth_align` instead of the raw
biosaur2 tsv) isoerror/isoerror2 are absent (lost in cross-run alignment) and the ensemble
correctly no-ops (232 < 260 deviation baseline; the `ml_n > best_n` guard leaves the
deviation result in place). So this only helps single-run biosaur2-native searches today —
getting isoerror-equivalent shape descriptors through `koth_align`'s consensus output is the
natural follow-up if this should also help multi-run/aligned searches.

NOT a fix for the larger gap: ms1searchpy's own honest ceiling on this exact run is ~1510
(native biosaur2, full pipeline). 249-252 is still ~6x below that.

## REAL LIGHTGBM BOUND IN, CLASSIFIER-QUALITY HYPOTHESIS FALSIFIED (2026-07-07)

Wrote `crates/uno-lightgbm-sys`: hand-written FFI over LightGBM's stable C API
(`LGBM_Dataset*`/`LGBM_Booster*`), linked against the real `lib_lightgbm.so` extracted from
the `lightgbm` PyPI wheel (no cmake needed — the C API is stable, validated the exact calling
convention with a Python ctypes probe before writing any Rust). `cv_train_predict_lgbm` in
`score.rs` replaced the pure-Rust `gbdt` crate in the qpred-ensemble path (same peptide-
disjoint CV / semi-supervised-label logic, LightGBM natively takes {0,1} labels so the old
gbdt {-1,+1} degenerate-gradient workaround is gone too).

Rationale for doing this: `dim_mean ~= DECOY_mean` (0.070 vs 0.068) on uno's pure-Rust gbdt
had been the standing excuse for weak per-PFM separation across TWO prior investigations
(2026-06-28 and this session) — "our classifier can't separate them, a real LightGBM
probably could." That hypothesis is now testable, not just assumed.

**Result on the SAME native-biosaur2 PXD015669 entrapment benchmark: real LightGBM gives
`dim_mean=0.068 DECOY_mean=0.067`, 258 confident proteins (257 real/1 yeast) vs the pure-Rust
gbdt ensemble's ~250 (249 real/1 yeast).** ~3% better, not the jump toward 1073 the
"classifier quality" framing implied.

**This FALSIFIES the classifier-quality hypothesis.** Real LightGBM, on uno's actual matched
PFM feature set, separates a dim-real peptide from a decoy just as poorly as uno's own gbdt
crate did. The bottleneck was never which gradient-boosting implementation uno uses — it's
that uno's per-PFM features (mass_z, rt_z, n_isotopes, n_scans, composition, group-consistency
stats) genuinely carry little signal beyond what the mass+RT deviation score already extracts,
on uno's own candidate/matching pipeline. The Python isolation test's 1073 number came from
ms1searchpy's OWN PFMs_ML.tsv running through ms1searchpy's OWN `final_iteration` — i.e. ms1's
OWN matching + biosaur2's full native feature richness + its two-stage RT calibration,
not "LightGBM vs gbdt" as an isolated variable. **Do not propose "swap the ML model" as a
lever again** — it has now been tried twice (uno's LDA, uno's gbdt, and now real LightGBM) and
never moved the needle beyond ~3-13%. If depth still matters, the next real candidates are
upstream of the classifier: (a) uno's own matching/candidate generation quality (why is
uno's tight-mass+RT candidate set so much noisier than ms1's PFMs_ML set to begin with?), or
(b) `final_iteration`'s iterative razor (peptide reassignment using PRIOR-round protein
scores across the 10 passes — this port's passes are independent, not chained).

Both the gbdt-ensemble and lightgbm-ensemble code paths are validated to safely no-op (never
regress) on koth-ALIGNED matrices where isoerror/isoerror2 are absent (242 < 260 deviation
baseline on the same run through `koth_align`, guard correctly left the baseline in place).
