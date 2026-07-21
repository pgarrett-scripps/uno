# uno vs. ms1searchpy — pipeline comparison

Stage-by-stage comparison of **uno** (this repo, pure-Rust MS1 protein search) against
**ms1searchpy** (`/home/patrick-garrett/Repos/ms1searchpy`, the Python DirectMS1 engine).

**Everything below is verified against source code**, cited as `file:line`. uno citations are
under `crates/`; ms1searchpy citations are under `ms1searchpy/ms1searchpy/`. Where a claim is
*not* from source (e.g. benchmark depth), it is marked **[empirical, not source]**.

> Provenance note: this doc reflects uno after the discriminant-path purge (commit `a91c6a9`) —
> uno is now **mass-anchor → proteins only**; the CV-gbdt/LDA PFM discriminant, `--out-peptides`,
> `--hills`, quant, and the `mass_anchor=off` fallback have been deleted.

---

## TL;DR — the essential differences

| Dimension | uno | ms1searchpy |
|---|---|---|
| Language | Rust | Python |
| Input | **single-run** biosaur2 feature file (`*.features.tsv`) | **single-run** biosaur2 feature file (`*.features.tsv`) — same input |
| Candidate set | **top-50** mass matches per feature, per target/decoy class | **all** matches inside the ppm window (no cap) |
| Decoys | **ON by default**, global-pool interior shuffle, ~1:1 | **OFF by default** (`-ad 0`); shuffle when enabled |
| RT model | learned **CNN** (`.urtc`), or DeepLC library file | **DeepLC** (always run `-deeplc 1`) |
| PFM ML rescorer | sequoia-boost (Rust XGBoost), ON | LightGBM, ON (`-ml 1`) |
| Protein score | binomial survival `-log10 P(X≥k)` | binomial survival `-log10 P(X≥k)` (same idea) |
| Determinism | byte-deterministic | non-deterministic (random search, `random.shuffle`) **[empirical, not source]** |

The two engines are the **same core algorithm** on the **same single-run input** (mass+RT match →
per-PFM ML filter → protein binomial with target/decoy FDR). They differ mainly in **candidate
breadth** (capped vs uncapped) and **ML gating shape**.

---

## Stage 1 — Input

Both engines consume the **same input: a single-run biosaur2 feature file** (`*.features.tsv`).
Neither reads raw data or finds features. This is deliberate — a single-run feature file carries
the per-feature quality columns (`cosine_score, isotope_score, ppm_error, neutron_offset,
isoerror, nScans, …`) that separate real peptides from random mass coincidences (and shuffled
decoys) in the PFM ML step. **Cross-run aligned/consensus matrices drop those columns and are
not a valid input** — uno rejects them and reads single runs only.

**uno** (`input.rs`) locates columns **by name**, order-independent (`input.rs:1-14`). Required:
`massCalib, mz, charge, rtApex`; optional `im, combined_score` (defaults 1.0), biosaur2 aliases
`nIsotopes`↔`n_isotopes` / `nScans`↔`n_scans`, quality columns `cosine_score, isotope_score,
ppm_error, neutron_offset, isoerror, isoerror2, rt_width`, and `intensityApex`/`intensitySum`
(`input.rs:160-210`). `Feature` carries neutral `mass, mz, charge, rt`, optional `im`,
`combined_score`, the quality columns, and a single `intensity` (`input.rs:20-62`).

**ms1searchpy** (`utils.py`) reads a single biosaur2 feature `.tsv` per run (or an mzML run
through biosaur2 inline, `utils.py:79-112`). Required: `nIsotopes, nScans, charge, massCalib,
rtApex, mz` (`utils.py:116-123`); optional quality columns `cosine_score, isotope_score,
combined_score, ppm_error, neutron_offset, isoerror, isoerror2` feed the ML if present
(`main.py:49-52, 577-585`). It loops over input files, processing each independently
(`main.py:511, 520-522`).

> **Same input, same evidence unit** — one run's biosaur2 features, with the quality columns
> intact. This is a prerequisite for a meaningful uno-vs-ms1 comparison: same run, same
> single-run features.

---

## Stage 2 — Digestion & peptide DB

| | uno (`database.rs`) | ms1searchpy (`search.py`) |
|---|---|---|
| Enzyme | `KR`, restrict `P`, C-terminal (`:198,200,293`) | `[RK]|{P}` trypsin (`:30`) |
| Missed cleavages | 0 (`:192`) | 0 (`-mc`, `:32`) |
| Length | 7–50 (`:193-194`) | 7–30 (`-lmin/-lmax`, `:28-29`) |
| Mass range | 500–5000 Da (`:199-200`) | (charge/tol driven) |
| Charges | matching `min_charge=1` (`matching.rs:32`) | 1–4 (`-cmin/-cmax`, `:33-34`) |
| Fixed mods | Carbamidomethyl@C `+57.02146` (`:185-186,201`) | Carbamidomethyl@C (`-fmods`, `:35`) |
| Variable mods | none; MS1 can't localize (`:188-191`) | none by default |

uno delegates digestion to `sage-core` (`database.rs:10, 285-313`); ms1searchpy uses
`pyteomics.parser.icleave`. Nearly identical digest, except **uno allows peptides up to length 50
vs ms1's 30**.

---

## Stage 3 — Decoy generation

Both use the **same method** — peptide interior-residue shuffle keeping terminal residues — but
differ in default and scope:

- **uno**: **ON by default** (`decoy_shuffle=true`, `database.rs:207-208`). Interior residues
  pooled **per protein**, deterministic Fisher-Yates (xorshift seeded by accession,
  `database.rs:138-146`), dealt back keeping first+last residue (`database.rs:152-158`).
  Composition-preserving but **not** mass-preserving — each decoy has a different mass, sampling
  the true random mass-coincidence background (`database.rs:59-65`). ~1:1 target:decoy
  (`database.rs:168-177`). Also offers reverse/Markov/DeBruijn protein-level generators
  (`decoy.rs:28-63`).
- **ms1searchpy**: **OFF by default** (`-ad 0`, `search.py:37`). When enabled, per-protein
  peptide interior `random.shuffle`, terminals fixed (`utils.py:246-263`), prefix `DECOY_`,
  writes a `_shuffled` DB (`utils.py:200, 292`).

> **Consequence:** run out-of-the-box, **ms1searchpy makes zero decoys** and reports FDR by a
> different route; you must pass `-ad 1` for target/decoy competition. uno always builds decoys.
> The shuffle recipe itself is essentially identical between the two.

---

## Stage 4 — Mass matching & candidate reduction

Both match on **neutral mass** (not m/z), same 10 ppm default tolerance:

- **uno** (`matching.rs`): `ppm_tolerance=10.0` (`:30`), isotope offsets `[-1,0,1]`
  (`:31`, koth-mono may be ±1 neutron off), `min_charge=1` (`:32`). Then **reduces to top-`k`=50
  per feature, split by target vs decoy class** (each bucket stable-sorted by `|mass_error_ppm|`,
  truncated, `matching.rs:150-156`; `top_k=50` at `score.rs:97`). A PFM = `(feature, peptide,
  offset)` with signed `mass_error_ppm` (`matching.rs:45-53`).
- **ms1searchpy** (`main.py`): symmetric ppm window from `-ptol=10` (`:591-592`), binary-search
  the mass-sorted features and **keep every match in `[m−dm, m+dm]`** — no top-k cap
  (`:615-624`). Also gates features by `-i` (min isotopes, default **2**, `search.py:23`) and
  `-sc` (min scans, default 2).

> **This is the deepest architectural difference.** ms1 keeps *all* mass matches (then leans on
> its downstream ML + isotope/scan gates to control the background). uno caps at 50/feature to
> bound memory/compute, then relies on the mass-anchor calibration loop. uno's set is broad
> (~17M PFMs on a 211k-feature run **[empirical, not source]**) because it needs breadth for its
> feedback-loop calibration.

---

## Stage 5 — Mass recalibration

Both estimate a systematic mass shift + width from a confident seed and re-express mass as a
z-score:

- **uno**: inside the mass-anchor bootstrap — 5 rounds, each recomputing `shift = median(errors)`
  and `sigma = clipped_sigma(...)` on confident-protein peptides, tolerance tightening from 10 ppm
  to `(3σ).min(10)` (`score.rs:341-411`). Final mass gate `(3σ).min(10 ppm)` (`score.rs:480`).
- **ms1searchpy**: stage-0 protein scoring → top-protein peptides → fit a Gaussian to the
  mass-diff histogram (`calibrate_mass`, `main.py:483-493`) giving `mass_shift_cor`,
  `mass_sigma_cor`; keep matches with `|md−shift|/sigma ≤ 3.0` (`main.py:981-984`). ML feature is
  the z-score `(md−shift)/sigma` (`main.py:1458`).

Same concept (recalibrate to ~0 ppm, gate at 3σ). uno folds it into an iterative 5-round
bootstrap; ms1 does a single-shot Gaussian fit on top-protein peptides.

---

## Stage 6 — RT prediction & calibration

**uno RT** (`rt.rs`, `rt_cnn.rs`): a learned **CNN** loaded from a `.urtc` file — AA embedding →
residual dilated Conv1d stack → masked mean+max pool → optional DeepLC-style **global
physicochemical features** (11-d in the v3 set) → MLP head → single scalar iRT
(`rt_cnn.rs:1-4, 242-361`). Dispatch by magic bytes (`URTC`→CNN, `UMLP`→MLP, else gbdt-JSON,
`rt.rs:151-166`). The CNN source header says it is **trained on experimental iRT
(ProteomeTools/Prosit)** (`rt_cnn.rs:1-4, 360`) — note: "DeepLC-distilled" is stated only for the
*gbdt/MLP* students (`rt.rs:1-2`, `rt_mlp.rs:5`), and "Chronologer" appears **nowhere in source**
(only in model filenames). uno can also take an **external DeepLC library** via `--rt-library`
(`peptide<TAB>rt`), which replaces the CNN as abstract RT (`main.rs:42-46, 154-170`). On-data
calibration (abstract→observed) is done downstream by the mass anchor via RANSAC + monotone PWL
(`calibrate.rs`, invoked in `score.rs`).

**ms1searchpy RT** (`main.py`): **DeepLC** — ms1searchpy is always run with `-deeplc 1` in this
project (a good RT predictor is the entire point of MS1-only search; the source's non-DeepLC
fallback is not a mode we use, benchmark, or consider). DeepLC is enabled at `:669-672` and
predicts in the second calibration stage. Two-stage calibration (`-ts 2`): fit `XRT_shift/XRT_sigma`
from `pred−obs` via a `noisygaus` Gaussian, reselect within ±3σ, refit DeepLC on the clean set
(`:1149-1157, 1319-1373`). Final RT gate: `(obs−pred−shift)²/σ² ≤ 9` i.e. ±3σ (`:1424-1430`).

> **RT models:** both engines use a learned sequence RT predictor — uno a self-contained CNN
> (no Python dependency), ms1searchpy **DeepLC**. Both calibrate on-data and gate at ~±3σ. All
> uno-vs-ms1 comparisons use ms1 with `-deeplc 1`.

---

## Stage 7 — PFM ML rescoring

Both train a **cross-validated gradient-boosted classifier** on per-PFM features, semi-supervised
(positives = peptides of confident proteins, negatives/decoys held out), fold-split by peptide to
avoid leakage, then use the score to gate which peptides pass to protein scoring:

| | uno (`score.rs`) | ms1searchpy (`main.py`) |
|---|---|---|
| Library | sequoia-boost (Rust XGBoost); LightGBM via `UNO_ML_LIGHTGBM` (`:1606-1622`) | LightGBM (`:28`) |
| Default | ON (`ml_rescore=true`, `:125`) | ON (`-ml 1`, `search.py:38`) |
| Positive labels | peptides of proteins at **q≤0.25** (`:963-967`) | peptides of **top-25%-FDR** targets (`:1899-1903`) |
| CV | 3-fold, peptide-disjoint (`:1097,1121`) | 3-fold via group `G` (`:1938-1953`) |
| Features | PFM (mass_z, rt, charge, isotopes, scans…) + AA composition + dipeptides + biosaur2 isoerror/rt_width (`:140-179, 1052-1094`) | same family: `mass_diff, rt_diff, Scans, Isotopes, mz, rt, ch, im`, koth-extra, AA counts, `c_DP/c_KP/c_RP` (`:1611-1668`) |
| HP search | fixed tuned params | **random search, 25 evals** (`:1928, 1553-1585`) |
| Gating | **ensemble**: rank by ML score, sweep cumulative inclusion bands, average each protein's binomial over bands 1–14 of 40 (`UNO_ENS_*`, `:1233-1303`) | bin preds into 50 quantiles, keep bins while decoy-fraction < `ml_correction`, re-bin to 10, pass `qpreds≤10` (`:2007-2051`) |

The classifiers and features are near-equivalent. **uno's distinctive piece is the ensemble
band-averaging** (average each protein's binomial score across ~top-35% inclusion thresholds
instead of one hard cutoff) — this is uno's main depth lever **[empirical, not source]**.
ms1searchpy uses a single decoy-calibrated quantile cutoff. uno's ML swaps in only if it strictly
beats the deviation-only result (`score.rs:1330-1333`).

---

## Stage 8 — Protein scoring & FDR

**Same statistical core:** protein score = binomial survival `−log10 P(X≥k | n, p)`, where `k` =
distinct matched peptides, `n` = theoretical peptides of the protein, `p` = random per-peptide hit
rate estimated from decoys.

- **uno** (`protein.rs`): exact binomial via regularized incomplete beta (`binom_sf_ge_k`,
  `:202-210`). `p` **dynamically re-estimated** as the decoy pool shrinks during greedy selection
  (`:889`). **Feature-level** competition (`score_proteins_feature_compete`): a winning protein
  claims all peaks under its peptides; a peptide dies only when its last peak is banned
  (`:940-996`), using lazy per-size max-heaps (`:854-881`). Observed peptides = a `FxHashSet` of
  distinct peptide indices (`:797-798`). q = running `(d+1)/t` (pyteomics correction=1),
  monotonized (`:1022-1036`).
- **ms1searchpy** (`main.py`, `utils.py`): `−log10 binom.sf(v−1, n, p)` (`utils.py:428-432`). `p`
  from decoys, recomputed dynamically in `worker` (`:231, 2266`). Greedy parsimony picks the
  best protein, assigns its peptides, decrements shared counts (`:285-301`). Threshold grid scan:
  `final_iteration` runs `worker` over **10 mass/RT-koef levels** and averages per-protein scores
  (`:357-402`). FDR during greedy: `decoy/target`, break at `≥12.5·fdr` (`:323-327`); final
  filter via `pyteomics.aux.filter(..., correction=1)` (`:424`).

Both do greedy parsimony + binomial + target/decoy protein FDR with a `+1` correction. uno
competes at the **feature/peak** level (one peak → one protein); ms1 competes at the peptide
level. ms1 averages the binomial over a 10-level mass/RT threshold grid; uno averages over its
ML-inclusion ensemble bands.

---

## Stage 9 — Output & headline metric

- **uno**: `--out-proteins` TSV — `accession, decoy, observed_peptides, theoretical_peptides,
  coverage, binomial_score, mean_rt_diff, q_value` (`main.rs`). Proteins at q≤1% logged.
- **ms1searchpy**: `_proteins.tsv` (FDR-filtered), plus `_PFMs.tsv`, `_protsN.tsv`,
  `_proteins_full.tsv` (`main.py:439-442`). Headline log metric **"Matched peptides for
  top-scored proteins"** = distinct-sequence PFMs on FDR-passing proteins (`:462-466`) — a
  lenient, non-peptide-FDR count. (uno had a matching `--out-matched-peptides` for like-for-like
  comparison; it was removed in the purge.)

---

## Where they genuinely diverge (summary)

1. **Candidate breadth** — uno caps at **top-50/feature**; ms1 keeps **all** in-window matches.
2. **ML gating shape** — uno averages the protein binomial over an **ensemble of inclusion
   bands**; ms1 uses a **single decoy-calibrated quantile cutoff** + a 10-level threshold grid on
   the protein side.
3. **Decoys** — uno **on** by default; ms1 **off** by default (`-ad 0`).
4. **Determinism** — uno is byte-deterministic; ms1 uses `random.shuffle` + 25-eval random HP
   search **[empirical, not source]**.

## Where they are essentially the same

- **Same single-run biosaur2 feature-file input**, with the per-feature quality columns intact.
- Neutral-mass matching at 10 ppm with isotope-offset handling.
- Learned sequence RT predictor (uno CNN / ms1 DeepLC), on-data calibrated, ~±3σ RT gate.
- Mass recalibration to ~0 ppm, 3σ mass gate.
- Interior-shuffle decoys keeping terminal residues.
- Semi-supervised 3-fold CV gradient-boosted PFM classifier, positives = confident-protein
  peptides.
- Protein score = binomial survival with decoy-estimated `p`, greedy parsimony, target/decoy
  FDR with `+1` correction.

**Bottom line:** uno is a faithful, pure-Rust, deterministic reimplementation of the DirectMS1 /
ms1searchpy algorithm on the **same single-run biosaur2 input**, with a self-contained CNN RT
model replacing DeepLC and an ensemble-of-bands ML gate replacing the single quantile cutoff. The
statistical core (binomial protein inference over mass+RT-matched peptides) is the same.
