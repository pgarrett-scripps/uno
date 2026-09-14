# TODO — what is left

Scoped to **what it takes to publish in JPR**. Anything that would merely make the
repo nicer has been cut; the full 11-lane pre-submission audit, with `file:line`
evidence for every item, is preserved at [`paper/AUDIT.md`](paper/AUDIT.md).

Three documents, no overlap:

| file | holds |
|---|---|
| **this file** | the code-side work, and the critical path |
| [`paper/SUBMISSION_TODO.md`](paper/SUBMISSION_TODO.md) | the manuscript work, item by item |
| [`paper/AUDIT.md`](paper/AUDIT.md) | the evidence: what was measured, what was verified correct, and every item cut from here |

**Critical path: the code fixes and the re-run are done (2026-09-11). What is left
is the author-blocked material in §3 and the submission bundle.**

---

## DONE 2026-09-03 — koth: split hills by FAIMS compensation voltage (`uno --raw`)

Found 2026-09-03 by the Figure 6 raw-mode comparison (`benchmark/raw_mode/`).
koth now reads `MS:1001581` (and Thermo's `FAIMS CV` trailer), keys hill
detection on -50/-65/-80 V with per-channel scan indices, and runs feature
assembly and scoring within each CV. The 05433 rerun reports 99,812 features and
2,175 proteins @1%, versus biosaur2's 100,949 features and 2,072 proteins. The
previous unsplit path reported 68,183 features and 867 proteins.

## DONE 2026-09-09 — `uno --seed`: the determinism claim was not a fair comparison

uno's byte-identical output came entirely from hardcoded constants: the booster
seed (42), the RANSAC sampling LCG constant, the decoy-shuffle hash salt, and
the CV fold-assignment hash. Nothing about the algorithm carries less randomness
than ms1searchpy's. Reporting one locked uno run against three unlocked
ms1searchpy runs measured the locking, not the engines.

`uno --seed N` now perturbs all four stages together. Seed 0 is the default and
leaves every constant untouched, so published results still reproduce: a seed-0
run on `hela5min_09766` is byte-identical to the committed benchmark output.
Unit tests in `decoy.rs` and `calibrate.rs` pin seed 0 to the historical stream
and assert that other seeds actually move the decoys and the sample sequence.

What the seed sweep (`benchmark/seed_sweep.sh`, `benchmark/spread_table.py`)
shows is that uno's seed-to-seed spread is the same order as ms1searchpy's
run-to-run spread, and that the published single value is not a typical draw on
every dataset. The manuscript's determinism section, Figure S6 and the cover
letter were rewritten around that; see `paper/SUBMISSION_TODO.md`.

The reciprocal also holds and is now demonstrated rather than asserted:
ms1searchpy's decoy database is byte-identical across processes once its
module-level generator is seeded. `benchmark/ms1_seeded.py` locks that plus the
LightGBM hyperparameter search, LightGBM training, and — the one that is easy to
miss — `PYTHONHASHSEED`, because ms1searchpy repeatedly applies a *seeded*
shuffle to a `list(set(peptides))` whose order randomized string hashing already
changed.

## 1. Code fixes that gate the re-run

**DONE 2026-07-31.** All five landed, plus two more found while verifying. Two runs
on one feature file are now byte-identical; `benchmark/determinism_check.sh` is the
check, and it passes on both `09744` and the FAIMS run.

- [x] **`est_rt_sig` residual ordering** (`score.rs`). Fixed at two levels: the
      `pep_pfms` walk now iterates sorted keys, and `clipped_sigma` sorts its sample
      before reducing, which covers its other three callers too. Confirmed real —
      with the sort removed, one sample gives `1.285371343411269` vs
      `1.2853713434112695` depending on order.
- [x] **Quant peptide-vector ordering** (`quant.rs`). `resolve_feature_competition`
      and `profile_reassign` both returned a `HashMap`'s values directly; both now
      emit in a canonical order. `group_abundance` is fixed transitively, since it
      sums in the order it is given.
- [x] **`pep.decoy` guard** (`score.rs`). The mass-anchor loop admitted decoy
      peptides of confident target proteins into the anchor set and the mass
      recalibration sample. `confident` is keyed `(is_decoy, accession)` and the
      loop asked only for the target entry, but uno's decoys carry their target's
      accessions, so the check did not exclude them. The RT path had always guarded
      this; the mass path had not.
- [x] **Two-run byte-diff determinism test.** Two forms: `benchmark/determinism_check.sh`
      for the end-to-end claim as the paper states it, and unit tests in `uno-core`
      (`clipped_sigma` order-invariance, feature-competition output stability). Both
      unit tests were verified to FAIL with their fix reverted.
- [x] **`benchmark/entrapment.py` tie-break.** Now `(-score, not is_yeast)`:
      deterministic, and **conservative** — a tied block is charged its entrapment
      hits before it is credited its true ones. The audit suggested target-first,
      which is equally deterministic but optimistic at every tie; that is the wrong
      direction for a paper whose claim is that its FDR is honest.

### Found while verifying — both were live

- [x] **`members_by_leading` collapsed target and decoy groups** (`group.rs`). It
      keyed the group table by leading accession alone, but uno's decoy peptides
      carry the *same* accessions as their targets, so a target group and a decoy
      group collide and `collect()` kept whichever came later in `groups` — an order
      set by the process hash seed. This was the last source of run-to-run
      difference, and it was invisible to every unit test: two runs agreed on all
      3325 proteins, every score, peptide count and q-value, and disagreed on the
      group members of 38 rows. Now keyed `(decoy, leading)`.
- [x] **`uno` did not build at HEAD.** `koth_ff::PipelineOptions` gained an
      `emit_ms2` field upstream, so `crates/uno/src/koth_ingest.rs` failed to
      compile — `cargo check -p uno-core` passed, which is why it went unnoticed.
      Set to `false`, which koth documents as byte-for-byte identical to the
      historic MS1-only path. **The re-run could not have been started without
      this.**

### Effect on depth

As predicted, the numbers move. On `hela5min_09746`, 1% goes 862 → **876** and 5%
goes 1209 → **1194** (old values from `benchmark/results/summary.tsv`). That figure
combines the decoy guard, the determinism fixes and the new conservative tie-break,
and it is a single dataset — it is a sanity check on magnitude and direction, not a
result. The re-run produces the real numbers.

Still open in this area, not a gate:

- [ ] `cargo test -p uno-lightgbm-sys` fails on a missing `lib_lightgbm.so`.
      Pre-existing and off the default path (the ML backend is sequoia-boost;
      LightGBM is opt-in via `UNO_ML_LIGHTGBM=1`), but it makes a bare `cargo test`
      red, which trains people to ignore it.
- [ ] `BroadFit::sigma_floor` (`broad_rt.rs:27`) is never read — bounded RT recovery sets
      the applied floor to zero, so the field records an input nothing consumes. It is the
      only warning `cargo build --release` emits from this workspace. Deliberately left
      alone until the manuscript's builds are done: removing it changes the binary, and the
      published results are pinned to `c376954d`.

## 2. The re-run — DONE 2026-09-11

The gate for every identification number in the paper, run once against one pinned
binary (commit `0710225`, binary SHA-256 `c376954d`) with `/usr/bin/time -v` on
every invocation. Repeat it with `benchmark/publication_rerun.py` followed by
`benchmark/refresh_summaries.py` whenever the engine changes again.

**2026-08-31 Technical Note scope:** the paper now reports only the six primary
DirectMS1-study files and three independent Astral files. Quantitation, the
run10 quality-column arms, and the separate three-replicate FAIMS analysis were
removed from the manuscript and no longer belong in the publication re-run.

**Preconditions**

- [x] §1 above landed and the determinism test passing.
- [x] Working tree committed. **Commit `0710225`**, binary SHA-256
      `c376954dbd7c725c8cd280b9c392ba3854fd9c831eb423601a3eb7ebc5539eda`. Both are in
      `paper/si/versions.typ`, along with the three out-of-tree path dependencies
      (sage-core, sequoia-boost, koth_ff), which a commit is the only way to pin.
- [x] The ms1searchpy patch captured as a diff: `benchmark/ms1searchpy_comparator.patch`
      (SHA-256 `f8d0c455`), taken from the installed checkout at commit `9df4a94`, which
      reports itself as 2.8.12 — a version that does not exist on PyPI. Without it, stock
      ms1searchpy reads none of the five quality columns the comparison depends on.
- [x] Log the uno build metadata alongside each run. Every search now writes a provenance
      stamp (command, input hashes, **binary hash**, commit, output hash) into
      `benchmark/results/publication_rerun/`, and the `/usr/bin/time -v` capture holds the
      engine's own log. Last round's `benchmark/results/*.log` were zero-length, which is
      why the build could not be identified afterwards.

**Scope — DONE 2026-09-11.** `benchmark/publication_rerun.py`, 156 uno searches on
the pinned binary: six primary files × biosaur2 × seeds 0-10, six primary files ×
`--raw` × seeds 0-2, twelve Astral files × both feature sources × seeds 0-2.
Seed-0 searches serialized (their `/usr/bin/time -v` output is what the resource
panels report), the rest two at a time.

- [x] **149 of the 156 outputs are byte-identical to the superseded build.** The
      seven that differ had all returned nothing before; all seven now report
      proteins (200 pg rep1 biosaur2 seeds 1-2, 200 pg rep2 raw seeds 1-2, 200 pg
      rep3 biosaur2/raw seed 1, 1 ng rep4 biosaur2 seed 1). One search still
      returns nothing: 200 pg rep1, `--raw`, seed 1. All 36 seed-0 searches are
      byte-identical, so no headline depth moved; only the re-measured timings did
      (wall +2 to +13%, peak RSS within ±12%).
- [x] ms1searchpy was not re-run: the comparator did not change. Its existing
      serialized searches carry through. The re-run did surface 20 comparator
      replicates that existed on disk but had never entered `spread.tsv`
      (reps 4-10 on 09746, 09766 and the FAIMS file); they are now in the
      inventory, which takes the comparator from 19 draws to 39 and leaves the
      headline unchanged (uno's median still wins 4 of 6 at 1% and 6 of 6 at 5%).
- [x] RT-model comparison: dropped from the manuscript 2026-09-01, not re-run.
      The RT-predictor accuracy benchmark WAS re-run on the pinned build in an
      isolated worktree: every accuracy metric and every per-peptide prediction
      came back byte-identical. Throughput (Figure S4) is a timing measurement and
      is re-measured separately on a quiet host.

**After — DONE 2026-09-11**

- [x] `benchmark/refresh_summaries.py` rebuilds every derived summary from the
      re-run tables: `spread.tsv`, the seed-sweep and variability-extension
      summaries, the uno rows of `summary.tsv`, `raw_mode/summary.tsv` and the
      Astral input-series summary.
- [x] Every referenced table is now written by a script.
      `analysis/scripts/gen_perdataset_table.py` was the last hand transcription;
      it reproduces the old cells exactly, which also validates the transcription
      it replaces.
- [x] Manuscript registry updated and rebuilt: `paper/si/versions.typ` pins the
      commit, the binary hash, the three out-of-tree path dependencies and the
      comparator patch; `just verify` is clean apart from build staleness that
      clears on commit.

---

## 3. Blocked on the author

Nothing can be drafted for these.

- [ ] **Funding block.** `paper.typ` carries the grant list copied verbatim from the
      d_noise manuscript. Confirm or correct, and file it in ACS's Funder Registry Tool.
- [ ] **Corresponding author's telephone number.** ACS requires it; absent today.
- [ ] **Release decision** (D5, "release all three, not yet"): the uno repo public, koth
      archived — at minimum a versioned Zenodo deposit, since it is load-bearing for all of
      §3.9 — and a Zenodo archive carrying the benchmark data, RT distillation labels and
      trained models. Every benchmark input and result file is gitignored today, so a
      reader cloning the repo gets none of the inputs behind any SI table. This unblocks
      the Data and Software Availability statement and `CITATION.cff`.

---

## 4. Small code items worth doing anyway

Not blockers. Each is minutes, and each is something a reviewer could hit.

- [x] `benchmark/quant/dilution_track.py` — `species()` treated any accession
      containing the substring `UPS` as UPS1, so `sp|Q9UPS8|ANR26_HUMAN` would count
      as a spike-in. **Done 2026-07-31**, in `peptide_fc_de.py` too, which had the
      same function. Now matches the deposit's `ups|` tag or a `_UPS` suffix.
      Verified behaviour-preserving on this dataset: 18/303/5 proteins, rho
      +0.71/+0.02/+0.25, 105 DE calls -> 22 UPS1 / 23 E. coli / 60 human, all
      unchanged. (`de_validate.py` has a third, different classifier that maps UPS
      accessions to HUMAN; it is not on the published path, but it is a trap.)
- [x] `benchmark/quant/peptide_fc_de.py` — the contrast is asymmetric (`>= 2`
      non-zero **high** replicates only), so proteins higher at the dilute end are
      structurally excluded. **Done 2026-07-31**: documented at the code site and
      stated in Methods, which now says the calls are not a two-sided result.
- [x] Three `#[test]`s silently `return`ed when hard-coded `/home/...` inputs were
      missing, which reports as PASSED. **Done 2026-07-31**: paths now resolve from
      `CARGO_MANIFEST_DIR` instead of one developer's home directory (one pointed
      into a sibling checkout entirely), the skip prints a marked SKIP line, and
      `UNO_REQUIRE_TEST_DATA=1` turns it into a failure for CI and the pre-re-run
      checklist. `matching.rs` now also asserts the target/decoy PFM ratio, index
      validity, the allowed isotope offsets, the charge floor, and that a
      substantial fraction of features match — the ppm bound alone held for an
      empty match set.

---

## What was cut, and why

Removed from this list as not required to publish. All of it survives in
[`paper/AUDIT.md`](paper/AUDIT.md) with full evidence if it ever matters:

- **Repo hygiene** — documenting the ~40 result-affecting `UNO_*` environment variables,
  removing the dead config surface (`discriminant`, `[protein]`, `calibration_fdr` and the
  rest), adding CI. Real, but no reviewer sees them.
- **Stale in-code documentation** — done in the 2026-07-30 pass (`rt.rs`, `rt_cnn.rs`,
  `database.rs`, `score.rs`, `CHANGELOG.md`, `docs/`).
- **The mechanical sweep** — em dashes, semicolons, >40-word sentences, `$m/z$` rendering,
  the `20{,}647` brace escape. All six house-style rules **now pass**; `just prose-check`
  keeps them passing.
- **Prose polish** — LLM tics, repeated framings, terminology harmonization, tense
  consistency. Judgment calls that do not gate acceptance.
- **Cosmetic layout** — vector figures over PNG, ACS horizontal-rule table styling, table
  justification, page-float placement. The 300 dpi PNGs already meet ACS's colour-art rule.
- **Optional citations** — the "shows command of the niche" additions. The load-bearing
  ones (entrapment methodology, AMT-tag lineage, multi-enzyme, the Percolator correction)
  are already in.

Two cuts that are **not** safe to forget, so they moved into `paper/SUBMISSION_TODO.md`
rather than disappearing: the SI's unrunnable S1 commands (D2) and the two
"GENERATED … do not edit by hand" labels that no generator produces (D1).

## uno low-input bootstrap floor (found 2026-09-04, PXD066701 input series)

On `CW0710_OldFAIMS_1ng__RA7_rep6_ch2_DIA200pg` (1 ng HeLa, Astral, 79,917
features) uno returned **0 proteins at 1%** while its replicates gave 187, 567,
688, 975. The mass-anchor bootstrap found only 23 confident proteins (205 RANSAC
inliers, RT residual 0.12 min, i.e. the calibration itself was fine), so the
stage-2 calibration-pair floor at `crates/uno-core/src/score.rs:998`
(`if cal_pairs.len() < 100 { break; }`) tripped and the run exited with an empty
list instead of falling back to the stage-1 (mass-anchor) protein set. Two other
1 ng replicates (rep4: 50 confident, 187 final) sit just above the floor, which
is why the 1 ng level is so variable. Fix candidates, in order: (1) when the
floor trips, score and report the stage-1 confident set rather than nothing;
(2) widen the anchor gate or lower the floor when the initial confident count is
small. Not changed during the preregistered input-series run (frozen binary).

**Fixed 2026-09-04 (uncommitted, `score.rs`):** the stage loop is now a
closure `run_stages(floor, seed)`; when stage 1 returns <20 confident proteins
it is retried with the mass-bootstrap proteins seeding stage 2, once at the
estimated RT sigma and once floored at 1.0 min, keeping the deeper result. The
two collapsed files recover (1 ng rep6 0 -> 523, 200 pg rep1 0 -> 389 at 1%);
1 ng rep1/rep4, 200 pg rep2, and 10 ng rep1 are byte-identical to the frozen
binary. `results/uno_fixed/` holds the reruns; `results/uno/` keeps the frozen
outputs. Needs a unit test for the collapse path and inclusion in the C1 re-run.

**Shipped 2026-09-04 as the default (`UNO_SMOOTH` level 2):** the stage loop
and the ML stage are `run_stages` / `run_ml` / `finalize` closures; the whole
pipeline runs at RT sigma ×1.0 and ×1.5 (plus the seeded retries on collapse)
and the candidate whose FINAL list has the larger soft depth Σ clamp(1 − q/0.05)
is kept. Against the current binary on the 14 Astral + 6 primary files: 1 ng
rep4 187 → 768, 1 ng rep2 688 → 925, 200 pg rep3 476 → 458 (5%: 786 → 929),
09950 1454 → 1466, 16 files byte-identical. Wall +50–80%. Levels 1 and 3 were
tested and rejected (see score.rs comment). `UNO_SMOOTH=0` restores the old
path. NOTE: the July primary reference outputs differ from the current binary
on 09734 (1035 → 946) and 05433 (2072 → 2327) independent of this change;
that is the C1 re-run.

## uno timsTOF collapse-to-zero (found 2026-09-11, HYE Ultra 2 `--raw`, fixed uncommitted in `score.rs`)

**Symptom.** Same `.d` file, 20 seeds: 5 runs report 0 proteins, the other 15
span 155–328. Nothing measurable against that floor (IM gate, multi-run, any
timsTOF number).

**Cause.** The mass-anchor bootstrap picks "confident" proteins with a walk over
the INSTANTANEOUS decoy/target ratio, so the set is essentially the targets
ranked above the first decoy protein. On timsTOF at background p≈0.22 the
mass-only binomial separates only ~20 proteins from the decoys; where the first
decoy lands is a property of the decoy draw. 21 targets above it with one seed,
7 with another. 7 proteins carry ~35 anchor peptides, below the 50 both
calibrators need, so mass stays at the 3.5 ppm default, RT calibration is
skipped and every stage sees an uncalibrated run. A proper suffix-minimum
q-value does not help (below rank ~20 targets and decoys interleave).

**Fix.** Post-loop fallback only: when the whole bootstrap ends under 50 anchor
peptides, take the top-scored target proteins regardless of q until ~250 unique
anchor peptides (contamination is absorbed by the median / clipped-sigma /
RANSAC calibrators), then estimate RT sigma on that path with the Gaussian-core
peak instead of `clipped_sigma` (the contaminated set inflates the clipped
estimate to 1.02 min vs a 0.38 peak; a healthy seed's clean set reads 0.50 /
0.38). Applying the rescue to any starved ROUND was tried first and rejected:
seeds whose round 0 starves but whose later rounds recover lost depth to the
contaminated set (302→124, 301→200, 277→137, 221→170).

**Result (20 seeds).** Collapses 5/20 → 0/20; the five rescued seeds report
205 / 197 / 265 / 144 / 187 (honest entrapment @1% 340 / 219 / 301 / 284 / 229),
inside the healthy seeds' range. All 15 healthy seeds byte-identical. Cannot
fire on any paper run: the minimum anchor count over every pinned Orbitrap and
Astral log is 165 (threshold 50); hela5min_09746 re-run byte-identical.

**Still open.** Three of the five rescued seeds (10, 16, 18) find 0 proteins at
stage 1 under every sigma multiplier and reach their result through the
bounded-RT recovery path, not the primary one. The rescued calibration is
sound (RANSAC 71 inliers, 0.17 min residual, same as the seeds that succeed);
stage 1 is simply at the q≤1% edge on this data. Keeping the RANSAC calibration
for the seeded retry (instead of `fit_robust` on the seed set) changed nothing.
The 155–328 spread among healthy seeds is the next stability item.
