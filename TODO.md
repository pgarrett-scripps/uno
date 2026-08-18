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

**Critical path: land the three code fixes → re-run → correct the numbers → submit.**

---

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

## 2. The re-run

The gate for every identification number in the paper. Do it **once**, correctly.
Budget half a day, serialized on one idle host, `/usr/bin/time -v` on every invocation.

**Preconditions**

- [x] §1 above landed and the determinism test passing.
- [ ] Working tree committed. **Record the commit hash** — it becomes the single value in
      `paper/si/versions.typ`.
- [ ] The ms1searchpy patch committed to a branch or captured as a diff, with **its**
      commit recorded. (The comparator is a local uncommitted +45-line patch, not released
      2.8.12 — which does not exist on PyPI. Without it, stock ms1searchpy reads none of
      the five quality columns that §3.5 and S6 depend on.)
- [ ] Log the uno build metadata alongside each run. Last round's `benchmark/results/*.log`
      are zero-length, which is why the build could not be identified afterwards.

**Scope**

- [ ] Six datasets: uno ×1 + ms1searchpy ×3.
- [ ] run10: uno ×1 + ms1searchpy ×3. run10 is a **tuning set**, not a benchmark.
- [ ] run10 column-stripped: uno ×1 + ms1searchpy ×3. Fills the one unauditable row of
      Table S5, whose caption currently admits the uno arm "was not retained".
- [ ] FAIMS: ms1searchpy **×3 per replicate**, not ×1. The previous round ran it once,
      which violates the paper's own stated standard.
- [ ] **Dump uno's PFM table** for run10, so Table S6's uno column can be re-derived.
- [ ] Quantitation: re-run `uno --matrix` at the pinned commit on the existing
      `uno_proc_newkoth/align_avg/` matrix, and **save the log**. Do *not* re-run koth
      alignment — that is the expensive step and its output is unchanged.

**After**

- [ ] Regenerate `benchmark/results/summary.tsv`, then `paper/figdata/*.tsv`, then
      `cd paper && just assets-all`.
- [ ] Re-derive every table from the new outputs **by script**, and keep the script.
- [ ] Then the manuscript work: `paper/SUBMISSION_TODO.md` C3, C4, C6, then C2 and C5.

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
