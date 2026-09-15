# KOTH next-release rerun plan

Status: planned; start after the new KOTH release is live. The release tag and
commit are intentionally unset. This plan covers KOTH-dependent evidence in
UNO, including feature generation, LFQ, downstream searches, and manuscript
outputs. No benchmark was launched while preparing this plan.

## 1. Freeze the release and prepare an isolated rerun

Create `benchmark/results/koth_<tag>_<date>/` with separate `raw_mode`, `matrix`,
`im`, `pilot`, and provenance directories. Preserve previous studies and the
current manuscript until the replacement evidence passes validation.

1. Resolve the published tag to its commit; verify release checksums. Build
   `uno`, `koth_ff`, and `koth_align` against that exact KOTH source. UNO uses a
   local path dependency (`Cargo.toml`), so installing standalone binaries does
   not update its embedded feature finder. Verify the actual build outputs
   rather than assuming a particular Cargo target directory.
2. Freeze the UNO revision and all local changes, dependency sources and locks,
   binary hashes, runtime versions, FASTAs, RT/IM models, raw-input checksums,
   run/condition manifests, seeds, thread counts, and effective configurations.
   The current UNO tree contains ongoing IM work: explicitly select and archive
   the intended state before building; a commit hash alone is insufficient.
3. Review the release's feature and LFQ configuration/schema changes. Record
   each required migration. Preserve the existing experimental settings where
   supported; declare any new default or estimator comparison as a separate
   arm before inspecting full-cohort results.
4. Adapt the existing drivers to accept the new study root and release identity.
   Remove implicit reuse of old generated features, hills, matrices, search
   receipts, and feature-dependent caches. Resume only when binary, complete
   input/configuration hashes, and output hashes match the new run.
5. Produce a dry-run inventory mapping every active KOTH-dependent claim to its
   inputs, jobs, evaluator, and final asset. Freeze this inventory before the
   full run, including any additional active KOTH studies landed since this plan.

**Driver changes needed:** `benchmark/matrix_release_study.py` hard-codes 0.5.0,
old source/configuration paths, and its output directory;
`benchmark/matrix_release_report.py` also asserts the old release tag.
`benchmark/publication_rerun.py` runs both input arms and writes canonical paths;
add scoped selection and isolated outputs before using it here. The IM drivers
inherit feature paths from old receipts. Hash every matrix bundle member used
by a search, not only `matrix_manifest.json`. Archive effective LFQ defaults,
including values not explicitly present in TOML.

## 2. Pass a small feature-generation and LFQ pilot

Use the prespecified four-run IonStar subset (two replicates from each extreme)
for LFQ checks and runtime estimates. Add representative low-input Astral,
multi-CV FAIMS, and Bruker timsTOF files for feature/interface coverage. Keep
pilot outputs separate; rerun them if the frozen settings subsequently change.

1. Compare standalone KOTH TSVs with embedded UNO `--raw` features under matched
   configurations and matching filters. Check feature membership, mass/RT/IM
   values, charge, quality columns, missing-value handling, and FAIMS separation.
   Account for TSV rounding; investigate differences rather than assuming
   byte-identical search results. Check Bruker-specific filters explicitly:
   `--features` does not automatically inherit `--raw` defaults.
2. Exercise long-matrix joint search and the supported ordinary matrix route.
   Validate schema, consensus/run identities, original observations, intensity
   cells, extraction q-values, direct/MBR labels, and complete sidecar joins.
   Verify integrated LFQ cells are counted once and missing cells remain missing.
3. Review feature counts, retained groups, detected/MBR support, missingness,
   intensity distributions, identification calibration, and IonStar ratios
   against the archived baseline. Changes are expected; unexplained failures,
   malformed values, and broken provenance stop promotion to the full run.
4. Measure extraction, alignment/LFQ, and search time, peak memory, and output
   size. Use those measurements to set worker limits, disk requirements, and
   a concrete full-run estimate. Avoid transferring old runtime estimates to
   the new algorithms.

**Gate:** all representative routes load and complete, identities and quantities
pass audit, configuration migrations are documented, and full-run settings are
frozen. Parameter exploration belongs in the pilot and must be recorded.

## 3. Regenerate every active KOTH-dependent experiment

The dependency order is **raw inputs → features/hills → alignment/LFQ → UNO
searches → evaluation → manuscript assets**. Reuse newly generated features
across arms only when raw inputs and full extraction settings match. A cached
feature-table search cannot stand in for a measured end-to-end `--raw` run.

| Evidence | Required rerun | Existing starting point |
| --- | --- | --- |
| Raw-mode comparison | Six primary files plus 12 Astral files, seeds 0–2: **54 `--raw` searches**, including fresh embedded feature generation. Refresh seed variability, overlap, entrapment, and raw-route resources. | `benchmark/publication_rerun.py`, `benchmark/variability_extension.py`, `benchmark/raw_mode/summarize.py` |
| Matrix/LFQ study | **32 fresh feature/hill extractions**, **eight independent alignments**, **144 searches**, and **four serialized timing controls**. Regenerate all long/wide LFQ exports, observations, q-values, details, and downstream proteins/precursors/abundances. | `benchmark/matrix_release_study.py`, `benchmark/matrix_release_report.py` |
| Current IM search comparison | Regenerate timsTOF **A2/B2/C2** KOTH features; rerun off/tree-soft/CNN-soft × three runs × seeds 0–2: **27 searches**. Refresh default versus explicit-model/disabled checks and their evidence. | `benchmark/im_cnn_study.py`, `benchmark/im_cnn_check.py`, `benchmark/im_default_20260914.json` |
| Other KOTH-dependent evidence | Rebuild KOTH-derived A1 observations/calibration/diagnostics if their lineage uses KOTH features; refresh any still-active quality, ablation, FAIMS, or legacy IM claims found by the inventory. | `benchmark/im_cnn_diagnostics.py`, `benchmark/im_published_comparison.py`, `experiments/im-timstof/uno/benchmark/im_timstof_study.py` |

The matrix design stays matched to the current release study: IonStar nested
5/10/20-run cohorts; Astral 200 pg/1 ng/10 ng/100 ng triplicates plus the 12-run
pool; separate searches, joint searches, and joint searches with cell-quality
summaries; seeds 0–2. Rebuild each subset's alignment independently. Do not slice
the full matrix to represent a smaller cohort.

The known core is **225 search jobs plus four matrix timing controls**, before
pilot checks, additional regressions, or new arms. It includes 32 matrix
extractions and three timsTOF search-file extractions; the 54 raw-mode jobs also
perform extraction internally. These counts are a planning baseline, not a
substitute for the frozen inventory.

Reuse immutable raw/mzML inputs, FASTAs, independent MS/MS evidence, and models
whose training data do not depend on KOTH. Public RT/IM training and public
holdout prediction benchmarks do not automatically require retraining. Audit
A1 observation lineage before deciding whether its prediction comparisons
also need regeneration.

Biosaur2-only and ms1searchpy results remain independent evidence. Verify matched
biosaur2 controls with the frozen UNO build; if non-KOTH changes alter its search
behavior, expand the inventory to the affected comparator searches. Historical,
retired experiments remain archived and explicitly versioned; any result used
as evidence for the new version must be regenerated.

## 4. Re-evaluate identification, LFQ, and resources

1. Recompute per-seed protein/precursor coverage, accession unions, overlap,
   internal-q lists, and empirical entrapment curves. Preserve whole tied-score
   blocks and independently audit the 1%, 3%, and 5% accession sets. Report all
   outcomes, including empty runs and regressions.
2. Separate cohort identification, original per-run feature detection, confident
   extraction, and MBR-only support. Recompute Astral MS/MS corroboration using
   the new assignments; an unmatched MS/MS record is not proof of a false ID.
3. Recompute IonStar expected ratios **1:1.5:2:2.5:3**, ratio error/compression,
   within-condition CV, human-background stability, and coverage/missingness.
   Preserve the normalization and rollup rules, but recalculate their numerical
   factors from the new data. Report common and arm-specific groups, without
   primary-analysis imputation. Map separate-search assignments to the same
   new LFQ intensities when comparing identification effects; apex intensities
   are a different estimator.
4. Re-measure serialized raw-route resources and the four warm-cache Astral
   matrix timing controls. Report extraction, alignment/LFQ, search-only, and
   end-to-end costs separately, with cache/concurrency policy explicit.
   Concurrent batch timings do not support serialized speed claims.
5. Produce an old-versus-new comparison with release/configuration identities,
   all seed outcomes, biological metrics, resources, and explained deviations.
   Audit every planned job and checksum before accepting completion markers.

## 5. Replace manuscript evidence after the audit passes

Update report readers to the new study root, regenerate summaries and statistics,
and refresh `fig5_raw_mode`, `figS7_koth_matrix`, the matrix SI tables, the IM search
table, and any other assets identified by the dependency inventory. Refresh A1
tables only if their underlying observations changed. Update Methods, Results,
captions, software versions, and limitations to reflect the actual outcome.

Run the existing manuscript pipeline from `paper/`: `just assets-all`, `just fmt`,
`just paper`, `just docx`, and `just verify`. Check provenance and stale-version
references, and visually inspect affected PDF and Word pages. Preserve coherent
multi-sentence manuscript paragraphs and remove repeated conclusions.

**Done means:** all planned KOTH-dependent jobs have successful audited receipts;
every current-version claim traces to new outputs; every intentionally retained
independent or historical result has explicit provenance; manuscript assets and
builds pass verification. Deliver a release-rerun report linking the frozen
manifest, comparisons, and rendered manuscript. Keep the previous version
available for reproduction.

## Runtime estimate and search-seed budget

Planning estimate from the archived job receipts: **4–6 hours of computation**
on the same host with approximately two concurrent workers, and **6–10 hours
total** including release integration, pilot review, reporting, and manuscript
checks. These are provisional allocations, not measurements of the unreleased
version; API/schema changes or substantially slower LFQ can extend them. Replace
the estimate with pilot measurements before the full launch.

| Archived workload | Sum of individual job wall times | Planning implication |
| --- | --- | --- |
| Matrix feature extraction and alignment | 30.1 minutes | This cost does not fall when search seeds are reduced. |
| Matrix searches | 341.5 minutes across 144 jobs | About 171 minutes at ideal two-worker utilization; allow scheduling overhead. |
| Matrix serialized timing controls | 5.1 minutes | Keep serialized. |
| Raw-mode searches | 14.6 minutes serialized plus 44.5 concurrent job-minutes | About 37 minutes with two workers for the repeated seeds. |
| IM searches | 41.7 minutes across 27 jobs | About 21 minutes with two workers, plus fresh Bruker extraction/checks. |

Retain **three search seeds (0, 1, 2)** for the final KOTH rerun. Every core
experiment in this plan already uses three, so reducing the unrelated historical
11-seed UNO sweep does not reduce the 225 planned searches. Three prespecified
seeds provide a practical descriptive check of search variability; they do not
establish precise uncertainty or replace independent acquisition replicates.
Keep every biological/acquisition run and report all seed outcomes.

The current `benchmark/results/spread.tsv` contains 11 UNO draws per primary
file and either three or ten ms1searchpy draws, depending on the file. Those
biosaur2-based historical comparisons are outside the KOTH-only rerun unless
the frozen UNO changes also affect their results. If that expansion is needed,
three prespecified UNO seeds can be the new primary comparison, with the larger
archived sweep retained as separately versioned supporting evidence.

The largest safe scheduling improvements are to extract each matrix/IM input
once per configuration, reuse its new features across seeds, share only validated
model-prediction caches, and run independent jobs concurrently within measured
memory limits. Additional workers may improve throughput; benchmark the new
version before assuming linear scaling.

IM caching is now implemented in the working tree. The new build precomputes
the complete target/decoy library at valid observed IM charges and retains raw
predictions by model/sequence/charge. Configure a shared `UNO_IM_CACHE_DIR` for
the rerun, alongside `UNO_RT_CACHE_DIR`; do not execute the frozen historical IM
driver unchanged to produce a new-version study. Warm the cache serially before
parallel timing runs, and measure cold precomputation separately. Model-unchanged
predictions remain reusable after KOTH changes; new decoy sequences and charges
are added as misses. The runtime estimates above predate this IM optimization
and must be updated from the new pilot rather than assuming a fixed speedup.

## Repository sources inspected for this plan

- `docs/matrix-release-20260913.md` and `benchmark/matrix_release_study.py`:
  current 32-run/eight-matrix/144-search design and configuration.
- `benchmark/publication_rerun.py`: 18 raw-mode files and three search seeds.
- `benchmark/im_cnn_study.py` and `docs/im-models.md`: KOTH-dependent paired IM
  searches and fixed public model training.
- `paper/FIGURE_MANIFEST.md` and `paper/analysis/scripts/gen_stats.py`: active
  figure/table/statistics dependencies and archived experiment distinctions.
- `docs/koth-matrix-benchmark-plan.md`: original pilot, evaluation, and
  normalization safeguards; the completed release study supplies the current
  joint-input design where the earlier protocol describes superseded arms.
