# UNO rerun with KOTH 0.7.0

Status: 23-file scope approved; matrix arm corrected to consensus-row search
following user feedback. Full rerun not started.
This replaces the earlier 32-file/eight-matrix/144-search matrix proposal.
Machine-readable selection: `benchmark/matrix_publication_koth070.json`.

Integration verification: all 93 UNO/uno-core tests pass, and release binaries
build offline against the released KOTH source. A two-file Astral compatibility
smoke regenerated features and schema-3 LFQ, preserved all 171,624 original
features on round-trip, and matched the wide export exactly for 23,333 groups.
Joint input retained 49,394 aligned observations; ownership-suppressed quantities
remained zero. One seed-zero joint search completed in 74.31 seconds with peak
RSS 1,234,556 KiB on four threads, using a fresh RT-cache directory. This is a
compatibility measurement, not a full-cohort runtime estimate or publication
result. Local logs, output hashes and receipts are in the ignored
`benchmark/results/koth070_compatibility/` directory. A consensus-row search
smoke (the completed search used the removed joint arm), Bruker and IonStar
pilots, standalone/embedded raw-route equivalence, and publication-driver
adaptation remain execution prerequisites below.

## 1. Use the release-compatible build on main

KOTH v0.7.0 is published. Its annotated tag resolves to
`88f9e7c54c9ad3d6ca5a56a906703a3d46982f43`, verified against GitHub. The release
page's automatically generated body mentions 0.6.0; the tagged source changelog
and tag commit establish the actual 0.7.0 changes. UNO retains its sibling path
build workflow, with an exact `=0.7.0` version requirement. Freeze source hashes
as well as the version; a local path is not an immutable source pin.

The release changes feature candidate ranking and best-prefix selection,
introduces bounded native MS1 streaming, replaces old LFQ admission settings
with group confidence, and applies exclusive signal ownership. Unlike the
0.6.0 LFQ-only update, **both feature/hill tables and LFQ must be regenerated**.
None of the old KOTH-derived separate-search results are reusable as new-version
results merely because they do not consume an LFQ matrix.

Use `benchmark/koth_uno.toml` for Orbitrap feature extraction and the embedded
Bruker profile for timsTOF. Both explicitly carry the release's isotope-evidence
parameters. Use `benchmark/koth_align.toml` for alignment/LFQ:

| Setting | Rerun value |
| --- | --- |
| Extraction RT/IM half-windows | 0.005 of run RT span / 0.015 absolute IM |
| Group coordinate limits | 20 ppm / 0.02 RT span / 0.05 IM |
| Group admission | `max_group_qvalue = 0.05`, at least two original runs |
| Quantitation | summed grid extraction, no upstream normalization, no averagine projection |
| Exports | schema-3 long bundle, wide TSV, details and decoys |

Remove `min_member_combined_score`, `min_seed_combined_score`, `min_group_size`,
and `allow_replicated_weak_seeds` from new LFQ configurations. UNO reads schemas
1–3, validates new confidence/ownership fields, and preserves zeroed shared or
ambiguous cells. Schema-1 bootstrap behavior stays reproducible; schema-2/3 input
uses the new group admission without silently adding the retired seed floor.
Original UNO feature/isotope anchor checks still apply. Group q-values, cell
extraction q-values and identification q-values remain distinct.

Parquet + Zstd IM caching is already committed on main (`3ba6be0`). It populates
automatically and shares model/sequence/charge predictions. Set an explicit
`UNO_IM_CACHE_DIR` for durable study storage: the default is an evictable OS-temp
cache. Set `UNO_RT_CACHE_DIR` too. Predictions remain reusable when only KOTH
changes; new shuffled decoys or charges cause additional misses. Do not reuse
feature-dependent outputs or calibrated residuals. Cold cache generation and
warm search timing must be reported separately.

## 2. Demonstrate direct identification from an LFQ matrix

| Cohort | Raw files | Alignments | Separate searches | Consensus-matrix searches |
| --- | ---: | ---: | ---: | ---: |
| Astral 100 ng | 3 | 1 | 9 | 3 |
| IonStar full concentration series | 20 | 1 | 60 | 3 |
| **Total** | **23** | **2** | **69** | **6** |

The question is whether UNO can identify precursors directly from an aligned
LFQ matrix. Compare two input routes:

1. **Individual files:** each original feature list is searched independently
   with `--features`; combine the resulting identification lists for the cohort
   comparison.
2. **Aligned matrix:** KOTH produces one consensus precursor per row, with
   aligned mass/RT/charge and one LFQ intensity per run. Search those rows once
   with `--matrix-long BUNDLE`, using the existing representative-consensus
   loader. Each row is one search feature; its intensity vector supplies
   quantification. Do not pass `--matrix-joint` or enable joint-cell scoring.

The consensus route uses one set of coordinates and representative original
feature-quality descriptors per row. It does not expand rows into their member
observations for multi-run scoring. LFQ quantities are used for quantification;
LFQ cell-quality summaries are not classifier features in this comparison.

Use seeds 0, 1 and 2 for both routes. Keep the database, digest, RT model, matching
filters and hardware policy matched. Compare identification depth and error
between the individual-file union and the consensus search. Assess the matrix
workflow's LFQ quantities against IonStar's known ratios; a second experiment
mapping individual-search assignments into LFQ is not required for this section.

This is **75 searches, 48% fewer than 144**, plus four serialized timing controls
(three Astral single searches and one Astral consensus-matrix search, seed zero).
The chosen scope keeps all four IonStar replicates at all five concentrations and all three
Astral 100 ng acquisition replicates. The input inventory is inherited by name
and checksum, not reselected by new results; all 23 source paths are present.
Recheck their checksums at execution freeze.

Drop the IonStar 5/10-run scaling matrices, the four other Astral matrices, and
the multi-observation joint-search arm and its cell-summary ablation from this
new publication run. Preserve the old full studies as versioned development evidence. They must not be relabeled
as 0.7.0 results. The narrower matrix section makes no low-input or cohort-scaling
claim. Astral 100 ng is a focused high-input demonstration, and IonStar provides
known-ratio quantitation. The previously unfavorable joint-search results remain
in the development archive; they are not results for the consensus-row route.
Both datasets have development history; this is not unseen validation.

## 3. Rerun the other active KOTH-dependent evidence

| Experiment | Required jobs | Why |
| --- | --- | --- |
| Raw-mode figure | 18 inputs × three seeds = **54 searches** with embedded extraction | New feature ranking affects every KOTH raw-mode result, including all 12 Astral input-series files. |
| Current timsTOF IM comparison | Fresh A2/B2/C2 features; three model settings × three files × three seeds = **27 searches** | Changed KOTH features affect candidate sets and calibrated IM evidence. |
| Diagnostics and controls | Route parity, model/default-disable checks, and any KOTH-derived A1 observations | Determine affected inputs by provenance; public training/holdout labels alone do not depend on KOTH. |

The known publication core is now **156 searches plus four matrix timing
controls**, rather than 225 plus four. Matrix extraction contributes 23 fresh
feature/hill tables; IM adds three; the 54 raw-mode searches also extract
internally. Reuse newly extracted tables across applicable arms only when raw
input and every effective extraction setting match. A cached TSV search cannot
replace an end-to-end `--raw` resource measurement.

Biosaur2-only and ms1searchpy evidence does not depend on KOTH. Keep it with its
recorded provenance and verify a matched biosaur2 control on the new UNO build.
Expand those reruns only if relevant UNO changes alter their behavior. Public
RT/IM model training does not need repeating solely because KOTH changed.

## 4. Freeze, pilot, and execute in dependency order

1. Archive the committed UNO tree, exact KOTH tag/source, other path dependencies,
   lockfiles, binaries, effective configurations, FASTA/model hashes, input
   hashes, run/condition mapping, seeds, worker limits and driver hashes in a
   new `benchmark/results/koth070_<date>/` directory. Historical scripts with
   hard-coded 0.5.0/0.6.0 paths are reproduction records, not launch commands.
2. Adapt the execution/report drivers to the approved JSON inventory and current
   TOMLs. Require complete input/output hashes for resume, including every long
   bundle file, and successful receipts for exactly 23 extractions, two matrices,
   75 searches and four timing controls. Keep raw-mode and IM outputs isolated.
3. Use small, prespecified compatibility pilots before the full run: representative
   Orbitrap/FAIMS and Bruker extraction; the established four-file IonStar extreme
   subset for LFQ. Check actual schema-3 bundle ingestion, original observation
   links, group-versus-cell confidence, ownership-suppressed zeros, and standalone
   versus embedded feature equivalence at TSV precision. No tuning on full-cohort
   outcomes. Archive pilot attempts independently.
4. Measure fresh extraction, ten-control group-confidence construction, exclusive
   LFQ, and search/cache time and memory. Set concurrency and the execution budget
   from these results. The old 3–3.5-hour matrix estimate covered a larger study
   and older LFQ; fewer searches do not imply a proportional 0.7.0 speedup.
5. Execute raw → features/hills → alignment/LFQ → matched individual/consensus-row searches
   → evaluations. Reuse predictions and completed matching jobs with valid hashes;
   serialize the resource controls and disclose cache state.

## 5. Publish only the audited replacement evidence

Recompute per-seed accession depth, overlap, additions/losses, internal-q lists
and empirical entrapment curves with whole tied-score blocks. Independently
check 1%, 3% and 5% accession sets. Keep failures, empty outputs and regressions.
The union of individually filtered lists is not automatically a calibrated
cohort list; group confidence is not calibrated identification FDR.

For IonStar, retain expected ratios **1:1.5:2:2.5:3**, E. coli ratio error and
compression, within-condition CV, human background stability, and coverage and
missingness. Fix normalization/rollup rules, rederive their factors from the new
outputs, compare common and arm-specific groups, and avoid primary-analysis
imputation. Count original detections, confident extractions and MBR-only cells
separately. Evaluate ownership and isotope-choice limitations rather than
assuming that exclusive extraction guarantees better ratios.

Keep the main-text matrix result to one concise developed paragraph and a compact
SI comparison for the two cohorts: identification, known-ratio accuracy/coverage,
and one clearly scoped search-resource control. Detailed per-seed records stay
in the reproducibility archive; remove scaling/low-input/cell-ablation claims
from the replacement section. Retain negative IonStar outcomes if they recur.

After all audits pass, repoint generators to the new reports and regenerate the
raw-mode figure, compact matrix figure/tables, affected IM tables and all derived
statistics. Update software versions, Methods, captions and limitations. Run
`just assets-all`, `just fmt`, `just paper`, `just docx`, `just verify`, and
`just check-stats-deep` from `paper/`; visually inspect changed PDF/Word pages.
Check manuscript paragraphs for repeated explanations and one-sentence fragments.

Current manuscript numbers remain explicitly tied to the archived 0.6.0 LFQ /
0.5.0 detection study until replacement results exist. Integration checks in
this task do not count as the 75-search publication rerun.
