# KOTH matrix input

## KOTH 0.7.0 compatibility

UNO accepts long-bundle schema versions 1, 2, and 3. KOTH 0.7.0 writes schema 3;
regenerate both features/hills and aligned LFQ bundles when upgrading because
feature ranking, group admission, and signal ownership have changed.

Use `benchmark/koth_align.toml` as the current alignment/LFQ profile:

```sh
koth_align BATCH_DIRECTORY --config benchmark/koth_align.toml --output BUNDLE
uno --matrix-long BUNDLE --matrix-joint --fasta proteins.fasta --out-proteins proteins.tsv
```

The extraction half-windows remain `rt_window_pct = 0.005` and
`im_tolerance = 0.015`. The new `[lfq.consensus] max_group_qvalue = 0.05`
replaces `min_member_combined_score`, `min_seed_combined_score`, `min_group_size`,
and `allow_replicated_weak_seeds`; the removed settings now fail configuration
loading. Groups need at least two original runs. Grouping coordinate limits
are separate from extraction windows.

Schema 2 adds `group_score` and `group_qvalue`; schema 3 adds signal-ownership
status, owned/excluded sample counts, competing-group identity, and preceding
signal fraction. UNO validates these fields, group consistency and admission,
and suppressed-cell zero intensity. It retains the original observations even
when shared or ambiguous LFQ signal has been zeroed. Group confidence is not
used as extraction confidence or peptide/protein q-values, nor added to the
classifier. Original feature quality and run-specific evidence remain distinct.

## Long-format evidence bundles

KOTH can export `matrix_manifest.json`, `lfq_matrix.long.tsv`, and
`feature_observations.tsv` with `[output] export_long = true` in its alignment
configuration. The observation table preserves all original features, including
alternatives and features excluded from the consensus. The cell table contains
each quantitative extraction once and links to its selected original observation.

```sh
uno --matrix-long aligned/matrix_manifest.json --fasta proteins.fasta \
    --out-proteins proteins.tsv --out-lfq precursors.tsv \
    --out-abundance abundance.tsv
```

This route restores the seed observation's original isotope count, scan count,
and shape descriptors for both calibration and classification. Native coordinates
remain in the observation table; consensus search uses aligned seed coordinates.
Original feature measurements are distinct from LFQ grid measurements. Missing
original observations in MBR cells do not acquire fabricated feature quality.

For independent per-run identification with shared computation:

```sh
UNO_RT_CACHE_DIR=./rt-cache uno --matrix-long aligned/matrix_manifest.json \
    --matrix-batch-out batch-results --fasta proteins.fasta
```

The batch command loads one database and each model once, retains RT predictions
in memory across runs, and calibrates/searches each original run independently.
It writes per-run protein, peptide, precursor and abundance tables plus an
accession union and batch manifest. The union is a union of per-run 1%-q calls;
it is not assigned a calibrated cohort q-value. Consensus minimum-run filtering
does not remove original observations from the batch route.

For controlled comparisons, `--matrix-no-feature-quality` disables restored
original descriptors and `--matrix-legacy-values` reproduces wide-TSV numerical
rounding. Neither is needed for normal long input. The optional diagnostic
`--matrix-prefer-isotope-support` selects the best actual member with at least
three isotopes, including that member's aligned coordinates and original quality;
it falls back to the seed when no such member exists. Setting
`[scoring] matrix_recurrence = false` suppresses only the contributing-run and
positive-cell-count classifier columns for multi-run input; it leaves input
rows, quantities, calibration and feature tie-breaking unchanged.

The manifest and table hashes, row identities, primary/seed links and distinct-run
support are checked on import. A changed bundle must have a new manifest.

## Existing wide matrices

```sh
uno --matrix intensity_matrix.tsv --fasta proteins.fasta \
    --out-proteins proteins.tsv --out-lfq precursors.tsv \
    --out-abundance abundance.tsv
```

This performs one search in the consensus reference RT space. Each assigned
precursor carries its intensity across the matrix's runs. A cohort identification
does not establish independent detection in every run.

Required metadata columns are `massCalib`, `mz`, `charge`, `rtApex`, and
`n_contributing_runs`; `im` and `combined_score` are optional. Run intensity
columns follow `n_contributing_runs`. Pass KOTH's `intensity_matrix.tsv`, not
`consensus_features.tsv` or `qvalue_matrix.tsv`. The default
`--matrix-min-runs 2` requires contributions from two original feature-finding
runs, not merely two positive LFQ cells. `--profile-reassign` is an independent,
opt-in operation on quantitative assignments after identification.

## Experimental extraction-quality evidence

Ordinary matrix input does not automatically consume KOTH's q-value or LFQ-detail
sidecars. The fresh matrix benchmark adds a strict preprocessing step that
checks feature identity, run names, target-cell coverage, intensities, and
contribution counts before constructing quality summaries. It rejects a sidecar
from a different matrix even when the dimensions match.

The optional schema consists of all nine columns below, placed before
`n_contributing_runs`. Partial schemas, nonfinite values, negative values, and
out-of-range fractions are rejected. The classifier receives these as distinct
LFQ features; they do not replace the single-run isotope-quality columns.

| Column | Meaning in the study preprocessing |
| --- | --- |
| `matrix_confident_fraction` | Fraction of all runs with positive intensity and extraction q <= 0.01 |
| `matrix_confident_direct_fraction` | Fraction of all runs with a confident original contribution |
| `matrix_confident_mbr_fraction` | Fraction of all runs with a confident transferred extraction |
| `matrix_median_log_confidence` | Median -log10 extraction q over runs, clipped at q = 1e-6 |
| `matrix_median_isotope_pattern` | Median isotope-pattern agreement among confident cells |
| `matrix_median_coelution` | Median isotope coelution among confident cells |
| `matrix_median_isotopes` | Median isotope count among confident cells |
| `matrix_median_abs_ppm` | Median absolute extraction mass residual among confident cells |
| `matrix_median_abs_rt` | Median absolute extraction RT residual in minutes among confident cells |

The denominator for each fraction is the total number of runs. Rows without
confident observations receive zero for the confident-cell medians and retain
their zero support fractions. KOTH extraction q-values assess extraction errors;
they are not UNO peptide or protein q-values.

`benchmark/matrix_evidence.py` writes separate confidence-masked, quality-enriched,
and all-zero-quality control matrices. The all-zero control holds the classifier
dimension constant when testing the added evidence. Ordinary matrices and all
single-run inputs keep the original feature-vector dimensions and order.

This is an experimental interface under evaluation in the
[fresh benchmark](koth-matrix-benchmark-plan.md), not a claim of improved accuracy.
# Experimental joint search

Use `uno --matrix-long BUNDLE --matrix-joint --fasta DATABASE --out-proteins proteins.tsv`
to perform one cohort search using original aligned observations. The matrix is
the sole experimental input: no individual-search protein lists are loaded.
FASTA digestion and RT prediction are shared, and overlapping member mass windows
are searched once per group.

Each candidate is assessed against actual mass/RT pairs and their original
isotope/shape measurements. At most one alternative per source run contributes
to agreement summaries, and duplicate source-table hashes share an evidence ID.
The classifier receives agreement as learned descriptors; the implementation
does not multiply replicate probabilities. Protein inference counts unique
peptides and estimates target/decoy error once for the cohort.

Joint input defaults to `--matrix-min-runs 1`, preserving singleton groups in
legacy schema-1 bundles; KOTH 0.7.0 itself requires two-run group support.
Missing runs and MBR-only extraction cells do not create direct peptide
observations. Optional extraction covariates are described below.
Quantification retains the original integrated LFQ cells exactly once;
it does not substitute original apex intensities for those quantities.

The experimental `[scoring] matrix_joint_agreement = false` ablation retains the
same candidate observations and disables agreement descriptors. Ordinary wide
and representative long-input searches retain their existing behavior and a
default two-run floor. The experimental route is not yet a recommended default;
see the pilot results for measured sensitivity and error.

Legacy schema-1 bundles can contain groups admitted by
`allow_replicated_weak_seeds`; UNO retains their original seed-quality floor for
initial calibration. Schema-2/3 bundles use the declared group-confidence gate
instead, without silently restoring the removed seed floor. UNO's existing
isotope/feature anchor criteria still apply. Later cohort recalibration uses
actual mass-compatible observations of confidently supported proteins.

## Optional joint extraction-quality scores

For `--matrix-long ... --matrix-joint`, set `[scoring]`
`matrix_joint_cell_quality = true` to append twelve learned row descriptors.
This experimental option defaults to false and leaves candidate generation,
original-observation support and integrated quantitative cell values unchanged.

The descriptors are MBR-filled count and fraction, confident MBR count, confident
cell fraction across the cohort, and mean q, worst q, q-value coverage and confident
fraction separately for direct and MBR cells. A filled cell has positive integrated
intensity; a confident cell also has an available extraction q-value <= 0.01.
Mean and worst q use only filled cells with known q-values. If none are known,
both use 1 with coverage 0, distinguishing unavailable values from perfect q=0.
Within-stratum coverage/confident fractions use that stratum's filled-cell count;
cohort fractions use the number of distinct source files. Counts are untransformed.

Original feature alternatives and duplicated source files cannot inflate these
counts. Synthetic extraction-decoy cells are excluded; the same target-cell row
summary is supplied to every target and decoy peptide candidate matching the row.
These are discriminant features, not multiplied probabilities or hard filters.
KOTH extraction q-values measure extraction confidence, not sequence identity.

The bundle preserves both measurement types. `feature_observations.tsv` contains
original feature measurements; `lfq_matrix.long.tsv` contains recomputed extraction
measurements such as `lfq_hybrid_score` (the geometric mean of RT, intensity,
isotope-pattern and coelution scores in Hybrid mode). Its plain `combined_score`
column is the original consensus seed score. This experiment uses extraction
q-value summaries and fill counts; it does not add raw hybrid-score summaries.

### Cell-score ablations and competition experiment

With `matrix_joint_cell_quality = true`, the optional scoring flags
`matrix_joint_cell_counts` and `matrix_joint_cell_q` each default to true.
Disabling counts zeros only MBR-filled count/fraction; disabling q zeros the
remaining ten q-derived descriptors, including confident counts and q coverage.
Column positions stay fixed, and original-observation agreement remains enabled.

`matrix_joint_alternative_groups = true` experimentally supplies the protein
competition engine with additional matrix rows having an actual member inside
the calibrated mass gate and combined mass/RT deviation <=9. It preserves the
original best row and counts each peptide once, even with several available rows.
The option defaults to false. Its pilot result was worse than the single-row
configuration, so it should not be promoted. Quantification still uses the chosen
best row and its integrated LFQ cells; this experiment changes identification
competition, not the quantitative estimator.

`diagnostic_peptides_path` optionally exports the selected calibration's peptide
candidates, ML rank, ensemble inclusion, and group-anchor observation IDs.
The export follows final candidate selection and does not overwrite itself with
a losing calibration's scores. It is a diagnostic candidate table, not a table
of accepted peptide identifications. It never becomes a search input.
