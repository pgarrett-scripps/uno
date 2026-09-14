# Fresh KOTH matrix-input benchmark

Prespecified planning record, 2026-09-12, written before the fresh searches.
Execution and verification are documented in `benchmark/matrix_study_README.md`.

## Question and scope

Can one UNO search of an aligned KOTH feature matrix replace separate searches
of its constituent runs, at comparable externally checked identification
confidence, while reducing repeated search work? Does cross-run evidence improve
identification, and does the resulting sample-level intensity matrix preserve
known abundance patterns?

Exclude all historical result values from the new paper analysis. Retain old
files as development records, but regenerate KOTH features, hills, alignment,
LFQ, UNO outputs, summaries, and figures in a new versioned output directory.
Reusing raw data is acceptable; reusing the old extracted matrix is not.

The primary contribution is joint identification across runs. Quantitative
accuracy and missingness validate the usefulness of the output; a new
differential-expression workflow is outside this experiment's scope.

## What the implementation currently supplies

UNO's matrix loader reads consensus mass, m/z, charge, reference RT, optional
IM, seed combined score, contributing-run count, and the intensity vector.
Its classifier already receives contributing-run count, summed intensity, and
the number of nonzero cells. Profile-guided reassignment can additionally use
cross-sample intensity profiles, but is opt-in and occurs downstream of the
identification stage. These are distinct mechanisms to test separately.

The matrix also loses individual-feature isotope counts, scan counts, isotope
quality and mass-error descriptors. More samples therefore do not automatically
mean more useful evidence under the current loader.

KOTH additionally exports `qvalue_matrix.tsv` and `lfq_details.tsv`. The latter
contains per-cell MBR status, isotope-pattern agreement, coelution, isotope
count, RT residual, peak width, mass residual, and IM residual. UNO currently
does not load these. An enriched-input experiment requires implementation;
it is not a description of existing support.

A positive matrix intensity is not an independent feature detection. Keep
original feature contributions, MBR extractions, and confidence-filtered
extractions distinct. KOTH extraction q-values and UNO identification q-values
address different errors and must not be substituted for one another.

Source anchors: `crates/uno-core/src/input.rs:97`,
`crates/uno-core/src/score.rs:214`, and
`../koth_rust/koth_ff/src/align_main.rs:401`.

## 1. Freeze inputs, versions, and a small route-equivalence check

Record UNO and KOTH revisions, any working-tree patch, binary hashes, FASTA and
model hashes, dependency locks, complete configurations, input checksums,
thread counts, and a run-to-condition manifest. UNO currently has uncommitted
changes; a commit identifier alone would not describe the build inspected here.

Generate KOTH single-run TSVs using exactly the configuration embedded by the
paired UNO `--raw` call. On a few prespecified representative files, compare
feature membership and numerical values, then search both inputs using matched
seeds. Report score/list agreement and any discrepancies. Check serialization
precision, feature order, missing IM handling, and instrument-specific defaults.
The TSV writer rounds values while the in-process path retains full precision,
so exact output identity is a hypothesis, not an established fact.

This is an interface check, not another large feature-finder benchmark. The
existing biosaur2-versus-raw comparison addresses feature-finder choice. If the
new experiment changes the KOTH revision used in that figure, refresh the
affected raw-mode searches on the same frozen build or clearly identify the
different experiment versions; never treat historical raw-mode numbers as a
matched baseline for a newly extracted matrix.

## 2. Use locally available data with complementary roles

**Primary recommendation: IonStar, PXD003881, all 20 runs.** Four replicates at
each of five E. coli spike levels in constant human background provide both
replicate agreement and known abundance changes. Raw files, converted mzMLs,
and processing recipes are present under
`/home/patrick-garrett/Data/PXD003881/`. Verify the complete file inventory and
condition mapping before execution. The deposit specifies fold changes of
1, 1.5, 2, 2.5, and 3; do not reuse the old script's illustrative 1:2:3:4:5
argument as quantitative ground truth.

Use a prespecified four-file subset (two replicates from each extreme) for
implementation checks and runtime estimation. Record every parameter explored
there. IonStar has prior development history and must not be called an unseen
validation dataset. Report its full cohort after settings are frozen.

**Validation extension: the paper's existing 12 Astral files, PXD066701.**
Use exactly the existing file manifest, with no new selection by results.
The deposited MS/MS evidence already used in the manuscript offers an
orthogonal check of added identifications, including at low input. Freeze all
matrix settings before this extension and disclose that the cohort already
appears in the single-run benchmark. Analyse matched-input replicate groups
as well as pooling input levels: a large sample supplying an identification
must not be presented as a new independent low-input detection.

PXD033350 remains an alternative if short-gradient MS1-only dilution data are
essential. Its old `/mnt/labdata/pgarrett/PXD033350_DirectMS1Quant` path is not
available on this host; prefer the verified local IonStar data over making that
download a prerequisite. No new dataset search is necessary for the core plan.

Dataset sources:

- https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD003881
- https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD066701

## 3. Separate workflow comparison from method development

Generate each file's KOTH features and hills once with the frozen current
version. Use those exact features for both separate searches and alignment.

| Arm | Input and behaviour | Question |
| --- | --- | --- |
| A: Separate searches | One UNO search per KOTH feature table | Baseline cohort coverage and repeated search cost |
| B: Current matrix search | One UNO search over the KOTH intensity matrix | What the shipped joint-search route achieves |
| C: Confidence-aware matrix search | Same consensus matrix plus validated per-cell confidence/quality inputs | Whether additional cross-run evidence improves B |

Run A and B before adding C. Keep FASTA, digest, RT model, seed set, and hardware
constant across arms. Use seeds 0, 1, and 2 as a prespecified minimum; show every
seed, including failures. Alignment is reused across search seeds. Record the
upstream randomness separately and do not claim search repeats assess alignment
uncertainty.

For C, first add strict row/run identity validation when joining matrix and
sidecar files; never rely on an unchecked row position. Separate genuine
contributions from transferred cells and add robust, label-independent
summaries of per-cell confidence, isotope/coelution quality, and alignment
residuals. Do not supply sample conditions, expected fold changes, or external
MS/MS identities to the classifier. Apply the same feature construction to
target and decoy peptide candidates.

Prespecified ablations should isolate: confidence masking of cells, added
quality summaries, and existing cross-run classifier features. For scoring
ablations, hold consensus rows and candidate sets fixed where possible. Test
profile reassignment separately, after identification, against its default-off
baseline. Report any simultaneous changes to filtering or candidate membership
instead of attributing their effect solely to a new score feature.

Use the shipped contributing-run floor of 2 as the initial B setting, and
prespecify sensitivity at 1 and 3; the historical value 6 is not a default or a
proven requirement for the new KOTH. Count true contributing runs for this
filter. Freeze the primary setting after the pilot, report the complete
sensitivity analysis, and never select it by best validation entrapment depth.

## 4. Measure identification, sample evidence, and computation separately

**Identification:** compare one pooled search with the deduplicated union of
separate searches, plus per-run coverage and intersection. Use a known-absent
entrapment proteome in the database for every arm; on IonStar both human and
E. coli are present, so neither can serve as absent entrapment. Freeze the
absent-proteome choice and treatment of shared peptides before searching.
Report empirical entrapment FDP for the final cohort lists and confidence
curves. A union of individually filtered 1% lists is not automatically a 1%
cohort list. Harmonize protein accessions/group membership for comparisons,
because separate searches may infer different indistinguishable groups.

**Per-sample support:** distinguish a cohort identification from evidence of
that precursor in each run. Report original detections, confident extractions,
MBR-only cells, missingness, and external MS/MS agreement for additional IDs.
An absent MS/MS match is lack of corroboration, not proof of a false MS1 ID.
KOTH q-value-filtered detections may still include MBR and must be labeled so.

**Quantification:** compare observed versus expected E. coli ratios, ratio
error/compression, within-condition CV, and human background stability.
Report coverage alongside accuracy so aggressive filtering cannot appear to
improve accuracy without revealing the lost observations. Use no floor
imputation in the primary assessment. Fix normalization and protein rollup
across quantification comparisons and report results on both common and
arm-specific quantified sets. Native apex intensities from single-run TSVs
and integrated LFQ intensities are different estimators; they cannot isolate
identification effects. For that comparison, map A's assignments onto the same
KOTH LFQ values used by B/C under one fixed assignment/rollup rule.

**Resources:** time feature finding, alignment/LFQ, and UNO searching
separately. Report search-only savings and full end-to-end wall time, CPU,
peak memory, and output size. Also report the total cost of A's separate
searches and their elapsed time under a stated equal-resource concurrency
policy. One search can be cheaper even if alignment makes total processing
slower; measure both rather than promising a speedup.

For scaling, preselect balanced IonStar cohorts of 5, 10, and 20 runs (one,
two, and four replicates per condition), using manifest order. Rebuild alignment
for each subset; slicing the 20-run matrix would leak the larger cohort's
consensus into smaller-cohort results. Compare each pooled search to A on its
own exact file set. Replicate subsets are related measurements, not independent
datasets.

## 5. Add only the new evidence to the paper

Add an SI section, **Joint identification from aligned KOTH feature matrices**,
with the frozen protocol, file/version manifest, route-equivalence check,
current-versus-enriched-input distinction, and complete seed/ablation results.

One figure can show cohort identification coverage and calibration, search-only
versus total resource cost, and missingness/ratio accuracy. A companion table
records file counts, true contributions, confident/MBR cells, and settings.
Add a short Methods paragraph describing matrix input and a Results paragraph
linking to the SI. Update the current statement that every run is searched
independently so it explicitly describes the single-run benchmark.

Do not restore the old quantitation section or its conclusions. Negative or
mixed matrix results remain useful evidence and must be retained. Numerical
claims enter the manuscript only through newly generated, provenance-tracked
tables, figures, and statistics. If scoring changes affect the default
single-run path, rerun the affected main-paper analyses as well.

Execution order: freeze manifest/builds, run the small pilot and equivalence
check, finish A/B, implement and test C only as a distinct extension, then
run the frozen validation cohort and generate the SI. Estimate full runtime
from the measured pilot before dispatching the complete analysis.
