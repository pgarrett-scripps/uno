# LFQ matrix fixes and long-format evidence

Status: grouping, long-format export and the original batch pilot are complete.
The subsequent one-matrix joint scoring implementation and KOTH replicated-group
retention experiment are recorded in [joint-search results](joint-matrix-search-results.md).
The joint pilot improves matrix sensitivity but remains below the separate-search
union at matched error; independent validation and paper promotion remain deferred.

## Objective and workspace

Improve KOTH's LFQ matrix and perform one UNO search using the original aligned
observations within that matrix, with shared FASTA digestion/RT predictions and
one cohort protein result. Recover more of the independent-search union at
comparable empirical error while preserving original measurements and integrated
LFQ cells. Independent searches serve as a comparator, not a source of rescued
matrix identifications.

Both repositories use branch `codex/lfq-matrix-fixes`:

- UNO: `/home/patrick-garrett/Repos/uno/experiments/lfq-matrix-fixes/uno`
- KOTH: `/home/patrick-garrett/Repos/uno/experiments/lfq-matrix-fixes/koth_rust`

The UNO worktree includes the source checkout's uncommitted tracked changes and
untracked study/audit files. KOTH starts from its clean current checkout.
`../worktree-baseline.json` records the starting commits and copied files;
`../uno-baseline.patch` records the inherited tracked changes. Original checkouts
were left unchanged. Sage and Sequoia resolve through sibling symlinks to their
existing source repositories; implementation work is confined to UNO and KOTH.

Existing frozen data remain at
`/home/patrick-garrett/Repos/uno/benchmark/results/matrix_20260912`.
New experiment outputs belong under this UNO worktree's
`benchmark/results/lfq_matrix_fixes/`, with a separate manifest and binary hashes.

## Design decision: long format is the primary evidence interface

Use an explicit versioned long-format bundle. A long table alone does not
prevent information loss: it must carry original feature measurements, not
only the LFQ grid measurements already present in `lfq_details.tsv`.

| Artifact | Row identity | Content |
| --- | --- | --- |
| `lfq_matrix.long.tsv` | `(consensus_id, run_id, extraction_kind)` | One target or decoy extraction cell; consensus coordinates, selected original observation ID, original feature quality, extraction intensity and quality, direct/MBR status |
| `feature_observations.tsv` | `(run_id, observation_id)` | Every original input feature, native and aligned coordinates, original descriptors, nullable consensus ID, primary/alternative status and any exclusion reason |
| `matrix_manifest.json` | One bundle | Schema version, run names/order, source hashes, reference run, alignment/configuration provenance, intensity estimator, units and table hashes |

`feature_observations.tsv` retains alternative and filtered original observations
because several candidates can map to one run/cell. Keeping only one long row
per cell would otherwise recreate the current loss. Quantitative intensity
occurs once per extraction cell; it is never summed repeatedly across observation
rows. Wide intensity/q-value matrices remain optional derived exports, and UNO
continues accepting the existing wide format.

Field contract:

- Stable IDs are scoped to the frozen input bundle; use input hashes and a
  deterministic feature index/canonical ordering, with explicit tie-breaking.
  Do not imply that consensus IDs remain identical across grouping algorithms.
- Preserve native mass/mz, RT, charge and optional IM; aligned mass/RT/IM occupy
  separate columns. Preserve original combined score, isotope count, scan count,
  cosine/isotope scores, ppm/neutron-offset descriptors, peak width and available
  isotope-fit errors, plus original apex/integrated intensity as separate fields.
- Prefix original measurements with `feature_` and extraction measurements with
  `lfq_`. In particular, `feature_n_isotopes` and `lfq_n_isotopes_found` are
  distinct. Original-feature fields are null for MBR-only cells and synthetic
  extraction decoys. Missing values are not zeros or fabricated quality scores.
- Keep `lfq_intensity`, extraction status, extraction score, extraction q-value,
  and observed/expected extraction coordinates. A positive intensity does not
  imply confident detection. Extraction q-values remain separate from UNO
  peptide/protein q-values; any paired target q-value attached to a decoy must
  have an explicit paired-value name, not masquerade as the decoy's own q-value.
- TSV floats round-trip without additional precision loss. Reject duplicate
  cell/observation keys, unknown runs, mismatched hashes, invalid coordinates,
  and broken observation links. Stream rows; avoid constructing a giant dense
  table merely to parse long input. Distinguish unattempted extraction from
  attempted extraction with zero signal.

## Implementation sequence

### 1. Correct consensus membership before tuning thresholds

KOTH files: `koth_ff/src/lfq/consensus.rs` and LFQ configuration.

Replace premature fragment filtering and unrestricted transitive merging with
deterministic bounded clustering. Process eligible observations in descending
quality order with stable ties; use mass-window candidate lookup and a fixed
representative. Require compatibility with the representative and bounds on the
whole group's mass/RT/IM span. When several groups qualify, choose by normalized
coordinate residual with stable ties. Retain every original member for inspection.

Apply minimum distinct-run count and minimum seed score only to completed groups.
Select each run's primary observation by that observation's actual score, with
documented ties; retain alternatives. Count distinct runs, never member rows or
positive extracted cells. A contradictory observation remains a separate
hypothesis rather than bridging incompatible groups.

Expose grouping mass/RT/IM limits independently of extraction windows. Initially
preserve the existing effective pairwise limits to isolate algorithm changes;
do not combine this correction with a broad tolerance sweep. Make the alignment
reference explicit in the experiment configuration.

Acceptance: convert the four audit reproductions into tests of corrected
behavior in KOTH. Add invariants for cluster bounds, distinct-run counts and
determinism under reordered run/feature input. Existing valid simple groups
remain stable. Instrument counts before member filtering, before/after grouping,
after group filtering and after UNO filtering.

### 2. Export and read the long-format bundle

KOTH files: consensus model, `lfq/mod.rs`, `align_main.rs`, output configuration.
UNO files: `uno-core/src/input.rs`, `uno/src/main.rs`, matrix input documentation.

Extend KOTH's existing long LFQ exporter into the bundle contract above. Populate
original descriptors by joining the retained observation IDs to the original
features. Keep complete observation provenance alongside selected per-cell
features. Reuse the existing extraction calculations; serialization changes
must not silently change quantification.

Add an explicit UNO long-input option (planned name `--matrix-long`) that reads
the bundle manifest and validates links before search. Separate observation data
from per-cell quantities internally so alternate observations cannot double-count
abundance. Continue to load legacy wide inputs with their documented missing
quality semantics.

Acceptance: fixed-consensus wide and long exports give identical intensities,
support counts and consensus coordinates; with original-quality use disabled,
the same reconstructed feature order yields identical search outputs. Validate
missing data, duplicates, inconsistent metadata, numeric round trips and multiple
observations per cell. Reconstruct each original single-file UNO input from the
observation table without losing any consumed descriptor.

### 3. Restore early quality evidence, then isolate consistency effects

UNO files: input representation, matching, `score.rs`.

First keep the consensus-coordinate search but restore representative original
quality in both mass-anchor selection and ML features. Use quality from an actual
observation; do not assemble an impossible feature from unrelated per-column
maxima. Start with the seed observation for a clearly attributable comparison.

Retain best-observation quality and absolute confident direct support as separate
signals. Use direct and transferred confidence separately. Do not use all positive
cells as a proxy for independent recurrence. Avoid a new hard requirement for
presence across all conditions or input amounts.

Test recurrence independently: disable recurrence-related classifier columns
while retaining exactly the same rows, quantities, quality, calibration and
tie-breaking. Test intensity effects separately if needed. Candidate mass/RT
consistency should summarize coherent matched observations rather than every
mass-compatible alternative; duplicated observations must not inflate distinct
peptide evidence or leak across training folds.

Acceptance: diagnostics show where original isotope quality enters calibration;
tests show MBR-only rows do not receive fabricated original quality. The
recurrence ablation changes only its intended classifier columns. Compare
calibration anchors, random-match background and target/decoy ranking as well
as final counts.

### 4. Add a batch route that retains per-run identification evidence

Use the same long bundle for one process with one digested database and a shared
in-memory sequence/model-keyed RT prediction map. Score original observations
with their native run coordinates and run-specific calibration. Alignment links
observations to quantitative cells after scoring; it does not restrict the search
to one seed's mass/RT. Preserve filtered/unassigned observations for this route.

Begin by reproducing independent per-run results in a single invocation with
shared computation. Report their accession union explicitly as a union of
per-run calls. It does not automatically have 1% cohort FDR. Any subsequent
joint protein evidence model must deduplicate peptide evidence, treat targets
and decoys symmetrically, and estimate cohort error rather than taking minimum
per-run q-values as calibrated cohort q-values.

Keep the disk RT cache enabled for development; point it at a writable,
study-specific directory. Load predictions once per batch and predict only
missing sequences. Time database loading, prediction, each calibration/scoring
stage, grouping and extraction separately. One command still contains per-run
scoring work; quantify the actual saving.

Acceptance: with fixed input ordering/seeds, each run's results match the
independent route; cached and uncached predictions/searches agree. Adding a run
cannot erase already available native observations. Cohort-level IDs and
per-sample confident detection remain distinct outputs.

### 5. Run a bounded comparison and gate paper changes on validation

Start with the existing five-run balanced IonStar cohort and its fixed reference,
same original feature/hill files, FASTA, RT model and seed zero. Reuse extraction
and baseline results after checking their manifests. No new dataset or feature
finder tuning is needed for the first decision.

| Arm | Purpose |
| --- | --- |
| A: unchanged grouping, long input, original-quality use off | Serialization/search parity with the existing wide baseline |
| B: corrected grouping, original-quality use off | Isolate grouping corrections |
| C: corrected grouping plus original quality | Measure restored early/ML evidence |
| D: C with recurrence columns disabled | Isolate the consistency hypothesis |
| E: original-observation batch route | Compare shared-computation search with independent runs |

This is five initial configurations, each with one seed. E performs five
run-specific searches inside one process; it is not one consensus calibration.
Run A–D first, then E when batch implementation passes parity checks. Add no
threshold grid or automatic three-seed expansion during development.

For every arm, record original features retained, independently assigned observations lost at
each stage, alternative coordinates retained, calibration anchors/background,
present/entrapment IDs and overlap, quantification bias/CV/completeness, time and
memory. Trace representative lost peptides from the independent calls through
the pipeline; label that as diagnostic evidence, not ground-truth identification.
Evaluate protein coverage at matched empirical error rates and show the curves;
the old 2,550-accession union at 5.78% FDP is not a 1%-FDR performance target.

After the implementation and settings are frozen, validate the selected corrected
route on the full 20-run IonStar and 12-run Astral cohorts, first at seed zero
(two configurations), then seeds one and two for that selected route only
(four additional configurations). Reuse existing independently searched reference
outputs only where matching seed/input/binary provenance supports the comparison;
report seed-zero union comparisons separately if repeated reference outputs are
missing. Do not restart the old unfinished 150-search queue as part of this plan.

The desired result is a smaller coverage gap without worse empirical error
control, while preserving useful LFQ ratios and reducing repeated computation.
No exact recovery percentage is assumed. If the gap remains, use the measured
retention/calibration diagnostics to select one next change. Integrate only the
validated final comparison into the paper and retain the old results as provenance.

## First implementation deliverable

Commit-sized KOTH correction: delayed final group filtering, bounded deterministic
grouping, actual per-run winner selection, preserved member provenance, and the
four corrected regression tests. Then build the long-format contract on that
representation. Implementation is authorized by the subsequent “ok go”. Initial validation uses
the existing five-run cohort before any expansion.
