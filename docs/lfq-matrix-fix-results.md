# Initial LFQ matrix implementation and pilot results

This records the initial grouping/representative pilot. The subsequent
[joint matrix search](joint-matrix-search-results.md) uses all original member
observations and includes the replicated-group retention experiment.

The long-format and batch routes are implemented in the paired
`codex/lfq-matrix-fixes` worktrees. The five-run batch reproduces every independent
protein table byte for byte and therefore preserves the full union. The corrected
representative-based consensus variants in this initial pilot did not improve
coverage at matched empirical error.

## Implemented behavior

KOTH forms bounded groups before applying final size/seed filters, chooses each
run's primary feature by its own quality, and retains alternative observations.
Grouping tolerances are independent of extraction windows; the alignment
reference can be fixed explicitly. Tests cover the four reported failure
mechanisms, missing IM, distinct-run counts and reordered inputs.

The versioned long bundle carries quantitative cells separately from original
feature observations. It includes all original features, even those excluded
from consensus identification. Original isotope counts come directly from the
source tables because KOTH's alignment reader reconstructs each feature as a
single-hill placeholder. Consumers validate hashes, cell identities, membership,
seed/primary links and numeric inputs. Wide input remains supported.

UNO restores original seed quality in calibration and classification. A separate
recurrence ablation changes only two classifier columns, preserving rows,
quantities and tie-breaking. An optional representative experiment selects the
best actual member with at least three isotopes and uses that member's aligned
coordinates and its own quality together, falling back to the seed.

The batch route loads one FASTA database and each model once, retains predictions
in memory between runs, and searches original observations with per-run
calibration. It produces independent outputs and an explicitly labeled union.
It does not assign a calibrated cohort q-value to that union.

## What was verified

- All 1,609,186 source observations from the five-run IonStar cohort round-trip
  with exact equality of every descriptor consumed by UNO.
- With original quality disabled and legacy numerical rounding requested,
  all 80,720 retained long-input features exactly match the corrected wide
  export. The old matrix retained 60,960 features at the same two-run floor.
- All five batch protein tables match the frozen independent-search outputs
  byte for byte: 2,099 present-species accessions and 15 yeast entrapment
  accessions in their union. Its estimated entrapment FDP is 2.92%, not 1%.
- Cache tests establish bit-identical predictions and no writes on a fully warm
  lookup. Later batch runs needed only 12,270, 4,275, 2,839 and 1,738 new
  predictions out of approximately 1.3 million requests per run.
- 213 tests pass: 159 KOTH, 48 UNO core and six CLI tests. Dataset-dependent unit fixtures remain optional;
  the round-trip and protein-parity statements above use the actual complete
  five-run study inputs and outputs.

## Consensus identification results

All cases use the same five-run cohort, fixed reference, seed zero, FASTA and RT
model. B–F use corrected group membership. B–D and F use the same legacy numeric
rounding to isolate the intended differences. The first two count columns expand
protein groups into accessions; the final column counts protein groups.

| Arm | Present accessions at internal q <= 1% | Estimated entrapment FDP | Present protein groups at empirical 1% entrapment q |
| --- | ---: | ---: | ---: |
| A: frozen matrix baseline | 1,170 | 0.35% | 1,323 |
| B: corrected grouping, original quality off | 1,038 | 0.79% | 1,133 |
| C: corrected grouping plus original seed quality | 1,105 | 0.74% | 1,107 |
| D: C with recurrence columns disabled | 1,116 | 0.73% | 1,241 |
| F: C with an isotope-supported representative | 1,265 | 0.97% | 1,288 |

F was one targeted follow-up after the planned initial ablations. The original
observation audit found 7,843 retained groups with a two-isotope seed and another
member having at least three isotopes. Selecting a real member that passes the
existing anchor isotope criterion improved internal-q depth, but did not exceed
the old baseline at empirical 1% entrapment q. The original seed remains the
default; this alternative is an explicit diagnostic flag.

These are development results on one cohort and seed. They establish the
implementation properties and identify limitations; they do not establish a
general improvement in sensitivity. The apparent gain from 1,170 to 1,265 at
internal q <= 1% accompanies a higher empirical error estimate.

## Decision and remaining research

The pilot did not pass the sensitivity gate for promoting a corrected joint
search into the paper, so the larger-cohort/repeated-seed validation and paper
updates were not triggered. The original paper results remain the baseline.
The corrected code is isolated for review; original working checkouts are
unchanged.

The working batch route supplies the union-preserving option now. A joint model
that retains multiple original observation hypotheses, and consistency statistics
computed from coherent matches rather than all mass-compatible candidates, remain
research steps. They should be evaluated separately from additional threshold
tuning. Quantitative-ratio/CV validation of the new grouping on the full 20-run
cohort also remains pending; this five-run cohort has only one run per condition.

Runtime records are diagnostic: a release build overlapped part of the batch
execution, the old reference timings used cold predictions, and batch output
includes an additional peptide table. No controlled wall-time speedup is claimed.

## Reproduction and artifacts

See [input documentation](matrix-input.md) for the CLI and
[implementation plan](lfq-matrix-fix-plan.md) for the design.

The driver is `benchmark/lfq_matrix_fixes.py`; its actions are `align5`, `audit5`,
`pilot`, `followup`, and `report`. It records commands, binary/config/input hashes,
exit codes, output hashes and resource logs. Completed job signatures are checked
before reuse. The initial executable is archived as `bin/uno-pilot-v1`; the
representative follow-up has a separate binary hash. Use the recorded executable
for an existing job rather than overwriting its signature with a later build.

All new study outputs are under `benchmark/results/lfq_matrix_fixes/` in this
worktree, including `pilot-report.json`, `representative-audit.json`,
`pilot-source-manifest.json`, the long bundle, per-arm output tables and logs.
The older frozen study under the original checkout was read only.

A compact [reviewable result bundle](../benchmark/lfq_matrix_fixes_pilot.json)
records the counts, parity checks, cache requests and final binary hashes. The
final executable reproduced the representative experiment's protein, precursor
and abundance output tables byte for byte after the input-validation changes.

Batch quantitative outputs currently retain each original run's apex intensities.
Mapping the independently identified union back onto integrated aligned LFQ cells
is a separate pending step; the consensus route already uses those integrated cells.
