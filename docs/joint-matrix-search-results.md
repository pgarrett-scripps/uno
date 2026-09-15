# Joint KOTH matrix search: implementation and pilot

Follow-up: [cell-quality scores](joint-cell-quality-results.md) improved the expanded
matrix to 1,533 accessions at empirical <=1% error; the original pilot below
remains the baseline comparison.

UNO now performs one joint search from the long KOTH matrix, scoring original
observations rather than merging independent-search results. On the five-run
IonStar pilot it recovered 1,458 present-species accessions at an empirical
1% error threshold, 8.9% above the frozen matrix baseline. It remains below
the 1,820-accession separate-search union
at the same error criterion.

| Search | Present accessions at empirical 1% |
| --- | ---: |
| Frozen original matrix | 1,339 |
| Previous isotope-supported representative | 1,303 |
| Joint original-observation scoring | 1,458 |
| Joint scoring, agreement descriptors disabled | 1,051 |
| Expanded KOTH retention, guarded calibration | 1,446 |
| Separate-search union, threshold adjusted | 1,820 |

Counts expand same-species protein groups into accessions. Whole score ties are
kept together, and thresholds are selected using entrapment for development
evaluation. These are one-cohort, one-seed results with only a few entrapment
accessions at the cutoff; they do not establish general sensitivity or calibrated
cohort error on independent data. Several targeted development variants were
examined and remain recorded in the JSON summary.

## Implemented changes

`--matrix-long BUNDLE --matrix-joint` searches the union of original member mass
windows once per aligned group, shares digestion and RT prediction, and estimates
one cohort protein result. Candidate-specific mass/RT consistency and original
isotope/shape descriptors feed the classifier. Alternatives contribute at most
one observation per source run; missing runs and extraction-only MBR cells do not
create identification evidence. Joint mode retains KOTH singleton groups by
default, and recalibration uses actual compatible original observations.

KOTH's optional `allow_replicated_weak_seeds = true` retains moderate-quality
groups supported by at least two runs under the existing coordinate bounds.
It increased searchable original observations from 524,060 to 853,370 and groups
from 99,325 to 182,889. The audit preserved every original measurement, all prior
group memberships, and all 496,625 prior integrated LFQ cells exactly. Newly
retained groups still require quantitative accuracy validation.

Expanded retention raised the initial random-match background from 0.325 to
0.442. Excluding groups below KOTH's original seed floor from initial calibration
restored 0.325 while keeping them available for identification. This did not
improve 1% depth over the smaller joint matrix, so expanded retention remains
opt-in. Its higher depth at looser error thresholds is recorded separately.

## Runtime and validation

The joint search took 57.1 seconds versus
390.4 seconds for five independent searches
(6.8 times faster). Both used eight threads, warm
prediction caches, the same FASTA/model/config and protein/LFQ/abundance outputs,
with serial execution. These are single timing measurements with different
identification coverage; they exclude feature finding and matrix alignment.
The expanded matrix separately required
45.4
seconds for alignment/LFQ construction and 76.9
seconds for search.

All 226 tests pass (160 KOTH, 59 UNO
core, 7 CLI). Golden comparisons preserve the previous
representative outputs, the unaffected strong-group joint outputs and all five
independent protein tables byte for byte. The new routes remain experimental in
the paired worktrees; no paper claims or default search behavior were promoted.

Reproduction is in `benchmark/joint_matrix_pilot.py`; full per-run outputs,
commands, hashes and timings are under `benchmark/results/lfq_matrix_fixes/`.
The compact result is `benchmark/joint_matrix_summary.json`. The remaining gap is
joint identification sensitivity at matched error, rather than preservation of
original observations or shared prediction caching.
