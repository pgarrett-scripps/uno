# Joint matrix cell-quality pilot

The [follow-up investigation](joint-gap-results.md) separates counts from q-derived
scores and traces the remaining losses. It also verifies a matrix advantage on
a separate Astral triplicate using fresh individual controls.

Adding MBR fill counts and extraction-q summaries improved the expanded KOTH
matrix from 1,446 to 1,533 present-species accessions at empirical entrapment
FDP <=1%. It reduced the standard joint matrix from 1,458 to 1,309, so the scores
remain opt-in. The best result is 5.1% above the previous best matrix result,
but still below the 1,820 separate-search union at the same error criterion.

| Matrix | Without cell scores | With cell scores |
| --- | ---: | ---: |
| Standard retained groups | 1,458 | 1,309 |
| Expanded replicated groups | 1,446 | 1,533 |

These are accession counts from one five-run IonStar cohort and seed zero,
selected using an empirical entrapment evaluation cutoff, not counts at newly
calibrated cohort q-values. At UNO's internal q<=0.01, the new standard and
expanded arms returned 1,290 and 1,441 present accessions, with empirical FDP
0.952% and 0.852%, respectively. The differing effects do not establish which
individual descriptor helps or hurts; count and q summaries were added together.

The twelve optional descriptors comprise MBR-filled count/fraction, confident
MBR count, cohort confident-cell fraction, and mean q, worst q, q coverage and
confident fraction separately for direct and MBR cells. Positive intensity
identifies filled cells; q<=0.01 defines confidence. Missing q is explicit,
duplicate source files cannot inflate support, and extraction-decoy cells are
excluded. Every peptide candidate for a row receives the same extraction
summary, including decoy peptide candidates. No individual-search output enters
the matrix search, and MBR cells do not create additional peptide observations.

Enable with `[scoring] matrix_joint_cell_quality = true` and `--matrix-joint`.
The comparison held inputs, candidate rules, seed, FASTA/model and eight-thread
execution fixed. Warm-cache search times were 56.9 seconds for baseline,
58.6 seconds with cell scores and 81.0 seconds for expanded plus cell scores.
The baseline protein, LFQ and abundance files reproduced archived outputs
byte for byte; all 63 core and 7 CLI tests passed. No KOTH code changed in this
follow-up, and no paper/default promotion was made.

KOTH's prior changes were deterministic grouping bounded in mass/RT/IM, final
retention after group formation, primary selection by each observation's own
quality, preservation of alternatives and original source measurements, explicit
alignment reference selection, a versioned long bundle, and optional retention
of weaker groups supported by at least two runs. The observation table preserves
original values; the cell table exports recomputed extraction measurements.
`lfq_hybrid_score` is the four-component geometric mean in Hybrid mode, whereas
plain `combined_score` in that table remains the original consensus seed score.

Reproduce with `python benchmark/joint_cell_quality_pilot.py`. Commands, binary,
source hashes, outputs and timing receipts are archived under
`benchmark/results/lfq_matrix_fixes/joint_cell_quality_v1/`; the compact result is
`benchmark/joint_cell_quality_summary.json`. Independent-cohort validation is
still required before promoting the expanded-plus-cell-score configuration.
