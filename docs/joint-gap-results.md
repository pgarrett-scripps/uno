# Matrix gap investigation

On the preselected Astral 100ng triplicate, the joint matrix search identified
1,723 present-species accessions versus 1,352 in the freshly rerun individual-search
union (27.4% higher), at empirical entrapment FDP <=1%. The matrix without
cell scores already identified 1,718; the new cell scores added only 5 here.
All three fresh individual protein tables matched the frozen outputs exactly.
This is one cross-dataset check, not evidence of a universal matrix advantage.

## IonStar: what accounts for the gap?

The 1,533-accession matrix result shared 1,362 accessions with the 1,820 union,
gained 171 and missed 458. Of those 458, 427 received a matrix protein score but
fell below the empirical cutoff. At least two retained peptide-associated
candidates entered protein inference for 443 of them, so simple loss of all
input evidence is uncommon. Partial loss still matters: missing proteins retained
a median 70% of their independent quant-associated sequences, versus 86% for
shared proteins. These quant assignments are diagnostic associations, not
validated peptide IDs or a causal partition of the missing proteins.

## Fixed ablations and targeted tests

| IonStar matrix | No cell scores | MBR counts only | Q-derived scores only | All cell scores |
| --- | ---: | ---: | ---: | ---: |
| Standard | 1,458 | 1,304 | 1,268 | 1,309 |
| Expanded | 1,446 | 1,404 | 1,530 | 1,533 |

Q-derived scores provide nearly all the expanded-matrix gain in this pilot.
The ablation keeps confident counts with q-derived scores, since computing them
requires q-values; the count-only arm uses MBR-filled count and fraction.
Allowing coherent alternative rows during protein competition reduced depth to
1,372. Lowering KOTH's member-quality floor from 0.5 to 0.3 increased groups from
182,889 to 242,623 but reduced depth to 1,096. Both changes remain unpromoted;
more candidate evidence did not translate into more reliable identifications.

## Validation and interpretation

All counts use whole tied-score cutoffs selected for empirical entrapment
evaluation, not a claim that internal q-values were recalibrated. Astral used
one current executable, matching FASTA/model/settings, and fresh individual
controls. An independent integer calculation reproduced the matrix and union
counts. The Astral files were part of the earlier manuscript benchmark; they
were selected for this follow-up before examining new scores, but are not a
pristine held-out dataset. Runtime observations are diagnostic because some
jobs overlapped; no new speed claim is made.

All 65 core and 7 CLI tests passed. The diagnostic control reproduced protein,
LFQ and abundance outputs byte for byte, and all 14 job output manifests were
verified. The selected-calibration diagnostic export and optional ablation
switches are in the isolated worktree. No paper or production defaults changed.
The matrix searches consume only KOTH bundles; independent-search outputs were
used solely as comparators and to trace losses.

The compact record is `benchmark/joint_gap_summary.json`. Reproduction drivers
are `joint_gap_pilot.py`, `joint_gap_audit.py`, `joint_gap_stage_audit.py`,
`joint_gap_alternatives.py`, `joint_gap_retention.py`, `joint_gap_validation.py`,
and `joint_gap_astral_controls.py` under `benchmark/`. Full receipts are under
`benchmark/results/lfq_matrix_fixes/joint_gap*/`.
