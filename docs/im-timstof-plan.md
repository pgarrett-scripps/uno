# timsTOF IM experiment

This work is isolated on `codex/im-timstof`, based on the frozen
`paper-2026-09-13` tree. It does not modify the submission manuscript.

The old `im_im2deep.gbdt` model declares 70 input features, while the current
encoder emits 99 RT descriptors plus charge (100). Its first tree expects
log-mass at index 61 and charge at index 69; the current encoder supplies
position-binned residue counts at those indices. The loader previously did
not validate the dimensions. The GBDT loader now rejects incompatible widths, and regression tests exercise
rejection and the charge position.
The missing historical encoder is not guessed or reconstructed.

A compatible Sequoia model uses the current 100 descriptors. It is trained on
the existing public Meier CCS dataset from ProteomicsML, not the HYE benchmark
identifications. CCS is transformed to a quantity proportional to inverse
mobility; each search estimates its multiplicative scale separately by charge.
The public validation split keeps complete sequences and all their charges
outside training. A separate HYE A1 acquisition calibrates relative error widths
only; its labels never train the model and it is excluded from search evaluation.

The paired experiment uses fresh KOTH 0.5.0 features from the 15-minute, 50 ng
Ultra 2 A2, B2 and C2 acquisitions, with seeds 0, 1 and 2. IM off is compared
with a primary three-sigma gate and a two-sigma sensitivity gate. All 27 searches
are retained. Search settings, model, features, and database are identical within
each pair; changing a seed changes decoys and scoring randomness together.

All three HYE species are present. Absent Arabidopsis supplies entrapment,
with its full proteome included in every search. Evaluation expands same-species
protein groups, excludes mixed groups, and includes entire tied-score blocks.
Depth at 1%, 3% and 5% empirical entrapment FDP is reported separately from
internally accepted q<=1% lists. These retrospective cutoffs are not deployable
q-values. Same-file Sage MS/MS corroboration is a secondary, incomplete check.

A consistent primary gain requires improved median 1% depth in all three files
and improvements in at least seven of nine paired searches. No threshold is
selected from the evaluation results. Interrupted incompatible-model and
unbounded-memory attempts are preserved outside the completed job inventory.
Prediction now deduplicates sequence/charge pairs and chunks feature matrices;
116,726 predictions, including repeated and reordered inputs, match the earlier
compatible implementation exactly.

Reproduction from this worktree uses `benchmark/im_timstof_study.py prepare`,
`benchmark/im_public_training.py`, and the Sequoia training command recorded in
`models/im_public_20260913.training.json`. Then run the study driver's `extract`, `calibrate`, and `run` actions,
followed by `benchmark/im_timstof_report.py` and `benchmark/im_timstof_check.py`.
The fitted model is retained as `models/im_public_20260913.sequoia`.
