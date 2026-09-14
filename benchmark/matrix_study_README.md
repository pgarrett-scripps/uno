# Fresh KOTH matrix study

Protocol: [`docs/koth-matrix-benchmark-plan.md`](../docs/koth-matrix-benchmark-plan.md).
Outputs: `benchmark/results/matrix_20260912/` (large files, gitignored).
Historical quantification outputs are never read by these scripts.

The frozen source trees include all local path dependencies and the original
working-tree changes. `frozen.json`, `enhanced-source.json`, `binaries.json`,
`inputs.tsv`, and `protocol-decisions.json` identify the experiment. The source
hashes are authoritative: building nested source snapshots caused KOTH's git
version discovery to pick up the enclosing UNO repository name. This limitation
is explicitly recorded rather than interpreting that string as the KOTH commit.

## Execution

Use Python with `matrix_study_requirements.txt`. The recorded environment is in
`python-runtime.json`. Rust builds use the saved locks and release profiles.

```sh
python3 benchmark/matrix_study.py freeze
# Build source/uno with Cargo; copy its binary into bin/uno. Build the KOTH
# binaries from that same workspace and copy into bin/. The saved study already
# contains these binaries; do not replace them during a running experiment.
python3 benchmark/matrix_study.py pilot-extract
python3 benchmark/matrix_study.py align --cohort pilot
python3 benchmark/matrix_study.py pilot-search
python3 benchmark/matrix_evidence.py --cohort pilot
# bin/uno_evidence is built from source_enhanced/uno, with the optional schema.
python3 benchmark/matrix_study.py pilot-evidence
python3 benchmark/matrix_study.py full-baseline
python3 benchmark/matrix_study.py full-evidence
python3 benchmark/matrix_astral.py run
python3 benchmark/matrix_report.py
python3 benchmark/matrix_astral_report.py
```

Do not start duplicate drivers while a search is active. Completed jobs are
reused only when command/config/binary signatures and output checksums agree.
Failed attempts retain their logs before retrying. The twenty primary individual
searches and three standard matrix searches at seed 0 run serially with eight
Rayon threads and RT prediction caching disabled. Extra seeds, sensitivity arms,
and validation searches may overlap and are excluded from resource claims. Pilot
resources are exploratory: some jobs overlap and some inherited a read-only
shared RT cache. Cold-cache measurements do not represent a warmed default
cross-process prediction cache.

The input manifest fixes the four-file pilot, all twenty IonStar files, nested
balanced 5/10/20-run cohorts, and seeds 0/1/2. Astral uses the manuscript's exact
twelve-file manifest, freshly processed with the same fixed settings.

## Interpretation

Compare joint identification with the accession-deduplicated union of separate
searches. Score that union's entrapment FDP explicitly; a union of independently
filtered 1% lists is not automatically a 1% cohort list. Mixed-species groups
are counted separately and excluded from species-specific metrics.

`masked` zeroes extraction q > 0.01 cells while preserving consensus membership.
`quality` adds nine extraction-quality summaries. `quality_null` adds the same
nine columns with zero values, controlling for classifier dimension and column
sampling. `run_evidence_null` retains the baseline floor-2 rows, sets counts and
intensities to constants, and is searched with floor 1: it is an ID-only ablation
that can also change intensity-based tie-breaking, never a quantitative output.
Profile reassignment is evaluated independently after identification.

Quantification uses known E. coli ratios and constant human background, with
label-independent median-ratio normalization and no missing-value imputation.
Native single-run apex intensities are not compared against integrated LFQ
intensities as if they were the same estimator. Astral MS/MS corroboration uses
the deposited DIA-NN report only for evaluation, never model inputs or tuning.

## Verification

```sh
python3 -m unittest discover -s benchmark -p test_matrix_study.py
cargo test --release --offline --locked -p uno-core --lib
```

`route_audit.rs`, compiled as an example against the frozen baseline source,
compares all UNO-consumed feature fields at KOTH TSV precision. Its pilot found
matching feature membership but different row order. The saved TSV and `--raw`
protein lists are not assumed identical. The unchanged-matrix regression checks
that adding the optional schema leaves ordinary matrix search output unchanged.

After every prespecified job has completed, run `python3 benchmark/matrix_finalize.py`.
It rejects missing cases, failed commands, changed outputs/configs/binaries/sources,
changed FASTA or RT model, and broken regression controls. It then regenerates the
reports and writes the complete case inventory and compact paper evidence bundle.
`--status` prints the completed/expected case counts without publishing results.
The prespecified inventory is 90 IonStar searches (including evidence arms) and
60 Astral searches, in addition to the exploratory pilot and route audit.

The long background execution sessions ended with SIGTERM before completion.
`matrix_queue.py --batch-size 3 --workers 3` resumes the fixed case inventory in
bounded batches. Run one batch at a time, and repeat until it reports no pending
cases. It preserves incomplete attempts before retrying and skips successful
completion stamps; finalization verifies their output checksums. The source,
binaries, inputs, search seeds, and cold RT-cache setting remain unchanged.
The original background drivers must be stopped before using this recovery queue.

## Published core scope

After the execution-time review, the paper uses the complete prespecified seed-0
comparison (30 IonStar and 20 Astral searches). `matrix_finalize.py --core`
verifies this scope and labels the bundle `complete-core`; it never claims all
150 planned searches completed. The bundle includes every completed repeat and
the explicit pending inventory. Default seed 0 was fixed before searching, not
selected by outcome. This core analysis makes no search-seed robustness claim.

At this checkpoint, 84 searches are complete and 66 repeat-extension cases remain.
The queue is stopped. Full-scope finalization still requires all 150 cases and
uses all three seeds; paper text must be revised if that extension is completed.

The RT prediction cache already exists, but all core full-cohort searches used
`UNO_RT_CACHE=0`. Only the serialized IonStar baseline timings support resource
claims. The current cache reparses the entire model-specific TSV per process;
a database-specific cache or persistent multi-file process could reduce overhead.
Across the 20 serialized IonStar individual searches, the logged prediction
block totals 721.9 of 2563.5 seconds (28.2%). Eliminating that entire block is
an upper bound, not a measured warm-cache speedup.
