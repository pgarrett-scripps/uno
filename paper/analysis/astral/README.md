# External Orbitrap Astral validation

Pilot data live outside the repository at
`/home/patrick-garrett/Data/PXD066701_Astral`.

## Preregistered dataset

- Deposit: PXD066701, Ryan Kelly laboratory, Brigham Young University.
- Sample: 10 ng Pierce HeLa digest.
- Instrument: Orbitrap Astral with FAIMS Pro at one fixed compensation voltage
  (-45 V).
- LC: dual-column short-cycle method (11.4 min total cycle, near-100% duty
  cycle), as described in the deposited project metadata.
- Selected complete matched set:
  - `CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng.raw`
  - `CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng.raw`
  - `CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng.raw`

Selection was made from file names and acquisition metadata before viewing any
uno or ms1searchpy result. Files will not be replaced based on search outcome.

## Pilot gate

The first file is processed through ThermoRawFileParser 1.4.5 and biosaur2
0.2.28. Production searches proceed only if conversion finishes, the feature
table contains the columns required by both engines, and both engines complete
RT calibration and protein inference. Failure is reported as an instrument or
feature-finder compatibility result rather than silently selecting another run.

## Comparison

Both engines receive the identical biosaur2 feature table, human-plus-yeast
target FASTA, and digest. uno uses the paper's frozen RT model. ms1searchpy uses
DeepLC and the same local quality-column patch as the main benchmark. All timed
runs are serialized and captured with `/usr/bin/time -v`. The common
`benchmark/entrapment.py` scorer reports human depth at 1% and 5% entrapment q.

The working tree is currently dirty, so this pilot records both the repository
HEAD and binary hashes. No Astral values enter the manuscript until the unified
benchmark is rerun from a committed, pinned build.

The completed production run and focused comparator controls are assessed in
`AUDIT.md`. `run_audit_controls.sh` tests stock-equivalent feature columns with
both DeepLC and ms1searchpy's default RT model; these are diagnostic controls,
not additional benchmark replicates.

## Data provenance

Official project and file API responses are preserved as
`project_metadata.json` and `file_metadata.json` in the external data root.
PRIDE reports no source checksum, so completed downloads are verified against
the deposited byte size and assigned a local SHA-256 checksum.

## RT-calibration controls (2026-09-03)

`RT_CALIBRATION.md` explains the production ms1searchpy collapse on these
files (contaminated DeepLC RT self-calibration) and bounds how much of the
gap is fixable on the comparator's side, using shipped options, an external
DIA-NN MS2 calibration list (`-pl`, `build_diann_calib_lists.py`), and a local
robust-calibration patch (`ms1searchpy_rt_calibration.patch`,
`-robust_calib 1`). `run_rt_calibration_controls.sh` reproduces the controls.
They are diagnostics, not benchmark replicates, and do not alter the
production numbers.

## Preregistered extension: the CW0710 input series (2026-09-03)

Declared before any of these files was searched. The deposit's CW0710
dual-column, fixed-CV (−45 V) series contains 36 Astral raw files over nine
HeLa inputs (200 pg, 500 pg, 1, 2, 5, 10, 20, 50, 100 ng). The three 10 ng
files above are the preregistered validation set; the other 33 files, listed
with deposited byte sizes in `input_series_files.tsv`, are taken in full as an
extension. Purpose: test whether the stock ms1searchpy calibration failure and
the uno lead depend on sample input, and characterise both engines across
input. No file will be dropped or replaced on the basis of search outcome; a
conversion or feature-detection failure is reported as such.

Protocol per file (`run_input_series.sh`, one job at a time, `nice`):
ThermoRawFileParser 1.4.5 → biosaur2 0.2.28 → uno ×1 (paper config, frozen
RT model) → stock ms1searchpy ×3 (`-deeplc 1 -ad 1`) → ms1searchpy `-pl` ×3
with the DIA-NN list for that file (`build_diann_calib_lists.py`). All timed
with `/usr/bin/time -v`. The 10 ng level keeps its existing ten stock runs and
three `-pl` runs. Same FASTA, digest, entrapment scorer, and binaries as the
10 ng set.

**Amendment 2026-09-04 (before any half-decade level was searched).** The full
33-file series was cut for compute time to the decade levels only: 200 pg, 1 ng,
10 ng, 100 ng, using the first three replicates by file name at each level
(1 ng keeps all five, which had already completed). The 500 pg, 2 ng, 5 ng,
20 ng, and 50 ng files were downloaded and converted but not searched. No
result from the dropped levels was seen before this decision; 200 pg rep2 and
rep3 were still queued when it was made and were run to completion.

**Amendment 2026-09-04 (after results were seen).** The 1 ng level is reported
with its first three files by name (rep1, rep2, rep4), the same rule as the
other levels, so the series is 12 files, three per level. 1 ng rep5 and rep6
had already been searched by every engine and arm; their rows are kept in
`results/input_series_dropped.tsv` (uno 567 / 523 at 1%, 1098 / 906 at 5%;
frozen uno 567 / 0; stock ms1searchpy 27–39; `-pl` 84–197) and they are
excluded from the figures, stats, SI table, and the raw-mode arm. This cut was
made after those results were seen, for symmetry only; the paper's Methods
states the removed values.

**uno fix (2026-09-04).** The frozen binary returned no proteins on 200 pg rep1
and 1 ng rep6. The retry described in `TODO.md` (stage loop re-seeded from the
mass-bootstrap proteins, RT sigma floored at 1.0 min, deeper result kept) was
applied and uno rerun on all 14 files into `results/uno_fixed/`; outputs are
byte-identical to the frozen binary on the other 12. The paper reports the
fixed binary. Summary: `summarize_input_series.py` →
`results/input_series_summary.tsv`.

**uno candidate selection (2026-09-04, shipped).** After the collapse fix, the
remaining low-input instability (1 ng rep4 at 187 with replicates at 567 to
975) was traced to every stage choosing by the count at 1% internal FDR. The
shipped default runs the pipeline at two RT sigmas and keeps the final list
with the larger soft depth (details in `TODO.md` and the `score.rs` comment).
All Astral uno outputs (`results/uno_fixed/`, `results/uno_raw_fixed/`) were
rerun with this binary; the previous fixed-binary outputs are in
`results/*_prev/`.

**Pinned publication re-run (2026-09-11).** Every uno search behind the
manuscript — both feature-source arms, seeds 0, 1 and 2 on all 12 reported
files — was re-executed on one pinned build (uno commit `0710225`, binary
SHA-256 `c376954d`), which adds seeded random streams, the evidence-retention
retry with discovery protection, and bounded RT recovery on top of the
2026-09-04 candidate selection. Outputs replace `results/uno_fixed/` and
`results/uno_raw_fixed/` in place; the superseded tables, their hashes and their
depths are archived under `benchmark/results/archive_574df1fd/`, and
`benchmark/results/publication_rerun/cases.tsv` records per search whether the
protein table changed and what its depth was before. Seed-zero searches were
serialized, so `results/input_series_summary.tsv` keeps reporting real
wall-clock and peak-RSS values for the uno rows.
