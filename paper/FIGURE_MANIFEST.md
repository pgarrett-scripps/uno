# Figure & build manifest — uno JPR manuscript

Single source of truth for figure filenames, panel contents, and data provenance.
Agents MUST follow the names here exactly so the manuscript's `#figure` blocks
and the generated PNGs line up.

## Hard rules

1. **Never invent a number.** Every value in every figure must trace to a file on
   disk. Allowed sources, in order of preference:
   - a results file under `benchmark/` or `/mnt/labdata/pgarrett/...`
   - a fresh run of `uno` / `ms1searchpy` you performed and saved under
     `benchmark/results/` (record the command in the script header)
   - a value transcribed from a checked-in `README.md`/`RESULTS*.md` into
     `paper/figdata/<name>.tsv` **with a provenance header comment naming the
     source file and line**. Transcription is the last resort.
2. If a panel's data does not exist and cannot be produced, **drop the panel and
   report it**. Do not substitute a plausible-looking value.
3. Figures are written to `paper/figures/<name>.png` at 300 dpi (also emit
   `.pdf` alongside where practical) by scripts in `paper/analysis/scripts/`.
   The script name IS the output name: `fig2_headtohead.py` -> `fig2_headtohead.png`.
4. All scripts import `paper/analysis/scripts/style.py` for fonts/colors/sizing so
   the set reads as one system. Colorblind-safe palette, no red/green pairs,
   legible at column width (single column ≈ 3.4 in, full width ≈ 7.0 in).
5. Scripts must be rerunnable end to end via `just assets` from `paper/`, which
   delegates to `analysis/`. Every live output is addressed by a semantic id in
   `assets.json`; the manifest records its output hash, generator, and declared
   input hashes. `just verify` fails on stale, unclaimed, or missing assets.

## Output format

Single combined PDF: `paper/paper.pdf` = main text, then Supporting Information
appended with S-numbered figures/tables.

"There is no separate supplementary.pdf" holds for the **draft build only**. The
SI source already lives in its own file (`si-body.typ`, `#include`d by
`paper.typ`); splitting the *output* is a submission-bundle task, blocked on the
11 main-text cross-references that point at SI labels. See `SUBMISSION_TODO.md`
item B4 and decision D12.

## Main-text figures

The prose references semantic labels and asset ids, so journal figure numbers
follow order of appearance rather than filenames.

| file | title | panels | data source |
|---|---|---|---|
| `fig1_pipeline.png` | uno pipeline & scoring cascade | schematic only | none (drawn) |
| `fig2_headtohead.png` | uno vs ms1searchpy across six DirectMS1 datasets | (a) depth @1% entrapment FDR, both engines as individual draws, medians and full ranges with sample sizes; (b) original serialized wall-clock versus peak RSS. Depth @5% is in the seed summary table and SI spread figure | `benchmark/results/spread.tsv` (depth), `benchmark/results/summary.tsv` (serialized resources) |
| `fig_astral_validation.png` | independent Orbitrap Astral validation | (a) seeds/reruns grouped within acquisition files, with per-file medians and ranges; lines join medians of the three file medians at each input; (b) original serialized resources | `benchmark/results/variability_extension/summary.tsv` and `~/Data/PXD066701_Astral/results/input_series_summary.tsv` |
| `fig4_overlap.png` | protein-level agreement | (a) every within-uno, between-engine and within-ms1searchpy pairwise Jaccard; (b) per-seed recovery of the intersection of available comparator lists; points, medians, full ranges | all primary protein outputs inventoried in `benchmark/results/spread.tsv`, loaded by `seed_data.py`; writes `figdata/overlap_analysis.tsv` |
| `fig5_raw_mode.png` | built-in feature detection vs biosaur2 tables | depth @1% entrapment FDR, three matched seeds per file/arm on all 18 files, shown as points with medians and ranges | `benchmark/results/variability_extension/summary.tsv`; original resource measurements remain in `benchmark/raw_mode/summary.tsv` |

## TOC / abstract graphic (SUBMISSION_TODO B7)

| file | contents | data source |
|---|---|---|
| `fig_toc_graphic.png` | ACS "For Table of Contents Only" graphic, exactly 3.25 in × 1.75 in at 300 dpi (cropping disabled). Left column schematic: an MS1 feature map for one run (the slide's own hand-placed feature blobs, not a random scatter) feeding a binomial protein score with an illustrative 1% FDR cutoff, three feature colours marking the three scored targets. Right column is the one part with numbers: per-dataset median depth at 1% entrapment FDR for both engines, the same medians Figure 2a reports. | `benchmark/results/spread.tsv` (bars); left column drawn |

For this draft it is embedded in `paper.typ` directly under the abstract, so it
is visible on page 1; at submission ACS takes it as a separate Paragon Plus
upload (and would place it on a "For Table of Contents Only" page).

`paper/toc_graphic.pptx` is the **design source**: the editable slide the
graphic was laid out on, drawn at 2x the ACS box (6.5 in x 3.5 in) with type at
2x printed size. `analysis/scripts/gen_toc_pptx.py` scaffolded the first version
of it and is one-shot -- rerunning it overwrites the hand layout.

The committed PNG is **not** an export of that slide: `fig_toc_graphic.py`
redraws the slide's layout, colours and type sizes in matplotlib, working in
slide units and halving every length on the way out, so the two stay diffable.
That keeps the graphic under rule 3 rather than outside it, and it is why the
script is no longer slated for deletion. The reason for redrawing rather than
exporting is the bar chart: the slide carried uno's single default-seed depth
against ms1searchpy's 3-replicate median, which overstated uno on 09746 and
09766; the script reads `spread.tsv` so the TOC graphic cannot drift from the
figure it advertises. Edit the slide for layout, then port the change into
`fig_toc_graphic.py` and rerun it -- and keep the slide's embedded worksheet in
step if the depths change.

## SI figures (S-numbered, appended to paper.pdf)

Six SI figures cover calibration, search variability, the RT model, prediction
accuracy, throughput, and the fresh matrix-input core comparison. Primary depth distributions and per-draw calibration
counts are generated tables.

| file | title | data source |
|---|---|---|
| `figS1_fdp_curves.png` | per-draw calibration curves for both engines, with pointwise medians | all primary outputs inventoried in `benchmark/results/spread.tsv`; writes `figdata/entrapment_curves.tsv` |
| `figS6_determinism.png` | primary 1%/5% depth distributions and archived default-seed ranks | `benchmark/results/spread.tsv`; excludes development run10 |
| `figS2_rt_model.png` | fixed RT model architecture, training corpus and holdout accuracy | shipped model and Chronologer corpus; no repeated search draws |
| `figS3_rt_predictor.png` | uno-training-disjoint CNN versus DeepLC accuracy, identical affine calibration for both: (a-c) residual vs observed RT per set (median + 10-90% band), (d-f) cumulative absolute-error distribution with the 95th percentile marked | `benchmark/rt_predictor/results/summary.tsv` and `results/predictions_<set>.tsv` (per-peptide test predictions), both written by `benchmark/rt_predictor/run.py` from the checked-in sets; errors normalized by each set's 5-95% RT range |
| `figS4_rt_throughput.png` | RT-prediction throughput vs batch size on a synthetic tryptic library: (a) wall time per batch, (b) peptides/s, medians of 5 repeats after warm-up, min-max band | `benchmark/rt_predictor/results/throughput.tsv` written by `benchmark/rt_predictor/throughput.py` (uno via `rt_throughput` Rust bin, DeepLC `make_preds(calibrate=False)`) on `sets/synthetic_tryptic.tsv` from `prepare_library.py` |

Retired quantitation, quality-column, duplicate-FAIMS, historical RT, and resource
graphics remain registered as archival assets but are not referenced by the
Technical Note. The analysis build skips the five retired quantitation, FAIMS,
resource, and quality-column generators; their scripts remain available for
manual archival reproduction when the historical inputs are present. The
determinism-detail graphic is now live as the seed-spread
figure. Schematics and fixed-model prediction plots do not summarize search
randomness; throughput already shows repeated measurement ranges.

## Removed: the pre-manifest PNGs

`uno_{mass_defect,mass_excess,td_candidates,td_winners,why_peaked}.png` predated
this manifest, had no generating script, sat at the repository root in `figures/`
rather than `paper/figures/`, were referenced by nothing, and were **deleted on
2026-07-31**. They are in git history if wanted (`git log -- figures/`), but a
figure with no script that reproduces it from on-disk data cannot be cited here
under rule 1, so they were dead weight.

## Variability tables

`gen_variability_table.py` generates `si/variability_summary.typ` (median and
full range at 1% and 5%, with n) and `si/entrap_reported.typ` (every primary
draw). `si/results_per_dataset.typ` remains the original serialized resource
cohort, explicitly labeled as such in the manuscript.

## Fresh koth matrix extension

`figS7_koth_matrix.py` writes `figS7_koth_matrix.{png,pdf}` from the audited
`figdata/koth_matrix_study.json` bundle, generated by
`benchmark/matrix_release_report.py`. The complete KOTH 0.5.0 study is archived
in `benchmark/results/matrix_20260913_koth050`: 32 fresh feature tables, eight
matrices, 144 searches over three seeds, plus four dedicated serialized timing
searches. All cohorts and seeds are included.

Panels show (a) Astral input-level and pooled accession coverage, (b) nested
IonStar cohort coverage, both at matched empirical entrapment cutoffs, (c) the
Astral 100 ng serialized warm-cache search comparison, and (d) known-ratio
response on common IonStar E. coli groups with shared normalization per seed.
No concurrent job timing is used as a performance claim.

`gen_koth_matrix_table.py` generates eight SI tables for frozen program identities,
route definitions, Astral identification and MS/MS corroboration, IonStar
identification, added/lost accessions, upstream per-sample extraction evidence,
quantification, and dedicated search timings. Historical numerical matrix
results are replaced rather than mixed into this study. The original single-run
engine comparisons retain their existing provenance.
