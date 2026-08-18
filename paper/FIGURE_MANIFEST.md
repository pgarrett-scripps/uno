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
   delegates to `analysis/` and then stamps `.assets-stamp`. That stamp is what
   `just check` compares, so a regenerated figure must be committed together with
   it: `just assets && git add figures .assets-stamp`. Use `just assets-all`
   after re-running a benchmark, since the data these scripts read lives outside
   the repo and its mtime is never consulted.

## Output format

Single combined PDF: `paper/paper.pdf` = main text, then Supporting Information
appended with S-numbered figures/tables.

"There is no separate supplementary.pdf" holds for the **draft build only**. The
SI source already lives in its own file (`si-body.typ`, `#include`d by
`paper.typ`); splitting the *output* is a submission-bundle task, blocked on the
11 main-text cross-references that point at SI labels. See `SUBMISSION_TODO.md`
item B4 and decision D12.

## Main-text figures

The 2026-08 restructure (main-text surgery) renumbered these: the old
`fig3_determinism` moved to the SI as `figS6_determinism` (minus its spread%
panel, which became fig2 panel c), old `fig4_entrapment` -> `fig3_entrapment`,
and old `fig5_quant` -> `fig4_quant` (minus its volcano panel, which became
figS5 panel c).

| file | title | panels | data source |
|---|---|---|---|
| `fig1_pipeline.png` | uno pipeline & scoring cascade | schematic only | none (drawn) |
| `fig2_headtohead.png` | uno vs ms1searchpy across six DirectMS1 datasets | (a) depth @1% FDR, (b) depth @5%, (c) ms1 replicate range as % of median (uno 0% baseline), (d) wall-clock s, (e) peak RSS GiB; uno bar + ms1 three replicates as points with median | `benchmark/results/summary.tsv` |
| `fig3_entrapment.png` | protein FDR calibration under two-proteome entrapment | (a) reported q vs entrapment-corrected FDP (uno vs ms1), (b) accepted human proteins vs entrapment FDP | `benchmark/results/uno_*.tsv`, `benchmark/results/ms1_*/**_proteins_full.tsv`, scoring logic from `benchmark/entrapment.py` |
| `fig4_quant.png` | spike-in dilution response (PXD033350) | (a) UPS1 protein abundance vs spike fmol, (b) Spearman rho distribution by species (UPS1 / E. coli / human probe) | `/mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc/align*/intensity_matrix.tsv` → `uno --matrix` → `benchmark/quant/dilution_track.py` |

## TOC / abstract graphic (SUBMISSION_TODO B7)

| file | contents | data source |
|---|---|---|
| `fig_toc_graphic.png` | ACS "For Table of Contents Only" graphic, exactly 3.25 in × 1.75 in at 300 dpi (cropping disabled). Conceptual, not a results plot: an MS1 feature map (seeded random scatter, "no MS/MS") matched by the engine against an in-silico peptide database on mass and predicted RT, out to an FDR-controlled protein list. No numbers shown. | none (drawn) |

For this draft it is embedded in `paper.typ` directly under the abstract, so it
is visible on page 1; at submission ACS takes it as a separate Paragon Plus
upload (and would place it on a "For Table of Contents Only" page).

The PNG currently in `figures/` is a drawn placeholder from
`fig_toc_graphic.py`. The **working source is `paper/toc_graphic.pptx`**, an
editable PowerPoint scaffolded by `analysis/scripts/gen_toc_pptx.py` (one-shot;
rerunning it overwrites hand edits). The finished slide gets exported over
`figures/fig_toc_graphic.png`, after which `fig_toc_graphic.py` should be
deleted so `just assets-all` cannot regenerate the placeholder over it. This is
a deliberate exception to rule 3 (script-generated figures): the graphic is
schematic, carries no numbers, and is hand-finished.

## SI figures (S-numbered, appended to paper.pdf)

Filenames match the S-number the figure receives from its order of appearance
in `si-body.typ`; the reproducibility section (S8) deliberately sits after the
quantitation detail so that `figS6_determinism.png` renders as Figure S6.

| file | title | data source |
|---|---|---|
| `figS1_fdp_curves.png` | per-dataset entrapment FDP curves, all six datasets | as fig3 |
| `figS2_faims.png` | FAIMS benchmark (PXD023977), uno vs ms1searchpy, 3 replicates | `benchmark/faims/faims_out/` run outputs, scored by `benchmark/faims/score_entrap.py`. (There is no `paper/figdata/faims.tsv` fallback; an earlier revision of this manifest named one that was never created.) |
| `figS3_resources.png` | wall-clock, user-CPU and peak RSS vs input feature count | `benchmark/results/summary.tsv`, `benchmark/data/*.features.tsv` line counts |
| `figS4_quality_columns.png` | biosaur2 quality-column ablation and per-column target/decoy AUC | `paper/figdata/quality_columns.tsv` transcribed from `benchmark/RESULTS.md` |
| `figS5_quant_detail.png` | per-species abundance trajectories, detection-loss trace, and the 50 vs 0.1 fmol volcano | as fig4, plus `benchmark/quant/peptide_fc_de.py` for the volcano |
| `figS6_determinism.png` | run-to-run reproducibility detail: ms1searchpy replicate spread vs uno's single deterministic value, + run10 close-up | `benchmark/results/summary.tsv`, run10 reps in `benchmark/results/ms1_rep{1,2,3}/`, `benchmark/RESULTS.md` |
| `fig_astral_validation.png` | independent Orbitrap Astral validation on three 10 ng HeLa runs | (a,b) entrapment depth at 1% and 5% FDR; uno once and all ten ms1searchpy repeats per raw file, (c) wall-clock versus peak RSS | `~/Data/PXD066701_Astral/results/astral_summary.tsv`, generated from primary protein tables by `analysis/astral/summarize.py` |

`figS6_rt_model.png` (predicted vs observed RT, calibration fit, forward vs
reverse RT shift) was listed here as optional and was **dropped**: reproducing it
needs a DeepLC install, and the reverse-sequence artifacts behind the third panel
did not survive (see `SUBMISSION_TODO.md`, decision D13). No script exists and
the manuscript does not cite it.

## Removed: the legacy PNGs

`uno_{mass_defect,mass_excess,td_candidates,td_winners,why_peaked}.png` predated
this manifest, had no generating script, sat at the repository root in `figures/`
rather than `paper/figures/`, were referenced by nothing, and were **deleted on
2026-07-31**. They are in git history if wanted (`git log -- figures/`), but a
figure with no script that reproduces it from on-disk data cannot be cited here
under rule 1, so they were dead weight.
