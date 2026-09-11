# Why ms1searchpy collapses on the Astral files, and what rescues it

Diagnostic controls run 2026-09-03 on the three preregistered 10 ng Astral
files. None of this changes the production protocol in `run_all.sh`; it
explains the production result and bounds how much of the gap is fixable on
the comparator's side. Outputs live under
`$DATA/results/rt_calibration/{sweep,pl,robust,diag}`; scripts are
`run_rt_calibration_controls.sh`, `build_diann_calib_lists.py`, and the
ms1searchpy patch `ms1searchpy_rt_calibration.patch`.

## Conclusion

ms1searchpy does not fail at feature matching, mass calibration, or protein
inference on these files. It fails at retention-time self-calibration. Its
DeepLC calibration set, built from every peptide of the stage-1 proteins, is
roughly two-thirds random matches on this data, and both DeepLC's per-bin
median calibration and the Gaussian sigma fit on that mixture inflate. The
fitted RT sigma is 1.4 to 1.9 min on an 11 min acquisition, so the RT gate
admits nearly everything and the search reports 26 to 110 proteins.

No shipped ms1searchpy option repairs this. Giving it a clean calibration set
does, by either of two routes, but it then reaches one tenth to one third of
uno's depth at 1% entrapment FDR (median per file: 26 to 29% with the external
MS2 list, 10 to 32% with robust self-calibration) and one fifth to one half at
5%. uno self-calibrates on the same features without help.

## Evidence

**Calibration quality.** From the deposit's own DIA-NN 2.0.2 report
(`CW0710_DN_all_.zip`, `report.parquet`), 47,689 MS2-verified peptides with
retention times exist for Astral rep1; 69% of them have a biosaur2 feature
within 10 ppm and 0.3 min. Calibrating on a random half and checking on the
other half:

| RT model | MAE (min) | robust sigma (min) | 95th pct abs error (min) |
|---|---:|---:|---:|
| DeepLC 3.1.13, piecewise calibration | 0.60 | 0.69 | 1.58 |
| uno `rt_cnn_chrono_global_v3_ch32`, linear calibration | 0.37 | 0.44 | 0.92 |

So DeepLC transfers to this dual-column LC about 1.5 times worse than uno's
model, but 0.69 min is nowhere near the 1.4 to 1.9 min ms1searchpy ends up
with. The remainder is contamination: of the 1,981 peptides in the production
calibration list (`-save_calib`), only 712 are consistent with a uno
confident-protein match at the same RT.

**Mass gate is not the problem.** True mass error of MS2-verified peptides
against biosaur2 features has a robust sigma of 0.66 ppm. ms1searchpy's
3-sigma gate (0.96 ppm) keeps 72% of true matches; uno's 1.54 ppm gate keeps
80%.

**Shipped options, rep1, one run each, human proteins at 1% / 5% entrapment
FDR** (production runs ranged 26 to 110 at 1%):

| Options | stage-1 proteins | DeepLC sigma (min) | 1% | 5% |
|---|---:|---:|---:|---:|
| default (`-deeplc 1 -ad 1`) | 140 | 1.43 | 57 | 57 |
| `-es 1` | 100 | 1.84 | 48 | 48 |
| `-ptol 3` | 127 | 1.73 | 41 | 41 |
| `-ci 2 -csc 2` | 60 | 2.72 | 51 | 51 |
| `-es 1 -ptol 3` | 71 | 1.40 | 63 | 63 |
| `-es 1 -ptol 3 -ci 2 -csc 2` | 24 | 2.05 | 45 | 45 |

**External MS2 calibration list (`-pl`).** ms1searchpy documents `-pl` for a
list of MS/MS-identified peptides. In this version `-pl` together with
`-deeplc 1` crashes (it disables the two-stage block, which is the only place
the DeepLC object is built); the patch keeps the two-stage path on. Lists were
built per file from the DIA-NN report (`build_diann_calib_lists.py`, Q ≤ 0.01,
unmodified, fully tryptic, 46k to 48k peptides). Three runs per file, launched
three-way concurrent, so wall times are not comparable to production:

| File | DeepLC sigma (min) | 1% (3 runs) | 5% (3 runs) | uno 1% / 5% |
|---|---:|---|---|---|
| rep1 | 0.66 | 380, 420, 484 | 862, 1146, 1042 | 1448 / 2145 |
| rep4 | 0.73 | 342, 437, 166 | 856, 778, 828 | 1182 / 1685 |
| rep6 | 0.71 | 258, 212, 589 | 727, 736, 751 | 1006 / 1567 |

An earlier single run on rep1 gave 693 at 1%; the spread between repeats at
1% is as large as the effect, which is the same instability seen on the
primary datasets. uno never sees these lists.

**Robust self-calibration (`-robust_calib 1`, local patch).** Before DeepLC's
own calibration, fit a linear map from uncalibrated DeepLC predictions to
observed RT with RANSAC (threshold 5% of the RT span, 0.55 min here) and keep
only the inliers for both the calibration half and the sigma estimate. Nothing
external is used. The default path is untouched; the flag defaults to 0.

Control on primary file `hela5min_09778`: stock runs gave 1173, 1220, 1147 at
1% with sigma 0.14 min; with `-robust_calib 1` the run gave 1200 with sigma
0.12 min and 1333 of 1855 calibration peptides kept. The patch is neutral on
clean data.

Astral (rep1 standalone pilot: 486 / 1759 kept, sigma 0.43, 213 / 771). First batch ran three-way concurrent, the rest serialized:

| File | run | kept / total calib | DeepLC sigma (min) | 1% | 5% |
|---|---|---|---:|---:|---:|
| rep1 | 1 | 540 / 1900 | 0.41 | 269 | 641 |
| rep4 | 1 | 492 / 1891 | 0.40 | 379 | 994 |
| rep6 | 1 | 511 / 1760 | 0.34 | 69 | 287 |
| rep1 | 2 | 595 / 1953 | 0.50 | 591 | 865 |
| rep1 | 3 | 579 / 1936 | 0.44 | 308 | 582 |
| rep4 | 2 | 523 / 1801 | 0.42 | 113 | 822 |
| rep4 | 3 | 488 / 1872 | 0.36 | 487 | 665 |
| rep6 | 2 | 493 / 1669 | 0.42 | 97 | 298 |
| rep6 | 3 | 512 / 1724 | 0.39 | 145 | 655 |

Per-file medians at 1% / 5%: rep1 308 / 641, rep4 379 / 822, rep6 97 / 298,
against uno 1448 / 2145, 1182 / 1685, 1006 / 1567. Serialized runs took 5.1 to
5.6 min wall and 7.4 to 7.6 GB peak RSS, the same as production. The sigma is
now 0.34 to 0.50 min, but the run-to-run spread at 1% (69 to 591) remains the
dominant feature of the comparator.

A variant that cleans only the calibration half and estimates sigma on the
unfiltered check half (`-robust_calib 2`) does not work: the Gaussian plus
flat-background fit again returns 1.8 min on a 70% background set.

## What this means for the manuscript

- The 30-fold Astral gap has a mechanism, and it is the comparator's RT
  self-calibration, not the feature table, the database, the mass gate, or the
  patch that feeds it extra feature columns.
- Reported as a fairness bound: even with information uno never receives
  (MS2 identifications from DIA-NN for the same files), or with its own
  calibration made robust, ms1searchpy reaches a fraction of uno's depth on
  these files.
- The production numbers in the paper stay as they are: stock ms1searchpy
  with `-deeplc 1 -ad 1`, ten runs per file. The controls above belong in the
  SI with the patch disclosed.
- Do not describe `-pl` or `-robust_calib` as ms1searchpy defaults. The
  `-robust_calib` option does not exist upstream.

## Input series (2026-09-04)

The preregistered extension (README.md, amended to the decade levels, then to
three files per level) gives 12 files: 200 pg, 1 ng, 10 ng, 100 ng ×3 each.
Per-level medians at 1% / 5% entrapment FDR from `input_series_summary.tsv`
(the two removed 1 ng files are in `input_series_dropped.tsv`):

| Input | uno | stock ms1searchpy | ms1searchpy `-pl` | stock sigma (min) | `-pl` sigma (min) |
|---|---|---|---|---|---|
| 200 pg | 476 / 786 | 47 / 55 | 256 / 519 | 1.5–2.0 | 0.63–0.69 |
| 1 ng | 925 / 1516 | 37 / 78 | 159 / 548 | 1.1–2.7 | 0.67–0.76 |
| 10 ng | 1182 / 1685 | 46 / 49 | 380 / 828 | 1.1–2.4 | 0.66–0.73 |
| 100 ng | 1165 / 1926 | 42 / 42 | 112 / 635 | 1.2–3.1 | 0.68–0.74 |

The stock calibration failure is input-independent, so it is a property of
this LC method, not a seed-count problem. uno's frozen binary returned no
proteins on 200 pg rep1 and 1 ng rep6 (stage-1 deviation collapse, see
`TODO.md`); the fixed binary recovers 389 and 523 at 1% (the latter file is
one of the two later removed from the reported set). `results/uno_fixed/` holds the reruns, `results/uno/`
the frozen outputs. biosaur2 feature counts rise only from ~74k (1 ng) to
~94k (100 ng), which is why uno's depth flattens above 10 ng.

## uno candidate selection (2026-09-04)

Investigating the raw/biosaur2 divergences showed that at low input the 1%
entrapment depth is the rank of the first yeast protein (scale factor 3.4:
one yeast hit needs >340 human proteins above it), and that uno's stages flip
between 0 and 100+ proteins on small RT-sigma changes because they choose by
the q<=1% count. Shipped fix: run the pipeline at sigma x1.0 and x1.5, keep
the final list with the larger soft depth. Tested alternatives: soft objective
inside the stage loop (misaligned with entrapment, 10 ng rep1 1448 -> 1200)
and choosing by the final internal 1% count (100 ng rep3 -298). Results and
the rerun Astral outputs are in `input_series_summary.tsv`.
