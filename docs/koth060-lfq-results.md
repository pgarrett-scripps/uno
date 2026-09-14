# KOTH 0.6.0 LFQ update

The released `koth_align 0.6.0 (7cc80fb84625)` rebuilt all eight LFQ bundles.
All 48 matrix searches and four serialized timing controls completed. Detection
and UNO were unchanged, so the audited 32 feature tables and 96 separate searches
were reused. Only the narrower RT extraction window affects these Orbitrap data;
none of their source features has measured ion mobility.

Three-seed median accessions at matched empirical entrapment FDP at most 1%:

| Cohort | Separate union | Joint | Joint + cell summaries |
|---|---:|---:|---:|
| Astral 100 ng | 1,352 | 1,718 | 1,618 |
| IonStar 20 runs | 1,911 | 782 | 883 |

All eight original-observation tables and all 24 joint protein tables without
cell summaries are byte-identical to the 0.5.0 study. Cell summaries still give
mixed identification results. Quantitation improves: the IonStar common-group
E. coli E/A ratios are 2.17 in both arms, versus 1.75 and 1.72 previously, with
3 expected. These common sets can change across versions. The joint all-group
ratio also improves from 1.76 to 2.15 while its protein identities stay fixed.
The new measurements remain compressed and do not establish cell-level FDR.

One warm-cache serialized control took 184.0 seconds for three separate searches
and 83.3 seconds for one joint cell-summary search (2.2-fold ratio). This excludes
feature finding, alignment, and LFQ construction; it is one timing repetition.

The manuscript keeps Astral first and retains every cohort, arm, and seed in the
SI. LFQ running prose fell from 1,549 to 923 words (40.4%), with Flesch–Kincaid
grade 12.27 to 11.85. `just paper` reports 4,857 main-text words and 2,950 SI words
(7,807 combined), with grades 12.1 and 11.6, respectively. Compact PDF tables
stay with their captions, and subsection labels stay with following text.

Validation passed: eight benchmark tests, exact independent accession-set
checks at all 216 cutoffs, reused-output/source/binary hashes, warm-cache output
parity, `just assets-all`, PDF and native Word builds, `just verify`, and
`just check-stats-deep`. All 129 non-matrix statistics are unchanged. All 51 Word
pages and the changed PDF pages were visually reviewed; the PDF overflow check
passed on every page. No one-sentence LFQ prose paragraphs were introduced.

Evidence is in `benchmark/results/matrix_20260914_koth060/`, including the frozen
protocol, release checksums, reports, and `verification/` logs. The committed
paper input is `paper/figdata/koth_matrix_study.json`. This revision remains on
`codex/koth060-lfq`, separate from the ongoing abstract edit on `main`; the older
`paper-2026-09-13` tag is preserved.
