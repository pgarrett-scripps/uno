# KOTH 0.5.0 and the UNO matrix manuscript rerun

The publication study and manuscript checks are complete in
`benchmark/results/matrix_20260913_koth050`: 32 extractions, eight alignments,
144 searches, and four serialized timing controls all passed their audits.
It replaces the old numerical matrix section, while preserving the primary
single-run manuscript evidence. The KOTH manuscript is maintained separately.

The fixed study includes 32 fresh feature extractions, eight matrices, three
search seeds, and 144 searches. Astral includes four input-specific triplicates
and the full 12-file pool; IonStar includes balanced nested 5-, 10-, and 20-file
cohorts. Every joint search reads the long bundle's original observations.
Independent-search protein outputs are used only to evaluate the comparator
union. The second matrix arm adds direct/MBR extraction-confidence summaries.

The release-studied KOTH options explicitly enable long export and replicated
weak-seed retention. These remain optional in KOTH 0.5.0. Group member/seed floors
are 0.5/0.75; minimum group size is one; bounded spans and automatic reference
selection are fixed before searching. Detection uses the existing UNO profile.

The evaluator compares expanded same-species accession sets at empirical
entrapment FDP <=1%, admitting entire tied-score blocks. Matrix ranking uses
protein score; union ranking uses minimum per-run internal q. A separately
implemented Decimal/integer audit checks counts and accession sets at 1%, 3%,
and 5%. These are retrospective known-entrapment operating points, not claimed
cohort q-value calibration. Internal-q lists are retained separately.

Dedicated warm-cache timing repeats the three Astral 100 ng single searches and
one joint cell-summary search serially. Execution-only amendments used concurrent
Astral and IonStar queues, followed by a four-worker IonStar resume, while
preserving every search setting. Interrupted attempts are archived separately;
every final job has a successful receipt. Main concurrent jobs are not used for
resource claims. Quantitative evaluation uses internal-q abundance outputs,
shared baseline median-ratio factors, common-group comparisons, and no imputation.
Upstream positive/confident cell counts do not imply peptide identification.

Reproduction from the worktree root:

```bash
python benchmark/matrix_release_study.py run
UV_CACHE_DIR=/tmp/uno-paper-uv-cache uv run --project paper/analysis --locked python benchmark/matrix_release_report.py
cd paper
UV_CACHE_DIR=/tmp/uno-paper-uv-cache just assets-all
UV_CACHE_DIR=/tmp/uno-paper-uv-cache just fmt
UV_CACHE_DIR=/tmp/uno-paper-uv-cache just paper
UV_CACHE_DIR=/tmp/uno-paper-uv-cache just docx
UV_CACHE_DIR=/tmp/uno-paper-uv-cache just verify
```

`frozen.json` archives source, binaries, settings, input inventory, and the driver
hash. Each job has input/output checksums and an exit-code receipt. The report
refuses incomplete inventories. `matrix_release_evidence_links.json` inventories
read-only links to existing primary-study evidence required to build the paper
in this isolated worktree; those results are not rerun or overwritten.

KOTH release source is PR 1, head `976a2297bc94974f4082f0f7bf8d7c0580ce7203`,
merged as `0908b50e5f3f856545e4b69e2b6f3cea6a957d4a`. The merge has the same
source tree as the tested head. Version `v0.5.0` triggers all four platform
archives. The release remains in the existing private repository; no crates.io
publication was requested. All CI checks and four real Bruker tests passed.
A .NET 8.0.27 Thermo build extracted 7,588 features and 48,160 hills from
`QEHF1_09744_JB.raw`. All nine tuned KOTH-paper configurations parse.

The existing GitHub plan does not permit branch protection (API 403), so all CI
was manually checked before merge/tag. No signing key is configured; the tag is
annotated, not cryptographically signed. Absolute QDA extraction-q calibration
remains an accepted limitation; paper citation metadata remains pending the
separate manuscript work. The release is published at
https://github.com/tacular-omics/koth_rust/releases/tag/v0.5.0. All four archives
match both GitHub digests and their accompanying SHA-256 files; both downloaded
Linux binaries report 0.5.0 at the tagged source revision.

The three-seed median Astral 100 ng counts are 1,352 for the separate union,
1,718 for joint searching, and 1,723 with cell summaries. Joint searching also
wins at 200 pg in all seeds, but loses at 1 ng; 10 ng and the 12-file pool vary
by seed. The union exceeds both matrix arms in every IonStar cohort and seed.
For all 20 IonStar runs, medians are 1,911 for the union, 782 for joint, and 897
for joint plus cells. These are matched empirical entrapment accession counts,
not internally reported protein-group counts. The paper leads with Astral while
reporting all cohorts, the IonStar losses, and the absence of a consistent gain
from cell summaries.

The warm-cache Astral 100 ng timing control took 209.4 seconds for three separate
searches and 97.1 seconds for one joint cell-summary search (2.156-fold).
Protein outputs matched the study repetitions byte for byte, and every timing
search had zero RT-cache misses. This is one serialized, search-only comparison.
On common IonStar groups, the three-seed median E. coli E/A ratios were 1.75
without cell summaries and 1.72 with them, against the expected 3.

All 129 non-matrix statistics are identical to the pre-edit manuscript. All 152
generated statistics passed full re-derivation; formatting, extractor
tests, prose rules, asset provenance, and build freshness passed `just verify`.
The active asset build now skips five explicitly retired figure generators,
records both pipeline-schematic formats, and supplies explicit Word widths for
the new tables. Archived adopted resource and quality-column PDFs were preserved.
Prose warnings were reviewed: repeated technical names in feature definitions
and existing long primary-study tables remain intentional; no rule was disabled.

The build reports 185 abstract words, 4,885 main-text words, and 3,589 SI words
(8,474 main plus SI, excluding references, floats, captions, math, and block code).
Flesch-Kincaid grades are 12.1 for the main text and 11.8 for SI; reading-ease
scores are 39 and 40, respectively. The 37-page PDF passed the overflow check,
and changed main-text and matrix pages were visually inspected. All 52 rendered
Word pages were inspected, including repeated table headers and the new matrix
tables. The final Word render was pixel-identical on every page to the inspected
render. The KOTH manuscript was not edited.
