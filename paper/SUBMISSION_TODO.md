# Minimum path to a JPR submission

Scope: what it takes to get **uno** submitted to and through *Journal of Proteome
Research*. This is the manuscript half of the work; the code half is in
[`../TODO.md`](../TODO.md), and the audit's full evidence base — `file:line` for
every finding, plus what was verified correct — is in [`AUDIT.md`](AUDIT.md).

Reference pipeline: `/home/patrick-garrett/Repos/d_noise-paper/paper`, a
near-complete JPR submission by the same author. Where uno is missing an artifact
that repo already has, it is called out below.

**Status legend:** `[ ]` to do · `[x]` done · `[!]` blocked on someone else.

---

## The honest summary

The manuscript is in good shape as a *document* and not yet defensible as a
*result*. Prose, structure, style and scholarship are largely fixed. The evidence
layer has not moved: every identification and calibration number still comes from
two different builds, one of which has a since-fixed ~20% order-sensitivity bug.

There is exactly one expensive thing on the critical path (the re-run, ~half a
day). Everything else is hours of editing plus tooling that already exists next
door in `d_noise-paper`.

**Minimum realistic path: the re-run, then §B, §C and §D below.** Roughly one
focused week, most of it not compute.

---

## A. Already done (do not redo)

- [x] Abstract trimmed 301 → **200 words** (JPR limit is 200).
- [x] Full ACS back matter: Associated Content, Author Information, Author
      Contributions, Notes (no competing financial interest), Ethics,
      Acknowledgments with funding **and** the required AI-use disclosure.
- [x] ORCID iDs removed from the manuscript body (ACS links them via Paragon Plus).
- [x] Section renamed to **Experimental Section** (JPR's required name).
- [x] **ACS bibliography style** — `#bibliography("paper.bib", title: [References],
      style: "american-chemical-society")`. Verified: renders ACS numbered format
      with a "References" heading. This was the fix `d_noise-paper/paper/paper.typ:1054`
      already had.
- [x] All six house-style rules pass; 39 tests pass; compiles with zero warnings.
- [x] Tables S2/S3 page overflow fixed (`#show figure: set block(breakable: true)`).
      Bbox scan across all 32 pages: **0 out-of-box words**.
- [x] Title updated to include quantitation.
- [x] Instrument corrected to Orbitrap Fusion Lumos with FAIMS Pro.
- [x] Methods now describes the shipped engine; 9 missing citations added.
- [x] **The directory now runs on `paper-scaffold` v1.0.0** (`just version`).
      `config.typ` holds the paper's identity; `paper.typ` carries the
      `BODY START`/`BODY END` markers the word counter, readability/density
      reports and narrator slice at; `paper.bib` → `references.bib`; the figure
      scripts moved to `analysis/scripts/` behind the `just assets` contract with
      a content stamp; `scripts/paper_stats.py` retired in favour of
      `readability.py` + `density.py` + `prose_check.py`, with its four unique
      checks ported. Added: `just density` (numerals, parentheticals, acronyms,
      nominalizations, passives per 1,000 words, flagged against the paper's own
      per-section median), `just prose-check`, `just test` (golden-file extractor
      tests, reflow-invariant), `just fmt`, and the audiobook pipeline.
- [x] Three British spellings ("behaviour" ×2, "analysed") and a false
      `@klammer2007rtprediction`-is-uncited report fixed. Found by
      `just prose-check` on its first run.
- [x] **2026-08-14 restructure: the main text was cut to five Results sections**
      (head-to-head incl. FAIMS / determinism / entrapment calibration /
      protein-level identity + LDA caution / quant dilution). Old §§3.6–3.7 (RT
      distillation + reverse-sequence benchmark), the LDA hardening detail, the
      quality-column ablation detail, and the quant DE post-mortem moved to the
      SI (new S8 reproducibility detail, S9 development-configuration
      experiments; superseded table now S10). The intro now claims three
      contributions, not four. Figures renumbered: old fig3_determinism → SI
      figS6 (its spread% panel promoted into fig2 as panel c), fig4_entrapment
      → fig3, fig5_quant → fig4 (volcano moved to figS5 panel c). Any section
      or figure number elsewhere in this file that predates this note refers to
      the OLD numbering (in particular C1's "Figures 2/3/4/S1/S3/S4" and C3's
      "§3.3", which is unchanged in content and still §3.3).

---

## B. Hard submission blockers

Cannot submit, or will bounce at editorial screening.

- [ ] **B1. Fill the Data and Software Availability statement.** Currently
      `[repository URL, to be assigned on release]` and
      `[10.5281/zenodo.XXXXXXX]`. JPR requires resolved deposition and access
      links. Blocked on the release decision (D5: "release all three, not yet").
      → uno public repo, koth archived (at minimum a versioned Zenodo deposit,
        since it is load-bearing for all of §3.9), and a Zenodo archive carrying
        the benchmark data, RT distillation labels and trained models.
      → Note **every benchmark input and result file is currently gitignored**
        (`benchmark/.gitignore`), so a reader cloning the repo today gets none of
        the inputs behind any SI table.
- [!] **B2. Confirm the funding block.** `paper.typ` currently carries the grant
      list copied verbatim from the d_noise manuscript (NIH R01 HL165168,
      R01 AG077046, R01 MH100175, R01 AG075862 to J.R.Y.; R01 MH132570 to L. Ye;
      U01 AG088679 to S. A. Lipton; Skaggs Graduate School for P.G.). That paper
      has a third author and a different scope. **Author must confirm or correct.**
      Also required separately in ACS's Funder Registry Tool at submission.
- [ ] **B3. Add the corresponding author's telephone number.** ACS requires it;
      currently absent. `d_noise-paper/paper/cover_letter.md` uses a `[PHONE]`
      placeholder for the same reason.
- [~] **B4. Split the SI into its own file.** **Source split DONE** — the SI body
      now lives in `si-body.typ` and is `#include`d by `paper.typ`, matching
      `d_noise-paper`. **The output split is still open**, and it is not a build
      flag: **10 main-text cross-references point at SI labels**
      (`@tbl:perdataset`, `@tbl:entrapreported`, `@tbl:faims`,
      `@tbl:qualityablation`, `@tbl:qualityauc`, `@tbl:quantdetail`,
      `@fig:fdpcurves`, `@fig:faims`, `@fig:resources`, `@fig:determinism`;
      the set changed in the 2026-08 restructure, see the note in §A), and
      Typst cannot resolve a reference across two documents. Note `d_noise-paper` has not solved this either — it also ships
      one combined PDF. Options at bundle time: (a) convert those 11 to plain text
      ("Table S5"), losing the hyperlinks; (b) compile the SI standalone and
      accept broken refs in the main text; (c) ask the editor whether a combined
      PDF is acceptable for review, with separate files at revision.
- [x] **B5. Produce `paper.docx`.** DONE. `typst2docx.py` ported from
      `d_noise-paper` and a `just docx` target added. `paper.typ` gained a
      `docx-mode` switch that bypasses arkheion (whose page rules HTML export
      ignores, dropping the title block and abstract) and wraps equations in
      `html.frame` so they survive as rasterized inline images. Verified output:
      4.7 MB, 176 paragraphs, **12 tables**, **34 numbered headings**, title block,
      abstract and keywords all present, 65 equations rasterized. `paper.docx` is
      gitignored; `just check` flags it stale by mtime.
- [x] **B6. Write the cover letter.** DONE — `paper/cover_letter.md`, following the
      `d_noise-paper` structure. Two placeholders remain: `[PHONE]` (see B3) and
      the repository/DOI line (see B1), plus a note to confirm the main-text word
      count against the current guidelines at submission.
- [x] **B7. TOC / abstract graphic.** DONE (2026-08-16) —
      `analysis/scripts/fig_toc_graphic.py` -> `figures/fig_toc_graphic.png`,
      exactly 3.25 in × 1.75 in at 300 dpi (the ACS "For Table of Contents
      Only" box). Pipeline schematic + wall-clock vs depth at 1% entrapment
      FDR, every number computed from `benchmark/results/summary.tsv`, so it
      For the draft it is embedded in `paper.typ` directly under the abstract;
      at submission move it to a separate Paragon Plus upload (ACS convention
      is a "For Table of Contents Only" page). Re-check the size spec against
      the current JPR guidelines at submission.
      **Working source is now `paper/toc_graphic.pptx`** (editable PowerPoint,
      6.5 in × 3.5 in = 2× the ACS box; scaffolded once by
      `analysis/scripts/gen_toc_pptx.py`, which OVERWRITES hand edits if
      rerun). Finish the slide there, export it as PNG at the ACS size over
      `figures/fig_toc_graphic.png`, and then delete
      `analysis/scripts/fig_toc_graphic.py` so `just assets-all` cannot
      clobber the export with the drawn placeholder.

---

## C. Integrity blockers

These will not stop submission. They will sink peer review, because each is a
claim a reviewer can falsify using the paper's own data or the cited deposits.
The audit found them; the decisions are settled; the edits mostly wait on the
re-run.

- [x] **C1 preconditions: the code fixes are DONE** (2026-07-31). All three D11
      fixes landed, plus a two-run byte-diff check (`benchmark/determinism_check.sh`)
      and unit tests, plus two defects found while verifying: a protein-group table
      keyed on accession alone that collapsed target and decoy groups, and a build
      break (`koth_ff::PipelineOptions::emit_ms2`) that meant **uno did not compile
      at HEAD**. §3.2's determinism paragraph and the Methods sentence were rewritten
      to describe what now holds. Depth moves: `hela5min_09746` 1% 862 -> 876,
      5% 1209 -> 1194 on a single dataset. Details in [`../TODO.md`](../TODO.md) §1.
- [ ] **C1. The re-run.** *This gates every number in Sections 3.1–3.5 and 3.9,
      Table 1, Tables S2/S3/S5/S6, and Figures 2/3/4/S1/S3/S4.*
      Preconditions and full scope are in [`../TODO.md`](../TODO.md) §1–2. In short:
      land the two determinism sorting fixes and the `score.rs:465` decoy guard
      first, add a two-run byte-diff test, commit and pin the tree, then run
      uno once and ms1searchpy three times over all six datasets plus run10, plus
      the run10 column-stripped arm and the FAIMS replicates. Budget half a day,
      serialized on an idle host. Do not write results prose before it completes.
      - Why it is not optional: `si/versions.typ` claims a single commit, but the
        Table 1 outputs carry an 8-column schema dated 07-22 and the FAIMS
        outputs a 10-column schema dated 07-28. Two builds, six commits apart.
        And commit `eccf0b2` reports run10 moving **1309 → 1426** after a
        row-order fix that landed *after* the six-dataset results were measured.
- [ ] **C2. Disclose the tuned defaults** (decision D2, settled). `top_k = 50`,
      the ensemble inclusion band, `min_child_weight = 80` and `cut_max = 8.0`
      were each selected by maximizing yeast-entrapment depth on run10/deephela —
      development data from the same family as the published benchmark, with the
      tuning recorded in the repo's own commit messages. Add the disclosure
      paragraph and relabel run10 a **tuning set**, not a benchmark.
      Include the counterweight: the shipped RT model was chosen *against* its
      run10 score, which is evidence of restraint.
- [ ] **C3. Re-verify the calibration paragraph after the re-run.** §3.3 is now
      correct for the *current* numbers (four of six conservative, 09778 at
      1.03%, 09746 at 2.16%, 12 of 18, 13 of 18). It is the most number-dense
      passage in the paper and every one of those values will move.
- [ ] **C4. Resolve the FDP-vs-q column** (decision D6: print the raw FDP).
      `si/entrap_reported.typ` is headed "Entrapment FDP (%)" but holds the
      running-minimum entrapment *q*, which Methods defines as a different
      quantity. Six sentences move with it. Fold into the post-re-run pass.
- [ ] **C5. Disclose the ms1searchpy patch.** The benchmarked comparator is a
      local uncommitted +45-line patch, not the released 2.8.12 — which **does not
      exist on PyPI** (latest is 2.8.9). Without the patch, stock ms1searchpy
      reads none of the five quality columns that §3.5 and S6 depend on. Pin the
      commit, ship the diff in the SI, and say plainly that the patch gives the
      comparator *more* information, so uno's margin is conservative.
- [ ] **C6. Fix `si/versions.typ`.** Per-results-block commits for uno; add the
      sage, koth and sequoia-boost commits (all out-of-tree path dependencies, so
      a commit is the only pin); add matplotlib, typst and OS/kernel; correct the
      NumPy/SciPy role line (the scoring scripts are stdlib-only).
- [x] **C7. `@tbl:rtdepth`'s "69-d" — resolved by removing the count.** The only
      surviving definition of the set is the initial commit's `rt.rs` header,
      which enumerates 20 AA + 20 N-term-window + 20 C-term-window + length +
      ln(1+mass) + intercept = **63**. "69" traces to a stale doc comment at
      `rt.rs:24` copied into `uno_jpr.md:213` and thence into the paper; no
      arithmetic anywhere supports it. Rather than swap one bare count for
      another, the prose now *describes* the encoding (AA counts, terminal-window
      counts, length, mass) and the table row reads "Composition-only". The
      dimensionality was decorative; the composition-vs-order-aware contrast is
      what carries the result. `rt.rs:24` corrected to 63.

---

## D. Required before the SI can be called honest

Cheap, but a reviewer checks these in seconds.

- [x] **D1. Fix the two "GENERATED … do not edit by hand" labels.** DONE — both
      relabelled "TRANSCRIBED BY HAND", with the reason and a pointer to the re-run's
      "re-derive every table by script" step, where the generator belongs (it can be
      validated against fresh output there, rather than against the numbers it was
      written from). Superseded text: No generator
      for `si/results_per_dataset.typ` or `si/entrap_reported.typ` exists anywhere
      on disk. Either write it or relabel them "Transcribed". The label actively
      discourages the only check that would catch drift.
- [x] **D2. Make S1's commands runnable as printed.** DONE — `si/commands.typ`
      rewritten. Every variable is now defined (`ROOT`, `DS`, `OUT`, `RT`, `MS1`,
      `DATA`, `R`, `QDATA`), the `cd`/copy that ms1searchpy needs is shown in both
      places, the omitted `score_entrap.py ms1` line and the
      `N_TARGET_DB`/`N_ENTRAP_DB` exports are back, `intensity_matrix.tsv` has a
      resolvable path, and the **koth_align step is documented for the first time**
      with its `[lfq] averagine_projection = true` config — so following S1 now
      reproduces the reported 326-protein averagine arm rather than the deprecated
      243-protein box-sum one. Header no longer claims to be verbatim, and names
      the two deliberate differences from the scripts. Superseded text: Undefined `${DS}`/`$DATA`,
      a missing `cd`/copy, the ms1searchpy FAIMS scoring line omitted entirely,
      `intensity_matrix.tsv` with no resolvable path, and no koth alignment
      command at all — so following S1 reproduces the deprecated 243-protein
      box-sum arm rather than the reported 326.
- [x] **D3. Add a dataset-name mapping table to S1.** DONE — `si/datasets.typ`,
      wired in as Table S2, mapping name in text -> deposit -> raw file -> feature
      file, with run10, QEHF1_09744_JB and deephela below a rule as the three runs
      that are NOT part of the six. Superseded text: A reader cannot map `09746`
      to a file in PXD015669. Note that **09744 ≠ 09746**.
- [~] **D4. Remaining figure-caption/units fixes.** Figures 3, S3 and S4 rewritten
      and the GB->GiB axis labels fixed in both scripts and the manifest (both
      figures re-rendered). Remaining: the Figure 4 caption never identifies its
      median curves or its two annotated points, and Figure 5 / S5 panel letters.
      Superseded text: (Figures 3, S3, S4; GB→GiB axis
      labels). Full list in [`AUDIT.md`](AUDIT.md) §8. Figures 1 and S1 are
      already done.
- [x] **D5. Re-render fig5 and figS5 under the venv.** DONE. They were the only
      two rendered by matplotlib 3.6.3 rather than 3.11.0. `just assets` flagged
      exactly those two as stale and re-rendered them under the
      `matplotlib==3.11.0` pin now recorded in `analysis/pyproject.toml`. The
      numbers they print are unchanged (22 UPS1 / 23 *E. coli* / 60 human;
      82% = 151/184; 1% = 1/82), so only the renderer moved.

---

## E. Build and release tooling

Not required for correctness, but it is what makes the submission bundle
reproducible. Most of this arrived wholesale when the directory moved onto
`paper-scaffold` (see section A); what is recorded below is what each item means
here now.

- [x] **E1. `typst2docx.py`** and the `docx` justfile target. DONE (see B5).
- [x] **E2. A `check` justfile target.** DONE, and now the scaffold's: `just
      check` runs `check-pdf`, `check-assets` (content stamp of `analysis/`
      against `.assets-stamp`), a `paper.docx` staleness check by mtime, an
      existence check on every embedded figure, and the PDF overflow scan.
      Currently reports `check: OK`. The house-style rules and the tooling's own
      tests moved to their own targets, `just prose-check` and `just test`.
- [x] **E3. `check-pdf`** staleness guard. DONE, over `paper.typ`, `config.typ`,
      `si-body.typ`, `references.bib`, `si/`, `figures/`, `figdata/` and
      `analysis/`. Caveat: it compares *commits*, so it cannot see uncommitted
      edits. Commit before trusting it.
- [x] **E4. Track the paper trees.** DONE — `paper/figures/`, `paper/si/`,
      `paper/analysis/`, `paper/figdata/`, `paper/tests/`, `paper/audio/`,
      `si-body.typ`, `config.typ`, `cover_letter.md` and both TODO files are
      tracked. `paper/.gitignore` covers `paper.docx`, the HTML intermediates,
      `__pycache__`, both venvs and every audio artifact. `paper.pdf` stays
      tracked deliberately, since `check-pdf` guards it.
      **Not yet committed** — so that the re-run's numbers and this tooling can
      land as separate, reviewable commits.
- [x] **E6. `check_overflow.py`** (new, not in `d_noise-paper`). Bbox scan
      via `pdftotext` for any glyph outside the printable area. This is the check
      that caught the Tables S2/S3 overflow, which `typst compile` reports as
      perfectly clean. Wired into `just check`.
- [~] **E5. Add `CITATION.cff`** to uno (and koth if released). **The
      environment half is DONE**: `paper/pyproject.toml` + `uv.lock` pin the
      manuscript toolchain and `paper/analysis/pyproject.toml` pins the figure
      dependencies (`matplotlib==3.11.0`), both built by `uv`, so a fresh clone no
      longer depends on the reader happening to have matplotlib. `CITATION.cff` is
      still to write, and is blocked on the same release decision as B1.

---

## F. Known-remaining, explicitly accepted

Not blockers. Recorded so they are not rediscovered as surprises.

- **Journal-name abbreviation.** The ACS style renders full journal names
  ("Mass Spectrometry Reviews") rather than CASSI abbreviations
  ("Mass Spectrom. Rev."). ACS copyediting normally handles this; if a reviewer
  objects, the fix is abbreviated `journal =` fields in `paper.bib`.
- **`arkheion` is a preprint template, not an ACS layout.** JPR accepts
  author-formatted manuscripts, and the figures are already sized to ACS column
  widths (7.0 in double-column). Not worth restyling before a decision.
- **Sections 3.6–3.8 are not reproducible** from the archived pipeline. The code
  behind them (the LDA rescorer, `--hills`, the on-data RT refit) was purged in
  `a91c6a9`. The SI now says so explicitly. This is disclosed, not hidden.
- **`@tbl:qualityauc` values are unauditable** (diagnostic scripts deleted) and
  contain an internal impossibility: `cosine_score` and `isotope_score` are given
  as exactly equal for both engines, which the data cannot produce. Decision D14
  is "re-derive from surviving data" — the ms1 side is on disk, the uno side needs
  a PFM dump from the C1 re-run.
- **§3.7's reverse-sequence benchmark has no surviving artifacts.** Decision D13
  is to reconstruct it (needs a DeepLC run). It is the paper's cleanest novelty
  claim and Lane G confirmed no prior art, so it is worth rebuilding rather than
  disclaiming.

---

## Suggested order

1. **C1 re-run preconditions** (determinism fixes + test), then the re-run itself.
   Start it early; it is the only long pole.
2. While it runs: **B5, B6, E1–E5** (docx, cover letter, tooling, commit the
   untracked trees) and **D1–D3** (SI honesty).
3. After it lands: **C3, C4, C6** (numbers, FDP column, versions), then **C2, C5**
   (the two disclosures), then re-verify every table by script.
4. **B1–B4, B7** last, at bundle time.
5. Final gate: [`AUDIT.md`](AUDIT.md) §10 — compile clean, bbox scan, all
   style rules, every number re-derived by script, clean clone builds.
