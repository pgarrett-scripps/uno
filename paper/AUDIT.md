# Pre-submission audit — the full record

> **This is the evidence base, not the work list.** For what is actually left to
> do, see [`SUBMISSION_TODO.md`](SUBMISSION_TODO.md) (the manuscript) and
> [`../TODO.md`](../TODO.md) (the code). Both were extracted from this file on
> 2026-07-31, scoped to what it takes to publish in JPR.
>
> Kept in full because it is what makes those lists actionable: every item cites
> `file:line` and the artifact it was checked against, the 16 settled decisions
> record *why* each call was made, and the appendix lists what was verified
> correct — which is what stops a finding being re-litigated later.
>
> **Checkbox states below are frozen at the audit date and were not maintained
> after it.** Several `[ ]` items have since been done (see §0b, and the
> 2026-07-30/31 passes). Trust the two work lists, not the boxes here.

Source: 11-lane pre-submission audit of `paper/paper.typ` (main text + SI), the analysis
scripts, and `crates/`, run 2026-07-29 against HEAD `7f2cece`.

Two kinds of entry:

- **`[ ]`** — obvious fix. The audit established the fact; no judgment call. Just do it.
- **`[?]`** — needs an author decision. Listed in **§0 Decisions** with the options.
  Do not start a `[?]` item until its decision is recorded there.

Ordering matters. See **§11 Ordering constraints** before picking up work.

Nothing in this list is a guess: every item cites `file:line` and the evidence file.

---

## 0. Decisions pending

Filled in as they are settled. Until a decision is recorded, every `[?]` item that
depends on it is blocked.

| id | question | decision | date |
|----|----------|----------|------|
| D1 | Re-run the benchmark at one pinned HEAD, or label each results block with its actual commit? | **RE-RUN EVERYTHING at one pinned HEAD.** uno once + ms1searchpy ×3 on all six datasets plus run10, serialized on one host. Every number in Table 1, Table S2, Figs 2/3/4/S1/S3, the abstract and the Conclusions is provisional until this lands. See §12 for the re-run protocol. | 2026-07-29 |
| D2 | Disclose that `top_k`, the ensemble band, `min_child_weight` and `cut_max` were tuned on run10/deephela entrapment depth? | **FULL DISCLOSURE PARAGRAPH.** Name each tuned default and its tuning set, state that the head-to-head is partly in-sample for uno, relabel run10 a tuning set, and include the RT-model counterweight (shipped model chosen *against* its run10 score). | 2026-07-29 |
| D3 | Describe the default-on GBDT peptide rescorer in Methods? | **DESCRIBE FULLY IN METHODS.** Add the stage, cite/declare sequoia-boost, and sharpen §3.8 to "learning helps at the peptide level, hurts at the protein level". "Deterministic binomial *protein* scorer" survives. | 2026-07-29 |
| D4 | Report the differential-test precision (21%, 60 human FPs), or drop the specificity claim and report nothing? | **REPORT 21% AND DROP THE CLAIM.** State 22 UPS1 / 23 *E. coli* / 60 human, precision 21%, with the imputation caveat; delete "high specificity" from the abstract, Limitations and Conclusions. Overrides the `benchmark/quant/README.md:70-72` decision not to publish a precision figure. | 2026-07-29 |
| D5 | Release koth, the uno repo, and a Zenodo archive before submission? | **RELEASE ALL THREE, BUT NOT YET.** Intent confirmed; timing deferred. Treat as a pre-submission gate, not current work: write Data and Software Availability last, once the URL and DOI exist. Until then koth stays private, so keep the §3.10 dependency disclosure accurate in the interim. Do not draft the availability statement early. | 2026-07-29 |
| D6 | `@tbl:entrapreported`: print the reported list's raw FDP, or rename the column "entrapment q"? | **PRINT THE RAW FDP.** Matches Methods:283 and the printed Human/Yeast columns. New values: uno 0.00–2.24%, median 0.744%; ms1 worst 2.11%, median 1.25%; at 5% 6.87/5.19% and the other four 2.47–3.64%. Six sentences move. | 2026-07-29 |
| D7 | On-data RT fine-tuning subsection: delete, or retire to S8 as a development experiment? | **DELETE THE SUBSECTION** (paper.typ:196-209), and delete the S8 sentence at :1041-1042 that cites the RT ablation as evidence for it. Also remove the three dead keys from the config listing. | 2026-07-29 |
| D8 | Keep the "Validation roadmap" section? | **CUT THE SECTION** (paper.typ:785-798). Fold the two informative items (the coarse entrapment bound, the trace-level sensitivity gap) into the last paragraph of Limitations. The release-hygiene admission becomes stale anyway once D5 and D11 land. | 2026-07-29 |
| D9 | Keep §3.9 (performance engineering) in the main text, or move to SI? | **CUT ENTIRELY** (paper.typ:642-667). The 40×/13× work survives only in `CHANGELOG.md:89-93` after this — accept that. **Knock-on: §3.10 renumbers to §3.9**, so every hand-written "Section 3.10" pointer must be updated (the audit confirmed all 30 section pointers currently resolve, so the set is known). Also drop §3.9 from the Methods scoping paragraph rather than adding it. | 2026-07-29 |
| D10 | Is 1.03% "conservative"? i.e. adopt a stated tolerance band applied to both engines, or count it as above nominal? | **COUNT BOTH AS ABOVE NOMINAL.** "Four of six conservative, one at the nominal level, one about twice optimistic." Keep ms1searchpy's 12-of-18 count as-is. Under D6's raw FDP the two uno values are 2.24% and 1.06%. | 2026-07-29 |
| D11 | Apply the two determinism sorting fixes and the `score.rs:465` decoy guard? (Changes output; gates D1.) | **FIX ALL THREE, THEN RE-RUN.** Land the `est_rt_sig` residual sort, the quant peptide-vector sort, and the `pep.decoy` guard at `score.rs:465`; add the two-run byte-diff test; then do the D1 re-run **once**. Note the decoy guard may move depth, so the re-run's numbers are genuinely new, not a reproduction. **Must land before the re-run.** | 2026-07-29 |
| D12 | Ship the SI as a separate file (JPR requirement) or keep the single combined PDF? | **SPLIT AT SUBMISSION TIME.** Keep the single-PDF build for drafting; split into `paper.pdf` + `si.pdf` only when assembling the bundle. The S-numbering show-rules already work, so leave them alone. Update `paper/FIGURE_MANIFEST.md:28-29`, whose "there is no separate supplementary.pdf" rule now holds only for the draft build. | 2026-07-29 |
| D13 | §3.7 reverse-sequence benchmark: reconstruct the artifacts, or disclose non-retention? | **RECONSTRUCT IT** as a retained, deposited script: 46k held-out tryptic peptides plus reverses, DeepLC labels, forward-minus-reverse delta correlation. Needs a DeepLC run. Makes "reusable" true and turns the paper's cleanest novelty claim into a shippable artifact. Deposit with the Zenodo archive (D5). | 2026-07-29 |
| D14 | `@tbl:qualityauc`: demote to qualitative, or re-derive the numbers? | **RE-DERIVE FROM SURVIVING DATA**, with the confident-target definition stated explicitly in the caption. ms1's side comes from the PFM tables already on disk; **uno's side needs a PFM dump, so add that to the D1 re-run's outputs** (see §12). | 2026-07-29 |
| D15 | Retitle to include quantitation? | **ADD "AND QUANTITATION"**: "uno: A Fast, Reproducible, Dependency-Light Rust Engine for MS1-Only Protein Identification and Quantitation". Matches the abstract's own framing. | 2026-07-29 |
| D16 | Back-matter content: funding sources, COI, author contributions, AI-use disclosure wording. | **ALL FOUR TO BE HANDLED NOW.** COI: standard no-competing-interest statement confirmed. Author contributions, AI-use disclosure: to be drafted for author review. Funding: **grant numbers and funder names still needed from the author** — this is the one item that cannot be drafted. | 2026-07-29 |

---

## 0b. Completed — 2026-07-30, low-compute pass

Verified after each batch: `typst compile` exits 0 with no warnings, and
`just test` passes (39 tests, up from 38).

**Manuscript structure**
- Title now ends "…Protein Identification and Quantitation" (D15).
- Typed `orcid:` fields removed from the author block; replaced with the ACS note
  that iDs are linked via Paragon Plus.
- Date "June 2026" → "July 2026".
- "Experimental Procedures (Methods)" → "Experimental Section".
- **On-data RT fine-tuning subsection deleted** (D7), and the S8 sentence that cited
  the RT ablation as its evidence rewritten to say the refit was evaluated during
  development and is not in the released engine.
- **§3.9 Performance engineering deleted** (D9). §3.10 renumbered to §3.9; the two
  live pointers updated (`:738` and the provenance comment at `:58`). Confirmed no
  other reference existed.
- **Validation roadmap deleted** (D8); its two informative items folded into the end
  of Limitations.
- **Full ACS back matter added**, adapted from `~/Repos/d_noise/paper/paper.typ`:
  Associated Content with an SI description covering all eight SI sections;
  Author Information with corresponding author, contributions, Notes (no competing
  financial interest) and an Ethics statement; Acknowledgments carrying the funding
  block **and the AI-use disclosure**; and Data and Software Availability naming all
  three deposits. Release URLs remain marked placeholders per D5 ("not yet").
- S1 preamble no longer claims "every measurement"; it now names what it covers and
  states plainly that the §3.6–3.8 experiments predate the archived pipeline and
  cannot be re-run from it.

**Config truth (§1)**
- `example_config.toml`: `missed_cleavages` 2 → **0** and `variable_mods` → **empty**,
  matching `DigestConfig::default()`; added the four decoy fields, `min_charge`, and
  the eleven `[scoring]` defaults that were missing, including `ml_rescore = true`;
  removed the three inert `ondata_rt*` keys. Added the tuning-provenance note on
  `top_k` and the ensemble band (D2).
- `paper/si/defaults.toml` regenerated from the config structs rather than copied,
  and the SI prose changed from "The complete default set is" to an accurate claim.
- `paper/si/uno.toml` header corrected: values verbatim, comments condensed, and only
  `max_len` actually differs from a uno default.

**Stale documentation (§1b)**
- `rt.rs` header: the 63-feature Klammer set is described as replaced, not current.
- `rt_cnn.rs` header: ProteomeTools/Prosit → **Chronologer**, with an explicit note
  that this model is not DeepLC-distilled.
- `database.rs`: "GLOBAL-POOL" → per-protein pool, "different composition" softened
  to "generally", and the absence of a decoy-vs-target collision filter stated.
- `score.rs`: `ml_rescore` doc comment "OFF by default" → ON, active in every
  reported benchmark; ensemble comment corrected from "top-30% (lo=1, hi=3)" of 10
  bands to the shipped 14-of-40, with the tuning provenance recorded.
- `CHANGELOG.md`: the on-data RT entry struck through and marked removed in `a91c6a9`.
- `docs/uno_vs_ms1searchpy.md`: dated note that `--out-peptides` and the quant path
  are back at HEAD, and that `ml_rescore` is on.
- `PLAN.md`: "Decisions (locked)" marked superseded, with what actually ships.
- `fig5_quant.py` / `figS5_quant_detail.py` provenance blocks corrected to the
  new-koth averagine matrix, the real row count, the gitignored status of the
  outputs, and the two different detection denominators.

**Tooling (`paper_stats.py`)**
- Taught the parser to recognize `#heading(numbering: none)[…]` as a level-1 heading.
  Without it every bold run-in label in the ACS back matter was miscounted as bold in
  running prose (rule 3 went from PASS to 8 false failures).
- Journal-mandated back matter is now exempt from the semicolon and sentence-length
  rules: the SI description is *required* to be a non-sentence semicolon-separated
  list, which was producing 12 unfixable false violations that would have masked real
  ones during the §9 sweep.
- Added `test_backmatter_sections_are_recognized`; scoped
  `test_no_table_cells_or_metadata_in_body` to exclude back matter, since the
  corresponding author's email now legitimately appears there.

**House-style counts are back to the pre-existing baseline**: 4 em-dashes,
6 semicolons, 4 sentences > 40 words — all of them the pre-existing ones in §2.7,
§3.9 (quantitation) and the Figure 5 caption, i.e. the §9 sweep. Rules 2, 3 and 11
PASS.

### Second batch — Methods↔code, citations, mechanical sweep

**Methods now describes the shipped engine (§3)**
- Decoys: "global peptide pool" → per-protein pool, with the code's own stronger
  rationale; "different composition and different mass" softened to "generally",
  with the mass shift stated and the absence of a collision guarantee admitted.
  Added the honest caveat that the anti-reversed-decoy argument is about
  composition-based predictors and is weaker for uno's own order-aware model.
- Matching: isotope offsets restated as ±1 carbon-13 neutron (1.00335 Da) on the
  neutral mass; top-$K$ corrected to per feature **per class**, so up to 100
  candidates.
- RT: "two students" → "two models" (only the first is distilled); the 13 scalar
  descriptors now enumerate 13; the 11-vs-13 global/scalar distinction made
  explicit; "shipped default" → "distributed with uno and passed explicitly".
- **Scoring rewritten.** The 2-D mass/RT grid description replaced by the shipped
  1-D combined-deviation sweep; the false claim that "decoys pass the identical
  filter, so the FDR estimate stays calibrated" replaced by an honest statement
  that the selected cutoff is an optimistically-biased order statistic, that the
  sweep cap bounds it, and that the entrapment probe measures the residual.
  The `(d+1)/t` convention now stated, with the note that ms1searchpy uses it too.
- **The default-on learned rescorer is now described (D3)**, with sequoia-boost
  cited, and explicitly distinguished from the protein-level LDA of §3.8.
- Quantitation: feature extraction claim corrected for the `--raw`/koth_ff path;
  profile reassignment marked opt-in; the estimators attributed to koth and the
  averagine projection marked non-default and explicitly enabled; the minimum-runs
  floor stated as a requirement; "stops there" corrected to precursor- and
  protein-level matrices; the DirectMS1Quant relationship softened to "in spirit".
- Input: required-vs-optional columns separated; `nScans`/`neutron_offset` mismatch
  resolved; the four omitted columns added; ms1searchpy's extra requirements noted.
- Limitations: **TMT reattributed** from DeepLC to uno's own modification-blind RT
  encoder; **IM/CCS corrected** to "trained and loadable but not wired into
  scoring"; the search-breadth claim scoped to one specificity and reconciled with
  the multi-enzyme result.

**Citations (§5)** — nine entries added and verified: `wen2025entrapment`
(*Nat Methods* 2025, the entrapment-methodology reference the estimator needed),
`madej2024entrapment`, `xiao2025swaps`, `zimmer2006amt`, `ivanov2017multienzyme`,
`yu2021ionquant`, `kall2007percolator`, `petritis2006sequence`, `moruz2017rtreview`,
plus `sequoiaboost`. Wired in at the FDP definition, the Introduction, §3.7's
premise, §3.8's safeguards (alongside, not replacing, `the2016percolator`) and the
match-between-runs discussion. `note`-field locators moved to `url`/`version`/
`howpublished`, which the renderer does not discard. The Introduction's false
"DeepLC and LightGBM stages are stochastic and unseeded" replaced with the decoy
shuffle and unseeded hyperparameter search, matching §3.2.

**FDP formula** promoted to a display equation with its symbols defined, and the
exchangeability assumption of the size correction now stated rather than implied.

**Mechanical sweep (§9) — all six house-style rules now PASS** (from 3 failing):
em-dashes 4→0, semicolons 6→0, sentences >40 words 4→0. Also fixed the
`$20{,}647$` brace escape (was rendering as `20{, }647`), the three `$m/z$`
built-up fractions, and removed the global `#set raw(lang: "rust")` that was
syntax-highlighting file names.

**Verified:** compiles with zero warnings, 39 tests pass, 0 uncited bib entries,
29 references, 37 cross-references.

**Still open and unchanged:** the abstract is still 301 words (the trim is §7, and it
must follow the §2 claim corrections); every results number is still provisional
pending the §12 re-run.

---

## 1. Code and docs: align with reality

Do this **before** rewriting Methods (§3), so Methods is written against a truthful
config listing.

- [ ] `example_config.toml:4` — `missed_cleavages = 2` contradicts the code default `0`
      (`crates/uno-core/src/database.rs:194`). Fix the value, or stop calling the file
      a listing of defaults. Note the code default **agrees** with the paper's
      Limitations claim; the documentation is what is wrong.
- [ ] `example_config.toml:16` — `variable_mods = { M = [15.994915] }` contradicts the
      code default (empty, `database.rs:191`). As published it implies the benchmark
      searched Met oxidation, which it did not.
- [ ] Write a generator that emits `paper/si/defaults.toml` from `DigestConfig::default()`,
      `MatchingConfig`, `ScoreConfig::default()` and `ProteinConfig::default()`.
      Current file lists 19 of 44 defaulted settings; missing include `ml_rescore`
      (default `true`), `decoy_mass_shift`, `decoy_shuffle`, `decoy_method`,
      `decoy_order`, `min_charge`, `n_folds`, `calibration_fdr`, `bootstrap_rounds`,
      `iterations`, `max_depth`, `shrinkage`, `max_train_rows`, `seed_min_runs_frac`,
      `mass_anchor_ppm`, and the whole `[protein]` table.
- [ ] Restore the `top_k` tuning comment that `paper/si/defaults.toml:24-27` elides
      relative to `example_config.toml:22-32` (or handle via D2).
- [ ] Remove or wire up the dead config surface: `ondata_rt`, `ondata_rt_finetune`,
      `ondata_rt_lambda` (read zero times), `discriminant`, `calibration_fdr`,
      `bootstrap_rounds`, `iterations`, `max_train_rows`, `seed_min_runs_frac`,
      `rt_calibrated`, and the entire `[protein]` table (deserialized at
      `crates/uno/src/main.rs:119`, never passed to `score()`).
- [ ] Document the ~40 result-affecting `UNO_*` environment variables in a `[dev]`
      section or the SI: `UNO_ENS_HI/LO/STEPS`, `UNO_CUT_MAX/STEP`, `UNO_RT_STAGES`,
      `UNO_RT_SIG`, `UNO_NO_FDR_PLUS1`, `UNO_ML_LIGHTGBM`, `UNO_DYNAMIC_P`,
      `UNO_CONTINUOUS`, `UNO_STAGE_MASS`, `UNO_ML_*`, `UNO_PROT_PERMUTE`, `UNO_RT_CACHE`.
      Each is a way for a reproduction attempt to differ silently from the paper.
- [ ] Add a determinism regression test: two runs on one feature file, byte-diff the
      output. Add minimal CI (`.github/` does not exist).
- [ ] Fix the three `#[test]`s that silently `return` when hard-coded `/home/...` inputs
      are missing (`database.rs:403-406`, `input.rs:380-383`, `matching.rs:175-178`),
      and make `matching.rs:168` assert something beyond the ppm bound.
- [ ] `benchmark/quant/dilution_track.py:45-47` — `species()` treats any accession
      containing the substring `UPS` as UPS1, so `sp|Q9UPS8|ANR26_HUMAN` would count as
      spike-in. Latent only (no published number affected); match `ups|` instead.

### 1b. Stale in-code documentation a reviewer will read

- [x] `crates/uno-core/src/rt.rs:5-8` — claims the Klammer set is "63 features";
      `rt.rs:22` says 69-d; `RT_N_FEATURES` is 99. Reconcile all three.
      **Done 2026-07-30:** 63 is the only count with arithmetic behind it (the
      header enumerates it); `rt.rs:24` now says 63, and the paper drops the count
      entirely in favour of describing the encoding. See SUBMISSION_TODO C7.
- [ ] `crates/uno-core/src/rt_cnn.rs:1-4` — says trained on "ProteomeTools/Prosit";
      it is Chronologer. The paper is right, the source header is wrong.
- [ ] `crates/uno-core/src/score.rs:117-121` — doc comment says `ml_rescore` is
      "OFF by default"; `:156` sets it `true`.
- [ ] `crates/uno-core/src/score.rs:1297-1307` — comment documents "top-30% (lo=1, hi=3)"
      of 10 bands; code ships `n_steps=40, lo=1, hi=14`. This is the only documentation
      of the ensemble's operating point.
- [ ] `crates/uno-core/src/database.rs:59-65` — function name and header say
      "GLOBAL-POOL … pooled, shuffled once"; the implementation at `:87-161` builds the
      pool **per protein** and says so. The paper copied the stale header.
- [ ] `CHANGELOG.md:78-86` — documents "On-data RT fine-tuning … on by default" under
      Unreleased; the feature was purged in `a91c6a9`.
- [ ] `docs/uno_vs_ms1searchpy.md:11-12` — asserts `--out-peptides` and quant "have been
      deleted"; both are present at HEAD (`main.rs:80-110,331-469`).
- [ ] `PLAN.md:14-25` "Decisions (locked)" — all four are superseded (RT model is a CNN
      not `gbdt`+Klammer; scorer is deviation+binomial not GBDT discriminant; decoys are
      shuffled not reversed; v1 scope now includes quant). Mark stale or delete.
- [ ] `paper/scripts/figures/fig5_quant.py:8-34` — provenance block names
      `uno_proc/align/intensity_matrix.tsv` (box-sum) and says "this figure plots the
      DEFAULT (box-sum) matrix throughout"; `ab.tsv`/`lfq.tsv` are the **averagine**
      (`align_avg`) matrices. Also `:31` claims files are "checked in" that are gitignored.
- [ ] `paper/scripts/figures/figS5_quant_detail.py:29-30` — docstring names
      `uno_proc/align/qvalue_matrix.tsv` and "11226 LFQ rows"; code reads
      `uno_proc_newkoth/align_avg/` and `lfq.tsv` has 11520 rows.
- [x] `paper/scripts/README.md:9,12,174` — three rows contradict the justfile
      (`just pdf` does NOT regenerate figures; `just figures` is timestamp-gated, not
      "every script"; the skip list also excludes `entrapment_data.py` and `conftest.py`).
      Table also omits `figures-all`, `all`, `report`.
      **Done 2026-07-30:** the file moved to `paper/analysis/README.md` and was
      rewritten against the new `analysis/justfile`. Its `paper_stats.py` half is
      gone with that script; the manuscript tooling is documented in
      `paper/README.md` and `paper/CLAUDE.md`.
- [x] `paper/FIGURE_MANIFEST.md:46,51-53,57` — the five legacy PNGs are at repo-root
      `figures/`, not `paper/figures/`; the `paper/figdata/faims.tsv` fallback does not
      exist; record that optional `figS6_rt_model.png` was dropped.
      **Done 2026-07-30**, all three, plus the script paths (now
      `paper/analysis/scripts/`) and the D12 caveat on "no separate
      supplementary.pdf".

---

## 2. Wrong claims and numbers

**Do all of §2 before any length trimming (§8), or the abstract gets trimmed twice.**

### 2a. FDR calibration

- [?] (D10) `paper.typ:35` (abstract), `:442`, `:471-472` (Fig 4 caption), `:812`
      (Conclusions) — "conservative on five of the six datasets and about twice
      optimistic on one". uno's own-1% values are 2.16, 1.03, 0.78, 0.67, 0.27, 0.00%:
      **two** exceed nominal 1%, not one. The same paragraph counts ms1searchpy's 1.04%
      as optimistic, which is a double standard against the comparator.
      Evidence: `paper/si/entrap_reported.typ`, `paper/figdata/entrapment_curves.tsv`.
- [ ] `paper.typ:445-447` — "the same two uno datasets are again optimistic, at 6.65%
      and 5.17%". Not the same two: at 5% they are 09746 and **09766** (09778 is 3.61%,
      conservative). Also "the same two" has no antecedent, since the preceding sentence
      says "one". Name the datasets.
- [ ] `paper.typ:38` (abstract), `:474` (Fig 4 caption) — "by a smaller margin each
      time" / "less severely in any single run". True only as an extreme-value statement
      (ms1 max 2.07% < uno max 2.16%). Per dataset, **13 of the 18** ms1searchpy runs
      have a higher empirical FDP than uno on the same file. Reword, and state plainly
      that uno's worst case is the worse of the two.
- [?] (D6) `paper/si/entrap_reported.typ:8` — the column headed "Entrapment FDP (%)"
      holds the **running-minimum entrapment q**, not the FDP that Methods `:283`
      defines. A reader recomputing from the printed Human/Yeast columns gets 2.24%
      where the table prints 2.16%. Knock-on if switched to raw FDP: uno range
      0.00–2.24%, median 0.744%; ms1 worst 2.11%, median 1.25%; at 5% 6.87/5.19% and
      the other four 2.47–3.64%. Six sentences move either way.
      Evidence: `paper/scripts/figures/entrapment_data.py:114-143`.
- [ ] `paper.typ:914-916` — "A single yeast hit therefore moves an individual panel by
      **more than** the entire difference between the two engines' medians". It is
      **less**: one hit is worth 0.17–0.47 pp, the median difference is 0.478 pp. No row
      reaches it. Change to "comparable to".
- [ ] `paper.typ:762-767` — "0.17 to 0.41 percentage points". The smallest accepted
      human count in `@tbl:entrapreported` is **723** (09746 ms1 rep 3), giving 0.47 pp,
      not 0.41. The sentence also takes the yeast maximum (8) from ms1 rows while taking
      the human range from uno rows only. State 0.17–0.47 over the whole table, or say
      "over uno's own lists". The provenance comment at `:766-767` also pairs the
      endpoints the wrong way round.
- [ ] `paper.typ:454-456` — "ms1searchpy's internal 1% count was 1319 proteins but only
      1220 **of those** survive". Three defects: the 1319 includes 5 yeast, so the
      human-comparable count is **1314**; 1220 is an independent cutoff, not a subset
      (one of the 1220 is not on the internal list); and both figures are **rep 1 only**,
      which the sentence does not say.
- [ ] `paper.typ:435-437`, `:452-454` — two claims say "on one 5-minute run" without
      naming it. Both check out: 830-vs-1035 is **09734**, the 0.39–1.62% replicate
      range is **09766**. Name them. (Note the Fig 4 caption already names 09746 while
      the body withholds it — pick one convention.)
- [ ] Add a 5% block to `paper/si/entrap_reported.typ` so the four percentages at
      `:445-447` are checkable. Data already on disk (`entrapment_curves.tsv`,
      `cutoff_type=reported_q`, cutoff 0.05).

### 2b. Quantitation

- [ ] `paper.typ:302` — "24 runs on a **Q Exactive HF**". It is an **Orbitrap Fusion
      Lumos with a FAIMS Pro interface** (CVs −50/−65/−80 V). Evidence: the mzML
      `instrumentConfiguration` under
      `/mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc/mzml_ab/`, the `LUM2_`
      prefix on all 24 run names, PRIDE's PXD033350 record, and the DirectMS1Quant
      Methods. Found independently by two lanes.
- [ ] Consequence of the above: PXD033350 is **also** FAIMS data, which bears on
      `paper.typ:750-752` ("the FAIMS run of Section 3.4 is read without its
      compensation-voltage structure, so that dimension is likewise unexploited").
- [?] (D4) `paper.typ:773-774`, `:793`, `:815-817` — "calls changed proteins with high
      specificity" / "high differential specificity". Contradicted by the paper's own
      Figure 5c, whose legend prints **60 known-absent human false positives**. Measured:
      22 UPS1 TP, 23 *E. coli* FP, 60 human FP, 105 calls, **precision 21%** — which
      `benchmark/quant/results/analysis_stdout.txt` itself records as `precision=21%`.
- [ ] `paper.typ:690-693` — the "42 changed UPS proteins" comparison is
      threshold- and database-mismatched. The value appears nowhere in the
      DirectMS1Quant text; it is readable only off **Figure 2a**, at **5% identification
      and 5% quantitation FDR**, against a **combined E. coli+UPS database with no human
      probe**. uno reports at BH 1% with a human probe. State both thresholds and both
      databases, and attribute the number to the figure.
- [ ] Related: the cited deposit's own `directms1quant_out_50_vs_0_1.tsv` is 42 rows,
      **all 42 UPS1, zero *E. coli*, zero human — 100% precision**. Any head-to-head
      framing must account for that.
- [ ] `paper.typ:693-694` — **delete** "The gap is not in the differential test, since
      the per-protein measurement count and statistic match". Neither matches: median
      uno n / DirectMS1Quant `n_arr` = **1.33** (uno's n is *larger*), median
      |log2 FC| ratio = **0.60**. Nothing on disk performs this comparison. Replace with
      the version that is true and supports the paragraph: measurement count cannot
      explain the gap, the ~40% fold-change compression can.
- [ ] `paper.typ:680-682` and Fig 5 caption `:732-733` — qualify "spike-in abundance
      tracks the dilution". The UPS1 median trace is **flat from 0.1 to 5 fmol**
      (−2.79, −3.08, −2.65, −2.61, −2.62 log2 vs 50 fmol) against an ideal −8.97 at
      0.1 fmol: a 50-fold concentration range with no response. The dose response exists
      only from 5 to 50 fmol. Evidence: `benchmark/quant/results/fig5_stdout.txt:2`.
- [ ] `paper.typ:697-700` and Fig S5 caption `:1019-1024` — the two detection traces use
      **different denominators** presented as one population: 82% is 151/**184**
      (precursors whose cell is filled at 50 fmol), ~1% is 1/**82** (precursors
      q-confident at 50 fmol). The figure legend prints both n; the prose does not.
      Evidence: `figS5_quant_detail.py:156-157`, `figS5_stdout.txt:5`.
- [ ] `paper.typ:700` — "falls to a few percent" → it falls to **zero** at three of the
      eight levels (2%, 0%, 0%, 0%, 1%). Say "zero or near zero".
- [ ] `paper.typ:707-713` — rescope the `min_scans = 1` collapse. The quoted numbers
      (p 0.15→0.51, no proteins accepted at 1%) hold at `--matrix-min-runs 3`. At
      **min-runs 6**, the floor used for every headline number in this section,
      `min_scans = 1` gives p = 0.118 and **584 proteins vs 119** — five times more, no
      collapse. `benchmark/quant/results/minscans_ab/README.md:179-186` drafts wording
      reporting **both** arms; the paper adopted only the collapse half. Report the
      strict-floor result too (precision 50%→29%, entrapment FDP 0.00%→0.23%).
- [ ] `paper.typ:707-708` — "Of the precursors lost, the signal is physically present
      but sub-threshold" is unqualified. **42%** of the 113 lost precursors reappear in
      the `min_scans = 1` consensus; 58% do not appear even there; only 2% recover as
      independent detections. Evidence: `minscans_ab/recovery.txt`.
- [ ] `paper.typ:683-684` — "the quantified human probe stay[s] flat at median rho of
      +0.25" rests on **n = 5**, and +0.25 is 12× the *E. coli* background it is
      bracketed with. Either add the argument that makes it defensible (+0.25 sits at the
      ~77th percentile of the 303 *E. coli* rho values, so it is consistent with the flat
      null) or state plainly that the absent-species check rests on five proteins.
- [ ] State `--matrix-min-runs 6` in Methods as a **requirement**, not an SI detail:
      without the floor the mass-anchor background rises and the search returns nothing.
      Currently it appears only in `paper/si/commands.typ:37`.
- [ ] State that `averagine_projection = true` was passed explicitly — it is **not**
      koth's default (`../koth_rust/koth_ff/src/lfq/mod.rs:221` has `false`).
- [ ] State the quant digest. The quant run passed **no `--config`**, so it used code
      defaults (`max_len = 50`, not the identification benchmark's 30). `paper.typ:268-271`
      is scoped to the ID benchmark and §3.10 never states its own digest.
- [ ] `paper/si/quant_detail.typ:7-10` — "the raw box-sum regressed the honest
      protein-ID count by 21% (306 -> 243) while averagine **held** it (305 -> 326)".
      All four numbers are correct (old/new koth × box-sum/averagine, verified against
      `uno_align.log:23`, `ab_boxsum.tsv`, `uno_align_avg.log:23`, `ab.tsv`), and the
      paper's 326 is the right one — but +6.9% is a rise, not "held", and the two
      baselines are never labelled.
- [ ] `paper.typ:1030-1031` — "n is the number of proteins quantified in at least 8 of
      the 24 runs" is inert: 18+303+5 = 326 = every row of `ab.tsv`, and the minimum
      observed is 23 of 24. As described it implies a selection that does not occur.

### 2c. Head-to-head framing

- [ ] `paper.typ:318-334` — add that uno's margin over the ms1searchpy median is
      **smaller than ms1searchpy's own replicate range on all six datasets**
      (25 v 30, 71 v 124, 155 v 181, 49 v 73, 161 v 300, 589 v 600). The paper already
      makes this point for the FAIMS side-analysis at `:496-497`. Soften the abstract
      (`:30-32`) and Conclusions (`:804-805`) accordingly.
      Evidence: `benchmark/results/summary.tsv`.
- [ ] `paper.typ:326-328` — "reproduces the 2000-protein figure the reference method
      reports for that acquisition". That figure is a **three-replicate combination**
      (2100 proteins, "after combining three technical replicates"); the reference's
      **single-replicate** FAIMS number at the matching min-3-scans setting is **1775**.
      uno's 2072 is one run, and at *entrapment* FDR rather than the reference's own
      target/decoy FDR. Compare like with like or drop "reproduces".
- [ ] `paper.typ:336-340` — report the **user-CPU factor (3.6–4.6×)** beside the wall
      factor (5.3–7.0×), state both engines' thread settings (ms1searchpy defaults to
      `-nproc 4`, `~/Repos/ms1searchpy/ms1searchpy/search.py:41`; uno uses rayon's
      default on a 20-thread host), and delete "user CPU time separated the two engines
      in the same way" — it is not true of the magnitude.
- [ ] `paper.typ:481-502`, `:806-808` — disclose that ms1searchpy was run **once per
      FAIMS replicate** (`benchmark/faims/run_faims.sh:28-38` has no replicate loop),
      against the paper's own standard at `:399-401` ("every ms1searchpy value in this
      paper is a median over three runs"). The three "replicates" are three acquisitions,
      not three engine runs. Also give uno's own spread (1460–1672, 13% of its mean) and
      drop "again" at `:495-496`.
- [ ] `paper.typ:490` — "The margin at 5% was 7.6%" is 7.54% from the exact means, or
      7.57% from the rounded means. Use 7.5%, or say it is computed from rounded means.
      (The 1% "8%" is exact.)
- [ ] `paper.typ:390-393`, abstract `:33-34` — "spanned ranges from 30 to 600 proteins,
      equal to 3.6% and 40.5%" never states the FDR threshold. These are the 1% figures;
      at 5% the spreads are 84–343 proteins, 6.1%–14.7%. Insert "at 1% FDR".
- [ ] `paper.typ:354-355` (Table 1 caption) — "its median and full replicate range are
      given" is true only of the 1% column; the 5% column is median-only.
- [ ] `paper.typ:259-262`, `:325-326` — "the deepest of the 5-minute runs" / "the
      sparsest 5-minute run" implicitly excludes the FAIMS run from the 5-minute class,
      but PXD023977's runs are also 5-minute gradients (per the cited reference's own
      title). Say "the five PXD015669 runs", or state the sixth run's gradient.
- [?] (D1) Label run10's numbers with their build. `eccf0b2` fixed a ~20% feature-row-
      order sensitivity and reports run10 **1309 → 1426**, landing after the six-dataset
      results were measured (2026-07-22). Both `1309` (`:341`) and `1306` (`:523`) come
      from the pre-fix build.
- [ ] `benchmark/results/run_console.log` records a second uno run10 measurement:
      wall 54.9 s / 3919 MiB, against `uno.time`'s 46.21 s / 3887 MiB. The paper quotes
      the faster one (`:342`). Either report the repeat or state that resource figures
      are single measurements. Note the two runs returned **identical depth**
      (1309/1887), which is independent corroboration of the determinism claim and could
      be cited as such.

### 2d. Ablations, RT, LDA, performance

- [ ] `paper.typ:524-526` — "The ms1searchpy median falls from 1155 to 1062 at 1%, a
      drop **within its own replicate spread**". 1062 is *below* the with-columns range
      1079–1220. The source argument is different: the with-reps (1079–1220) **overlap**
      the without-reps (1023–1134) (`benchmark/RESULTS.md:62-64`). Reword.
- [ ] `paper.typ:175-179` — "The remaining **13** are scalar descriptors" then
      enumerates **10**. The code has 13 (`crates/uno-core/src/rt.rs:113-125`); the three
      missing are the mean, maximum and minimum Kyte–Doolittle values.
- [x] `paper.typ:541-543`, `:562` — "the prior **69-d** Klammer set". The code documents
      63 (`rt.rs:6-8`) and says 69 three lines later (`rt.rs:22`); Klammer et al. 2007
      states no feature count at all. Reconcile, and cite Klammer for the
      SVR-composition *approach* rather than for a dimensionality it never gives.
      **Done 2026-07-30:** count removed from both prose and `@tbl:rtdepth`; the
      prose now describes the encoding and cites Klammer for the approach.
- [ ] `paper.typ:579-582`, `@tbl:reverse` — "captured 72% (Pearson r = 0.68)" conflates
      two different quantities (the table has separate "order-Δ corr." = 0.68 and
      "order magnitude" = 72% columns; r² = 0.46, so no stated transform links them).
      Define "captured" explicitly as the magnitude ratio and stop presenting r as its
      parenthetical. Rename the "order magnitude" column, which reads as "order of
      magnitude".
- [ ] `paper.typ:571-572` vs `@tbl:reverse:593` — prose says a composition model is
      "mathematically forced to predict the same value for both"; the table gives
      composition-only R²_rev = 0.95. Reconcilable (forward and reverse DeepLC RTs are
      themselves correlated) but unexplained, so the row reads as contradicting the
      prose. Add one clause.
- [ ] `paper.typ:547-548`, `@tbl:rtdepth` — "restored identification depth to DeepLC
      parity" rests on "~180 / ~152" against DeepLC's exact "176 / 152". An approximate
      value cannot support a parity claim, and `paper/uno_jpr.md:239-241` says these
      counts are "noise-dominated by the on-data calibration step; reference-overlap
      (~152) is the robust metric". Report exact counts or make the claim on overlap only.
- [ ] `paper.typ:604-613` vs `@tbl:entrap:635-638` — prose says 218 → 242 → 223, the
      table says 201 / 208 / 209. **Not an error**: they are two different searches
      (human-only DB vs human+yeast entrapment DB, per `paper/uno_jpr.md:191,296,328`).
      Add one distinguishing clause; as printed a reader sees a contradiction on one page.
- [ ] `paper/si/superseded.typ` row 2 vs `paper.typ:604` — S8 retires "uno protein count
      on the same feature set: 175–179" while §3.8 uses **218** as the
      same-configuration baseline. Per `uno_jpr.md:191` these are different koth feature
      sets. Distinguish them in one of the two places.
- [ ] `paper.typ:646-647` — "about **99** (RT, mass) cells". The only retained grid
      definition (`git show e16693b:crates/uno-core/src/protein.rs:1527-1541`) gives
      12 × 8 = **96** at shipped defaults; the code's own comment says "~100".
- [ ] `paper.typ:648-650` — "up to 17.7 million … per cell … roughly 1 billion across
      the grid". 17.7e6 × 99 = **1.75 billion**. Use "of order a billion" or give the
      measured total. Also `:648` "17.7 million" and `:655` "about 18 million" are the
      same measurement rounded two ways in adjacent paragraphs.
- [ ] `paper.typ:644` vs `:662` — "profiled … using 20 threads" vs "fell from about 40 s
      to about 1 s **single-threaded**". The pre-optimization grid was serial
      (parallelism was added *by* this work), so "using 20 threads" mis-describes the
      43 s baseline. State which figure is which.
- [ ] `paper.typ:518-519` — "the value stays at 0.50 even in the full unreduced
      41-million-match pool". `benchmark/RESULTS.md:79-80` records that check for
      **`cosine_score` only** (0.5005). Name the column.
- [ ] `paper.typ:506-509` — "target and decoy peptides occur at approximately 50:50"
      among matches with RT deviation < 0.01 min: no count, no interval, and the
      diagnostic scripts were deleted. Directionally supported for ms1 on run10
      (decoy fraction 41.6% against a 45.7% pool base rate, odds ratio 1.18), but
      unauditable for uno as written.
- [ ] `paper.typ:506-507` vs `:527-528` — "mass and RT cannot separate target from decoy
      peptides" and, 20 lines later, "What does separate target from decoy is retention
      time". Both defensible; reads as self-contradiction. Qualify the first as
      "cannot separate them decisively".
- [?] (D14) `paper/si/quality_auc.typ` — two problems. The headline rows (RT distance
      0.136 uno / 0.080 ms1) match **neither** plausible definition of "confident target"
      when recomputed from the surviving `run10.features_PFMs_ML.tsv`: 0.026–0.030 using
      all targets, 0.171–0.188 using the `qpreds==0` bin. And `cosine_score` /
      `isotope_score` are given as *exactly equal* for both engines, though the two
      columns correlate r = 0.31 with **zero** row-wise ties. The qualitative conclusion
      ("the columns are a non-lever") is robust; these specific values are not.
- [ ] `paper/si/quality_auc.typ` — rows 1–5 are raw AUCs and rows 6–7 are |AUC − 0.5|,
      in the same two unlabelled columns. Add a metric column or express all seven the
      same way (0.021, 0.021, 0.026, 0.003, 0.000, 0.136, 0.171).
- [ ] `paper/si/quality_ablation.typ` — the "ms1searchpy median, columns stripped" row
      mixes replicates (1062 @1% is rep 1, 1447 @5% is rep 3). Standard for per-threshold
      medians; note it in the caption.
- [?] (D9) `paper.typ:642-667` — §3.9 optimizes a 2-D (RT, mass) threshold grid that
      **no longer exists in the shipped engine**. `ProteinConfig::{rt_window,
      mass_z_window, optimize_thresholds}` are declared but consumed nowhere; the shipped
      path is a 1-D deviation sweep (`score.rs:766-790`). So "The wall-clock figures in
      Section 3.1 are measured after these changes" (`:666-667`) is true only
      chronologically. The section also never says its profiled input is a smaller,
      different 45,155-feature koth-aligned matrix, so its 3.4 s collides with §3.1's
      32–45 s. Found independently by two lanes.
- [ ] `paper.typ:309-314` — Methods promises "Their scope is stated where they are
      reported" but lists only §3.6, 3.7, 3.8. Add §3.9 and its koth-aligned-matrix
      configuration.
- [ ] `paper.typ:533`, `:184-186` — §3.6's title claim is about the distilled MLP, which
      Methods states is **not** the model used in §3.1–3.5. Add one sentence so a
      reviewer cannot read it as validating the shipped CNN.
- [ ] `paper.typ:746-747` — "peptide-centric extraction from raw hills does not add
      identifications" is a finding from a deleted code path (`--hills`, removed per
      `docs/uno_vs_ms1searchpy.md:11-12`) and is not scoped as such, unlike §3.6–3.8.
- [?] (D13) `paper.typ:568-596`, `@tbl:reverse` — §3.7 has **no surviving script, data
      or log** (no reverse-sequence tool in `crates/uno-rt-train/src/bin/`, no artifact
      matching 46,304). Unlike `@tbl:qualityauc`, the caption does not disclose
      non-retention, and the prose calls the benchmark "novel and reusable".

---

## 3. Methods ↔ implementation

Gated on §1.

- [?] (D7) `paper.typ:196-209` — the whole "On-data RT fine-tuning" subsection describes
      code that **does not exist at HEAD**. `cfg.ondata_rt`, `ondata_rt_finetune`,
      `ondata_rt_lambda` are declared and defaulted true/true/1.0
      (`crates/uno-core/src/score.rs:88-96,144-146`) but read **zero times**; the
      implementation was deleted in `a91c6a9` ("Purge discriminant path… -2810 lines"),
      and `score.rs:990-997` is the tombstone ("was tried here … not a validated win
      anywhere"). At minimum, "On-data readout fine-tuning is enabled by default and was
      enabled in every search reported here" must go. Also remove the S8 sentence at
      `:1041-1042` that cites the RT ablation as evidence for this design.
- [ ] `paper.typ:201-203` — if the section survives in any form: "the penultimate-layer
      activations" is wrong. `RtCnn::embed` returns the masked mean+max pooled conv
      vector plus an intercept (`rt_cnn.rs:417-421`), i.e. the *input* to the MLP head,
      and it excludes the 11 global descriptors the head receives.
- [ ] `paper.typ:228-231`, `:620-621`, `@tbl:entrap` caption `:629-631` — there is **no
      LDA implementation at HEAD** (grep finds only comments). The only handle is
      `ProteinConfig::ml_rescore`, and `ProteinConfig` is dead. "Provided but disabled by
      default" implies a reviewer can switch it on; they cannot. Reword to "evaluated on
      an earlier development configuration, not present in the released engine", and say
      in the caption that `@tbl:entrap` is not reproducible from the cited artifact.
- [?] (D3) `paper.typ:211-231` — Methods describes only a deviation discriminant plus a
      binomial, but the shipped default **also runs a learned peptide-level rescorer**:
      `ScoreConfig::ml_rescore = true` (`score.rs:121,156`), consumed at `:1031`, a
      3-fold cross-validated gradient-boosted-tree classifier over ~45 PFM features
      (`:1128-1163`) with a 14-of-40-band cumulative inclusion ensemble
      (`:1305-1407`) that replaces the deviation result whenever `ml_n > best_n`.
      Default backend is **sequoia-boost** (`:1688-1708`), cited nowhere in `paper/`.
      Active in every reported search. `docs/uno_vs_ms1searchpy.md:25` states it.
      Note: this does not weaken §3.8 — the paper's own point is that learning helps at
      the peptide level and hurts at the protein level.
- [ ] `paper.typ:153-158` — "interior-residue shuffling against a **global peptide
      pool**". The pool is built **per protein**, deliberately
      (`crates/uno-core/src/database.rs:87-90`: "so the interior pool is built PER
      PROTEIN … A global pool would average organisms and make the binomial p
      under-estimate the in-organism random rate"). The code's rationale is stronger than
      the paper's; use it.
- [ ] `paper.typ:154-155` — "yields a decoy of different composition and different mass
      from any target". Nothing enforces this: `global_shuffle_decoy_db` dedups decoys
      against each other only (`database.rs:168-171`), with no decoy-vs-target sequence
      or mass collision filter, and interiors are dealt cyclically from the same
      protein's pool. Soften to "in general", or add the filter and report the rate.
      This sentence is load-bearing for the argument against reversed decoys.
- [ ] Also `database.rs:77-83` — decoy masses use `mass::residue(aa).unwrap_or(0.0)`, so
      a non-standard residue silently contributes **0 Da** to a decoy, while the target
      path falls back to Sage's f32 mass.
- [ ] `paper.typ:221-224` — correct the grid description **and** its FDR defence.
      Shipped: a 1-D sweep over a single combined deviation cutoff
      `dev = mass_z² + (rt_diff/rt_sig)²`, grid 1.0–8.0 step 0.5 = 15 cells
      (`score.rs:766-775`); the mass gate and `rt_sig` come from monotone recalibration,
      not a grid. And "mirroring ms1searchpy's `worker`" — `worker` is a **function** at
      `~/Repos/ms1searchpy/ms1searchpy/main.py:2160`, not a file, it is not a mass/RT
      grid, and it **averages** its 10 passes rather than taking an argmax.
- [?] (D2) `paper.typ:221-224` — replace the defence "Because decoys pass the identical
      filter at every grid point, the FDR estimate stays calibrated regardless of which
      cutoff maximizes yield". The objective is `nt = count(!decoy && q <= 0.01)`
      (`score.rs:839`) and the winner is `if nt > best_n` (`:865`), i.e. an argmax over
      cells on a noisy decoy-derived estimate. The code's own comment says what the paper
      denies — `score.rs:756-765`: "the raw best-cut-by-target-count selection **over-fits
      to a loose, decoy-heavy high cut**. Capping at 8 drops them." `cut_max = 8.0` is an
      undisclosed regularizer against exactly this bias, and `9064d0b` shows it was
      chosen by measuring entrapment depth on deephela and run10.
- [ ] `paper.typ:217-224` — state the FDR convention. Shipped protein FDR is
      **(d+1)/t** with a running minimum (`protein.rs:1022-1036`); `UNO_NO_FDR_PLUS1`
      switches to d/t; `score_proteins_greedy` uses plain d/t (a silent inconsistency
      between scorers). ms1searchpy's reported list also uses (d+1)/t, so the answer
      favours the paper.
- [ ] `paper.typ:757-760` — reattribute the TMT limitation. It blames DeepLC's
      r ≈ 0.71, but uno does not use DeepLC at search time and its shipped RT models
      **cannot represent a modification at all**: `rt_features` takes a stripped sequence
      (`rt.rs:64-66`) and the CNN tokenizes only the 20 standard residues
      (`rt_cnn.rs:25`). A TMT peptide gets the unmodified peptide's predicted RT.
- [ ] `paper.typ:748-750` — "uno predicts collisional cross-section via an analogous
      distilled model where mobility is present". The model and loader exist
      (`crates/uno-core/src/im.rs`, `models/im_im2deep.gbdt`, `--im-model`) but
      `score()` names the parameter **`_im_model` and never uses it**
      (`score.rs:1911`), `featurize`'s `im_diff` is called with a hard `None`
      (`:1129`), and `ProteinConfig::im_window` is dead. Reword to "trained and loadable
      but not yet wired into scoring".
- [?] (D11) `paper.typ:226-228` — "No result-affecting path uses process-randomized
      iteration order or an unseeded RNG". Two surviving channels:
      (a) `est_rt_sig` (`score.rs:547-586`) builds its residual vector by iterating a
      `std::collections::HashMap` (`:369`) and reduces it with `clipped_sigma`, whose
      loop does **unsorted sequential f64 summation** (`:1788-1789`); `rt_sig` feeds the
      `dev <= cut` acceptance. The same bug class was fixed for the calibration pairs
      and for RANSAC but not here.
      (b) Quantitation: `resolve_feature_competition` collects over a `HashMap`
      (`quant.rs:52,64`), so `--out-lfq` row order is process-random and
      `group_abundance` sums f64 in that order (`:247-260`), so `--out-abundance` values
      differ in the low bits between runs. `group.rs:89` has an order-dependent fallback.
      Fix requires re-running (gates D1); or soften the claim.
- [ ] `paper.typ:405-409` — "sorting the calibration seed on **feature index**". Both
      fixes sort on **value**: `fit_ransac` sorts by `(pairs[k].0, pairs[k].1, k)`
      (`calibrate.rs:164-174`) and `build_cal_pairs` sorts by `(abstract_rt,
      observed_rt)` (`score.rs:955`); the index is only a final tiebreak. The second
      fix ("resolving ties on accession") is accurate.
- [ ] `paper.typ:163-165` — "reduced to the top-K best-mass candidates **per feature**".
      `k` is applied per feature **per class**: targets and decoys are truncated
      separately (`matching.rs:150-155`), so up to 100 candidates per feature at
      `top_k = 50`, as `example_config.toml:23` itself says.
- [ ] `paper.typ:162-163` — "isotope-error offsets (0, ±1 **Da/charge**)". The offset is
      applied to the **neutral mass** in units of the C13–C12 spacing:
      `center = feature.mass - offset * NEUTRON`, `NEUTRON = 1.003354835`
      (`matching.rs:71,123`, `mass.rs:13`).
- [ ] `paper.typ:118-123` vs `:978-980` — two different, both-incomplete column lists.
      Methods has `nScans` and no `neutron_offset`; SI S6 has `neutron_offset` and no
      `nScans`. **`neutron_offset` is the correct member of the quality set**
      (`input.rs:9-10,54-57`, `benchmark/RESULTS.md:49-50`,
      `paper/figdata/quality_columns.tsv`). uno requires only four columns
      (`massCalib, mz, charge, rtApex`); both lists omit `nIsotopes`, `isoerror`,
      `isoerror2`, `rt_width`, which uno reads and at least two of which are
      load-bearing. ms1searchpy additionally **requires** `nIsotopes` and `nScans`.
      Make the two lists identical and mark required vs optional per engine.
- [ ] `paper.typ:150-152` — "fixed/**variable** modifications" delegated to sage-core.
      With variable mods enabled the decoy path breaks: decoy mass is computed from the
      stripped sequence plus **static** mods only (`database.rs:74-83`), and one decoy is
      keyed per unique *sequence* (`:149,159`), so N modforms of a target yield 1 decoy,
      unbalancing the target/decoy ratio. Either note that variable mods are unsupported
      on the decoy path, or drop `max_variable_mods` from the sentence.
- [ ] `paper.typ:237-238` — "uno does not perform its own feature extraction". At HEAD it
      can: `--raw` runs koth's finder **in-process** (`main.rs:26-32,140-147`,
      `crates/uno/src/koth_ingest.rs`, `koth_ff` path dep), which is what `7f2cece`
      reworked. Reword to "does not implement its own; where it reads raw data it embeds
      koth's finder, and every result here came from precomputed biosaur2 files".
- [ ] `paper.typ:239-241` — profile reassignment is **opt-in** (`--profile-reassign`),
      not part of an unconditional "also". It was passed for the reported run.
- [ ] `paper.typ:241-248` — the two per-cell intensity estimators live in **koth**, not
      uno. Attribute them, and show the `[lfq]` config line in S1.
- [ ] `paper.typ:248-249` — "uno emits the protein-group × sample abundance matrix and
      stops there". It also emits a precursor × sample LFQ matrix (`--out-lfq`) and a
      per-peptide table (`--out-peptides`) — and the paper's own differential test
      consumes the LFQ file.
- [ ] `paper.typ:185-186` — "This is the shipped default". The CLI has **no** default
      `--rt-model` (`main.rs:52-54`; without it RT is unused). The model was passed
      explicitly in all three benchmark invocations. One word would make this exact.
- [ ] `paper.typ:270` — note that the digest mass bounds `[500, 5000]` are applied
      **after** modifications (`../sage/crates/sage/src/database.rs:197-204`). Correct
      behaviour for an MS1 engine, worth one clause.

---

## 4. Disclosures to add

- [?] (D2) Tuning-provenance paragraph. Four shipped defaults were selected by
      maximizing honest yeast-entrapment depth on datasets the paper then evaluates on:
      - `top_k = 50` — `example_config.toml:24-26` ("25→105 IDs, 50→208, 100→275, all at
        ~1% entrapment FDP"), `CHANGELOG.md:3-15`. Same probe, sample type and metric as
        the paper. The manuscript's only stated rationale is cost (`:163-165`).
      - the ensemble inclusion band (`UNO_ENS_STEPS=40, LO=1, HI=14`,
        `score.rs:1305-1307`) — commit `683333b` ("run10 775 → 1218 @1%; deephela
        1585 → 1752"), `184d8b5` ("run10 @1% 1218 → 1448").
        `docs/uno_vs_ms1searchpy.md:187-189` calls it "uno's main depth lever".
      - `min_child_weight = 80` (`score.rs:1574`) — commit `79725a6` ("run10 1218;
        deephela 1954 (+11.5%)… chose the moderate 80, not the run10 peak at 200").
      - `cut_max = 8.0` (`score.rs:767`) — commit `9064d0b`, chosen on deephela and run10.
      `paper.typ` contains no occurrence of `top_k`, `sequoia`, `min_child_weight`,
      `mass_anchor_ppm` or "tuned". run10 must be labelled a **tuning set**, not a
      benchmark. Counterweight worth including: `models/rt_cnn_README.md:14-20,74-93`
      records that the shipped RT model was chosen *against* its run10 score (248 vs 401
      for the 5-global variant), which is evidence of restraint.
- [ ] ms1searchpy patch disclosure. The comparator is a **local uncommitted patch**
      (+45 lines in `main.py`, +6 in `utils.py`), not released 2.8.12 — which does not
      exist on PyPI (latest is 2.8.9) and has no GitHub tag past v2.7.3. Stock 2.8.12
      reads none of the five quality columns (`utils.py:116-128`); `KOTH_EXTRA_COLUMNS`
      (`main.py:49-52,577-580`) is what makes §3.5/S6 possible. The patch's own comment
      claims those columns are "Absent for plain biosaur2 input → behavior is unchanged",
      which is **false for the benchmark's own files**. Pin the commit, ship the diff in
      the SI, reword `:980`, and note the patch gives ms1searchpy *more* features, so
      uno's margin is conservative.
- [ ] Entrapment-estimator assumptions (`:278-289`, `:762-765`, `:786-789`): state
      exchangeability explicitly, and report the scale-factor sensitivity. Entry-count
      scale is 3.404 (what `benchmark/entrapment.py:22` uses); peptide-count scale is
      3.711, so the estimator understates false human acceptances by ~9%. In the paper's
      favour: peptide mass distributions match closely (median 1354 vs 1359 Da).
- [ ] Equivalence bounds or explicit downgrades for the untested null claims: `:523`
      (1309→1306), `:525-526`, `:516-521`, `:968-973` (Fig S3 "scale smoothly … constant
      in ratio"), `:683-685`, `:685-687` ("stable across koth versions and both
      estimators" — no cross-version numbers appear anywhere), `:693-695`, `:685`
      ("cleanly separated").
- [ ] Spearman inference: no p-value or CI is attached to any rho, and the resampling
      unit is never defined (24 runs over 8 distinct x-levels, per-protein rho then
      medianed across proteins with peptide-within-protein sharing ignored). Measured
      bootstrap CIs: UPS +0.712 [+0.630, +0.766]; *E. coli* +0.024 [+0.003, +0.068];
      human +0.247 [−0.218, +0.405]. "89% above 0.5" is 16/18, Wilson CI [67%, 97%].
- [?] (D4) Differential test: state the decision rule (BH at **1%** plus a 2-fold
      `|mean log2 FC| >= 1` gate, `benchmark/quant/peptide_fc_de.py:98-107`), and
      disclose that BH runs **after** the gate over 852 of 2228 proteins — filtering on
      the numerator of the statistic whose p-value is then BH-adjusted. Also that the
      test runs on the unfiltered precursor matrix: **1939 of 2265 accessions (86%) are
      not in uno's 326-protein 1% FDR set**, including 1172 of 1177 human entrapment
      accessions, which is why 60 human proteins are called changed. Two further
      defects: missing low-condition cells are imputed to a single global 1st-percentile
      floor (`:99`), and peptide × replicate-pairs are treated as iid.
- [ ] Note that the DE contrast is non-monotone across levels (TP counts 22, 22, 12, 22,
      15, 4, 0 for lo = 0.1, 0.25, 1, 2.5, 5, 10, 25) — direct evidence against "the gap
      is not in the differential test". Only 50-vs-0.1 is reported, which is the correct
      pre-registered extreme contrast, so this is not cherry-picking; but no multiplicity
      family is stated.
- [ ] `benchmark/entrapment.py:57` — score ties are broken by file row order, and the
      two engines' files have different row orders. uno's binomial score is quantized so
      ties are common. Make it deterministic (`key=lambda r: (-r[0], r[1])`) and re-run.
      Low expected magnitude; affects both engines.
- [ ] `benchmark/quant/peptide_fc_de.py:88-90` — the contrast is asymmetric
      (`>= 2` non-zero **high** replicates only), so proteins higher at the dilute end
      are structurally excluded. Correct for a spike-down design; report it.

---

## 5. Citations and prior art

Additions (all verified against Crossref / Europe PMC / PubMed by the audit):

- [ ] **Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U**, "Assessment of
      false discovery rate control in tandem mass spectrometry analysis using
      entrapment", *Nature Methods* 2025;22(7):1454-1463,
      doi:10.1038/s41592-025-02719-x, PMID 40524023. **Highest-priority omission**: the
      entrapment estimator carries the paper's entire depth axis and is currently
      uncited. State which of that paper's estimator families the formula at `:283`
      corresponds to. Cite at `:278-289`, `:106-109`, `:786-789`.
- [ ] **Madej D, Lam H**, "On the use of tandem mass spectra acquired from samples of
      evolutionarily distant organisms to validate methods for false discovery rate
      estimation", *PROTEOMICS* 2024;24(15):2300398, doi:10.1002/pmic.202300398 — a
      viewpoint directly on the two-proteome design used here. Cite when justifying the
      yeast / *E. coli* choice.
- [ ] **Madej D, Lam H**, "Query Mix-Max Method for FDR Estimation Supported by
      Entrapment Queries", *J Proteome Res* 2025;24(3):1135-1147,
      doi:10.1021/acs.jproteome.4c00744 — same journal, same operation.
- [ ] **Xiao Z, Tüshaus J, Kuster B, The M, Wilhelm M**, "SWAPS: A Modular Deep-Learning
      Empowered Peptide Identity Propagation Framework Beyond Match-Between-Run",
      *J Proteome Res* 2025;24(4):1926-1940, doi:10.1021/acs.jproteome.4c00972,
      PMID 40052690. MS1-driven identification in the target journal that **critiques
      DirectMS1 by name**. Cite in the Introduction.
- [ ] **Zimmer JS, Monroe ME, Qian W-J, Smith RD**, "Advances in proteomics data analysis
      and display using an accurate mass and time tag approach", *Mass Spectrom Rev*
      2006;25(3):450-482, doi:10.1002/mas.20071, PMID 16429408 (and optionally
      **Pasa-Tolić L et al.**, *BioTechniques* 2004;37(4):621-639,
      doi:10.2144/04374rv01). The AMT-tag lineage is the intellectual precursor of
      MS1-only ID; the current framing at `:69-81` reads as claiming DirectMS1 invented
      mass+RT identification.
- [ ] **Ivanov MV et al.**, "MS/MS-Free Protein Identification in Complex Mixtures Using
      Multiple Enzymes with Complementary Specificity", *J Proteome Res*
      2017;16(11):3989-3999, doi:10.1021/acs.jproteome.7b00365, PMID 28905631. Cite in
      the Introduction **and** at `:754-758`, where the paper's "expanding search breadth
      is counterproductive" claim makes the opposite move to this group's published
      finding. Rewrite that Limitation to engage with it or scope it to single-enzyme
      single-run search.
- [ ] **Ivanov MV et al.**, "Improving the Protein Inference from Bottom-Up Proteomic Data
      Using Identifications from MS1 Spectra", *J Am Soc Mass Spectrom* 2021;32(8),
      doi:10.1021/jasms.1c00061, PMID 33900766. Cite at `:217-224`.
- [ ] **Ivanov MV et al.**, "Reanalysis of DIA Data Demonstrates the Capabilities of
      MS/MS-Free Proteomics…", *J Am Soc Mass Spectrom* 2024,
      doi:10.1021/jasms.4c00134, PMID 38938158. Cite so the reference engine is described
      at its current state.
- [ ] **Yu F, Haynes SE, Nesvizhskii AI**, "IonQuant Enables Accurate and Sensitive
      Label-Free Quantification With FDR-Controlled Match-Between-Runs", *Mol Cell
      Proteomics* 2021;20:100077, doi:10.1016/j.mcpro.2021.100077, PMID 33813065. The
      existing FDR-controlled treatment of the exact problem §3.10 rediscovers. Cite at
      `:698` and `:779`.
- [ ] **Cox J et al.** (MaxLFQ), *Mol Cell Proteomics* 2014;13(9):2513-2526,
      doi:10.1074/mcp.M113.031591, and **Ammar C et al.** (directLFQ), *Mol Cell
      Proteomics* 2023;22(7):100581, doi:10.1016/j.mcpro.2023.100581. Position the
      feature-competition / profile-reassignment method at `:233-254` against standard
      LFQ, and say what differs.
- [ ] **Petritis K et al.**, "Improved peptide elution time prediction for reversed-phase
      liquid chromatography-MS by incorporating peptide sequence information",
      *Anal Chem* 2006;78(14):5026-5039, doi:10.1021/ac060143p, PMID 16841926, and
      **Moruz L, Käll L**, "Peptide retention time prediction", *Mass Spectrom Rev*
      2017;36(5):615-623, doi:10.1002/mas.21488 — the latter states the paper's `:571-572`
      argument verbatim and traces it to Sanger 1952. Attribute §3.7's *premise* to them
      and narrow the novelty claim at `:583-584` to the *metric*. Also acknowledge
      AutoRT's `--add_reverse` flag.
- [ ] **Bouwmeester R et al.**, "Transfer learning in DeepLC improves LC retention time
      prediction across substantially different modifications and setups",
      *Nature Communications* 2026;17:2601, doi:10.1038/s41467-026-68981-5. Natural
      comparison at `:195-209`, which contrasts uno's adaptation with a frozen predictor.

Corrections:

- [ ] `paper.typ:609-613` — **replace `@the2016percolator`**. Percolator 3.0 is about
      throughput and picked-protein inference and introduces none of the three named
      mechanisms. Use **Käll L, Canterbury JD, Weston J, Noble WS, MacCoss MJ**,
      *Nat Methods* 2007;4:923-925, doi:10.1038/nmeth1113 for semi-supervised iterated
      positives, and **Granholm, Noble & Käll**, *J Proteome Res* 2012 for target/decoy
      pair co-folding. Keep The 2016 only if protein-level FDR is meant.
- [ ] `paper.typ:84-86` — "The DeepLC and LightGBM stages are stochastic and unseeded" is
      falsified by the source: `SEED = 50` (`main.py:29`), LightGBM training explicitly
      seeded (`:1593`, `:1603`), CV permutations and ID shuffles seeded, and DeepLC is
      only called for `make_preds`/`calibrate_preds` with no retraining. The genuinely
      unseeded paths are the **decoy shuffle** (`utils.py:247`, bare `random.shuffle`)
      and the **hyperparameter random search** (`main.py:1567`, bare
      `np.random.RandomState()`) — which is exactly what the paper's own §3.2 says at
      `:399-400`. The Introduction currently contradicts the Results.
- [ ] Abstract `:23-24` and `:86-88` — soften "depends on the deep-learning RT predictor
      DeepLC". DeepLC is absent from `requirements.txt`/`install_requires`, imported
      lazily, and `-deeplc` **defaults to 0**; the README calls it "recommended". Use
      "recommends, and in the configuration benchmarked here requires". (The
      *benchmarking* decision to always use `-deeplc 1 -ad 1` stays.)
- [ ] `paper.typ:86-88` — "its RT-prediction step dominates wall-clock time for a single
      run" is the authors' own measurement, not a claim in the cited paper. Point at
      their own timing data instead of leaving it inside a cited clause.
- [ ] `paper.typ:243-248` — `@senko1995averagine` supports the averagine pattern and ion
      population estimation, but contains nothing about chromatography or co-eluting
      interference (it is an FT-ICR paper). Keep the citation for the pattern; mark the
      interference-rejection property as this work's or koth's.
- [ ] `paper.typ:250-252` — "following the DirectMS1Quant approach" overstates it.
      DirectMS1Quant explicitly applies **no correction** to peptide p-values, gates at
      raw p<0.05 plus a fold-change threshold, counts co-directed significant peptides
      k of n, computes a **binomial** protein probability, applies BH to *that*, and
      accepts at **p<0.05**. uno does a per-protein one-sample t-test on aggregated
      peptide fold changes with BH at **1%**. The machinery and the threshold both
      differ, and the 5× stricter threshold plausibly accounts for much of the
      sensitivity gap reported at `:690-705`. (BH itself is correctly attributed, and
      `@benjamini1995fdr` is the right primary reference.)
- [ ] `paper.typ:199-201` — add `@ivanov2021ms1searchpy` for "the DirectMS1 on-data
      trick" (described there: coefficients "calculated on-the-fly for specific data and
      separation conditions", plus the two-stage refit).
- [ ] `paper.typ:756` — add `@ivanov2021ms1searchpy` for "consistent with DirectMS1's
      minimal digest" (that paper states "no missed cleavage sites").
- [ ] `paper.bib:94-104` — `abdrakhimov2021biosaur` has `pages = {e9045}` with no volume
      or issue; Crossref now assigns **volume 39, issue S1**. Add them. Also consider a
      parenthetical noting that biosaur2 has no separate publication, so the "Biosaur"
      reference is intentional.
- [ ] `paper.bib:119-126` — `ke2017lightgbm` has no `pages`; the printed NIPS 2017
      proceedings pagination is 3146-3154.
- [ ] Move locators out of `note` fields into real `url` / `publisher` fields (or use
      `@software`). The renderer **silently discards `note`**, so `@koth` renders with no
      URL and no version, `@gbdtcrate` with no version and no crates.io URL,
      `ke2017lightgbm` loses its proceedings URL, and `wilburn2023chronologer` loses
      "Preprint, posted 2 June 2023". Two of three software citations currently carry no
      locator at all in the output. Also remove the LaTeX `\url{}` / `\texttt{}` macros,
      which would render literally under any style that emits `note`.
- [ ] Switch to ACS style and rename the heading: `#bibliography("paper.bib",
      style: "american-chemical-society", title: "References")`. JPR requires
      **superscript numbers without parentheses** and **CASSI journal abbreviations**;
      the current IEEE default renders bracketed `[1]`, truncates author lists to
      "et al." after three, and prints full journal names with vol./no./pp. labels.
      Reference *ordering* is already correct (order of first mention).
- [ ] Check with institutional access: whether **Krokhin & Spicer**
      (*Expert Rev Proteomics* 2011, doi:10.1586/epr.11.79, or *Anal Chem* 2007,
      PMID 17939646) claim SSRCalc predicts "mirror sequences". A secondary source quotes
      that wording; if it is in the primary text it is the closest competing prior art on
      the §3.7 priority claim and must be cited. This is the one live threat to that claim.
- [ ] Optional, shows command of the niche: **Penanes PA et al.**, "Potential of
      Negative-Ion-Mode Proteomics: An MS1-Only Approach", *J Proteome Res* 2023,
      doi:10.1021/acs.jproteome.3c00307; **Postoenko VI et al.**, "IQMMA…",
      *J Proteome Res* 2023, doi:10.1021/acs.jproteome.3c00075 (relevant to the
      feature-finder choice at `:263-266`).
- [ ] Also cite or at least declare **sequoia-boost** (public at
      github.com/pgarrett-scripps/sequoia-boost, not on crates.io) — it is the default ML
      backend and appears nowhere in `paper/`.
- [ ] `paper.typ:153-158` vs `:568-584` — make the qualifier "under composition-based RT
      models" explicit at `:155`. uno's shipped model is **not** composition-based, and
      §3.7 reports it capturing 72% of the forward−reverse signal, so a reviewer will ask
      why reversed decoys were still rejected.

---

## 6. Provenance, reproducibility, release

- [ ] Commit the working tree, and commit `paper/figures/`, `paper/si/`,
      `paper/scripts/`, `paper/figdata/` — all currently untracked, so **a checkout of
      HEAD does not compile** (every `image(...)` and every `include "si/*.typ"` is
      missing).
- [?] (D1) Rewrite `paper/si/versions.typ:8`. "uno 0.1.0, commit `7c80199`" cannot
      describe every measurement: the Table 1 / Fig 4 outputs carry an **8-column**
      protein schema dated 2026-07-22, the FAIMS outputs a **10-column** group-aware
      schema dated 2026-07-28 — two builds, six commits apart, spanning
      `mass_anchor_ppm`, group-aware ID, feature competition, and the row-order fix.
      Also HEAD is `7f2cece`, one commit later, with a dirty tree.
- [ ] `paper/si/versions.typ` — add rows for **sage-core** (path dep, unpublished
      `0.15.0-beta.2`, `../sage` at `0639176`) and **sequoia-boost** (`../Sequoia` at
      `5c2772e`, the default ML backend). Path deps carry no `Cargo.lock` hash, so a
      commit is the only pin. Also add **matplotlib** (all ten figures need it),
      **typst 0.14.2** (the manuscript's own build tool), and **OS/kernel**
      (Ubuntu 24.04.4, 7.0.0-28-generic).
- [ ] `paper/si/versions.typ:1-4` — the provenance line names `benchmark/README.md`,
      which contains no version or hardware information.
- [ ] `paper/si/versions.typ:20-21` — "NumPy / SciPy … scoring scripts" is wrong: all
      four scoring scripts are stdlib-only. The scripts that use numpy/scipy are the
      figure scripts, which run under `uno/.venv` (**scipy 1.17.1**, not the listed
      1.16.2). Two biosaur2 versions exist on the host (0.2.28 in `FeatureFinder/.venv`,
      0.2.35 in `uno/.venv`) — disambiguate which produced the feature files.
- [ ] `paper/si/versions.typ:12` — DeepLC's role line ("ms1searchpy retention-time
      model") hides that DeepLC is also the **distillation teacher** and the label source
      for `@tbl:rtdepth` and `@tbl:reverse`.
- [ ] `paper/si/versions.typ:26-29` — "all timings" cannot cover §3.6–3.9, which are an
      earlier development configuration. Scope the table.
- [ ] Relabel `paper/si/results_per_dataset.typ:1-3` and
      `paper/si/entrap_reported.typ:1-4` — both say "GENERATED … Do not edit by hand"
      but **no generator exists anywhere on disk**, so both are hand-maintained and the
      label discourages the only check that would catch drift. Either relabel
      "Transcribed" or write the generator. (Both were re-derived by the audit and are
      currently numerically exact.)
- [ ] `paper/si/uno.toml:1-4`, `paper/si/defaults.toml:1-3` — "Copy of" is not accurate:
      values are byte-identical but comments were rewritten and abridged. Say "values
      copied verbatim; comments condensed".
- [ ] `paper/si/defaults.toml:2-3` — "overrides only the three `[digest]` fields"
      overstates by one: `min_len = 7` equals the default, so only `max_len` and
      `missed_cleavages` change anything.
- [ ] `paper/si/commands.typ:1-2` — "Verbatim from…" is inaccurate and the block is not
      verbatim. Fixes needed: block 0 comes from `benchmark/README.md:57-61`, which the
      header does not name; every path was rewritten root-relative; `run.sh:43-49` scores
      with an inline `python3 -c` importing `entrapment` as a module, not the CLI form
      shown; the `score_entrap.py ms1 …` line (`run_faims.sh:39`) is **omitted**, so the
      ms1searchpy row of Table S4 has no command; `export N_TARGET_DB/N_ENTRAP_DB`
      (`run_faims.sh:24-26`) is omitted (harmless only because the defaults happen to
      match); `ms1searchpy LUM2_${R}_VG.features.tsv` **fails as printed** (the file is
      under `$DATA`; the script copies it into `$OUT` and `cd`s there first); and
      `${DS}`, `${OUT}`, `${R}`, `$DATA` are never defined.
- [ ] `paper/si/commands.typ:37` — `intensity_matrix.tsv` has no resolvable path (the
      real one is `/mnt/labdata/.../uno_proc_newkoth/align_avg/intensity_matrix.tsv`),
      and **the koth alignment step that produces it is documented nowhere in S1**,
      although `:304-305` says the runs were "cross-run-aligned by koth". Add the
      `koth_align` command and its `[lfq] averagine_projection = true` config; without
      them §3.10 is not reproducible from the SI. Following S1 as printed reproduces the
      **deprecated 243-protein box-sum arm**, not the reported 326.
- [ ] Add a **dataset-name mapping table** to S1: manuscript name → deposit raw file →
      local feature file. The paper says `09746`, the file is
      `benchmark/data/hela5min_09746.features.tsv`, and `si/commands.typ:12`'s `${DS}`
      only works as `hela5min_09746`. Note that **09744 ≠ 09746** — `QEHF1_09744_JB`
      (`:311`) is a different run, given in a different format, and is not one of the six.
- [ ] `paper.typ:851-854` and `:881` — narrow "the **exact** commands, versions and
      configuration behind **every** measurement" and "The **complete** default set".
      S1 covers only the six-dataset, FAIMS and quant-post-processing steps. Not covered:
      the §3.5 column-stripping, §3.6 RT comparison, §3.7 reverse benchmark, §3.8 LDA,
      §3.9 profiling, the run10 numbers quoted in §3.1, and the koth alignment.
      Preferred fix: narrow the promise **and** add one line stating plainly that the
      development-configuration experiments predate the archived pipeline and cannot be
      re-run — this converts a checkable overclaim into a disclosed limitation.
- [ ] `paper.typ:851-854` — the stated SI contents omit **S4 (FAIMS), S5 (Resource
      scaling) and S6 (quality columns)**, 3 of 8 sections. Extend the sentence.
- [?] (D5) Deposit the results data. `benchmark/data/`, `benchmark/results/`,
      `benchmark/faims/faims_out/` and the quant matrices are **all gitignored**, so a
      reader cloning the repo gets none of the inputs behind any SI table. Enumerate what
      a reviewer would actually receive today, and what the Zenodo archive must contain.
- [?] (D5) `paper.typ:820-824` — fill the two placeholders. `#link("https://github.com/")
      [\[URL\]]` renders literally as "[URL]"; the repo 404s; the Zenodo DOI is
      "[DOI/Zenodo, TODO]". **Hard submission blocker.**
- [?] (D5) `@koth` — private (404), DOI-less, no `CITATION.cff`, and load-bearing for all
      of §3.10, Fig 5, Fig S5 and Table S7. Either release it (a versioned Zenodo archive
      suffices) or state explicitly in Data and Software Availability that the
      quantitation path depends on unreleased software.
- [ ] `paper.typ:796-798` — name **all three** blocking path dependencies, not just
      sage: `../sage/crates/sage`, `../koth_rust/koth_ff`, `../Sequoia/crates/sequoia-boost`.
- [ ] Add a `CITATION.cff` to uno (and koth if released), and a `requirements.txt` or
      `pyproject.toml` for the figure/scoring scripts — a fresh clone has no venv, so
      `just figures` silently depends on the reader having matplotlib and scipy.
- [ ] Consider citing **Perez-Riverol Y et al.**, "Open-Source and FAIR Research Software
      for Proteomics", *J Proteome Res* 2025;24:2222-2234,
      doi:10.1021/acs.jproteome.4c01079 — same journal, directly on the open-source claim
      at `:822`, and it sets the FAIR expectations the placeholder statement fails.

---

## 7. Compliance and back matter

- [ ] **Trim the abstract from 301 to ≤200 words** (JPR hard limit). Do this **after
      §2**, since four of its sentences are being corrected. Natural cuts: the
      determinism-range sentence (`:33-34`), the 12-of-18 sentence (`:37-38`), and part
      of the quant block (`:39-43`).
- [ ] Add `ASSOCIATED CONTENT` / `AUTHOR INFORMATION` / `ACKNOWLEDGMENTS` back matter.
      Every required element is currently absent:
      - [ ] Corresponding author marked with an asterisk, **with telephone number**
      - [ ] Supporting Information description **in the main text**, labelled
            "Supporting Information", non-sentence form, with the file-type extension
            ("(PDF)"), placed **before** Acknowledgments and References, plus the
            "available free of charge" statement — and covering **all eight** SI sections
            plus Figures S1–S5 and Tables S1–S8
      - [?] (D16) Author Contributions
      - [?] (D16) Conflict of Interest / Notes (if none: "The authors declare no
            competing financial interest.")
      - [?] (D16) Funding Sources (required in the manuscript **and** the Funder Registry
            Tool)
      - [?] (D16) Acknowledgments — **including the AI-use disclosure**, which ACS
            requires inside Acknowledgments "with a description of when and how the tools
            were used". This is owed: `paper/STYLE_PROMPT.md` is an LLM copy-edit prompt
            applied to this manuscript.
- [ ] `paper.typ:5-17` — remove the typed `orcid:` fields. ACS: "ORCID iDs should **not**
      be typed into the manuscript"; link them in the submission system instead.
- [ ] `paper.typ:114` — rename "Experimental Procedures (Methods)" → **"Experimental
      Section"** (JPR's section name). This changes the three `(Methods, …)` pointers at
      `:538`, `:674`, `:864`.
- [ ] Add a one-line statement that no new human material was collected (all three
      datasets are public re-analyses of HeLa, an immortalized line), so an editorial
      ethics query never arises.
- [ ] Check whether JPR requires a **TOC / abstract graphic** for Articles (the
      guidelines reference a separate specification without reproducing it) and produce
      one if so. None exists.
- [?] (D12) Split the SI into a separate file for submission. JPR requires SI "submitted
      as a separate file(s)"; the single combined PDF is a documented decision
      (`paper/FIGURE_MANIFEST.md:28-29`, `paper.typ:830-833`).
- [?] (D8) "Validation roadmap" (`:785-798`). JPR does not forbid extra sections, but
      this one is neither results nor interpretation — it is four proposed experiments
      plus an admission that the release is incomplete, in a paper whose availability
      statement is already two placeholders. The two items carrying real information (the
      coarse entrapment bound, the trace-level sensitivity gap) are already in
      Limitations.
- [ ] Conclusions (`:800-818`) restate §3.1–3.4 and §3.10 numerically with little new
      interpretation; ACS asks that Conclusions "not summarize information already
      presented".
- [?] (D15) Title says "Protein Identification" only, while quantitation occupies §2.7,
      §3.10, Figure 5, Table S7, Figure S5 and four abstract sentences.
- [ ] `paper.typ:49` — `date: "June 2026"` is stale (measurements are 2026-07-28) and is
      not an ACS element.
- [ ] Populate PDF metadata: `#set document(keywords: (...))` — Keywords and Subject are
      currently empty although the 7 keywords render on page 1.

---

## 8. Figures, tables, layout

- [ ] **Make Tables S2 and S3 breakable** (`breakable: true`, repeated header row, or
      split per-dataset). Both overflow their page: the last data row prints **on top of
      the caption and the page folio** and is unreadable (60 and 20 out-of-box words at
      y=769.8 pt against a 756.9 pt margin). Present in the committed PDF and a fresh
      compile.
- [ ] Set an explicit width on the Dataset column of Tables S2 and S3 —
      `LUM2_05433_VG` overflows into the neighbour (renders as `LUM2_054331000949` and
      `LUM2_05433_VGuno`). Also shorten "ms1searchpy rep N", which wraps in all 18 rows.
- [ ] `paper.typ:922-928` — **rewrite the Figure S1 caption.** It describes Figure 4b's
      axes. The figure plots x = engine-reported q, y = entrapment FDP, and the vertical
      line is at **reported q = 1%**, not the entrapment operating point
      (`figS1_fdp_curves.py:56-80`).
- [ ] `paper.typ:134-145` — **rewrite the Figure 1 caption.** It describes an
      identification path "(left)" and a quantitation path "(right)"; the schematic is
      one vertical cascade whose left lane is the **database** branch and right lane the
      **feature-file** branch (`fig1_pipeline.py:116,140`). No aligned-matrix entry
      point, cell mapping or profile competition is drawn. Either recaption or add the
      quantitation path to the figure.
- [ ] `paper/scripts/figures/fig1_pipeline.py:176-177` — the schematic labels the shipped
      RT model "a pure-Rust CNN **distilled from DeepLC** (chrono_global_v3_ch32)",
      contradicting Methods `:183-186` (trained directly on Chronologer experimental RT)
      and undercutting the paper's own "no deep-learning dependency" framing. Relabel;
      the DeepLC-distilled student is the MLP of §3.6.
- [ ] `fig1_pipeline.py:125` — the schematic states a "+11 Da decoy mass shift" that
      Methods never mentions (it is real: `database.rs` `decoy_mass_shift 11.0`). State
      it in Methods or drop it from the figure.
- [ ] `paper.typ:417-423` — the Figure 3 caption describes only panel (a); the figure
      draws four labelled panels (depth @1%, depth @5%, replicate range as % of median,
      run10 close-up). Add `*(a)*`–`*(d)*`.
- [ ] `paper.typ:464-477` — the Figure 4 caption never identifies its two most prominent
      elements: the thick across-dataset median curves and the two annotated points
      (0.7% / 1.2%) that the §3.3 claim rests on.
- [ ] `paper.typ:966-973` — the Figure S3 caption says "the three ms1searchpy
      replicates"; the figure plots the **median with a min–max error bar**
      (`figS3_resources.py:96-104`). It also silently includes **run10 as a seventh
      dataset** — the point that carries the scaling claim and is not one of the six.
      And soften "scale smoothly … close to constant in ratio": uno's plotted peak RSS is
      non-monotone (4.10 GiB at 101k features → 3.80 at 211k) and the wall-clock ratio
      spans 5.26–7.79× across the plotted points.
- [ ] `paper.typ:988-993` — the Figure S4 caption says "area under the ROC curve … An
      AUC of 0.5 is no separation"; panel (b) plots **|AUC − 0.5|**, where **0** is no
      separation (`figS4_quality_columns.py:96-101`). The caption also omits the two
      non-column rows the panel shows (RT, full classifier).
- [ ] `figS4_quality_columns.png` — panel (b) row label reads "full ML classifier" vs
      "full match classifier" in Table S6 and at `:528-529`. Harmonize.
- [ ] `fig2_headtohead.py:57` and `figS3_resources.py:107` — axis label "peak RSS (GB)"
      for a quantity that is `peak_rss_mb / 1024` where `peak_rss_mb` is already KB/1024,
      i.e. **GiB** (7.4% below the GB value). The main text at `:342` correctly says GiB,
      so the axes contradict the prose. Also update `paper/FIGURE_MANIFEST.md`, which
      says "peak RSS GB".
- [ ] **Re-render `fig5_quant` and `figS5_quant_detail` under the venv.** They were
      rendered by **matplotlib 3.6.3** (bare `python3`) while the other eight used
      3.11.0, so they are the only two that do not regenerate byte-identically. Numeric
      content is identical (verified against `fig5_stdout.txt` and `figS5_stdout.txt`
      line for line), so this is rendering consistency — but `just figures-all` will
      silently change both. **Do this before any figure regeneration.**
- [ ] `fig5_quant.py:228` — the "BH 1%" annotation is clipped by the y-spine and
      overplotted by scatter markers (renders as "PH"). Move it inside the axes and raise
      its zorder.
- [ ] `figS2_faims.png` — panels (a) and (b) use different y-limits (0–2000 vs 0–2500),
      so equal bar heights mean different counts across panels.
- [ ] `paper.typ:1017-1025` — the Figure S5 caption has no panel letters although the
      figure has them and the Figure 5 caption uses them. Same for Figures 3, S2, S4.
- [ ] Add a per-replicate **5% column** to `paper/si/faims.typ`. The prose at `:493-494`
      rests on all six replicate-by-threshold cells and on the "1878 against 1863" pair;
      the table gives 1% only. The values exist (`faims_out/figS2_stdout.txt`, annotated
      on Fig S2b) and were reproduced exactly by the audit.
- [ ] Renumber the SI items into citation order. Main-text first-mention order is
      **Table S5 → Figure S3 → Table S2** → S3 → S4/Figure S2 → S6 → Figure S4 → S7 →
      Figure S5. S numbering is assigned by definition order in the SI.
- [ ] **Figure S1 is never cited in the main text** (only inside the SI at `:915` and in
      a `//` comment at `:449`).
- [ ] `paper.typ:963-975` — SI section **S5 has no body text**; the heading is
      immediately followed by S6, and its figure floats to the next page. **S7** likewise
      has zero body text. Add one sentence each, or fold S5 into S2.
- [ ] Page 11 is ~40% blank because Table 2 floated to page 12; Figure 5 is first
      mentioned on page 13 and printed on page 15, inside Limitations. One `placement`
      pass.
- [ ] Consider embedding the `.pdf` figures instead of `.png` — a vector
      `fig1_pipeline.pdf` already exists, and Figure 1 is text-heavy line art for which
      ACS asks 1200 dpi. Switching all ten satisfies every ACS resolution tier and
      shrinks the PDF. Note the current 300 dpi PNGs **do** meet the colour-art rule.
- [ ] Tables render with a **full grid including vertical rules**; ACS uses horizontal
      rules only. Add `stroke: (x: none, y: 0.5pt)` or equivalent.
- [ ] SI table cells are justified, producing large inter-word gaps (Table S1's Role
      column renders `MS1-only␣␣␣␣␣search␣␣␣␣␣and`). Set `par(justify: false)` inside
      tables and widen the affected columns.

---

## 9. Mechanical sweep

Land as **its own commit**, last, so it does not obscure substantive diffs.

- [ ] Run the STYLE_PROMPT copy-edit over **§2.7, §3.10 and the Figure 5 caption only**.
      All 14 house-style violations live there — the block revised at 20:05, after the
      11:27 style pass. Fixing that one section clears three failing rules at once.
      - [ ] 4 em-dashes, all in prose: `:244`, `:703` (×2), `:712`. Zero in markup,
            comments, code spans or n/a cells.
      - [ ] 6 semicolons, none defensible under the rule's exemptions: `:242`, `:676`,
            `:682`, `:692`, `:699`, `:733` (Fig 5 caption).
      - [ ] 4 sentences > 40 words: `:709` (59 w), `:697` (47 w), `:681` (42 w),
            `:702` (41 w).
      - [ ] 4 rhetorical colons: `:246`, `:646`, `:709`, `:981`.
- [ ] `paper.typ:913` — `$20{,}647 \/ 6{,}066 = 3.40$` renders as
      **`20{, }647/6{, }066`**. Typst math treats `{}` as literal glyphs, not grouping.
      Take the numbers out of math; the main text already writes them as plain prose at
      `:268`.
- [ ] `paper.typ:242`, `:717`, `:747` — `$m/z$` renders as a **built-up fraction**
      mid-sentence, three times. Use `_m/z_` or `$m slash z$`.
- [ ] `paper.typ:283` — the FDP formula, the load-bearing definition of the paper, is set
      inline, wraps mid-formula, mixes a built-up fraction with a literal slash, and has
      **no operator between `n_yeast` and the parenthesis**. Make it a display equation
      with one notation, and define the symbols in words. Suggested:
      `$ "FDP"(k) = (n_"yeast"(k)) / (n_"human"(k)) dot (N_"human") / (N_"yeast") $`.
- [ ] `paper.typ:52` — remove the global `#set raw(lang: "rust")`. Nothing in the
      manuscript is inline Rust, and it puts spurious highlighting inside file names and
      flags: a **red `*`** in `*.features.tsv`, a red `=` in `min_scans = 2`, a red `--`
      in `cargo build --release` (which also line-breaks mid-flag), a red `/` in
      `benchmark/entrapment.py`.
- [ ] `paper.typ:756` — `DirectMS1's` renders with a **prime** (U+2032), not an
      apostrophe: Typst converts `'` after a digit. Recast as "the minimal digest of
      DirectMS1", or insert the character literally.
- [ ] `paper.typ:663-664` — `$40 times$` / `$13 times$` render with a leading space
      ("roughly 40 × on its own"). Use `40×` / `13×` as literal text, matching `:175`.
- [ ] Expand abbreviations. **In the abstract**: LDA (`:44`), FAIMS (`:29`), UPS1
      (`:39`). **In the body**: TMT (`:757`), XIC (`:237`), ROC / AUC (`:517`), RANSAC
      (`:190`), MLP (`:184`, expanded at `:180` but never linked), RSS
      (`si/versions.typ:24`, expanded three pages later at `:896`). Introduce **FDR** and
      **FDP** properly — both are used abbreviated without a parenthetical introduction
      in either the abstract or the body, and ACS requires each to stand alone. RT is
      defined only in the abstract; the body never redefines it.
- [ ] `paper.typ:500` and `@tbl:entrap:638` — **CV means two different things**
      (compensation voltage; cross-validation), and a third sense (coefficient of
      variation) is spelled out at `:795`. Write "compensation-voltage fractions" at
      `:500` and "hardened cross-validation" in the table.
- [ ] Cut redundancy (keep the instance named):
      - [ ] "depth here is entrapment depth, not the engine's reported q" — **5×**
            (`:285-289`, `:431-433`, `:475-477`, `:894-896`, `:906-907`). Keep Methods +
            the Table S2 caption clause.
      - [ ] "same single-run input to both engines" — **7×** (`:26-27`, `:91-93`,
            `:118-129`, `:143-144`, `:320-321`, `:376-377`, `:803`). Keep Methods + the
            Fig 2 caption.
      - [ ] "uno run once because deterministic, ms1searchpy three times" — **5×**
            (`:271-273`, `:354-355`, `:379-381`, `:400-401`, `:859-861`). Keep Methods.
      - [ ] the quant "detection, not the estimator" argument — **5×** (`:689-721`,
            `:773-783`, abstract, Conclusions, roadmap, Fig 5 caption). Keep §3.10;
            reduce the fifth Limitation to one sentence plus a pointer.
      - [ ] single-dataset FDP quantization — **4×** (`:457-459`, `:910-917`, `:926-927`,
            `:762-765`). Keep S3, which has the arithmetic.
      - [ ] mass and RT cannot separate at the peptide level — **3×** (`:215-217`,
            `:506-512`, `:741-744`). Keep §3.5.
      - [ ] the FAIMS database description — **3×** near-verbatim (`:293-295`,
            `:487-489`, `:939-943`). Keep Methods.
      - [ ] the 830-vs-1035 example — **2×** (`:436-437`, `:908`). Keep S3, which names
            the dataset.
      - [ ] "decoys match real features by mass coincidence, so the columns describe the
            feature" — **2×** near-verbatim (`:521-523`, `:981-983`). Keep §3.5; S6 should
            defer.
      - [ ] Captions restating the adjacent paragraph: Fig 2 (`:382-383`), Fig 3
            (`:420-423`), Fig 4 (`:471-474`), Fig 5 (`:733-735`), Fig S2 (`:950-953`).
- [ ] LLM tics that survived the copy-edit (the big ones are already gone — zero "it is
      worth noting", "importantly", "notably", "crucially", "strikingly"):
      - [ ] "which is why" **×5** (`:45`, `:108`, `:400`, `:703`, `:814`) — and `:108`
            and `:814` are the *same sentence* in the Introduction and Conclusions.
      - [ ] "rather than" ×13 — the X-not-Y frame is the default rhythm.
      - [ ] emphasis-by-adverb: "is **precisely** the signal" (`:573`), "is **itself**
            the result" (`:697`), "robust **precisely** because" (`:619`),
            "**by construction** rather than **by luck**" (`:403`).
      - [ ] signposting: "That is the subject of Section 3.2." (`:334`), "What transfers
            is the ordering" (`:498`), "Note that the horizontal axis of (b)…" (`:474`).
      - [ ] aphorisms: "A single reported number from such an engine is a draw, not a
            measurement" (`:400`), "This is a transferable cautionary result." (`:621`),
            "exposed the gain as illusory" (`:605`), "no room to separate" (`:158`).
- [ ] Terminology:
      - [ ] "honest" (`:248`, `:672`) — internal provenance jargon from
            `benchmark/RESULTS.md`, undefined in the paper, and it implies the *other*
            numbers are not. Use "entrapment-scored" or nothing.
      - [ ] "mass-anchor background probability" (`:712`) renames the "random-match
            background rate p estimated from decoys" of `:218-219`. Use one name.
      - [ ] "box-sum" (`:242`, `:248`) vs "**raw** box-sum" (`:676`).
      - [ ] "two students" (`:170-186`) — only the first is distilled; the second is
            trained directly on Chronologer. Rename the second "the convolutional model".
      - [ ] "deviation discriminant" (`:212`) is used once, never defined, and its
            relation to "the full match classifier" (`:528`) is never stated.
      - [ ] "over-fill" (`:1024`) and "transferred fills" (`:702`) are each coined once
            and never defined.
- [ ] Grammar and wrong-word:
      - [ ] `:442` "calibrated **to** conservative" — ungrammatical; "calibrated **or**".
      - [ ] `:707` "**Of the precursors lost**, the signal is…" — dangling modifier;
            "**For** the precursors that are lost…".
      - [ ] `:337-339` "That is a factor of 5.3 to 7.0. **It** also ran in roughly
            half…" — nearest antecedent is "factor"; use "uno".
      - [ ] `:391-392` "spanned **ranges** from 30 to 600 proteins, **equal to** 3.6%
            **and** 40.5%" — a range equated to two point values; the abstract says
            "3.6 **to** 40.5%". Harmonize.
      - [ ] `:528` "further from chance **again**" — redundant.
      - [ ] `:268` "Both engines searched **one** human plus *S. cerevisiae* FASTA" →
            "a single combined … FASTA".
      - [ ] `:541-543` "A naive composition-feature model … **trailed** at R² = 0.972" —
            trailed what? The comparator arrives in the next sentence.
      - [ ] `:544` "**fewer of them** overlapped the reference search" — ambiguous
            between absolute (119 vs 152) and proportional (84% vs 86%).
      - [ ] `:429` "behaviour" — the only Britishism; ACS uses American spelling.
      - [ ] `:473` "five **of six**" vs `:35`, `:442`, `:812` "five **of the six**".
      - [ ] Tense: §3.1–3.5 are past ("uno **identified**"), §3.10 and the abstract's
            quant sentence are present ("uno **identifies**"). ACS asks for consistency.
- [ ] Amend `paper/STYLE_PROMPT.md`: the "~76-column hard wrap" is not held (414 of 1051
      lines exceed 76; the actual wrap is ~82). **Amend the rule, do not reflow 414
      lines.** Genuine outliers worth a manual reflow: `:183` (142 chars), `:669` (106),
      `:561` (105), `:592` (98), `:563` (91), plus ~12 un-reflowed short mid-paragraph
      lines marking previous edit sites (`:179, 191, 252, 265, 300, 485, 525, 542, 545,
      610, 776, 817`).
- [ ] Also codify in STYLE_PROMPT: ranges are written "to" in prose and with an en-dash
      in tables (currently coherent but unstated, so a future copy-edit could "fix" it
      the wrong way); the thousands-separator convention (counts get commas, protein
      depths never do); and the en-dash source form (`paper.typ` uses literal `–`,
      `si/*.typ` uses Typst `--`).
- [ ] Remove `paper/scripts/__pycache__/` (7 `.pyc`) and the superseded
      `paper/uno_jpr.md` from the submission bundle. Consider untracking the built
      `paper/paper.pdf` (2.4 MB binary).
- [ ] `paper/scripts/paper_stats.py` — the checklist covers rules 1–6 only. Consider
      adding rule 11 (species italics) and a §3.10/Table S7 guard; the most recently
      revised section is the one with no assertions. Relabel the "Journal word count"
      header, which both excludes captions (ACS counts them) and includes the
      Experimental Section (ACS excludes it).
- [ ] `paper/scripts/README.md` and the `--report` header use em-dashes; the house-style
      checker is not itself house-style clean.

---

## 10. Verification gate

Run all of these before declaring the manuscript ready.

- [ ] `typst compile` — zero warnings, zero errors.
- [ ] **Bbox scan every page for out-of-box text**, not just S2/S3. This is the check
      that would have caught the table overflow.
- [ ] `paper_stats.py --report` — all six house-style rules PASS; abstract ≤ 200 words.
- [ ] Every number in the paper re-derived from disk **by a script**, including the new
      5% table, the DE figures and the min-runs-6 arm. Keep the script.
- [ ] All ten figures regenerate byte-identically under one interpreter.
- [ ] Label / citation sets unchanged except intended additions; 0 uncited bib entries.
- [ ] Rendered reference list proofread **in ACS style** (the switch will expose
      author-list truncation, journal abbreviations and discarded fields).
- [ ] A clean clone builds the paper, and every command in S1 runs as printed.
- [ ] Two-run byte-diff determinism test passes in CI.

---

## 11. Ordering constraints

All 16 decisions are settled (§0), so the gate is now the re-run.

1. **§1 (config truth) + the D11 code fixes first.** Methods must be written against a
   truthful default listing, and the determinism/decoy fixes must land before the re-run.
2. **§12 re-run second.** Every number in Table 1, Table S2, Table S3, Table S5, Table S6,
   Figs 2/3/4/S1/S3/S4, the abstract and the Conclusions is **provisional** until it
   completes. Do not write results prose before then.
3. **§2 and §3 (claim corrections) after the re-run, before §7's abstract trim** — or the
   abstract gets trimmed twice. The abstract must lose ~100 words *and* four of its
   sentences are being corrected.
4. **D9's renumbering (§3.10 → §3.9) before any section-pointer edits**, so pointers are
   updated once.
5. **Commit the tree and pin the re-run commit before touching `si/versions.typ`.**
6. **Write Data and Software Availability last**, once D5's URL and DOI exist. Do not
   draft it early.
7. **Switch the bibliography style before proofreading the reference list** — the switch
   exposes author-list truncation, journal abbreviations and the discarded `note` fields.
8. **Fix table breakability before the final layout and page-count pass.**
9. **Re-render fig5/figS5 under the venv before any `just figures-all`.**
10. **§9 mechanical sweep lands last, as its own commit**, so it does not obscure
    substantive diffs.

Critical path: **§1 + D11 → §12 re-run → §2 → §3 → §7 (abstract) → §10.**
§5 (citations), §6 (provenance), §8 (figures) and §9 (mechanical) hang off that spine and
can run in parallel with the re-run.

---

## 12. Re-run protocol (D1)

The re-run is the gate for most of §2. Doing it once, correctly, is cheaper than doing it
twice — so fold in everything that needs fresh output.

**Preconditions** (all must be true before the first run):

- [ ] D11 fixes landed: `est_rt_sig` residual sort, quant peptide-vector sort,
      `pep.decoy` guard at `score.rs:465`.
- [ ] Two-run byte-diff determinism test written and **passing**.
- [ ] `benchmark/entrapment.py:57` tie-break made deterministic
      (`key=lambda r: (-r[0], r[1])`).
- [ ] `example_config.toml` defaults corrected (§1), since the re-run publishes them.
- [ ] Working tree committed. **Record the commit hash** — it becomes the single value in
      `si/versions.typ`.
- [ ] ms1searchpy patch committed to a branch or captured as a diff, with **its** commit
      recorded (§4).
- [ ] Log the uno build metadata alongside each run. The previous round's
      `benchmark/results/*.log` are zero-length, which is why the build could not be
      identified afterwards.

**Scope.** Serialize everything on one idle host; capture `/usr/bin/time -v` for every
invocation.

- [ ] Six datasets: uno ×1 + ms1searchpy ×3 = 18 ms1 runs.
- [ ] run10: uno ×1 + ms1searchpy ×3. Note run10 is now labelled a **tuning set** (D2),
      not a benchmark.
- [ ] **run10 column-stripped: uno ×1 + ms1searchpy ×3.** This finally fills the one
      unauditable row of Table S5, whose caption currently admits the uno
      column-stripped run "was not retained". Opportunity created by D1.
- [ ] **FAIMS: ms1searchpy ×3 per replicate**, not ×1. The previous round ran it once per
      replicate, which violates the paper's own stated standard at `:399-401`. Fixes that
      finding rather than disclosing it.
- [ ] **Dump uno's PFM table** for run10 so Table S6's uno column can be re-derived (D14).
- [ ] Quantitation: re-run `uno --matrix` at the pinned commit on the existing
      `uno_proc_newkoth/align_avg/` matrix, and **save the log** — the 20:05 run's log was
      overwritten, which is why the 326-protein result currently has no run record.
      Do **not** re-run koth alignment; that is the expensive step and its output is
      unchanged.
- [ ] Separately (D13): rebuild the reverse-sequence benchmark with DeepLC labels.

**Rough cost.** ~33 ms1searchpy runs at 200–400 s ≈ 3–3.5 h; ~11 uno runs at 30–50 s
≈ 10 min; plus the quant re-run and the DeepLC pass. Budget half a day, serialized.

**After the re-run:**

- [ ] Regenerate `benchmark/results/summary.tsv`.
- [ ] Regenerate `paper/figdata/*.tsv`, then all ten figures **under the venv**
      (`just figures-all`), so fig5/figS5 stop being matplotlib-3.6.3 outliers.
- [ ] Re-derive every table from the new outputs by script, and keep the script (§10).
- [ ] Expect the decoy-guard fix to move depth. The new numbers are genuinely new, not a
      reproduction of the old ones — do not assume the direction or the magnitude.

---

## Appendix: what the audit verified as correct

Do not re-litigate these; they were recomputed from primary data, not read off a
generated table.

- All **24 benchmark depth values** re-derive exactly from the protein tables on disk,
  by two independent implementations. Table 1's 36 cells and Table S2's 24×8 cells are
  exact, as are all six feature counts.
- **5.3–7.0× wall, 1.9–2.2× memory, 3.0–39.7% margins, "30 to 600 proteins",
  "3.6% and 40.5%"** — all confirmed from the `/usr/bin/time -v` records directly,
  including every MiB/GiB conversion. The paper's correction of `RESULTS.md`'s MiB/1000
  error is right.
- **"deeper at 1% on all six"** against the median, **"wins three of the six"** against
  the best replicate (09746, 09734, 09778), and **"tied at 5% on the deepest 5-minute
  run"** (09950, 1740 = 1740) are all correct and honestly stated.
- The **entire FAIMS section** reproduces exactly by re-running the scorer on the nine
  protein tables: 1642/1460/1672 and 1662/1205/1554 at 1%, 2112/1863/2204 and
  2007/1878/1861 at 5%, means 1591/1474 and 2060/1915, the 8% margin, and both
  ms1searchpy wins. The target-only FASTA is exactly 20,193 human + 4,431 *E. coli*.
- **`entrapment_curves.tsv` and `quality_columns.tsv` regenerate byte-for-byte**;
  **8 of 10 figures regenerate byte-identically**.
- **All nine cells of Table S7** reproduce from scratch, and the table is correctly scoped
  to the new-koth averagine arm (the other three matrices do not reproduce it). The 326
  count and the 0.34% FDP both re-derive independently.
- `@tbl:qualityablation` reproduces row by row; the ablation was genuinely performed
  (`ms1_noqual_rep1` contains none of the five columns).
- **Zero cross-reference defects**: 22 labels, 19 bib keys, 37 references, no duplicates,
  no unresolved keys, no orphans, and **all 30 hand-written section pointers resolve**.
  The SI renumbering mechanism works exactly as designed.
- **All 16 journal DOIs resolve** with matching metadata; no encoding damage; "John R.
  Yates III" renders correctly; reference ordering is already ACS-correct.
- **The comparison is not rigged**: the entrapment scorer reads both engines' output
  identically, the suspected protein-group-splitting asymmetry **does not exist** (zero
  species-mixed groups), ms1searchpy ran with `-deeplc 1 -ad 1` everywhere,
  `benchmark/run.sh` sets no `UNO_*` variable, and the sage digest is stock.
- **Novelty survives three independent adversarial searches**: no Rust MS1-only engine
  exists; no published third-party same-input ms1searchpy benchmark exists; and across
  five search angles no prior work reports forward−reverse RT Δ correlation as an
  evaluation metric. **No undisclosed author self-overlap** in the published record.
- `@tbl:superseded`'s provenance claim is correct — a previous draft does exist in git,
  and every retired value appears there verbatim.
- The FAIMS `run_faims.sh` path-bug correction is complete, with an explicit retraction
  in `benchmark/faims/README.md`; the paper quotes only the new values.
- House-style rules 2, 3 and 11 are fully clean (0 bullets, 0 bold in running prose,
  19/19 species binomials italicised with no false positives). "Significant" appears
  nowhere in the manuscript.
