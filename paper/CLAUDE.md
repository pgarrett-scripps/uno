# Working in this directory

Instructions for an agent editing the uno manuscript. Read [README.md](README.md)
for what the pipeline does and [STYLE.md](STYLE.md) for prose conventions. This
directory is built on paper-scaffold 3.15.1 (`just version`); HISTORY.md records
the migration and the upstream history.

**Read [`SUBMISSION_TODO.md`](SUBMISSION_TODO.md) before editing the paper**, and
[`AUDIT.md`](AUDIT.md) when you need the evidence behind an item. They hold a
pre-submission audit's fix list and 16 settled decisions, each citing `file:line`.
Most of what looks like an error in this manuscript is already recorded there with
a decision attached. Main single-run identification numbers come from the pinned
re-run of 2026-09-11 (`benchmark/publication_rerun.py`); re-run it and then
`benchmark/refresh_summaries.py` if the engine changes again.
The fresh matrix extension is separate: `figdata/koth_matrix_study.json` records
the complete KOTH 0.5.0 three-seed rerun: 144 searches, separately serialized
timing controls, verified outputs, source snapshots, and runtime-library hashes.
[`../TODO.md`](../TODO.md) is the code-side half.

## The contract for generated things

Figures and tables are declared in `assets.json` and referenced by id —
`fig("fig.fig2-headtohead", width: 100%)`, `tbl("tbl.versions")` — never by
filename; an undeclared id fails the compile, and `prose-check` errors on a
declared file named directly. Numbers the analysis computes are declared in
`analysis/scripts/gen_stats.py` and read back as `#s("id")`; guards there fail
the build when a re-run changes a number's meaning under a sentence.

Every asset the Technical Note references now has a generator that calls
`record()`. The assets still marked **adopted** (`origin.by = "adopted"` in
assets.json) are the retired figures and tables kept for the archive, plus the
three hand-maintained SI tables (`tbl.commands`, `tbl.datasets`,
`tbl.versions`) — edit those by hand, then `just adopt note="..."`. Result values in the
prose already use `#s()` reads: `gen_stats.py` derives reproducible values, and
hand-entered registry entries carry provenance notes when the generating
analysis is external or no longer retained. Use `#lit("...")` only for deliberate
protocol constants and identifiers; `just prose-check` reports anything else.

## Never do these

**Never hand-edit files under `figures/` or `si/`.** Both are declared assets;
`just check-assets` compares their hashes against assets.json and fails on a
silent edit. The SI tables are hand-transcribed today (SUBMISSION_TODO D1) — a
deliberate correction is: fix the source the number came from, edit the file,
then re-record with `just adopt note="..."` so the change is accepted on the
record. Once generators exist, an edit survives only until the next
`just assets`.

**Never type a result into the prose that no generated table, figure, or
stats.json entry backs.** If a number is worth stating, it is worth being
traceable; the ways out are `#s("id")`, a hand entry in stats.json with an
origin note, or `#lit("...")` for genuine prose.

**Never delete the `// >>> BODY START` / `// <<< BODY END` markers** in
`paper.typ`. The word counter, the readability report, and the narrator all
slice the prose at them, and each hard-fails without them.

## Before saying the work is done

Run `just verify`. It runs every gate — formatting, extractor tests, prose
rules, the stats/assets declarations, and staleness (which includes
`tools/check_overflow.py`, the uno-specific bbox scan that caught Tables S2/S3
printing over the page number while `typst compile` reported the build clean).
"I edited the text" is not done; "`just verify` is clean" is done. It rebuilds
nothing, so build first: `just paper && just verify`.

If you changed prose, quote the word count and readability numbers `just paper`
prints rather than estimating, and read the prose-check warnings rather than
silencing them (`prose-check.toml` takes an exception only with a written
reason). `just density` shows which section is densest relative to the paper.

If you changed inline markup, math, links, or cross-references, `just test`
asserts the extractors still handle every construct; they are regexes in
`tools/typst_prose.py` that a reflow can break silently, which has happened
three times. If you taught an extractor a new construct, add a case to
`tests/fixture.typ` and regenerate goldens with `just test-update`, reading the
diff before committing. Do not add coverage via `paper.typ`; that prose gets
deleted.

Before submitting: `just preflight` (fresh builds + verify + re-derived stats +
the Crossref DOI audit; needs network).

## When adding a table or figure

Read [FIGURE_MANIFEST.md](FIGURE_MANIFEST.md) first: it is the single source of
truth for figure content, panel contents and data provenance, and its first
hard rule is **never invent a number**.

A figure: add `analysis/scripts/figN_name.py`, importing `style.py` before
creating any figure (its rcParams apply at import). The script name is the
output name — `fig2_headtohead.py` writes `../figures/fig2_headtohead.{png,pdf}`
— and the analysis justfile picks up any `fig*.py` with no wiring. The script
must also call `record(...)` (see `analysis/scripts/_assets.py`) declaring the
png and its data inputs, and the prose references it as `fig("fig.<id>")`. Set
`metadata={"Software": None}` and seed any RNG, or every regeneration churns
the PNG bytes. The `.pdf` twin is declared under `<id>-pdf` and allowed as an
orphan in prose-check.toml — it is for the submission bundle, not the compile.

A table: generators go in `analysis/scripts/` as `gen_*_table.py`, write a bare
`#table(...)` into `../../si/` with no caption or label, call `record(...)`,
and are picked up with no further wiring; si-body.typ wraps `tbl("id")` in a
`#figure` where the caption and label live.

`just assets` regenerates (mtime-gated); `just assets-all` forces every figure
— use it after re-running a benchmark, because the data these scripts read
lives outside the repo and its mtime is not consulted.

## Editing the Typst preamble

The primary Word export is `just docx`: `tools/resolve_typst.py` resolves every
project helper (`#s()`, `fig()`, `tbl()`, `refn()`, `@fig:x`) into plain Typst
(`paper.resolved.typ`, gitignored — read it when the export looks wrong), and
pandoc converts it with native equations and citeproc references. The
`docx-mode` block in `paper.typ` serves only the `just docx-html` fallback and
is inert on the PDF path; if you change it, verify with `just docx-html` and
confirm the headings and abstract survive, not just that the command exits 0.

Typst gotcha worth knowing: a method chain broken across lines after `#let x =`
or inside `[...]` ends at the first newline, and the continuation is parsed as
literal text. The error points at a closure parameter and reads
`unknown variable: a`. Wrap the chain in a code block `{ ... }`.

## Python

One environment per concern, both managed by uv. The manuscript toolchain is
`pyproject.toml` at the root (the scripts live in `tools/`); the analysis has
its own in `analysis/`. Run things with `uv run`, never a bare `python3` that
picks up whatever is on PATH, and never `uv run --with X` inline: add the
dependency to the right pyproject so the lock stays honest.

## Scope

Do not restructure the pipeline to fix a one-off problem. The declaration
checks, the build stamps, the adopted-asset contract, and the resolver-based
docx path each exist because a specific failure happened — HISTORY.md records
which. Ask before removing one.
