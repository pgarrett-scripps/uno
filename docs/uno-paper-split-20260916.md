# Standalone uno manuscript — 2026-09-16

The identification manuscript now lives in the private repository
[uno-paper](https://github.com/pgarrett-scripps/uno-paper), beside the uno software
checkout. It retains the 15 commits affecting the original paper subtree and
all manuscript working-tree edits, including the LFQ split and timing corrections.
`provenance/migration.json` records the exact source state and SHA-256 hashes;
`provenance/source/` preserves the complete pre-split manuscript snapshot.

The `paper` entry in this repository is a relative compatibility symlink to
`../uno-paper`. Existing benchmark commands can continue to use `paper/analysis/`.
Clone the manuscript beside uno to use these commands. The Rust software remains
independent of the manuscript checkout; neither the code nor raw data was moved.

```sh
git clone git@github.com:pgarrett-scripps/uno-paper.git ../uno-paper
cd ../uno-paper
just setup
just paper
just docx
just verify
```

The manuscript's scientific generators read software and benchmark inputs from
the sibling `../uno/` by default; `UNO_ROOT=/path/to/uno` selects another checkout.
The generated figures, tables, and statistics are committed, so compiling the
manuscript does not require raw data. Heavy benchmark outputs and their frozen
scientific provenance remain in this repository and the original external data
locations. The LFQ manuscript remains separate in `../uno-lfq-paper/`.

The migration rebuilt all 11 active figure generators and the active table and
statistics generators. All 147 statistic values match the original manuscript.
The main text, abstract/configuration, SI text, and bibliography are byte-identical
to the pre-split version. PDF and Word builds, `just verify`, and
`just check-stats-deep` pass. Main text remains 4,659 words, with Flesch–Kincaid
grade 11.8; Supporting Information contains 2,362 words.
