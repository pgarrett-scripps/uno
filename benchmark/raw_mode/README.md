# Built-in feature detection (`uno --raw`) vs biosaur2 tables

Figure 6 of the manuscript. Every benchmark and Astral file is searched a
second time straight from its converted mzML with `uno --raw` (koth feature
finder in-process, bundled biosaur2-matching profile), with the FASTA, RT
model and `benchmark/uno.toml` unchanged. The two protein lists per file are
scored by `benchmark/entrapment.py` and compared.

- `run_primary.sh` — downloads the six primary `.raw` files from PRIDE
  (`urls.txt`), converts with ThermoRawFileParser 1.4.5 (docker, as in
  `benchmark/README.md`), runs `uno --raw`, and stages the protein tables as
  `benchmark/results/uno_raw_<ds>.tsv` with `<ds>.uno_raw.time` beside them.
  Working copies live in `~/Data/primary_raw/`.
- `paper/analysis/astral/run_raw_mode.sh` — the same for the three Astral
  mzML files already on disk (`~/Data/PXD066701_Astral/results/uno_raw/`).
- `summarize.py` — writes `summary.tsv` (tracked): per file and arm, feature
  count, depth at 1%/5% entrapment FDR, shared/union of the two 1% lists,
  wall time and peak RSS. The raw arm's wall time includes feature detection.
- `paper/analysis/scripts/fig6_raw_mode.py` reads `summary.tsv`.

Thermo `.raw` needs uno built with `--features thermo` (.NET 8); the runs here
went through mzML, which needs no extra runtime.

Multi-CV FAIMS is supported: koth partitions MS1 hills, feature assembly, and
scoring by compensation voltage. On the three-CV Lumos file
(`faims_LUM2_05433_VG`), the raw arm reports 99,812 features and reaches 2,175
proteins at 1% entrapment FDR, versus 100,949 features and 2,072 proteins for
the biosaur2-table arm. This rerun used koth_ff 0.2.0 at commit `f64fe17`.
