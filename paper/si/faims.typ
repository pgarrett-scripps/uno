// Measured 2026-07-28, not transcribed. Source:
// benchmark/faims/faims_out/figS2_stdout.txt, the console output of
// paper/scripts/figures/figS2_faims.py, which scores both engines' protein
// tables in benchmark/faims/faims_out/ with benchmark/faims/score_entrap.py.
// Replicate order is LUM2_05433 / 05434 / 05435. Database sizes for the
// entrapment scale factor come from the target-only FASTA that
// benchmark/faims/run_faims.sh builds: 20,193 human + 4,431 E. coli.
// The ms1searchpy row supersedes the 1432 / 1883 / "1659 / 1136 / 1503"
// recorded in benchmark/faims/README.md. Those did not reproduce, because
// run_faims.sh used relative MS1/ROOT/OUT paths inside a `cd "$OUT"`, so the
// committed script never actually invoked ms1searchpy. The bug is fixed.
#table(
  columns: 4,
  align: (left, right, right, left),
  table.header([Engine], [1% FDR], [5% FDR], [Per-replicate, 1%]),
  [uno], [1591], [2060], [1642 / 1460 / 1672],
  [ms1searchpy], [1474], [1915], [1662 / 1205 / 1554],
)
