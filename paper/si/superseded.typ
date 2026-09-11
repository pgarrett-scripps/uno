// Retired measurements, kept for the record. All were made on an earlier
// development configuration of uno (koth-processed features, koth-aligned
// human + E. coli benchmark matrix, distilled MLP retention-time model) and
// are superseded by the six-dataset benchmark in the main text. Values
// transcribed from the corresponding sections of the previous draft of
// paper/paper.typ (git history), which is their only source on disk.
#table(
  columns: 3,
  align: (left, left, left),
  table.header([Measurement], [Value], [Superseded by]),
  [ms1searchpy protein count, six repeats on one koth feature set],
  [152, 153, 155, 188, 203, 303],
  [Current reproducibility benchmark],

  [uno protein count on the same feature set],
  [175--179],
  [Current reproducibility benchmark],

  [Frozen RT model with monotone calibration, koth-aligned matrix],
  [39 proteins at 1%],
  [Current head-to-head benchmark],

  [Fresh composition-feature ridge refit on the same matrix],
  [9 proteins at 1%],
  [Current head-to-head benchmark],

  [Readout fine-tune over the frozen embedding, same matrix],
  [59 proteins at 1%],
  [Current head-to-head benchmark],

  [Ridge-penalty sweep of the readout fine-tune],
  [$lambda = 1 -> 59$, $10 -> 56$, $100 -> 19$, $1000 -> 4$],
  [n/a (retired)],

  [Entrapment of the frozen and fine-tuned RT paths],
  [34 and 38 real proteins, 0 yeast in both],
  [Current entrapment-calibration analysis],
)
