// PXD033350, 8 spike levels x 3 replicates, uno --matrix --profile-reassign on
// the averagine-projection (koth_align) matrix -- the standard uno quant path.
// rho = Spearman correlation of protein abundance against log spike amount,
// computed by benchmark/quant/dilution_track.py (true statistics.median), n =
// proteins quantified in >= 8 of the 24 runs. Measured 2026-07-28 on the current
// koth feature finder, taken from benchmark/quant/results/fig5_stdout.txt (panel
// b) / analysis_stdout.txt. The averagine estimator is used because, under the
// current koth, the raw box-sum regressed the honest protein-ID count by 21%
// (306 -> 243) while averagine held it (305 -> 326); dilution tracking itself is
// stable across both estimators.
#table(
  columns: 5,
  align: (left, left, right, right, right),
  table.header([Species], [Role in the design], [n], [Median $rho$], [Fraction $rho > 0.5$]),
  [UPS1], [spike, 0.1--50 fmol], [18], [+0.71], [89%],
  [_E. coli_], [constant background], [303], [+0.02], [4%],
  [_Homo sapiens_], [known absent, entrapment probe], [5], [+0.25], [0%],
)
