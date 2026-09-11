// Transcribed from benchmark/RESULTS.md, "Depth with vs. without the quality
// columns (honest entrapment)" table, plus the replicate ranges given in the
// paragraph beneath it. run10 features (211,064), human+yeast entrapment.
#table(
  columns: 4,
  align: (left, right, right, left),
  table.header([Engine and input], [1% FDR], [5% FDR], [ms1searchpy replicate range, 1%]),
  [uno, quality columns present], [1309], [1887], [n/a (single seed)],
  [uno, quality columns stripped], [1306], [1772], [n/a (single seed)],
  [ms1searchpy median, columns present], [1155], [1651], [1079--1220],
  [ms1searchpy median, columns stripped], [1062], [1447], [1023--1134],
)
