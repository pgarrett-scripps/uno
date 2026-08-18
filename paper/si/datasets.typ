// TRANSCRIBED BY HAND. Manuscript name -> deposit -> local feature file, read off
// benchmark/data/ and the two PRIDE records. Written because the manuscript names
// runs by number ("09746") while the files on disk are prefixed, so a reader could
// not map one to the other; and because 09744 is a DIFFERENT run from 09746, given
// in a different format, and is not one of the six.
#table(
  columns: 4,
  align: (left, left, left, left),
  stroke: (x: none, y: 0.5pt),
  table.header(
    [Name in text], [Deposit], [Raw file], [Feature file (`benchmark/data/`)],
  ),

  [09746], [PXD015669], [`..._09746_...raw`], [`hela5min_09746.features.tsv`],
  [09766], [PXD015669], [`..._09766_...raw`], [`hela5min_09766.features.tsv`],
  [09734], [PXD015669], [`..._09734_...raw`], [`hela5min_09734.features.tsv`],
  [09778], [PXD015669], [`..._09778_...raw`], [`hela5min_09778.features.tsv`],
  [09950], [PXD015669], [`..._09950_...raw`], [`hela5min_09950.features.tsv`],
  [LUM2_05433_VG], [PXD023977], [`LUM2_05433_VG.raw`],
  [`faims_LUM2_05433_VG.features.tsv`],

  table.hline(stroke: 0.5pt),
  [run10 (tuning set)], [PXD015669], [`..._run10_...raw`], [`run10.features.tsv`],
  [QEHF1_09744_JB], [PXD015669], [`QEHF1_09744_JB.raw`], [`09744.features.tsv`],
  [deephela], [PXD004732], [240-min HeLa], [`deephela.features.tsv`],
)
