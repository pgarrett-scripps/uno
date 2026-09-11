// TRANSCRIBED BY HAND. Manuscript name -> deposit -> local feature file, read off
// benchmark/data/, the PRIDE records, and the Astral validation directory
// (12-file input series; Astral feature files live in the external data root).
// Written because the manuscript uses compact run names while the files on disk
// carry longer acquisition names.
#set text(size: 8pt)
#show raw: it => text(
  font: "DejaVu Sans Mono", size: 7pt,
  it.text.replace("_", "_\u{200b}").replace(".", ".\u{200b}"),
)
#table(
  columns: (18%, 15%, 34%, 33%),
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

  [Astral 200 pg rep1], [PXD066701], [`CW0710_OldFAIMS_200pg__RA9_rep1_ch1_DIA200pg.raw`],
  [`...RA9_rep1_ch1_DIA200pg.features.tsv`],
  [Astral 200 pg rep2], [PXD066701], [`CW0710_OldFAIMS_200pg__RA9_rep2_ch2_DIA200pg.raw`],
  [`...RA9_rep2_ch2_DIA200pg.features.tsv`],
  [Astral 200 pg rep3], [PXD066701], [`CW0710_OldFAIMS_200pg__RA9_rep3_ch1_DIA200pg.raw`],
  [`...RA9_rep3_ch1_DIA200pg.features.tsv`],
  [Astral 1 ng rep1], [PXD066701], [`CW0710_OldFAIMS_1ng__RA7_rep1_ch1_DIA200pg.raw`],
  [`...RA7_rep1_ch1_DIA200pg.features.tsv`],
  [Astral 1 ng rep2], [PXD066701], [`CW0710_OldFAIMS_1ng__RA7_rep2_ch2_DIA200pg.raw`],
  [`...RA7_rep2_ch2_DIA200pg.features.tsv`],
  [Astral 1 ng rep4], [PXD066701], [`CW0710_OldFAIMS_1ng__RA7_rep4_ch2_DIA200pg.raw`],
  [`...RA7_rep4_ch2_DIA200pg.features.tsv`],
  [Astral 10 ng rep1], [PXD066701], [`CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng.raw`],
  [`...RA4_rep1_ch1_DIA10ng.features.tsv`],
  [Astral 10 ng rep4], [PXD066701], [`CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng.raw`],
  [`...RA4_rep4_ch2_DIA10ng.features.tsv`],
  [Astral 10 ng rep6], [PXD066701], [`CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng.raw`],
  [`...RA4_rep6_ch2_DIA10ng.features.tsv`],
  [Astral 100 ng rep2], [PXD066701], [`CW0710_OldFAIMS_100ng__RA1_rep2_ch2_DIA10ng.raw`],
  [`...RA1_rep2_ch2_DIA10ng.features.tsv`],
  [Astral 100 ng rep3], [PXD066701], [`CW0710_OldFAIMS_100ng__RA1_rep3_ch1_DIA10ng.raw`],
  [`...RA1_rep3_ch1_DIA10ng.features.tsv`],
  [Astral 100 ng rep4], [PXD066701], [`CW0710_OldFAIMS_100ng__RA1_rep4_ch2_DIA10ng.raw`],
  [`...RA1_rep4_ch2_DIA10ng.features.tsv`],

  table.hline(stroke: 0.5pt),
  [run10 (tuning set)], [PXD015669], [`..._run10_...raw`], [`run10.features.tsv`],
  [deephela], [PXD004732], [240-min HeLa], [`deephela.features.tsv`],
)
