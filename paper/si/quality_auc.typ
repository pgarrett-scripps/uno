// Transcribed from benchmark/RESULTS.md, the two AUC tables in "Feature
// analysis". These values exist only in that record: the diagnostic scripts
// that produced them were removed after use (RESULTS.md, final paragraph), so
// unlike the ablation depths they could not be recomputed from disk. AUC is
// for separating confident-target from decoy matches, where 0.5 is no
// separation. The last two rows report the distance from 0.5 rather than the
// AUC itself, as in the source.
#table(
  columns: 3,
  align: (left, right, right),
  table.header([Per-match feature], [uno], [ms1searchpy]),
  [`cosine_score`], [0.521], [0.514],
  [`isotope_score`], [0.521], [0.514],
  [`combined_score`], [0.526], [0.517],
  [`ppm_error`], [0.503], [0.502],
  [`neutron_offset`], [0.500], [0.500],
  table.hline(),
  [retention time (|AUC $-$ 0.5|)], [0.136], [0.080],
  [full match classifier (|AUC $-$ 0.5|)], [0.171], [0.134],
)
