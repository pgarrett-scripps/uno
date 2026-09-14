// =============================================================================

#import "stats.typ": lit, s
// Everything project-specific about the manuscript's identity lives here.
//
// Single source of truth: paper.typ imports it for the PDF and for the Word
// front matter, wordcount.typ counts the abstract out of it, and audio/config.py
// reads the title straight out of this file so the narration can never announce
// a title the paper no longer has.
// =============================================================================

#let paper-title = "uno: A Fast Rust Engine for MS1-Only Protein Identification"

// Short form used on the audiobook cover art. Keep it to a couple of words.
#let paper-wordmark = "uno"

// Shown under the wordmark on the cover; \n breaks the line.
#let paper-cover-subtitle = "MS1-Only Protein Identification\nin Rust"

// ORCID iDs are linked via the corresponding author's ACS Paragon Plus profile
// (per ACS policy, they are not typed into the manuscript).
#let paper-authors = (
  (
    name: "Patrick Garrett",
    email: "pgarrett@scripps.edu",
    affiliation: "Department of Chemical Physiology, The Scripps Research Institute, 10550 North Torrey Pines Road, La Jolla, CA 92037, USA",
  ),
  (
    name: "John R. Yates III",
    email: "jyates@scripps.edu",
    affiliation: "Department of Chemical Physiology, The Scripps Research Institute, 10550 North Torrey Pines Road, La Jolla, CA 92037, USA",
  ),
)

#let paper-keywords = (
  "MS1-only proteomics",
  "DirectMS1",
  "retention time prediction",
  "reproducibility",
  "Rust",
  "open-source software",
)

#let paper-date = "August 2026"

// Shown on the audiobook cover under the author line.
#let paper-institution = "The Scripps Research Institute"

// The bibliography style. Typst ships CSL styles by name; JPR is an ACS journal.
#let paper-bib-style = "american-chemical-society"

// The abstract. Kept as its own binding, rather than inline in the template
// call, because three separate consumers slice it out of this file by name: the
// Word export path, the word counter (JPR caps the abstract at 200 words
// separately), and the audiobook narrator.
#let paper-abstract = [
  MS1-only protein identification uses precursor mass and retention time without
  MS/MS fragmentation. We present uno, an open-source Rust engine that accepts
  the same feature tables as ms1searchpy. Across #s("bench.dataset_count")
  DirectMS1-study datasets, uno had greater median depth on #s(
    "bench.uno_median_wins",
  ) datasets at 1% two-proteome entrapment false discovery rate (FDR) and #s(
    "bench.uno_median_wins_5pct",
  ) at 5%. Searches run one at a time were #s("bench.speed_factor_min") to #s(
    "bench.speed_factor_max",
  ) times faster and used about half the peak memory. Both engines showed search
  variability. Across #s("bench.uno_seed_count") seeds per primary file, uno's
  full depth range was #s("bench.uno_seed_spread_1pct_pct_min") to #s(
    "bench.uno_seed_spread_1pct_pct_max",
  )% of its median at 1% entrapment FDR. Fixed seeds provided repeatability
  rather than evidence of lower variability. On #s("astral.files") independent
  HeLa Orbitrap Astral files, the median of file-median uno depths ranged from
  #s("astral.uno_seed_level_median_min") to #s(
    "astral.uno_seed_level_median_max",
  ) across input levels. ms1searchpy's own retention-time calibration failed on
  that liquid chromatography method. External MS/MS calibration improved its
  depth. One low-input search using uno's built-in feature detection returned no
  proteins. Protein-list overlap and empirical false discovery proportion (FDP)
  varied in both engines. uno provides a fast alternative for short-gradient
  MS1-only identification, with external entrapment calibration defining its
  present scope.
]

// -----------------------------------------------------------------------------
// Derived values. Nothing to edit below here.
// -----------------------------------------------------------------------------

// Unique affiliations in first-appearance order, so the Word front matter can
// number them the way the PDF template does. Deriving this rather than typing a
// second author line by hand removes the drift that a duplicated list invites.
#let paper-affiliations = {
  let seen = ()
  for a in paper-authors {
    if a.affiliation not in seen { seen.push(a.affiliation) }
  }
  seen
}

#let affiliation-number(affil) = (
  paper-affiliations.position(x => x == affil) + 1
)

// "Patrick Garrett^1, John R. Yates III^1" for the Word front matter, which has
// no template to build an author line for it. Derived rather than retyped, so the
// superscript markers cannot drift out of step with the PDF.
// Wrapped in a code block because a method chain broken across lines after
// `#let x =` would otherwise end at the first newline.
#let paper-author-line = {
  paper-authors
    .map(a => a.name + super(str(affiliation-number(a.affiliation))))
    .join(", ")
}

// "Garrett, Yates" -- used as the audiobook artist tag and on the cover.
//
// The last whitespace-separated token is NOT the surname when the name carries a
// generational suffix: "John R. Yates III" would otherwise give "III". Drop
// trailing suffix tokens first. audio/config.py has the same list, and must.
#let name-suffixes = ("Jr", "Jr.", "Sr", "Sr.", "II", "III", "IV")
#let surname-of(name) = {
  let parts = name.split(" ")
  while parts.len() > 1 and parts.last() in name-suffixes {
    let _ = parts.pop()
  }
  parts.last()
}
#let paper-surnames = paper-authors.map(a => surname-of(a.name))

// Restated on the Supporting Information title page.
#let si-authors = paper-authors.map(a => a.name).join(", ")
