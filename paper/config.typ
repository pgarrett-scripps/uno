// =============================================================================
// Everything project-specific about the manuscript's identity lives here.
//
// Single source of truth: paper.typ imports it for the PDF and for the Word
// front matter, wordcount.typ counts the abstract out of it, and audio/config.py
// reads the title straight out of this file so the narration can never announce
// a title the paper no longer has.
// =============================================================================

#let paper-title = "uno: A Fast, Reproducible, Dependency-Light Rust Engine for MS1-Only Protein Identification"

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
  "label-free quantification",
)

#let paper-date = "July 2026"

// Shown on the audiobook cover under the author line.
#let paper-institution = "The Scripps Research Institute"

// The bibliography style. Typst ships CSL styles by name; JPR is an ACS journal.
#let paper-bib-style = "american-chemical-society"

// The abstract. Kept as its own binding, rather than inline in the template
// call, because three separate consumers slice it out of this file by name: the
// Word export path, the word counter (JPR caps the abstract at 200 words
// separately), and the audiobook narrator.
#let paper-abstract = [
  MS1-only protein identification (DirectMS1) infers peptide and protein
  identities from precursor mass and retention time (RT) alone, without MS2
  fragmentation, enabling short gradients and high throughput. The reference
  implementation, ms1searchpy, is written in Python with a deep-learning RT
  predictor. We present uno, an open-source MS1-only search and identification
  engine in Rust. uno reads the same single-run feature files, so both engines
  see identical input. On six datasets from the DirectMS1 papers, uno identified
  more proteins at 1% entrapment false-discovery rate (FDR) than the median of
  three ms1searchpy replicates, by 3.0 to 39.7%. It ran 5.3 to 7.0 times faster
  in half the peak memory. On three independent 10 ng HeLa Orbitrap Astral
  runs, uno returned 1006 to 1448 proteins at 1% entrapment FDR, against
  ms1searchpy medians of 45 to 47, while running approximately 10 times faster.
  uno is byte-deterministic. ms1searchpy's depth on
  identical input varies by up to 40.5% of its median. Two-proteome entrapment
  finds uno's reported 1% protein FDR conservative on four of six datasets and
  about twice optimistic on one. On the DirectMS1Quant UPS1 dilution
  series, uno-derived abundances tracked the spike design (median Spearman
  $rho = +0.71$); trace-level differential testing remains outside the validated
  scope.
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
