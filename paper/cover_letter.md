# Cover letter — Journal of Proteome Research

**To:** The Editors, *Journal of Proteome Research*

**Re:** Submission of a Technical Note

---

**Manuscript title**

uno: A Fast Rust Engine for MS1-Only Protein Identification

**Corresponding author**

Patrick Garrett

Department of Chemical Physiology, The Scripps Research Institute

10550 North Torrey Pines Road, La Jolla, CA 92037, USA

Telephone: [PHONE]

Email: pgarrett@scripps.edu

**Co-author**

John R. Yates III, Department of Chemical Physiology, The Scripps Research
Institute

---

Dear Editors,

Please consider the enclosed manuscript for publication as a Technical Note in
the *Journal of Proteome Research*.

MS1-only proteomics identifies proteins from precursor mass and retention time
without MS2 fragmentation, supporting short gradients and high throughput. We
present uno, an open-source Rust engine for this task. uno reads the same
single-run feature tables as the reference implementation, ms1searchpy, allowing
the engines to be compared on identical feature input, database, and digest.

On six datasets from the DirectMS1 studies, uno had higher median protein
depth on four files at 1% two-proteome entrapment FDR and all six at 5%. The
serialized searches ran several-fold faster in approximately half the peak
memory. Twelve independent Orbitrap Astral files from another laboratory test
transfer beyond the primary instrument and deposits. Both engines are
randomized, and the revised figures show individual draws and medians. Protein
list agreement and empirical FDP are also evaluated across seeds. Fixed-seed
repeatability is an operational benefit, not evidence of lower variability.

The manuscript connects these comparisons to the engine's mass-search,
retention-time calibration, scoring, and protein-inference architecture. It
separates identification depth from internal FDR calibration by scoring each
engine's protein lists against the same absent-species probe. The depth
advantage depends on the dataset and threshold; no universal superiority or
uniform internal FDR calibration is claimed. The Astral and built-in feature
finder comparisons include a smaller three-seed extension, and acquisition
variation is distinguished from search variation. These repeats also reveal
remaining seed-dependent zero-depth failures at low input, which are retained
and described as a limitation.

This work fits the Technical Note format as a concise software advance with a
directly demonstrated advantage over the existing implementation. The source,
trained retention-time model, benchmark scripts, exact comparator patch, and
result archive will be released at the repository and DOI below. The manuscript
also discloses which defaults were tuned on related development data and keeps
those data outside the validation claim.

**Supporting Information.** One PDF containing the software environment and
benchmark inputs, complete per-run measurements, entrapment counts and
calibration curves, and seeding and repeatability checks.

**Length.** The abstract is within the journal's 200-word limit, and the text is
within the Technical Note recommendation. The Experimental Section is reported
separately in accordance with the journal's counting guidance.

**Data and software.** [Repository URL and Zenodo DOI to be inserted on
release.] The mass spectrometry data are available through PRIDE under
PXD015669, PXD023977, and PXD066701.

**Declarations.** This manuscript is original, has not been published
previously, and is not under consideration elsewhere. Both authors approved the
submission and declare no competing financial interest. The work is a secondary
computational analysis of public data and required no additional ethics review.
Use of AI tools during manuscript preparation is disclosed in the
Acknowledgments.

Because John R. Yates III is a co-author and the journal's Editor-in-Chief, we
respectfully request conflict-free handling by an independent editor, with no
editorial involvement by the co-author.

Thank you for your consideration.

Sincerely,

Patrick Garrett

On behalf of both authors
