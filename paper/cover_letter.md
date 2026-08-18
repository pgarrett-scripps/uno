# Cover letter — Journal of Proteome Research

**To:** The Editors, *Journal of Proteome Research*

**Re:** Submission of an original Article

---

**Manuscript title**

uno: A Fast, Reproducible, Dependency-Light Rust Engine for MS1-Only Protein
Identification and Quantitation

**Corresponding author**

Patrick Garrett
Department of Chemical Physiology, The Scripps Research Institute
10550 North Torrey Pines Road, La Jolla, CA 92037, USA
Telephone: [PHONE]
Email: pgarrett@scripps.edu

**Co-author**

John R. Yates III (Department of Chemical Physiology, The Scripps Research Institute)

---

Dear Editors,

Please consider the enclosed manuscript for publication as an Article in the
*Journal of Proteome Research*.

MS1-only proteomics identifies peptides and proteins from precursor mass and
retention time alone, without MS2 fragmentation. It is the route to very short
gradients and high sample throughput, and the open-source reference implementation
for it, ms1searchpy, has been developed and published in this journal. We present
uno, an MS1-only search and quantitation engine written in Rust. uno reads the same
single-run feature files that ms1searchpy reads, which lets the two engines be
compared on genuinely identical input rather than on separately prepared pipelines.

Two properties motivated the work. The first is reproducibility. The reference
engine's decoy shuffle and hyperparameter search draw on an unseeded random number
generator, so repeating a search of the same feature file returns a materially
different protein list; across our benchmark that spread reached 40% of the
engine's own median. A single reported number from such a pipeline is a draw rather
than a measurement. uno is deterministic by construction and returns byte-identical
output across repeated runs, which is what makes its single column a fair
comparator for a median of three. The second is dependency weight. uno predicts
retention time in-process with a pure-Rust convolutional model trained on
experimental retention times, so the Python and deep-learning stack that dominates
single-run wall-clock time in the reference pipeline is not required at search time.

On six datasets taken from the DirectMS1 papers' own deposits, uno identified more
proteins at a 1% entrapment false-discovery rate than the median of three
ms1searchpy replicates, while running several times faster in roughly half the peak
memory. We report depth by a two-proteome entrapment probe scored identically for
both engines rather than by either engine's own target/decoy estimate, and we use
the same probe to ask what each engine's reported 1% list is actually worth. We
think that second question is the more useful contribution: the probe shows a
learned protein-level rescorer inflating the empirical false-discovery rate while
raising apparent counts, which is why uno's protein score is a transparent binomial
model. We report that as a cautionary result rather than as a feature.

We have tried to be explicit about the limits. Several of uno's defaults were
selected on development data drawn from the same family as the benchmark, and the
manuscript says so and identifies which. The quantitation path is validated for
dilution tracking, which is robust, and not for trace-level differential testing,
where it trails the reference method in both sensitivity and precision. We report
that gap with its measured numbers rather than omitting it.

This work fits *Journal of Proteome Research* in both subject and form. The journal
publishes the MS1-only method papers this tool is measured against, the entrapment
and false-discovery-rate methodology we rely on, and the open-source proteomics
software this manuscript sits among.

**Supporting Information.** A single PDF containing Sections S1 to S8, Figures S1
to S5, and Tables S1 to S8. It covers the benchmark protocol with software
versions, hardware and the exact commands; the full per-dataset and per-replicate
results; entrapment calibration of each engine's own reported list with
per-dataset false-discovery-proportion curves; the FAIMS benchmark across all three
replicate runs; resource scaling against input size; the per-feature quality-column
ablation with per-column target/decoy separation; extended quantitation detail; and
the retired measurements this revision supersedes.

**Length.** The manuscript complies with the journal's limits. The abstract is 200
words. [Confirm main-text word count and figure/table counts against the current
guidelines at submission.]

**Data and software.** [Repository URL and Zenodo DOI to be inserted on release.]
The mass spectrometry data reanalyzed here are publicly available from
ProteomeXchange via PRIDE under accessions PXD015669, PXD023977 and PXD033350.

**Declarations.** This manuscript is original, has not been published previously,
and is not under consideration elsewhere. All authors have approved the submission.
The authors declare no competing financial interest. The study is a secondary
computational analysis of publicly deposited data and required no ethics review.
Use of AI tools during preparation is disclosed in the Acknowledgments.

We suggest the manuscript is best handled by an editor with an interest in
proteomics software, false-discovery-rate methodology, or high-throughput
bottom-up workflows.

Thank you for your consideration.

Sincerely,

Patrick Garrett
On behalf of both authors
