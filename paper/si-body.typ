// Included from paper.typ; an included file has its own scope. BODY-ONLY by
// the scaffold's convention: the visual SI title block lives in paper.typ at
// the include site, because `just docx` resolves this file verbatim and
// synthesizes its own SI title -- a literal #paper-title here reaches pandoc
// as an unknown identifier and kills the export.
#import "stats.typ": lit, s
#import "assets.typ": fig, tbl

// Supporting Information uses independent S-prefixed counters.
#show figure.where(kind: image): set figure(numbering: n => "S" + str(n))
#show figure.where(kind: table): set figure(numbering: n => "S" + str(n))
#counter(figure.where(kind: image)).update(0)
#counter(figure.where(kind: table)).update(0)
#set heading(numbering: "S1.1")
#counter(heading).update(0)

= Benchmark protocol and software environment <sec:si-protocol>

The main single-run protocol and its entrapment procedure are given in
@sec:methods-benchmark. Three additional ms1searchpy searches per Astral file
used its `-pl` option with an external MS/MS calibration list. The 10 ng runs of
that arm ran concurrently, so only their depth is reported. The Astral
single-run uno values come from the same pinned build as the primary engine
comparison (@sec:methods-benchmark). The outputs of the superseded builds are
retained in the archive, as are the searches of the two 1 ng files removed from
the reported set. The built-in feature-detection arm (@sec:results-raw) searched
the same mzML files with `uno --raw`, the bundled koth profile, and otherwise
identical options. Its wall time includes feature detection, and its per-file
counts, list agreement, and resources are in `benchmark/raw_mode/summary.tsv`.

@tbl:datasets maps the compact dataset names used in the text to the deposited
raw files and local feature tables, and @tbl:versions lists the software
environment and benchmark host.

#figure(
  caption: [Dataset names, deposits, and local feature files. The primary and
    Astral rows support the claims in the main text. run10 and deephela are
    tuning records and are not counted as independent validation.],
  tbl("tbl.datasets"),
) <tbl:datasets>

#figure(
  caption: [Software and hardware used for the main single-run comparisons.],
  tbl("tbl.versions"),
) <tbl:versions>

The primary benchmark used peptides of length 7 to 30, no missed cleavages, and
fixed carbamidomethyl cysteine. The patched ms1searchpy #lit("2.8")-series state
and every uno search option are pinned in the benchmark archive. The only
nondefault uno digest setting aligns its peptide bounds with the comparator:

#raw(read("figdata/uno.toml"), lang: "toml", block: true)

The remaining uno defaults are generated from the configuration structs:

#raw(read("figdata/defaults.toml"), lang: "toml", block: true)

*Retention-time calibration controls on the Astral files.*

Self-calibrated ms1searchpy reported a second-stage DeepLC RT standard deviation
($sigma$) of #s(
  "astral.ms1_stock_rt_sigma_min",
) to #s("astral.ms1_stock_rt_sigma_max") min on the #s("astral.files") Astral
files, with no dependence on input. Its calibration set is every peptide of its
first-pass proteins. On these files most of that set was random matches, and
both DeepLC's per-bin median calibration and the Gaussian $sigma$ fit on the
mixture widened accordingly. On the first Astral file, we ran one search with
each of the shipped options `-es 1`, `-ptol 3`, and `-ci 2 -csc 2`, one with
`-es 1 -ptol 3`, and one with all three together. Every one left $sigma$ above
#lit("1.4") min and the depth between #lit("41") and #lit("63") proteins at 1%
entrapment FDR, inside the self-calibrated range. The mass gate was not the
cause: against the MS/MS-verified peptides, the true feature mass error had a
robust standard deviation of #lit("0.66") ppm, and ms1searchpy's three-$sigma$
gate kept #lit("72")% of true matches against #lit("80")% for uno's gate.

The external list for `-pl` was built from the deposit's DIA-NN report as
described in @sec:methods-benchmark. On the same MS/MS-verified peptides, each
model was calibrated on a random half and checked on the other half. DeepLC
reached a robust standard deviation of #lit("0.69") min, and uno's shipped RT
model reached #lit("0.44") min. The remainder of the self-calibrated $sigma$ is
therefore attributable to the contaminated calibration set, not to the RT model.
The controls, the list-building script, and the one-line `-pl` fix are in the
benchmark archive under `paper/analysis/astral/`.

= Retention-time model <sec:si-rt-model>

The RT network combines residue embeddings, an initial convolution, and three
residual blocks with dilations 1, 2, and 4 (@fig:rtmodel a). Padding is masked
after each convolution. Masked mean and maximum pooling produce a fixed-length
vector, which is joined to eleven global descriptors before a two-layer
regression head. The descriptors supply length, composition, and hydropathy
@kyte1982hydropathy information, including terminal residues (@tbl:rtfeatures).

Training used the Chronologer library's eleven datasets on a harmonized
hydrophobic-index scale. Peptides with modifications other than fixed
carbamidomethyl cysteine were excluded, and repeated sequence observations were
averaged. A sequence hash assigned peptides to disjoint training, validation,
and holdout sets. Training minimized smooth L1 loss with Adam and a
plateau-triggered learning-rate reduction, selecting the epoch with the lowest
validation error. Predictions were mapped to each run's gradient by the
calibration described in the main text.

The corpus filters and holdout errors are shown in @fig:rtmodel b–d. Error was
lowest for peptides of 9 to 15 residues and increased toward the length
extremes, where training examples were less abundant. Accuracy was similar
across source datasets with enough holdout peptides to evaluate.

#figure(
  placement: auto,
  fig("fig.figs2-rt-model", width: 100%),
  caption: [uno's retention-time model. *(a)* The network, with layer shapes
    read from the shipped model file. *(b)* The training corpus by source
    dataset: peptides the amino-acid-only filter keeps, and those it drops for
    carrying another modification. *(c)* Holdout error by peptide length. Bar
    labels are the number of holdout peptides in each bin. *(d)* Holdout error
    by source dataset for sources with at least 200 holdout peptides. The dashed
    line is the overall holdout error. Errors are in hydrophobic-index units,
    the corpus's harmonized scale.],
) <fig:rtmodel>

#figure(
  caption: [The eleven global physicochemical descriptors, computed for the
    holdout peptide #lit("DFEAVANCINAK") and shown as the model receives them.],
  tbl("tbl.rt-features"),
) <tbl:rtfeatures>

= Retention-time predictor benchmark <sec:si-rt-predictor>

The benchmark used the Chronologer holdout, unmodified ProteomeTools/Prosit
holdout peptides, and DeepLC's public #lit("unmod.csv") example set. Peptides
contained 7 to 30 canonical residues and were absent from uno's training and
validation partitions. Fixed carbamidomethyl cysteine was supplied to DeepLC;
other modifications were excluded. Repeated sequence observations were reduced
to their median RT. DeepLC's training membership could not be verified, so
training exclusion is established only for uno.

A sequence hash assigned 20% of each set to calibration and 80% to testing. Both
models used a least-squares affine map fitted only on calibration observations.
DeepLC's piecewise calibration was replaced with this common mapping to compare
predictors under the same calibration procedure. Test errors were normalized by
each set's experimental RT 5th to 95th percentile range (@fig:rtpredictor).
Per-peptide predictions, raw errors, calibration times, model hashes, and
software versions are retained in the benchmark output.

#figure(
  placement: auto,
  fig("fig.figs3-rt-predictor", width: 100%),
  caption: [Retention-time predictor accuracy on the held-out 80% of each test
    set, as a percentage of that set's experimental RT 5th to 95th percentile
    range. *(a, b, c)* Residual (predicted minus observed) against observed RT:
    the line is the median and the band the 10th to 90th percentile range in
    twelve equal-count bins. *(d, e, f)* Cumulative distribution of absolute
    error. The open circle marks the 95th percentile, labeled with its value.],
) <fig:rtpredictor>

Throughput was measured on a hash-ordered in silico tryptic digest of the
benchmark's human FASTA entries, restricted to unique canonical sequences of 7
to 30 residues. Both models predicted the first $N$ peptides for batch sizes
from 100 to #lit("100,000"), with five repeats after a warm-up at the largest
size. Only inference was timed, without calibration, on the same workstation
with each engine's default threading (@fig:rtthroughput).

#figure(
  placement: auto,
  fig("fig.figs4-rt-throughput", width: 100%),
  caption: [Retention-time prediction throughput on a synthetic tryptic peptide
    library. *(a)* Wall time per batch against batch size. *(b)* Peptides
    predicted per second, with the steady-state value at 100,000 peptides
    labeled. Lines are medians of five repeats, and bands span the minimum to
    maximum range.],
) <fig:rtthroughput>

#pagebreak()
#show <tbl:imcohorts>: set block(breakable: false)
#show <tbl:impublished>: set block(breakable: false)
#show <tbl:imsearch>: set block(breakable: false)

= Ion-mobility model and evaluation <sec:si-im-model>

The ion-mobility (IM) model uses the RT sequence CNN and training engine with
separate weights. Charge divided by four is appended to the eleven global
descriptors, giving one prediction per sequence and charge. The shared workflow
supports full-network adaptation from independently assigned sequences and
observed mobility; usage and lineage checks are documented in #lit(
  "docs/im-models.md",
).

Training used the Meier trapped ion mobility spectrometry (TIMS) collisional
cross-section (CCS) corpus @meier2021ccs, distributed through ProteomicsML.
Peptides were restricted to 7 to 30 canonical residues, charges 2 to 4, and
measured neutral mass within #lit("0.02") Da of fixed carbamidomethyl-cysteine
chemistry. Sequence/charge replicates were reduced to their median. CCS was
converted to proportional inverse mobility as
$"CCS" sqrt(m m_g / (m + m_g)) / (1000 z)$, where $m$ is neutral peptide mass in
Da, $m_g = 28.0134$ Da, and $z$ is charge. This output requires
acquisition-specific calibration.

Deterministic sequence hashes kept all charges together and selected #s(
  "im.train_seq",
), #s("im.validation_seq"), and #s("im.test_seq") disjoint training, validation,
and test sequences. These contained #s("im.train_obs"), #s("im.validation_obs"),
and #s("im.test_obs") observations, respectively. Training from scratch used
Adam, smooth L1 loss, learning rate #lit("0.001"), batch size 512, seed 0, and
labels standardized using training statistics. The run allowed 20 epochs and
five epochs without improvement in validation mean absolute error (MAE); epoch
#s("im.best_epoch") was selected.

Without test calibration, holdout MAE was #s("im.mae"), root mean squared error
was #s("im.rmse"), and $R^2$ was #s("im.r2"). Reversing internal residue order
while preserving composition and termini increased MAE from #s(
  "im.order_original",
) to #s("im.order_reversed") on a hash-selected #s("im.order_n")-observation
subset, confirming sensitivity to sequence order. Neither test labels nor search
outcomes selected the model.

== Comparison with IM2Deep

We compared uno with IM2Deep's single-output TIMS predictor
@declercq2025tims2rescore, using #lit("IM2DeepUni.ckpt") from package #lit(
  "2.0.2",
). Fixed carbamidomethyl cysteine was supplied explicitly, and CCS outputs
underwent the same conversion as uno's labels. The cohorts were the reserved
public test and acquisition #lit("A1"), restricted to the training chemistry,
lengths, and charges and excluding all uno training/validation sequences. A1
labels were sequence/charge medians from Sage target identifications with
peptide and protein reported q-values at or below 1%.

A sequence hash assigned 20% of each cohort to calibration, keeping all charges
together (@tbl:imcohorts). For both models, per-charge median
observed-to-predicted ratios were fitted on those observations alone. Errors
were measured on the remaining sequences, with no test-label calibration.

#figure(
  caption: [Partitions used for the matched IM2Deep comparison. Observations are
    unique sequence/charge pairs; all charges of a sequence remain together.
    Every sequence is absent from uno's training and validation lineage.],
  tbl("tbl.im-cohorts"),
) <tbl:imcohorts>

IM2Deep had lower median and tail error on both cohorts (@tbl:impublished). Its
training membership could not be excluded, including possible overlap with the
public Meier corpus. The comparison therefore does not establish accuracy on
data unseen by both models. The uno model's shared RT/IM training workflow
provides consistency, but no accuracy advantage was observed here.

#figure(
  caption: [Held-out prediction error after identical per-charge calibration.
    Median and 95th-percentile (P95) absolute relative errors are percentages of
    each observed mobility value, pooled across charges. $R^2$ is computed on
    the cohort's observed scale. Public training overlap cannot be excluded for
    IM2Deep.],
  tbl("tbl.im-published"),
) <tbl:impublished>

== Use after RT shortlisting

IM is enabled by default after mass/RT shortlisting, with no IM rejection
threshold. Unique sequence/charge pairs with valid observed IM are predicted
once per scoring stage. A per-charge median observed-to-predicted ratio
calibrates all target and decoy candidates equally, requiring at least 20 valid
candidates. The absolute relative residual and an availability flag enter the
classifier. Missing IM skips prediction; #lit("--no-im") disables the feature.
Joint aligned-member evidence does not use it.

The follow-up search experiment paired IM off, an internal tree predictor, and
the CNN across three seeds in 15 min, 50 ng timsTOF Ultra 2 acquisitions #lit(
  "A2",
), #lit("B2"), and #lit("C2"). Feature tables, FASTA, RT model, and settings
were held fixed. The tree was trained on a larger public split and is not
IM2Deep. All draws are reported at 1% empirical entrapment FDP (@tbl:imsearch),
with same-species groups expanded to accessions, mixed-species groups excluded,
and tied scores kept intact. These searches are separate from the pinned primary
and matrix benchmarks.

#figure(
  caption: [Present-species accessions at 1% empirical entrapment FDP in the
    paired IM experiment. Human, yeast, and _E. coli_ were present; absent
    _Arabidopsis thaliana_ supplied entrapment accessions, scaled by the FASTA
    present-to-entrapment ratio. Each row retains all three seeds, including
    failed bootstraps with zero depth. The internal tree arm is not IM2Deep.],
  tbl("tbl.im-search"),
) <tbl:imsearch>

CNN IM improved #s("im.improved") of nine pairs, worsened #s("im.worse"), and
tied #s("im.tied"). All acquisition medians increased, but the predefined
requirement of at least seven improved pairs was not met. The identification
benefit was therefore variable, despite the higher medians.

= Complete primary benchmark results <sec:si-results>

@tbl:variability summarizes all primary random draws. @tbl:perdataset gives the
serialized resource-measurement set: the default-seed uno search, re-measured on
the pinned build, and the three initial ms1searchpy searches per file. The
additional comparator reruns contribute to the depth distribution but not to
this resource table. Individual draw depths are archived in
`benchmark/results/spread.tsv`.

#figure(
  caption: [Primary identification depth across random draws. Each depth cell
    gives the median followed by the full observed range in parentheses. Draws
    are uno seeds or independent ms1searchpy reruns on the same feature file.
    They are not separate acquisitions.],
  tbl("tbl.variability-summary"),
) <tbl:variability>

#figure(
  caption: [Original serialized searches: default-seed uno and the first three
    ms1searchpy runs per file. Identification depth, CPU time, wall time, and
    peak resident memory. Memory is maximum resident set size in MiB reported by
    `/usr/bin/time -v`.],
  tbl("tbl.results-per-dataset"),
) <tbl:perdataset>

= Seeding and repeatability checks <sec:si-reproducibility>

Four stages of the search consume randomness: the decoy shuffle, the RANSAC
retention-time calibration's pair sampling, the cross-validation fold
assignment, and the classifier's row and column subsampling. All four derive
their stream from a single command-line seed. The default value reproduces the
streams of the archived default runs. Any other value perturbs all four stages
together, which is how the seed-to-seed spread in the main text was measured.

The primary comparison used #s("bench.uno_seed_count") uno seeds per file and
#s("bench.ms1_draw_min") to #s("bench.ms1_draw_max") comparator searches. The
smaller Astral and raw-mode extension used seeds #lit("0"), #lit("1"), and #lit(
  "2",
) in every single-run arm. All #s("rerun.searches") searches in that inventory,
were executed on the single pinned build of @tbl:versions rather than carried
over from an earlier one. Re-executing that inventory left #s(
  "rerun.unchanged",
) protein tables byte-identical to those of the build used in earlier drafts.
The #s("rerun.changed") that differ are searches that had previously returned
nothing, #s("rerun.recovered") of which now report proteins (#s(
  "astral.uno_seed_recovered",
) on Astral biosaur2 tables, the rest in the built-in feature-detection arm).
The remaining #s("rerun.empty") search still returns nothing. The comparator was
not rebuilt, so its outputs and timings are carried over unchanged. Seed-zero
searches ran serialized, because their wall time and peak memory are what the
resource panels report. The remaining seeds ran two at a time, and their timings
are recorded but not reported. Each search stores its command, its input hashes,
its output hash and the build's hash beside its protein table, under
`benchmark/results/publication_rerun/`. The per-file, per-feature-source,
per-seed summary is `benchmark/results/variability_extension/summary.tsv`.
Matching seeds between feature sources provides a paired comparison of each
file.

#figure(
  placement: auto,
  fig("fig.figs6-determinism", width: 100%),
  caption: [Primary search variability. *(a, b)* Depth at 1% and 5% entrapment
    FDR across all available draws, with open circles, median ticks, and
    minimum-to-maximum connectors. *(c)* Percentage of the other sampled uno
    seeds whose 1% depth was strictly below the archived default seed's depth.
    The dashed line marks 50%. These finite-sample ranks are not population
    percentiles or evidence of deliberate seed selection.],
) <fig:seedspread>

Fixing the seed is not by itself sufficient for a byte-identical result.
Development-time differences between repeated uno searches traced to hash-map
iteration order, which is randomized per process, reaching floating-point
aggregation steps or first-seen tie-breaks. The affected RT-calibration,
mass-recalibration, feature-competition, and protein-group paths now sort on
value or accession before reduction. Unit tests cover the individual reductions,
and `benchmark/determinism_check.sh` runs two complete searches and
byte-compares their protein, peptide, competition, and group outputs. This
full-file check is part of the release gate for the benchmarked version.

The comparator is seedable on the same terms. ms1searchpy's decoy database is
byte-identical across processes once its module-level generator is seeded, and
its hyperparameter search and classifier accept seeds as well. One further step
is needed and is easy to overlook. Several stages apply a seeded shuffle to a
list built from a Python set of peptide strings, and the order of that set
already depends on randomized string hashing. The interpreter's hash seed must
therefore be fixed before the process starts. We reported the comparator as
shipped, unseeded, because that is how it is used in practice, and we made the
same measurement for uno by varying its seed rather than by holding it fixed.

= Entrapment calibration <sec:si-calibration>

For each engine, @tbl:entrapreported gives the human and yeast proteins on its
internally reported 1% list, the scale-corrected empirical FDP, and the depth at
1% entrapment $q$. The correction factor is #s("bench.human_db") divided by #s(
  "bench.yeast_db",
), or #s("cal.db_scale"). Because the accepted lists contain few yeast proteins,
one additional yeast hit can materially move a single-dataset estimate. The
across-dataset pattern is therefore more stable than any one value.

The per-dataset curves in @fig:fdpcurves show where each reported list departs
from the external entrapment proteome across the whole threshold range. Circles
mark the two operating points the main text quotes: the entrapment FDP each
engine's own 1% and 5% list actually carries. A curve above the identity line is
optimistic, and thin curves show individual seeds or reruns in both engines.
Bold curves are pointwise medians where every search has a nonempty list.
@tbl:entrapreported gives the corresponding counts at the 1% point.

#figure(
  placement: auto,
  fig("fig.figs1-fdp-curves", width: 100%),
  caption: [Each engine's internally reported $q$ against the entrapment FDP of
    that list, one panel per primary dataset. Vertical lines mark reported $q$
    of 1% and 5%, and open circles mark every search's FDP there. Thin curves
    show individual uno seeds (blue) or ms1searchpy reruns (orange), and bold
    curves show pointwise medians. Each panel labels the uno and ms1searchpy
    sample sizes. Shading marks the region where the reported FDR is optimistic.
    Neither the curves nor their ranges are confidence intervals.],
) <fig:fdpcurves>

#figure(
  caption: [Each engine's internally reported 1% protein list scored against the
    common yeast entrapment proteome.],
  tbl("tbl.entrap-reported"),
) <tbl:entrapreported>

= Identification from koth LFQ matrices <sec:si-matrix>

#show table: set text(size: 9pt)
#show table: set par(justify: false)
#show figure: set block(breakable: false)

We compared individual-file searches with direct searches of aligned consensus
rows in two cohorts. The Astral HeLa cohort (#lit("PXD066701")) contains the
#lit("100") ng triplicate. IonStar (#lit("PXD003881")) @shen2018ionstar contains
#lit("4") replicates at each of #lit("5") _E. coli_ concentrations in a constant
human background, with expected ratios #lit("1"), #lit("1.5"), #lit("2"), #lit(
  "2.5",
), and #lit("3"). All #lit("23") feature/hill tables and both LFQ matrices were
regenerated with koth #lit("0.7.0"). Across seeds #lit("0"), #lit("1"), and
#lit("2"), the study comprises #lit("69") individual searches and #lit("6")
consensus-matrix searches. The cohorts were selected before these rerun
outcomes, but both had informed development. Seed ranges describe search
variability conditional on fixed upstream outputs, not biological replication.
Program identities are listed in @tbl:matrixversions, and @fig:matrixcohorts
summarizes identification and ratio recovery.

#figure(
  caption: [Programs used for the complete matrix rerun. Full source snapshots,
    local patches, binary hashes, configurations, and driver amendments are
    recorded in the study archive.],
  tbl("tbl.koth-matrix-versions"),
) <tbl:matrixversions>

== Feature rows and quantities

Feature finding uses uno's bundled Orbitrap profile. Alignment uses automatic
reference selection, and consensus grouping requires at least #lit("2") original
runs with group q at most #lit("0.05"). Group coordinate limits are #lit("20")
ppm, #lit("0.02") of the RT span, and #lit("0.05") inverse-mobility units. LFQ
uses summed grid extraction with half-windows of #lit("0.005") of each run's RT
span and #lit("0.015") inverse mobility units, extraction target-decoy
competition, no upstream normalization, and no averagine projection. IM settings
are inactive on these Orbitrap inputs. The release's exclusive-signal ownership
can suppress shared or ambiguous cells to zero. Group confidence, extraction
confidence, and identification q-values refer to different quantities.

An individual search receives one original feature table. A matrix search uses
`--matrix-long` to receive one mass/RT/charge coordinate set and representative
original quality descriptors per consensus row, together with its LFQ intensity
vector. It does not expand a row into original member observations or add LFQ
cell-quality summaries to the classifier. The two matrices contain #s(
  "matrix.astral_groups",
) Astral and #s("matrix.ionstar_groups") IonStar rows. Databases, digest
settings, RT model, seeds, and matching filters are matched between routes
within each cohort.

== Identification comparison

Same-species protein groups were expanded to accessions and deduplicated;
mixed-species groups were excluded. The individual-file union was ranked by each
accession's minimum reported q-value across files, and matrix accessions by
binomial score. Whole tied-score blocks entered together. Yeast provided
entrapment sequences, with false-discovery proportions scaled by the
present-to-yeast database size. We evaluated the largest accession set meeting
each empirical #lit("1"), #lit("3"), or #lit("5")% cutoff and independently
reproduced its counts and accession set using decimal ranks and integer
arithmetic. These retrospective cutoffs are not deployable cohort q-values. The
FDP of each route's reported #lit("1")% list is therefore shown separately
(@tbl:matrixidentification).

#figure(
  fig("fig.figs7-koth-matrix", width: 100%),
  caption: [Direct consensus-row searches using koth #lit("0.7.0"). *(a)*
    Present-species accessions at empirical entrapment FDP at most #lit("1")%
    for the union of individual searches and the consensus-matrix search. *(b)*
    IonStar protein-group ratios versus condition A after median-ratio
    normalization. Markers summarize three seeds and error bars show their full
    range; dotted lines show expected ratios.],
) <fig:matrixcohorts>

#figure(
  caption: [Identification depth at empirical entrapment cutoffs, with medians
    and full ranges across three seeds. Astral counts human accessions; IonStar
    counts human and _E. coli_ accessions. The final column reports entrapment
    FDP for the independently defined internal #lit("1")% lists.],
  tbl("tbl.koth-matrix-identification"),
) <tbl:matrixidentification>

== Quantitative checks

IonStar quantitation uses summed protein-group abundances accepted at internal
protein q at most #lit("0.01"), independently of the empirical accession
cutoffs. Within each seed, complete protein groups determine label-independent
median-ratio normalization factors. No values are imputed. Endpoint ratios and
within-condition CV require at least #lit("2") positive replicates; the table
reports eligible-group counts and positive protein-group cell coverage. Before
identification and protein rollup, the IonStar matrix contains #s(
  "matrix.positive_cells_millions",
) million positive precursor cells, of which #s(
  "matrix.confident_positive_percent",
)% pass extraction q at most #lit("0.01"). Positive or extracted cells do not by
themselves establish peptide identity in that sample.

The median E/A ratio across seeds is #s("matrix.ecoli_ratio_consensus") for _E.
coli_, compared with the expected #lit("3"), and #s(
  "matrix.human_ratio_consensus",
) for the constant human background. Ratio compression remains despite
exclusive-signal extraction. The high protein-group coverage is conditional on
identified groups and sums across precursors; it does not describe precursor
completeness or sample-level identification confidence
(@tbl:matrixquantification).

#figure(
  caption: [IonStar consensus-matrix quantities, summarized by three-seed
    medians and full ranges. Absolute log2 error compares each group's E/A ratio
    with its species expectation. CV summarizes eligible group-condition pairs;
    positive coverage includes all protein-group sample columns.],
  tbl("tbl.koth-matrix-quantification"),
) <tbl:matrixquantification>

The complete commands, source and output hashes, accepted accession sets, and
per-seed metrics are archived in `benchmark/results/matrix_20260915_koth070/`
and `paper/figdata/koth_matrix_study.json`. A separate four-file IonStar pilot
preceded the publication queue. Four seed-zero timing repetitions reproduced the
corresponding protein outputs, but overlapped other work and do not support a
speed comparison. The earlier joint-member, cell-summary, low-input matrix, and
cohort-scaling experiments remain development records rather than part of this
focused comparison.
