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
#set heading(numbering: n => "S" + str(n))
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

@fig:rtmodel a shows the network architecture. A peptide is a sequence of
residue tokens, and a learned embedding maps each token to a vector. An initial
convolution layer expands that to the working number of channels. Three residual
blocks with dilations 1, 2, and 4 then widen the context each position sees
without adding layers. Padding positions are masked after every convolution, so
a position never reads beyond the end of the peptide and inference runs over the
real length rather than a fixed maximum. Masked mean and maximum pooling over
positions give a fixed-length vector for a peptide of any length. The eleven
global descriptors are concatenated at that point, and a two-layer fully
connected output returns a single predicted value. The whole model is 41,665
parameters.

The pooled vector is an average and a maximum over positions, which destroys two
kinds of information that retention depends on. A sum is not recoverable from a
mean without the length, and position-weighted effects vanish when every
position is weighted equally. The descriptors restore exactly those: total
hydropathy alongside its mean, and the hydropathy @kyte1982hydropathy of the
first and last residue alongside whole-peptide composition. @tbl:rtfeatures
lists the vector for one holdout peptide, with each value as the model receives
it. The table is generated by importing the training code's own feature
function, and the generator fails if the feature count there disagrees with the
Rust inference path, so the two implementations cannot drift apart silently.

The corpus is the Chronologer harmonized library, whose eleven source datasets
were acquired on different gradients and mapped onto a single hydrophobic-index
scale. uno's model encodes amino acids only, with cysteine implicitly
carbamidomethylated, so peptides carrying any other modification were dropped.
@fig:rtmodel b shows what that filter costs per source. Two sources are removed
almost entirely, one phosphorylation library and one isobaric-label library, and
the tryptic and human leukocyte antigen peptide libraries survive largely
intact. Repeated observations of one peptide were averaged, and a hash of its
sequence assigned it to the training, validation, or holdout partition. No
peptide appears in more than one partition. Training minimized a smooth L1
(Huber) loss with the Adam optimizer and a plateau-triggered learning-rate
reduction, selecting the epoch with the best validation error.

@fig:rtmodel c and d show where the model is weak. Error is lowest for peptides
of 9 to 15 residues, which dominate a tryptic digest, and rises at both ends:
short peptides carry little sequence to condition on, and peptides beyond 30
residues are rare in the corpus. Accuracy is similar across the source datasets
that contribute enough holdout peptides to measure, so no single contributor
carries the result. The model does not predict minutes. It predicts a position
on the harmonized scale, and the per-run calibration described in the main text
supplies the mapping into the observed gradient.

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

The predictor benchmark used the peptide-disjoint Chronologer holdout,
unmodified ProteomeTools/Prosit holdout peptides, and DeepLC's public #lit(
  "unmod.csv",
) example set. Every peptide contained only canonical residues, was 7 to 30
residues long, carried no variable modification, and was absent from the
Chronologer train and validation partitions. Fixed carbamidomethyl cysteine was
supplied to DeepLC and is implicit in uno. Repeated sequence observations were
represented by their median experimental RT. The exact peptide membership of the
packaged DeepLC training corpus is not available here, so the exclusion is
established for uno but cannot be established for DeepLC.

A stable sequence hash assigned 20% of each set to calibration and the remainder
to testing. Every error in @fig:rtpredictor is divided by the 5th to 95th
percentile range of the corresponding experimental labels, allowing comparison
across the three RT scales. Both models received the same calibration: a
least-squares affine map from the model's native output scale to the observed
scale, fitted on the 20% partition. DeepLC's own piecewise calibration, which
ms1searchpy applies in-search, was not used here because with about 200
calibration peptides it extrapolates poorly at the ends of the elution range and
would measure the calibrator rather than the model. Prediction throughput was
measured separately (@fig:rtthroughput). The complete accuracy output also
records per-peptide test predictions, calibration time, raw MAE, 95th-percentile
absolute error, $R^2$, Spearman correlation, package versions, model hash, and
host environment.

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

Throughput depends on the peptides, not on which experimental set they came
from, so it was measured on one synthetic library. The library is an in silico
tryptic digest of the human entries in the benchmark FASTA, restricted to
canonical residues, 7 to 30 residues, and unique sequences, deterministically
ordered by sequence hash. Both models predicted the first $N$ peptides of that
library for $N$ from 100 to 100,000, with five repeats each after a warm-up call
at the largest size. Only inference was timed, with no calibration, on the same
workstation with each engine's default threading. The wall time per batch in
@fig:rtthroughput is flat at small batches, where each engine's fixed per-call
overhead dominates, and rises in proportion to batch size once that overhead is
amortized. The steady-state throughput is the value at the largest batch.

#figure(
  placement: auto,
  fig("fig.figs4-rt-throughput", width: 100%),
  caption: [Retention-time prediction throughput on a synthetic tryptic peptide
    library. *(a)* Wall time per batch against batch size. *(b)* Peptides
    predicted per second, with the steady-state value at 100,000 peptides
    labeled. Lines are medians of five repeats, and bands span the minimum to
    maximum range.],
) <fig:rtthroughput>

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

= Joint identification from aligned koth feature matrices <sec:si-matrix>

#show table: set text(size: 9pt)
#show table: set par(justify: false)

This study evaluates one uno search over a koth LFQ bundle against the union of
separate searches over its constituent feature tables. All features, alignments,
LFQ extractions, and searches were rerun with koth #lit("0.5.0") and a frozen
uno build. The Astral HeLa input series (#lit("PXD066701")) contains #lit("3")
files at each of #lit("4") input levels: #lit("200") pg, #lit("1") ng, #lit(
  "10",
) ng, and #lit("100") ng. Each level forms a matrix, and a further matrix pools
all #lit("12") files. IonStar (#lit("PXD003881")) @shen2018ionstar supplies
#lit("4") replicates at each of #lit("5") _Escherichia coli_ levels in a
constant human background, with expected _E. coli_ factors of #lit("1"), #lit(
  "1.5",
), #lit("2"), #lit("2.5"), and #lit("3"). Nested IonStar cohorts contain #lit(
  "5",
), #lit(
  "10",
), and #lit("20") files, balanced across conditions. These datasets were
examined during method development; the rerun measures reproducibility across
the complete fixed cohorts rather than performance on previously unseen data.

The fixed inventory comprises #lit("32") fresh feature tables, #lit("8") LFQ
bundles, and #lit("144") identification searches across seeds #lit("0"), #lit(
  "1",
), and #lit("2"). Every planned seed and cohort is included. Features and
alignment were generated once per cohort, so seed ranges measure search
variability conditional on those upstream outputs. Nested cohorts reuse files,
and different seeds reuse the same measurements; neither supplies independent
biological replication. Commands, input hashes, source snapshots, binaries, and
successful-output hashes are archived in
`benchmark/results/matrix_20260913_koth050/`. The generated evidence bundle is
`paper/figdata/koth_matrix_study.json`. These builds are separate from the
primary single-file comparison (@tbl:matrixversions).

#figure(
  caption: [Frozen matrix-study programs. Full source manifests and binary
    hashes are archived; the uno source snapshot includes local changes, so HEAD
    alone does not identify that build.],
  tbl("tbl.koth-matrix-versions"),
) <tbl:matrixversions>

*Measurements and search routes.*

Feature finding uses the profile bundled with uno. koth groups observations
using bounded mass, retention-time, and ion-mobility spans, then applies final
group filters. This study permits single-run groups and replicated groups
without a strong individual seed, retaining the fixed member and strong-seed
score floors of #lit("0.5") and #lit("0.75"). Alignment selects its reference
automatically. LFQ uses summed uniform-grid intensities, extraction target-decoy
competition, and no upstream normalization or averagine projection. These
settings were fixed before the rerun and were identical across the two datasets.
Replicated weak-seed retention is an explicit experimental setting rather than a
new release default.

The long bundle separates original feature observations from recomputed LFQ
cells. Matching uses the aligned coordinates of original observations while
retaining their own isotope descriptors, feature scores, and apex intensities.
Cell intensities and extraction scores describe integration at the aligned
location, including match-between-runs (MBR) fills. The exported hybrid
extraction score is the geometric mean of RT, intensity, isotope-pattern, and
coelution scores; it does not replace the original feature score. MBR fills have
no original feature identity. A pooled protein identification therefore does not
establish detection in every sample.

Single-file searches consume the original koth feature tables. Both joint arms
consume the long bundle with `--matrix-long` and `--matrix-joint`, sharing
digestion and RT prediction across the cohort (@tbl:matrixroutes). Candidate
matching searches the union of original member mass windows within each aligned
group. No separate-search protein scores or accepted lists enter either joint
search. Groups below the strong-seed threshold are excluded from initial mass/RT
calibration but remain available for identification.

For each peptide candidate, the joint scorer selects at most one original
observation per source run, preferring a match inside the mass gate and then the
smallest joint mass/RT deviation. The deviation is the sum of squared
standardized mass and RT residuals. An observation is compatible when it passes
the mass gate and its deviation is at most #lit("9"). The selected observation
supplies the original-quality columns. A further #lit("14") descriptors
summarize direct and compatible run counts, compatible fraction, second and
median deviation, median and median absolute deviation of mass and RT residuals,
median feature quality, isotope and scan counts, selected quality, and a
second-run indicator. Counts use log-one-plus transforms; deviations are capped
at #lit("100") for summary features, with a missing second observation assigned
that cap and distinguished by its indicator. These are learned classifier inputs
rather than independent replicate likelihoods. Alternative observations cannot
add repeated support from the same source run, and an MBR fill does not create
another peptide observation.

The second arm additionally supplies #lit("12") extraction-cell summaries:
filled MBR count and fraction, confident MBR count and overall confident
fraction, and mean q, worst q, q availability, and confident fraction separately
for direct and MBR cells. Filled cells have positive intensity; confidence means
extraction q at most #lit("0.01"). Cohort-wide fractions divide by the number of
distinct source tables. Within direct or MBR strata, q availability and
confident fractions divide by the number of positive cells in that stratum.

Mean and worst q use available values only. When none is available in a stratum,
both summaries are #lit("1"), distinguished by zero availability; q equal to
zero is retained. The same feature-derived summaries are available to matched
target and sequence-decoy candidates. This comparison measures the combined
effect of added covariates, including changed classifier feature sampling; it
does not isolate a causal benefit of any individual column.

#figure(
  caption: [Fresh input routes. The separate-search union is an evaluation
    comparator; its outputs are never supplied to a joint search. Ordinary raw
    mode was not rerun as an additional arm of this study.],
  tbl("tbl.koth-matrix-routes"),
) <tbl:matrixroutes>

*Matched-error identification comparison.*

Protein-group membership can change between searches. We expand same-species
protein groups to accessions, deduplicate them, and exclude mixed-species groups
from species-specific counts. Human is present in Astral; human and _E. coli_
are present in IonStar. Yeast supplies absent-species entrapment. The correction
is
$
  "FDP" = n_"yeast" / n_"present" dot N_"present" / N_"yeast",
$
where lowercase counts refer to accepted accessions and uppercase counts to
entries in the searched FASTA. Joint lists are ranked by protein score; the
separate-search union is ranked by the minimum internal q for each accession
across its constituent files. Whole tied-score blocks enter together. For each
route we report the deepest cutoff with empirical entrapment FDP at most #lit(
  "1",
)%, together with its yeast count and achieved FDP. This retrospective cutoff
uses the known absent species for evaluation and is not a deployable cohort q
value. Simply combining individually filtered lists would not establish cohort
FDR control. Counts and FDP at internal protein q at most #lit("0.01") are
retained separately in the evidence bundle.

@fig:matrixcohorts presents Astral first, with all input levels retained,
followed by the nested IonStar cohorts. Identification summaries show the median
and full range across the three seeds, not confidence intervals. At #lit("100")
ng, joint searching recovered #s("matrix.astral_joint_median") human accessions
(range #s("matrix.astral_joint_min")–#s("matrix.astral_joint_max")), versus
#s("matrix.astral_union_median")
(#s("matrix.astral_union_min")–#s("matrix.astral_union_max")) for the separate
union. The cell-summary arm recovered #s("matrix.astral_joint_cell_median") (#s(
  "matrix.astral_joint_cell_min",
)–#s("matrix.astral_joint_cell_max")). Joint searching also exceeded the union
in all three #lit("200") pg comparisons, with median counts of #s(
  "matrix.astral_low_joint_median",
) and #s("matrix.astral_low_union_median"), respectively. The #lit("1") ng
comparison favored separate searching in every seed, whereas the #lit("10") ng
result varied by seed (@tbl:matrixastral). Pooling all files also gave gains in
some seeds and losses in another. These results support an input-dependent
identification gain, with no consistent additional benefit from the cell
summaries at the two input levels where the joint arm won across all seeds.
Added and lost accessions are reported in @tbl:matrixadditional.

#figure(
  fig("fig.figs7-koth-matrix", width: 100%),
  caption: [Fresh matrix searches with koth #lit("0.5.0"). *(a)* Astral and
    *(b)* IonStar present-species accessions at matched empirical entrapment FDP
    at most #lit("1")%. Points, medians, and full ranges show all three seeds.
    *(c)* One dedicated serialized timing repetition for the Astral #lit("100")
    ng cohort, with a shared warm RT cache. *(d)* _E. coli_ abundance ratios on
    common protein groups in the #lit("20")-run IonStar cohort, using shared
    normalization within each seed; points and ranges summarize the three
    seed-specific group medians. Nested cohorts and repeated search seeds are
    not independent biological replicates.],
) <fig:matrixcohorts>

#figure(
  caption: [Astral accession coverage at matched empirical entrapment cutoffs.
    Entries are three-seed medians with full ranges. MS/MS is cohort-level
    corroboration of accepted human accessions, not sample-specific detection.
    NA denotes an undefined estimate.],
  tbl("tbl.koth-matrix-astral"),
) <tbl:matrixastral>

Astral corroboration uses the deposited DIA-NN report solely for evaluation.
Target precursor and protein-group q must both be at most #lit("0.01"). An
accession is corroborated if it appears anywhere in the same searched cohort;
absence from that report is not proof of a false identification. The evidence
bundle also retains corroboration specifically among accessions added by joint
searching. These comparisons do not supply peptide-level error control or
establish detection in every matrix cell.

#figure(
  caption: [Accessions added and lost relative to the separate-search union at
    each route's matched empirical cutoff. Jaccard compares the accepted
    present-species sets; entries summarize all three seeds.],
  tbl("tbl.koth-matrix-additional"),
) <tbl:matrixadditional>

IonStar provides a complementary mixture-series assessment
(@tbl:matrixidentification). The separate-search union exceeded both joint arms
in every IonStar cohort and seed. For all #lit("20") runs, joint searching with
cell summaries recovered #s("matrix.ionstar_joint_cell_median") present-species
accessions, compared with #s("matrix.ionstar_union_median") for the separate
union and #s("matrix.ionstar_joint_median") for joint searching without cell
summaries. The substantial loss in this larger matrix shows that the Astral gain
does not establish a general advantage across instruments and cohort sizes. The
small accepted yeast counts also make empirical FDP sensitive to individual
errors, even when routes meet the same nominal evaluation threshold.

#figure(
  caption: [IonStar matched-cutoff accession coverage and calibration. Present,
    yeast, and FDP columns refer to the empirical cutoff at most #lit("1")%. The
    final column separately evaluates each internally reported #lit("1")% list,
    including the deduplicated union. Entries are medians and full ranges across
    three seeds.],
  tbl("tbl.koth-matrix-identification"),
) <tbl:matrixidentification>

*Extraction confidence and quantitative evaluation.*

Before identification, the #lit("20")-run IonStar bundle contains #s(
  "matrix.positive_cells_millions",
) million positive target extraction cells; #s(
  "matrix.confident_positive_percent",
)% pass extraction q at most #lit("0.01"). The all-file Astral bundle likewise
distinguishes positive and confident direct and MBR cells (@tbl:matrixsamples).
These are upstream extraction counts, not identified peptides or
protein-supported assignments. Absolute calibration of koth extraction q has not
been independently established, so the exported threshold is used as an evidence
descriptor rather than a verified cell FDR.

#figure(
  caption: [Upstream target extraction cells in the all-file Astral bundle.
    Positive means intensity greater than zero; confident columns additionally
    require koth extraction q at most #lit("0.01"). Counts precede UNO
    identification and do not claim peptide identity or cell-level FDR
    control.],
  tbl("tbl.koth-matrix-sample-support"),
) <tbl:matrixsamples>

Quantitative evaluation uses integrated LFQ abundances for protein groups at
internal protein q at most #lit("0.01"), separately from the empirical accession
lists above. Median-ratio normalization factors are estimated from complete
protein groups in the joint arm and reused for the cell-summary arm within each
seed. This assumes that most groups are stable. We report both all accepted
groups and the exact group-membership intersection between arms, without
imputation. E/A ratios require at least #lit("2") positive replicates at each
endpoint; within-condition CV likewise requires at least #lit("2") positive
replicates. Ratio-group counts and positive-cell coverage accompany these
conditional metrics (@tbl:matrixquantification).

The common _E. coli_ groups yield median E/A ratios of #s(
  "matrix.ecoli_ratio_joint",
) without cell summaries and #s("matrix.ecoli_ratio_joint_cell") with them,
against the expected #lit("3"). Both arms use the same integrated cell values,
so changes in quantitative summaries can reflect identification and group
selection. These scoring covariates do not themselves correct integration bias
or establish improved quantitative accuracy.

#figure(
  caption: [IonStar #lit("20")-run quantitative evaluation, with three-seed
    medians and full ranges. Absolute log2 error compares each group's observed
    E/A ratio with the expected species ratio. CV summarizes eligible
    group-condition pairs; positive coverage includes all sample columns.
    Common-group comparisons share groups and normalization within each seed.],
  tbl("tbl.koth-matrix-quantification"),
) <tbl:matrixquantification>

*Resource measurements.*

The main rerun used up to four concurrent jobs with #lit("8") threads each and a
shared RT-prediction cache; those timings are not used for speed claims. A
separate serialized repetition searched the three Astral #lit("100") ng feature
tables and their joint cell-summary bundle at seed #lit("0"), with the RT cache
already warm. The identification outputs were checked against the corresponding
study outputs. Search wall time was #s("matrix.search_speedup")-fold lower for
the single joint search (@tbl:matrixresources). This is one measured repetition,
without a timing uncertainty estimate, and excludes feature extraction,
alignment, and LFQ construction. It demonstrates the observed search cost with
prediction caching already available to both routes, rather than an end-to-end
speedup or a hardware-independent performance guarantee.

#figure(
  caption: [Dedicated serialized Astral #lit("100") ng search timings with #lit(
      "8",
    ) threads and a warm shared RT cache. Separate-route wall and CPU times sum
    the three searches; peak RSS is their maximum process RSS. The joint route
    performs one search. Upstream processing is excluded.],
  tbl("tbl.koth-matrix-resources"),
) <tbl:matrixresources>
