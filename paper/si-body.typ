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

The run protocol and the entrapment procedure behind every depth number are
given in @sec:methods-benchmark. Three additional ms1searchpy searches per
Astral file used its `-pl` option with an external MS/MS calibration list; the
10 ng runs of that arm ran concurrently, so only their depth is reported. The
Astral uno values come from the same pinned build as every other uno number here
(@sec:methods-benchmark). The outputs of the superseded builds are retained in
the archive, as are the searches of the two 1 ng files removed from the reported
set. The built-in feature-detection arm (@sec:results-raw) searched the same
mzML files with `uno --raw`, the bundled koth profile, and otherwise identical
options; its wall time includes feature detection, and its per-file counts, list
agreement, and resources are in `benchmark/raw_mode/summary.tsv`.

@tbl:versions lists the software environment and benchmark host, and
@tbl:datasets maps the compact dataset names used in the text to the deposited
raw files and local feature tables.

#figure(
  caption: [Software and hardware used for the reported measurements.],
  tbl("tbl.versions"),
) <tbl:versions>

#figure(
  caption: [Dataset names, deposits, and local feature files. The primary and
    Astral rows support the claims in the main text. run10 and deephela are
    tuning records and are not counted as independent validation.],
  tbl("tbl.datasets"),
) <tbl:datasets>

The primary benchmark used peptides of length 7 to 30, no missed cleavages, and
fixed carbamidomethyl cysteine. The patched ms1searchpy #lit("2.8")-series state
and every uno search option are pinned in the benchmark archive. The only
nondefault uno digest setting aligns its peptide bounds with the comparator:

#raw(read("figdata/uno.toml"), lang: "toml", block: true)

The remaining uno defaults are generated from the configuration structs:

#raw(read("figdata/defaults.toml"), lang: "toml", block: true)

*Retention-time calibration controls on the Astral files.*

Stock ms1searchpy reported a second-stage DeepLC sigma of #s(
  "astral.ms1_stock_rt_sigma_min",
) to #s("astral.ms1_stock_rt_sigma_max") min on the #s("astral.files") Astral
files, with no dependence on input. Its calibration set is every peptide of its
first-pass proteins; on these files most of that set was random matches, and
both DeepLC's per-bin median calibration and the Gaussian sigma fit on the
mixture widened accordingly. On the first Astral file, one search each of the
shipped options `-es 1`, `-ptol 3`, `-ci 2 -csc 2`, `-es 1 -ptol 3`, and all
three together left the sigma above #lit("1.4") min and the depth between #lit(
  "41",
) and #lit("63") proteins at 1% entrapment FDR, inside the stock range. The mass
gate was not the cause: against the MS/MS-verified peptides, the true feature
mass error had a robust sigma of #lit(
  "0.66",
) ppm, and ms1searchpy's three-sigma gate kept #lit("72")% of true matches
against #lit("80")% for uno's gate.

The external list for `-pl` was built from the deposit's DIA-NN report as
described in @sec:methods-benchmark. On the same MS/MS-verified peptides,
calibrated on a random half and checked on the other half, DeepLC reached a
robust sigma of #lit("0.69") min and uno's shipped RT model #lit("0.44") min.
The remainder of the stock sigma is therefore attributable to the contaminated
calibration set, not to the RT model. The controls, the list-building script,
and the one-line `-pl` fix are in the benchmark archive under
`paper/analysis/astral/`.

= Complete primary benchmark results <sec:si-results>

@tbl:variability summarizes all primary random draws. @tbl:perdataset gives the
serialized resource cohort: the default-seed uno search, re-measured on the
pinned build, and the three initial ms1searchpy searches per file. The
additional comparator reruns contribute to the depth distribution but not to
this resource table. Individual draw depths are archived in
`benchmark/results/spread.tsv`.

#figure(
  caption: [Primary identification depth across random draws. Each depth cell
    gives the median followed by the full observed range in parentheses. Draws
    are uno seeds or independent ms1searchpy reruns on the same feature file;
    they are not separate acquisitions.],
  tbl("tbl.variability-summary"),
) <tbl:variability>

#figure(
  caption: [Original serialized searches: default-seed uno and the first three
    ms1searchpy runs per file. Identification depth, CPU time, wall time, and
    peak resident memory. Memory is maximum resident set size in MiB reported by
    `/usr/bin/time -v`.],
  tbl("tbl.results-per-dataset"),
) <tbl:perdataset>

= Entrapment calibration <sec:si-calibration>

For each engine, @tbl:entrapreported gives the human and yeast proteins on its
internally reported 1% list, the scale-corrected empirical FDP, and the depth at
1% entrapment $q$. The correction factor is #s("bench.human_db") divided by #s(
  "bench.yeast_db",
), or #s("cal.db_scale"). Because the accepted lists contain few yeast proteins,
one additional yeast hit can materially move a single-dataset estimate. The
across-dataset pattern is therefore more stable than any one value.

The per-dataset curves in @fig:fdpcurves show where each reported list departs
from the external probe across the whole threshold range. Circles mark the two
operating points the main text quotes: the entrapment FDP each engine's own 1%
and 5% list actually carries. A curve above the identity line is optimistic, and
thin curves show individual seeds or reruns in both engines. Bold curves are
pointwise medians where every search has a nonempty list. @tbl:entrapreported
gives the corresponding counts at the 1% point.

#figure(
  placement: auto,
  fig("fig.figs1-fdp-curves", width: 100%),
  caption: [Each engine's internally reported $q$ against the entrapment FDP of
    that list, one panel per primary dataset. Vertical lines mark reported $q$
    of 1% and 5%, and open circles mark every search's FDP there. Thin curves
    show individual uno seeds (blue) or ms1searchpy reruns (orange), and bold
    curves show pointwise medians. Each panel labels uno/ms1searchpy sample
    sizes. Shading marks the region where the reported FDR is optimistic.
    Neither the curves nor their ranges are confidence intervals.],
) <fig:fdpcurves>

#figure(
  caption: [Each engine's internally reported 1% protein list scored against the
    common yeast entrapment probe.],
  tbl("tbl.entrap-reported"),
) <tbl:entrapreported>

= Seeding and repeatability checks <sec:si-reproducibility>

Four stages of the search consume randomness: the decoy shuffle, the RANSAC
retention-time calibration's pair sampling, the cross-validation fold
assignment, and the classifier's row and column subsampling. All four derive
their stream from a single command-line seed. The default value reproduces the
streams of the archived default runs; any other value perturbs all four stages
together, which is how the seed-to-seed spread in the main text was measured.

The primary comparison uses #s("bench.uno_seed_count") uno seeds per file and
#s("bench.ms1_draw_min") to #s("bench.ms1_draw_max") comparator searches. The
smaller Astral and raw-mode extension uses seeds #lit("0"), #lit("1"), and #lit(
  "2",
) in every arm. All #s("rerun.searches") uno searches, in every arm and at every
seed, were executed on the single pinned build of @tbl:versions rather than
carried over from an earlier one. Seed-zero searches ran serialized, because
their wall time and peak memory are what the resource panels report; the
remaining seeds ran two at a time, and their timings are recorded but not
reported. Each search stores its command, its input hashes, its output hash and
the build's hash beside its protein table, under
`benchmark/results/publication_rerun/`; the per-file, per-feature-source,
per-seed summary is `benchmark/results/variability_extension/summary.tsv`.
Matching seeds between feature sources provides a paired comparison of each
file.

#figure(
  placement: auto,
  fig("fig.figs6-determinism", width: 100%),
  caption: [Primary search variability. *(a,b)* Depth at 1% and 5% entrapment
    FDR across all available draws, with open circles, median ticks, and
    minimum-to-maximum connectors. *(c)* Percentage of the other sampled uno
    seeds whose 1% depth was strictly below the archived default seed's depth.
    The dashed line marks 50%; these finite-sample ranks are not population
    percentiles or evidence of deliberate seed selection.],
) <fig:seedspread>

Fixing the seed is not by itself sufficient for a byte-identical result.
Development-time differences between repeated uno searches traced to
process-randomized hash-map iteration reaching floating-point reductions or
first-seen tie-breaks. The affected RT-calibration, mass-recalibration,
feature-competition, and protein-group paths now sort on value or accession
before reduction. Unit tests cover the individual reductions, and
`benchmark/determinism_check.sh` runs two complete searches and byte-compares
their protein, peptide, competition, and group outputs. This full-file check is
part of the release gate for the benchmarked version.

The comparator is seedable on the same terms. ms1searchpy's decoy database is
byte-identical across processes once its module-level generator is seeded, and
its hyperparameter search and classifier accept seeds as well. One further step
is needed and is easy to overlook: several stages apply a seeded shuffle to a
list built from a Python set of peptide strings, whose order randomized string
hashing has already changed, so the interpreter's hash seed must be fixed before
the process starts. We report the comparator as shipped, unseeded, because that
is how it is used in practice, and we make the same measurement for uno by
varying its seed rather than by holding it fixed.

= Retention-time model <sec:si-rt-model>

@fig:rtmodel a gives the network. A peptide is a sequence of residue tokens, and
a learned embedding maps each token to a vector. A convolution stem lifts that
to the working channel width. Three residual blocks with dilations 1, 2 and 4
then widen the context each position sees without adding layers. Padding
positions are masked after every convolution, so a position never reads beyond
the end of the peptide and inference runs over the real length rather than a
fixed maximum. Masked mean and maximum pooling over positions give a
fixed-length vector for a peptide of any length. The eleven global descriptors
are concatenated at that point, and a two-layer head returns one value. The
whole model is 41,665 parameters.

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
carbamidomethylated, so peptides carrying any other modification were dropped;
@fig:rtmodel b shows what that filter costs per source. Two sources are removed
almost entirely, one phosphorylation library and one isobaric-label library, and
the tryptic and HLA peptide libraries survive largely intact. Repeated
observations of one peptide are averaged, and a hash of its sequence assigns it
to the training, validation, or holdout partition. No peptide appears in more
than one partition. Training minimized a smooth L1 loss with the Adam optimizer
and a plateau-triggered learning-rate reduction, selecting the epoch with the
best validation error.

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
    carrying another modification. *(c)* Holdout error by peptide length; bar
    labels are the number of holdout peptides in each bin. *(d)* Holdout error
    by source dataset for sources with at least 200 holdout peptides; the dashed
    line is the overall holdout error. Errors are in hydrophobic-index units,
    the corpus's harmonized scale.],
) <fig:rtmodel>

#figure(
  caption: [The eleven global physicochemical descriptors, computed for the
    holdout peptide #lit("DFEAVANCINAK") and shown as the model receives them.],
  tbl("tbl.rt-features"),
) <tbl:rtfeatures>

= Retention-time predictor benchmark <sec:si-rt-predictor>

The predictor benchmark uses the peptide-disjoint Chronologer holdout,
unmodified ProteomeTools/Prosit holdout peptides, and DeepLC's public #lit(
  "unmod.csv",
) example set. Every peptide contains only canonical residues, is 7--30 residues
long, carries no variable modification, and is absent from the Chronologer train
and validation partitions. Fixed carbamidomethyl cysteine is supplied to DeepLC
and is implicit in uno. Repeated sequence observations are represented by their
median experimental RT. The exact peptide membership of the packaged DeepLC
training corpus is not available here, so the exclusion is established for uno
but cannot be established for DeepLC.

A stable sequence hash assigns 20% of each set to calibration and the remainder
to testing. Every error in @fig:rtpredictor is divided by the 5th--95th
percentile range of the corresponding experimental labels, allowing comparison
across the three RT scales. Both models receive the same calibration: a
least-squares affine map from the model's native output scale to the observed
scale, fitted on the 20% partition. DeepLC's own piecewise calibration, which
ms1searchpy applies in-search, is not used here because with about 200
calibration peptides it extrapolates poorly at the ends of the elution range and
would measure the calibrator rather than the model. Prediction throughput is
measured separately (@fig:rtthroughput). The complete accuracy output also
records per-peptide test predictions, calibration time, raw MAE, 95th-percentile
absolute error, $R^2$, Spearman correlation, package versions, model hash, and
host environment.

#figure(
  placement: auto,
  fig("fig.figs3-rt-predictor", width: 100%),
  caption: [Retention-time predictor accuracy on the held-out 80% of each test
    set, as a percentage of that set's experimental RT 5th--95th percentile
    range. *(a--c)* Residual (predicted minus observed) against observed RT: the
    line is the median and the band the 10th--90th percentile range in twelve
    equal-count bins. *(d--f)* Cumulative distribution of absolute error; the
    open circle marks the 95th percentile, labelled with its value.],
) <fig:rtpredictor>

Throughput depends on the peptides, not on which experimental set they came
from, so it is measured on one synthetic library: an in-silico tryptic digest of
the human entries in the benchmark FASTA, canonical residues, 7--30 residues,
unique sequences, deterministically ordered by sequence hash. Both models
predict the first $N$ peptides of that library for $N$ from 100 to 100,000, five
repeats each after a warm-up call at the largest size; only inference is timed,
with no calibration, on the same workstation with each engine's default
threading. The wall time per batch in @fig:rtthroughput is flat at small
batches, where each engine's fixed per-call overhead dominates, and rises in
proportion to batch size once that overhead is amortized. The steady-state
throughput is the value at the largest batch.

#figure(
  placement: auto,
  fig("fig.figs4-rt-throughput", width: 100%),
  caption: [Retention-time prediction throughput on a synthetic tryptic peptide
    library. *(a)* Wall time per batch against batch size. *(b)* Peptides
    predicted per second, with the steady-state value at 100,000 peptides
    labelled. Lines are medians of five repeats; bands are the min--max range.],
) <fig:rtthroughput>
