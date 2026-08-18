// ============================ Supporting Information =======================
// Appended to the same document, as a single combined PDF. S-prefixed
// numbering for both figures and tables, following the mechanism used in
// d_noise/paper/supplementary.typ.
#show figure.where(kind: image): set figure(numbering: n => "S" + str(n))
#show figure.where(kind: table): set figure(numbering: n => "S" + str(n))
#counter(figure.where(kind: image)).update(0)
#counter(figure.where(kind: table)).update(0)
#set heading(numbering: none)

#align(center)[
  #text(15pt, weight: "bold")[Supporting Information]
  #v(2pt)
  #text(13pt)[uno: A Fast, Reproducible, Dependency-Light Rust Engine for
    MS1-Only Protein Identification]
  #v(2pt)
  #text(10pt, style: "italic")[Patrick Garrett, John R. Yates III]
]

#v(1em)

This section gives the commands, versions and configuration behind the
identification, FAIMS and quantitation benchmarks of Sections 3.1 to 3.5. It
then gives the full per-dataset and per-replicate results, extended entrapment
and quantitation detail, the FAIMS and resource-scaling comparisons, the
quality-column ablation and the run-to-run reproducibility detail. The
retention-time and rescorer experiments of Section S9 were made on an earlier
development configuration that predates the archived pipeline and cannot be
re-run from it. Figures and tables are numbered S1, S2, and so on.

= S1: Benchmark protocol, software versions and commands

Both engines received the identical single-run feature file, the identical FASTA
and an identical digest specification on every dataset. uno was run once per
dataset because its output is byte-identical between runs, ms1searchpy three
times because its output is not. Runs were serialized on one workstation so that
concurrent load could not corrupt the timings, and resource use was captured
with `/usr/bin/time -v`. Depth was scored for both engines by the same
entrapment script, described in the main text (Methods, Benchmark data and
comparison). Software versions and the benchmark hardware are given in
@tbl:versions, and @tbl:datasets maps each run's name in the main text onto its
deposit and its feature file on disk.

#figure(
  caption: [Software and hardware used for every measurement reported here.],
  include "si/versions.typ",
) <tbl:versions>

#figure(
  caption: [Dataset names, deposits and local feature files. The six benchmark
    runs are above the rule. Below it are the three runs used elsewhere in this
    work and not part of the head-to-head: run10, which Section 3.2 reports and
    which is a tuning set rather than a benchmark; QEHF1_09744_JB, used for the
    development-configuration experiments of Sections 3.4 and S9, which is a
    different run from 09746 despite the similar number; and deephela, a
    240-minute gradient used only where stated.],
  include "si/datasets.typ",
) <tbl:datasets>

The commands below regenerate the identification, FAIMS and quantitation
results.

#include "si/commands.typ"

The only uno configuration used in the identification benchmark is the digest
alignment to ms1searchpy's defaults.

#raw(read("si/uno.toml"), lang: "toml", block: true)

Every other uno setting took its default value. The defaults, read from the
configuration structs rather than transcribed, are:

#raw(read("si/defaults.toml"), lang: "toml", block: true)

= S2: Per-dataset and per-replicate results

@tbl:perdataset gives every exact value behind @fig:headtohead, including the
three ms1searchpy replicates individually and the resource measurements
summarized in the main text.

#figure(
  caption: [Per-dataset and per-replicate identification depth, CPU time,
    wall-clock time and peak resident memory. Depth is the number of human
    proteins accepted at the stated two-proteome entrapment FDR, not a count
    from either engine's own reported $q$. Memory is peak resident set size in
    MiB, as reported by `/usr/bin/time -v`.],
  include "si/results_per_dataset.typ",
) <tbl:perdataset>

= S3: Entrapment calibration

@tbl:entrapreported holds the per-dataset counts behind the calibration
statement in Section 3.3. For each engine it gives the human and yeast proteins
on its own reported 1% list, the resulting scale-corrected empirical
false-discovery proportion, and the entrapment depth at 1% used everywhere else
in this paper. The last two columns differ because the entrapment depth is the
deepest cutoff whose entrapment $q$ meets 1%, which need not coincide with the
engine's own 1% list. On 09734 the gap is 830 against 1035 proteins.

These per-dataset proportions are quantized and should be read as such. The
correction factor is 20,647 divided by 6,066, or 3.40. One accepted yeast
protein is therefore worth roughly 0.34 percentage points of estimated FDP at
1000 accepted human proteins, and the accepted sets in @tbl:entrapreported
contain between 0 and 8 yeast proteins. A single yeast hit therefore moves an
individual panel of @fig:fdpcurves by about as much as the entire difference
between the two engines' medians. The pooled median over datasets is the
trustworthy statistic, not any one panel.

#figure(
  placement: auto,
  image("figures/figS1_fdp_curves.png", width: 100%),
  caption: [
    Each engine's own reported $q$ against the entrapment false-discovery
    proportion of the list it reports at that $q$, one panel per dataset, for
    uno and for each ms1searchpy replicate. This is the per-dataset version of
    @fig:entrapment(a), and the vertical line marks a reported $q$ of 1%, not
    the entrapment operating point. Individual panels are quantized by the small
    number of accepted yeast proteins behind each curve.
  ],
) <fig:fdpcurves>

#figure(
  caption: [Each engine's own reported 1% protein list, scored against the
    entrapment probe.],
  include "si/entrap_reported.typ",
) <tbl:entrapreported>

= S4: FAIMS benchmark

The FAIMS comparison of Section 3.1 uses the three replicate 200 ng HeLa
Lumos-with-FAIMS runs of PXD023977 and the deposit's own human plus #emph[E.
  coli] database, with the premade decoy entries stripped so that each engine
generates its own. #emph[E. coli] is the absent entrapment species there, so
this comparison is independent of the yeast probe used elsewhere.

#figure(
  placement: auto,
  image("figures/figS2_faims.png", width: 100%),
  caption: [
    FAIMS benchmark on PXD023977, three replicate runs, #emph[E. coli]
    entrapment. uno against ms1searchpy at 1% and 5% FDR. uno is deeper on
    average at both thresholds and on four of the six replicate-by-threshold
    comparisons, and ms1searchpy's replicate spread exceeds the difference
    between the engines.
  ],
) <fig:faims>

#figure(
  caption: [FAIMS benchmark, mean over the three replicate runs, with the
    per-replicate 1% counts in run order LUM2_05433 / 05434 / 05435.],
  include "si/faims.typ",
) <tbl:faims>

= S5: Resource scaling

#figure(
  placement: auto,
  image("figures/figS3_resources.png", width: 100%),
  caption: [
    Wall-clock time, user CPU time and peak resident memory against input
    feature count. uno is plotted as its single deterministic value; ms1searchpy
    as the median of its three replicates with a bar spanning their minimum and
    maximum. Seven points are shown per panel, not six: run10 is included, which
    is the tuning set rather than one of the six benchmark datasets. Both
    engines grow with input size over the range benchmarked. The separation
    between them is not constant: the wall-clock ratio spans 5.3 to 7.8 across
    these points, and uno's peak memory is not monotone in feature count.
  ],
) <fig:resources>

= S6: The per-feature quality columns are not a depth lever

The feature finder emits per-feature quality columns (`cosine_score`,
`isotope_score`, `combined_score`, `ppm_error`, `neutron_offset`) that a
machine-learning rescoring stage can read. Both engines do read them. Neither
depends on them, because they carry almost no target/decoy information: the
columns describe the feature, and both engines' decoys match real, well-formed
features by mass coincidence.

On run10, separating confident-target matches from decoy matches, the columns
give an area under the ROC curve of 0.500 to 0.526 in uno and a comparable range
in ms1searchpy (@tbl:qualityauc). For uno the value stays at 0.50 even in the
full unreduced 41-million-match pool, so candidate reduction is not stripping a
signal. There is none at any pool size. Stripping all of the columns from the
input changes uno's depth on run10 from 1309 to 1306 at 1% FDR, and from 1887 to
1772 at 5%. The ms1searchpy median falls from 1155 to 1062 at 1%, a drop within
its own replicate spread (@tbl:qualityablation, @fig:qualitycols). // benchmark/RESULTS.md
What does separate target from decoy is retention time, with a distance from
chance of 0.136 for uno against 0.080 for ms1searchpy, and the full match
classifier is further from chance again for both engines (@tbl:qualityauc).

#figure(
  placement: auto,
  image("figures/figS4_quality_columns.png", width: 100%),
  caption: [
    The per-feature quality columns as a depth lever, on the run10 dataset.
    *(a)* Identification depth with and without the columns, for both engines.
    *(b)* Per-column separation of confident-target from decoy matches, plotted
    as the absolute deviation of the ROC area from chance, so that *zero* is no
    separation and larger is more. Two rows are not quality columns: retention
    time, and the full match classifier, both included for scale.
  ],
) <fig:qualitycols>

#figure(
  caption: [Per-feature target/decoy separation, run10, as recorded from a
    one-off diagnostic run whose scripts were not retained. The last two rows
    give the distance from chance for retention time and for the full match
    classifier, in the form the source records them, rather than a raw AUC.],
  include "si/quality_auc.typ",
) <tbl:qualityauc>

#figure(
  caption: [Depth with and without the per-feature quality columns, run10, human
    plus yeast entrapment. Every row except the column-stripped uno row was
    recomputed from protein tables still on disk. That row survives only in the
    benchmark record, because the run behind it was not retained.],
  include "si/quality_ablation.typ",
) <tbl:qualityablation>

= S7: Quantitation detail

uno quantifies from koth's averagine-projected intensities (Experimental
Section, Quantitation), which project each cell onto its mass-derived isotope
pattern and reject the co-eluting off-pattern interference a raw box-sum would
absorb. It is the more robust of the two estimators for this task.

The differential contrast behind Section 3.5 is the extreme 50 against 0.1 fmol
pair. DirectMS1Quant @ivanov2022directms1quant reports 42 changed UPS proteins
there, read from its Figure 2a and accepted at a 5% quantitation FDR against a
database without a human entrapment probe. Our own test is stricter on both
counts, and the comparison should be read with that in mind. A
peptide-fold-change test on uno's matrix recovers 22 of the UPS1 proteins, about
half as many. It also accepts 23 #emph[E. coli] and 60 human proteins, and since
human is the known-absent probe those are false positives outright, so precision
on this contrast is 21%, against a reference output that contains no false
positives at all (@fig:quantdetail c). Both numbers are sensitive to how missing
cells at the dilute end are imputed, so the contrast is indicative rather than a
performance figure. The measurement count is not what limits the test. uno
carries slightly more peptide measurements per protein than the reference does,
and the shortfall is in the quantitative values feeding the test, whose origin
is detection rather than the estimator.

Tracing UPS1 precursors that are confident at 50 fmol down to 0.1 fmol gives two
divergent answers whose distance is itself the result. Of the precursors whose
aligned-matrix cell is filled at 50 fmol, 82% still have a filled cell at 0.1
fmol, because the aligner transfers a match into runs where the peptide was
never independently found. Of the precursors independently detected at 50 fmol,
scored at koth's own per-cell $q$ of 0.01 or below, essentially none survive to
the dilute end (@fig:quantdetail b). The two percentages have different
denominators and are not two views of one population. Transfer of this kind can
be brought under explicit FDR control, as IonQuant does for match-between-runs
@yu2021ionquant, which koth's aligner does not currently do. A fold change built
on transferred fills rather than genuine detections is fragile, which is why the
dilution trend, robust to over-fill, is the primary quantitative result of
Section 3.5.
// benchmark/quant/results/figS5_stdout.txt, panel b

For a substantial minority of the precursors lost, the signal is physically
present but sub-threshold. Dilute peptides elute in only 1 to 2 scans, and
koth's hill-builder (`min_scans = 2`) discards their traces. Lowering the
threshold to `min_scans = 1` does not help. It multiplies the feature count by
6.4 with noise, and an end-to-end test on the six extreme-contrast runs confirms
a net loss. At the looser minimum-runs floor used for that test, the
random-match background rate $p$ rises from 0.15 to 0.51, so half of all
features match something by chance and uno accepts no proteins at all at 1% FDR.
At the stricter floor used for the reported results the search still returns
proteins, but precision falls and the entrapment FDP rises, so the conclusion is
the same in both arms. Untargeted detection cannot separate a one-scan peptide
from one-scan noise.
// benchmark/quant/results/minscans_ab/README.md

#figure(
  placement: auto,
  image("figures/figS5_quant_detail.png", width: 100%),
  caption: [
    Quantitation detail on PXD033350. *(a)* per-species abundance trajectories
    across the eight spike levels. The trajectories separate the spike from the
    constant and the absent species. *(b)* detection-loss traces that follow
    spike-in precursors confident at 50 fmol down to 0.1 fmol. The two traces
    are the aligned matrix cell staying filled and the precursor being found
    independently at a per-cell $q$ of 0.01 or below, and the gap between them
    is match-between-runs over-fill. *(c)* the volcano of the 50 against 0.1
    fmol contrast, limited by detection at the dilute end and reported as
    indicative rather than as a precision figure.
  ],
) <fig:quantdetail>

#figure(
  caption: [Dilution tracking by species on PXD033350. $rho$ is the Spearman
    correlation of protein abundance against log spike amount, and n is the
    number of proteins quantified in at least 8 of the 24 runs.],
  include "si/quant_detail.typ",
) <tbl:quantdetail>

= S8: Run-to-run reproducibility detail

Every source of run-to-run nondeterminism found in uno's development traced to
Rust's process-randomized hash-map iteration order (Section 3.2). The affected
paths were the RANSAC retention-time calibration seed, the retention-time width
estimate, the mass-anchor recalibration sample, the peptide list published to
the quantitation rollup, the feature-competition and profile-reassignment
outputs, and the protein-group table. Each was made order-independent by sorting
on value before reducing, or on accession where a tie genuinely had to be
broken. @fig:determinism plots the resulting contrast dataset by dataset.

#figure(
  placement: auto,
  image("figures/figS6_determinism.png", width: 100%),
  caption: [
    Run-to-run reproducibility on fixed input. *(a)* and *(b)* plot the three
    ms1searchpy replicates against uno's single deterministic value at 1% and at
    5% entrapment FDR. *(c)* is a close-up of run10, the tuning set. Repeating
    the identical search moves ms1searchpy's answer by 30 to 600 proteins at 1%
    FDR depending on the dataset, while uno's output is byte-identical between
    runs. The scale-free form of the spread, the replicate range as a percentage
    of each dataset's own median, is panel (c) of Figure 2 in the main text.
  ],
) <fig:determinism>

= S9: Development-configuration experiments: retention-time model and protein rescorer

The experiments in this section were made on the earlier development
configuration (run QEHF1_09744_JB of PXD015669, koth-called features, the
distilled MLP retention-time model; Methods, Benchmark data). Their values are
not comparable to the six-dataset benchmark.

The distilled student model encodes a peptide by a 99-dimensional order-aware
feature set: amino-acid composition (20), N- and C-terminal residue identities
(first three and last three, 6), position-binned composition (three sequence
bins × 20 = 60), and 13 scalar descriptors, namely length, ln mass,
Kyte–Doolittle @kyte1982hydropathy, Krokhin-SSRCalc @krokhin2006ssrcalc and
Eisenberg @eisenberg1984hydrophobicity hydrophobicity sums, the mean, maximum
and minimum Kyte–Doolittle value, an N-terminus-weighted SSRCalc term, proline,
aromatic and basic counts, and a hydropathy moment. A small multilayer
perceptron over these features, trained on DeepLC predictions over a large,
diverse tryptic digest, is the distilled student.

A naive composition-feature model trailed at $R^2 = 0.972$ against DeepLC
predictions. Its encoding followed the composition-only approach of Klammer et
al. @klammer2007rtprediction: amino-acid counts, counts within fixed N- and
C-terminal windows, peptide length and mass. It gave 141 proteins against the
DeepLC library's 176, and fewer of them overlapped the reference search
(@tbl:rtdepth). Replacing the encoding with the 99-d order-aware feature set
raised the distilled model to $R^2 = 0.980$, matching the LightGBM ceiling, and
restored identification depth to DeepLC parity (@tbl:rtdepth), in-process and in
pure Rust. The feature encoding, not the boosting backend, was the lever. A
still-richer 215-d set gave only marginal improvement ($R^2 = 0.9824$),
indicating that the 99-d set is near the practical ceiling for hand-crafted
features.

#figure(
  caption: [RT model against identification depth at 1% FDR, measured on the
    earlier development configuration (run QEHF1_09744_JB of PXD015669,
    koth-called features). Reference overlap is the number of accepted proteins
    also found by the DeepLC-driven reference search.],
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    table.header(
      [RT model],
      [$R^2$ vs DeepLC],
      [Proteins / ref-overlap],
      [RT-stage cost],
      [Dependency],
    ),
    [Composition-only], [0.972], [141 / 119], [in-process], [none (Rust)],
    [Order-aware (99-d, this work)],
    [0.980],
    [\~180 / \~152],
    [in-process],
    [none (Rust)],

    [DeepLC library], [(teacher)], [176 / 152], [\~130 s preproc.], [Python],
  ),
) <tbl:rtdepth>

A reverse-sequence benchmark quantifies how much sequence-order information a
model captures. A peptide and its reverse have identical amino-acid composition
but different true retention time. That sequence order carries retention
information beyond composition is long established @petritis2006sequence, and a
retention-coefficient model is by construction unable to represent it, since it
predicts one value for any permutation of a given composition
@moruz2017rtreview. An order-aware predictor must assign a peptide and its
reverse different RTs, while a composition-only model is forced to predict the
same value for both. This forward minus reverse RT difference is the signal that
lets a search engine distinguish a real peptide from a reversed decoy, and the
benchmark reports how much of it a given model actually captures.

The benchmark holds 46,304 held-out tryptic peptides and their reverses, which
are non-tryptic and out of distribution. DeepLC RT was predicted for both, and
each model was measured on how well it reproduced the forward minus reverse RT
shift. A composition-only model captured 0% of this signal, by construction.
uno's order-aware model captured 72% (Pearson $r = 0.68$) while predicting
reverse-sequence RT nearly as well as forward ($R^2_"rev" = 0.96$,
$R^2_"fwd" = 0.98$, @tbl:reverse). The residual 28% is the genuine gap between
hand-crafted features and DeepLC's sequence convolution.

#figure(
  caption: [Reverse-sequence RT benchmark, 46,304 held-out tryptic peptides and
    their reverses, on the earlier development configuration.],
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    table.header(
      [Model],
      [$R^2_"fwd"$],
      [$R^2_"rev"$],
      [order-$Delta$ corr.],
      [order magnitude],
    ),
    [Composition-only], [0.96], [0.95], [0.00], [0%],
    [uno order-aware (this work)], [0.98], [0.96], [0.68], [72%],
  ),
) <tbl:reverse>

The hardened variant of the LDA protein rescorer of Section 3.4 applied
high-confidence iterated positives ($q <= 0.01$, refreshed), target/decoy pairs
co-folded, and stronger covariance shrinkage. It removed the inflated apparent
count on the human-only search (242 to 223) but did not restore FDR control,
which stayed at 11.4% empirical under the entrapment probe (@tbl:entrap). These
safeguards follow established semi-supervised rescoring practice
@kall2007percolator @the2016percolator.

#figure(
  caption: [Two-proteome entrapment of three protein scorers, measured on the
    earlier development configuration. Empirical FDR is the
    database-size-corrected entrapment estimate at a reported 1%. The binomial
    is the shipped default; both LDA variants are disabled.],
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    table.header(
      [Protein scorer], [Real accepted], [Yeast (false)], [Empirical FDR]
    ),
    [Binomial (default)], [201], [1], [1.7%],
    [LDA rescorer (naive)], [208], [7], [11.5%],
    [LDA rescorer (hardened CV)], [209], [7], [11.4%],
  ),
) <tbl:entrap>
