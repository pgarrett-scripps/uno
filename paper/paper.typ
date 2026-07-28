#import "@preview/arkheion:0.1.0": arkheion, arkheion-appendices

#show: arkheion.with(
  title: "uno: A Fast, Reproducible, Dependency-Light Rust Engine for MS1-Only Protein Identification",
  authors: (
    (
      name: "Patrick Garrett",
      email: "pgarrett@scripps.edu",
      affiliation: "Department of Chemical Physiology, The Scripps Research Institute, 10550 North Torrey Pines Road, La Jolla, CA 92037, USA",
      orcid: "0000-0002-8434-9693",
    ),
    (
      name: "John R. Yates III",
      email: "jyates@scripps.edu",
      affiliation: "Department of Chemical Physiology, The Scripps Research Institute, 10550 North Torrey Pines Road, La Jolla, CA 92037, USA",
      orcid: "0000-0001-5267-1672",
    ),
  ),
  abstract: [
    MS1-only protein identification (the DirectMS1 approach) infers peptide and
    protein identities from precursor mass and retention time (RT) alone, without
    MS2 fragmentation, enabling very short LC gradients and high sample
    throughput. The reference implementation, ms1searchpy, is written in Python and
    depends on the deep-learning RT predictor DeepLC. We present *uno*, an
    open-source MS1-only search engine written in Rust that consumes aligned
    feature matrices from the koth feature finder. uno delivers four practical
    advantages over existing tooling at equal-or-greater identification depth:
    (i) *reproducibility* — on identical input, ms1searchpy returned 152–303
    proteins across repeated runs ($sigma approx 55$) whereas uno is deterministic
    ($sigma approx 1$), byte-identical by construction; (ii) *a distilled,
    pure-Rust RT model* that reproduces DeepLC's MS1 identification depth
    ($R^2 = 0.980$ vs DeepLC predictions) while running in-process and removing the
    Python/DeepLC dependency entirely; (iii) *on-data RT fine-tuning* — refitting
    only the pretrained RT network's linear readout over its frozen embedding on a
    run's own confident hits raises identification depth (+51% on the benchmark
    here) while remaining FDR-honest under two-proteome entrapment (zero
    entrapment false positives); and (iv) *speed* — algorithmic and parallel
    optimization of the protein-inference grid reduces an end-to-end search from
    43 s to 3.4 s ($approx 13 times$) with byte-identical output. A two-proteome
    entrapment harness confirms the protein-level FDR is approximately honest and,
    using the same harness, shows that a learned (LDA) protein rescorer inflates
    the empirical FDR despite raising apparent counts — underscoring why a
    deterministic binomial scorer is retained, and why the RT fine-tune (which
    passes the same test) is adopted instead.
  ],
  keywords: ("MS1-only proteomics", "DirectMS1", "retention time prediction", "reproducibility", "Rust", "open-source software", "label-free quantification"),
  date: "June 2026",
)

#set raw(lang: "rust")

// Draft note: §3.1–3.6 report measurements from the prior development branch
// (HeLa PXD015669 operating point). §3.7 (RT fine-tune) and §3.8 (performance)
// are measured on the current pipeline against the koth-aligned human+E. coli
// benchmark feature matrix used in this revision; their setup is stated inline.

= Introduction

Bottom-up proteomics conventionally identifies peptides from MS2 fragment spectra.
MS1-only identification (DirectMS1 @ivanov2020directms1) takes a different route:
it matches the accurate precursor mass of an LC–MS feature to a theoretical
peptide and confirms the assignment with predicted retention time, omitting MS2
entirely. Because no time is spent acquiring fragment spectra, MS1-only workflows
support very short gradients (minutes) and correspondingly high throughput, at the
cost of an information-poor identification problem: mass and RT alone are weak
discriminators, and identity must be concentrated at the protein level across
multiple consistent peptides.

The established open-source tool for this approach, *ms1searchpy*
@ivanov2021ms1searchpy, couples a precursor mass search to RT prediction by
*DeepLC* @bouwmeester2021deeplc and a LightGBM-based rescoring stage, with
protein-level inference by a binomial model. It is effective but has two practical
costs that matter at scale. First, its identification counts are *not reproducible
run to run*: the DeepLC and LightGBM stages are stochastic and unseeded, so
repeated searches of the same features yield materially different protein lists
(§3.1). Second, it carries a *heavyweight Python/deep-learning dependency* (DeepLC,
TensorFlow) whose RT-prediction step dominates wall-clock time for a single run.

We present *uno*, an MS1-only identification and label-free quantitation engine
written in Rust. uno is one component of an all-Rust pipeline: the *koth* feature
finder @koth performs alignment, deisotoping and XIC/LFQ; uno consumes koth's
aligned feature matrix and assigns peptide/protein identities, controls FDR, and
rolls intensities up into quant matrices. uno's design goals are determinism,
speed, and minimal external dependencies, while meeting or exceeding the
identification depth of ms1searchpy. The contributions of this work are:

+ A clean, open-source Rust reimplementation of DirectMS1-style MS1 search with
  *deterministic* output (§2).
+ A quantitative demonstration that *MS1-only identification counts are strongly
  implementation-dependent*: ms1searchpy spans 152–303 proteins on one fixed
  feature set, and uno sits at its median with $approx 50 times$ lower variance
  (§3.1).
+ A *distilled, pure-Rust RT model* with engineered order-aware features that
  reproduces DeepLC's MS1 identification depth ($R^2 = 0.980$) without a Python
  dependency (§3.3).
+ An *on-data RT fine-tuning* step that refits only the pretrained network's
  readout over its frozen embedding on a run's confident hits, raising depth while
  passing two-proteome entrapment (§3.7).
+ *Performance engineering* of the protein-inference stage (memoized binomial,
  best-first greedy parsimony with lazy invalidation, parallel threshold grid)
  that reduces a search from 43 s to 3.4 s with byte-identical output (§3.8).

= Experimental Procedures (Methods)

== Input and pipeline position

uno consumes the aligned feature matrix produced by koth: per-feature `massCalib`,
`mz`, `charge`, `rtApex`, ion mobility (`im`, where available), `combined_score`,
and per-sample intensities. RT is in koth's reference-run space, so a single
calibration suffices rather than per-file calibration.

== Peptide database and decoys

In-silico tryptic digestion is performed with sage-core @lazear2023sage (enzyme
rules, missed cleavages, length/mass bounds, fixed/variable modifications), and
masses are recomputed in `f64`. Candidate lookup is a binary search over a
mass-sorted peptide vector. Decoys are generated by sequence reversal with an
optional mass shift; we use mass-shifted decoys because they provide a clean
protein-level FDR null, whereas same-mass reversed decoys are RT-indistinguishable
"twins" of their targets under composition-based RT models (§3.4).

== Feature–peptide matching

Each feature is matched to all peptides within a ppm tolerance, enumerating
isotope-error offsets (0, ±1 Da/charge). The resulting peptide–feature matches
(PFMs) are reduced to the top-$K$ best-mass candidates per feature to bound the
downstream cost without discarding the RT-correct assignment.

== Retention-time model (this work, §3.3)

uno predicts RT with a model *distilled from DeepLC*: DeepLC predictions over a
large, diverse tryptic digest serve as training labels. Sequences are encoded by a
99-dimensional *order-aware feature set*: amino-acid composition (20), N- and
C-terminal residue identities (first three, last three; 6), *position-binned
composition* (three sequence bins × 20 = 60), and 13 scalar descriptors (length,
ln mass, Kyte–Doolittle @kyte1982hydropathy / Krokhin-SSRCalc @krokhin2006ssrcalc /
Eisenberg @eisenberg1984hydrophobicity hydrophobicity sums, an N-terminus-weighted
SSRCalc term, proline/aromatic/basic counts, and a hydropathy moment). The shipped
student is a small multilayer perceptron over these features (a gradient-boosted
@gbdtcrate variant is also supported); it is distilled offline once and predicts
RT for all candidate peptides in-process at search time.

Because the frozen predictor lives in DeepLC's abstract RT space, uno calibrates it
to the run's observed gradient using a RANSAC-seeded, monotone piecewise-linear
bootstrap on confident identifications, progressively tightened across rounds and
vouched by protein-level multiplicity.

== On-data RT fine-tuning (this work, §3.7)

A frozen predictor plus monotone calibration can only *rescale* the prediction
axis; it cannot correct per-peptide ranking errors the network makes on an unseen
gradient. Refitting a fresh linear model on raw composition features (the
DirectMS1 on-data trick) discards the pretrained representation and underperforms
(§3.7). Instead, uno *fine-tunes only the readout*: it freezes the MLP's learned
nonlinear embedding (the penultimate-layer activations) and refits a ridge-
regularized linear map from that embedding to the run's *observed* RT, using the
confident-hit set with iterative outlier rejection. This keeps the pretrained
representation while adapting offset, scale, and shape to the run. The fit is
closed-form and deterministic, so it preserves byte-identical output. The same
confident-hit selection and decoy treatment apply, so target/decoy FDR remains
calibrated (validated by entrapment, §3.7).

== Scoring and FDR

PFMs are scored by a deviation discriminant combining recalibrated mass error (in
$sigma$ units, via sigma-clipped re-estimation of the mass-error distribution) and
calibrated RT deviation. Because mass and RT alone do not separate target from
decoy at the single-peptide level in this regime (§3.2), identity is concentrated
at the protein level: a binomial survival score $-log_10 P(X >= k | n, p)$
evaluates whether a protein's $k$ observed peptides exceed the random-match
background rate $p$ estimated from decoys, followed by greedy parsimony,
target/decoy pair competition, and a target/decoy FDR walk. Mass/RT acceptance
thresholds are optimized on a grid to maximize target proteins at the chosen FDR
(mirroring ms1searchpy's `worker`); because decoys pass the identical filter at
every grid point, the FDR estimate stays calibrated regardless of which cutoff
maximizes yield.

All stages are implemented to be *deterministic*: there is no use of
process-randomized iteration order or unseeded RNGs in any result-affecting path,
so repeated runs on identical input produce byte-identical output (§3.1). An
optional learned (cross-validated LDA) protein rescorer is provided but *disabled
by default*, because entrapment validation shows it inflates the empirical FDR
(§3.6); the binomial model with on-data-fine-tuned RT is the validated default.

== Quantitation

Identified peptides/proteins are mapped onto koth's per-sample aligned intensities
to emit peptide- and protein-level quant matrices; uno does not perform its own
feature extraction (koth owns alignment, XIC integration, and match-between-runs).
Two per-cell intensity estimators are available at the alignment stage: the default
box-sum (all signal in the ($m/z$, RT) extraction window), and an optional
#emph[averagine projection] in which the reported intensity is the observed
isotopologue vector projected onto the mass-derived averagine isotope pattern — a
matched filter that retains on-pattern signal and rejects co-eluting off-pattern
interference. uno emits the protein-group × sample abundance matrix and stops there;
downstream differential testing (peptide fold-change, $t$-test, Benjamini–Hochberg
FDR, following the DirectMS1Quant approach) is left to standard tools, keeping uno an
identification-and-quantitation engine rather than a statistics framework.

== Benchmark data and comparison

Variance and RT-distillation results (§3.1, §3.3, §3.4) use the HeLa dataset
PXD015669 (run `QEHF1_09744_JB`, Q Exactive HF), with features called by koth.
The on-data fine-tune and performance results (§3.7, §3.8) are measured on a
koth-aligned benchmark feature matrix (20 samples, 45,155 aligned features)
searched against a human + #emph[E. coli] database (24,819 proteins) with the
distilled MLP RT model; the entrapment experiments add a #emph[S. cerevisiae]
proteome (6,068 proteins) as a known-absent probe. ms1searchpy and uno were run
on the identical koth feature set where compared.

Quantitation is evaluated on the DirectMS1Quant UPS1-into-#emph[E. coli]
dilution series (PXD033350: eight spike levels, 0.1–50 fmol, × 3 replicates, 24 runs;
Q Exactive HF), searched against the 48 UPS1 spike proteins, the #emph[E. coli]
proteome (constant background), and the human proteome as a known-absent entrapment
probe. Runs were converted, feature-found and cross-run-aligned by koth; uno searched
the aligned matrix, and peptide fold-change differential testing was scored against
the published DirectMS1Quant results and Spectronaut ground truth.

= Results and Discussion

== MS1-only identification counts are strongly implementation-dependent

Running ms1searchpy on a single fixed koth feature set six times yielded
*152, 153, 155, 188, 203, and 303 proteins* at a nominal 1% FDR ($sigma approx
55$; the logged random-match background $p$ ranged 0.016–0.037). The high outlier
(303) is consistent with the stochastic LightGBM rescoring occasionally
overfitting the target/decoy boundary. uno, run on the same features through its
deterministic path, returned *175–179 proteins* ($sigma approx 1$) — sitting at the
*median* of ms1searchpy's distribution with $approx 50 times$ lower variance.

uno's determinism is by design. Two sources of run-to-run nondeterminism in early
development were traced to Rust's process-randomized hash-map iteration order:
(i) the RANSAC RT-calibration seed was assembled from a hash map and sampled
positionally, and (ii) the protein greedy-parsimony and FDR steps broke score ties
by iteration order. Both were made order-independent (sort the calibration seed by
feature index; resolve ties on accession). After these fixes, independent runs on
identical input are *byte-identical*. Determinism is thus a guaranteed property of
the implementation, not a statistical observation.

== Identity is concentrated at the protein level, not the peptide level

We confirmed directly that in this regime mass+RT cannot separate target from decoy
peptides: among matches with near-zero RT deviation ($< 0.01$ min), target and
decoy peptides occur at $approx$ 50:50. Most target peptides in a proteome are not
actually present in a given sample, so a random target match is no more likely than
a decoy match. Discrimination therefore arises only from *multiplicity* — real,
present proteins accumulate several RT-consistent peptides while decoy proteins
remain singletons — which the binomial protein model exploits. This property is
intrinsic to MS1-only search and bounds what any single-peptide discriminant
(including ML rescoring) can achieve.

== A distilled pure-Rust RT model reaches DeepLC's identification depth

The RT predictor is the heaviest dependency in MS1-only search. A naïve
composition-feature model (the prior 69-d Klammer set) trailed: $R^2 = 0.972$ vs
DeepLC predictions, and 141 proteins / 119 reference-overlap on the benchmark (vs
the DeepLC library's 176 / 152). Replacing the encoding with the 99-d *order-aware*
feature set (§2.4) raised the distilled model to *$R^2 = 0.980$, matching the
LightGBM ceiling*, and restored identification depth to *DeepLC parity*
(@tbl:rtdepth), achieved *in-process in pure Rust* with no Python dependency. The
feature encoding, not the boosting backend, was the lever; a still-richer 215-d set
gave only marginal improvement ($R^2 = 0.9824$), indicating the 99-d set is near
the practical ceiling for hand-crafted features.

#figure(
  caption: [RT model vs identification depth on the HeLa benchmark (1% FDR).],
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    table.header([RT model], [$R^2$ vs DeepLC], [Proteins / ref-overlap], [RT-stage cost], [Dependency]),
    [Composition (69-d)], [0.972], [141 / 119], [in-process], [none (Rust)],
    [*Order-aware (99-d, this work)*], [*0.980*], [*\~180 / \~152*], [*in-process*], [*none (Rust)*],
    [DeepLC library], [(teacher)], [176 / 152], [\~130 s preproc.], [Python/TF],
  ),
) <tbl:rtdepth>

== A reverse-sequence benchmark quantifies sequence-order capture

A peptide and its reverse have identical amino-acid composition but different true
retention time; an order-aware predictor must assign them different RTs, while a
composition-only model is mathematically forced to predict the same value for both.
This forward−reverse RT difference is precisely the signal that lets a search engine
distinguish a real peptide from a reversed decoy.

We constructed a benchmark of 46,304 held-out tryptic peptides and their reverses
(non-tryptic, out-of-distribution), predicted DeepLC RT for both, and measured how
well each model reproduced the forward−reverse RT shift. A composition-only model
captured *0%* of this signal (by construction). uno's order-aware model captured
*72%* (Pearson $r = 0.68$) while predicting reverse-sequence RT nearly as well as
forward ($R^2_"rev" = 0.96$, $R^2_"fwd" = 0.98$; @tbl:reverse). The residual ~28%
is the genuine gap between hand-crafted features and DeepLC's sequence convolution.
This benchmark is, to our knowledge, a novel and reusable way to report RT-model
sequence-order sensitivity.

#figure(
  caption: [Reverse-sequence RT benchmark (46,304 held-out peptides + reverses).],
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, center, center, center, center),
    table.header([Model], [$R^2_"fwd"$], [$R^2_"rev"$], [order-$Delta$ corr.], [order magnitude]),
    [Composition-only], [0.96], [0.95], [0.00], [0%],
    [*uno order-aware (this work)*], [*0.98*], [*0.96*], [*0.68*], [*72%*],
  ),
) <tbl:reverse>

== Two-proteome entrapment confirms the protein-level FDR is honest

Target/decoy FDR can be optimistic if the decoys do not faithfully model the
random-coincidence behaviour of real-but-absent peptides. We tested uno's reported
FDR directly with a two-proteome entrapment: the sample was searched against the
real database plus an entrapment proteome — #emph[S. cerevisiae], absent from the
sample — together with decoys for both. Any accepted yeast protein is a known false
positive. We estimate the empirical false-discovery proportion as
$n_"yeast" times (N_"real") / (N_"yeast") \/ n_"real"$, the DB-size correction for
the smaller entrapment proteome. At a reported 1% protein FDR, uno's binomial
scorer is approximately honest (@tbl:entrap).

== On-data RT fine-tuning raises depth while preserving FDR

We evaluated the on-data RT options of §2.5 on the koth-aligned human+#emph[E.
coli] benchmark matrix. The frozen MLP with monotone calibration identified *39
proteins* at 1% FDR. Refitting a fresh composition-feature ridge on the confident
hits — the literal DirectMS1 on-data trick — *degraded* depth to *9 proteins*: a
from-scratch linear model is a far weaker RT predictor (residual $sigma approx
0.56$ min vs the frozen model's $approx 0.2$ min) and discards the pretrained
representation. Fine-tuning only the readout over the frozen MLP embedding (§2.5)
instead raised depth to *59 proteins* (+51% over frozen), the configuration now
shipped as default.

To distinguish real gain from FDR leakage we swept the ridge penalty: depth was
stable across $lambda = 1 -> 10$ (59 → 56) and collapsed only when over-regularized
into under-fitting ($lambda = 100 -> 19$, $lambda = 1000 -> 4$). A genuine overfit
would peak sharply at low $lambda$ and crumble with any regularization; the
observed plateau is the signature of real signal. We then ran the entrapment harness
of §3.5 (real DB + #emph[S. cerevisiae] probe). Both the frozen and fine-tuned
configurations accepted *zero* yeast proteins at 1% FDR — the nearest yeast protein
sat at $q = 0.15$, $15 times$ past threshold — while the fine-tune added real
proteins (34 → 38 on the entrapment DB; @tbl:entrap). The fine-tune is therefore
FDR-honest: it is the mirror image of the LDA rescorer (§3.6), gaining depth
*without* admitting entrapment false positives.

#figure(
  caption: [Two-proteome entrapment. Binomial and LDA rows: HeLa human sample vs
  human + yeast (prior operating point). RT-fine-tune rows: koth human+#emph[E.
  coli] benchmark vs real DB + yeast probe (this revision). "Empirical FDR" is the
  DB-size-corrected entrapment estimate at a reported 1%.],
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    table.header([Scorer / RT mode], [Real accepted], [Yeast (false)], [Empirical FDR]),
    [Binomial (default)], [201], [1], [1.7%],
    [LDA rescorer (naïve)], [208], [7], [11.5%],
    [LDA rescorer (hardened CV)], [209], [7], [11.4%],
    table.cell(colspan: 4, align: center)[#emph[RT fine-tune ablation (this revision)]],
    [Frozen MLP + calibration], [34], [0], [\~0%],
    [*MLP readout fine-tune (default)*], [*38*], [*0*], [*\~0%*],
  ),
) <tbl:entrap>

A limitation of the probe used here is resolution: with zero entrapment hits a
single yeast acceptance would scale (×4.09 DB-size factor) to $approx 10%$
estimated FDR, so "0%" means "below this test's floor," not "provably zero." The
*comparative* claim — that the fine-tune adds no entrapment false positives
relative to the validated-honest frozen path — is robust; tightening the absolute
estimate requires a size-matched probe proteome (§5).

== Learned protein rescoring raises apparent counts but inflates FDR

We separately evaluated a learned protein-level rescorer combining the binomial
score with log peptide count, coverage, and RT/mass consistency, implemented as a
deterministic closed-form Fisher LDA with 3-fold cross-validation by accession hash
(no RNG). It raised apparent identifications (218 → 242 proteins) and was
byte-identical across runs. However, the *same entrapment test exposed the gain as
illusory*: the LDA accepted *7 yeast (false) proteins versus the binomial's 1*, an
estimated empirical FDR of *11.5%* at a nominal 1% (@tbl:entrap).

Hardening the rescorer with textbook semi-supervised safeguards — high-confidence
iterated positives ($q <= 0.01$, refreshed), target/decoy pairs co-folded, and
stronger covariance shrinkage — *removed the inflated apparent count
(242 → 223) but did not restore FDR control* (still 11.4% empirical). The inflation
is therefore *not classical overfitting* but an *objective mismatch*: a learned
scorer is trained to separate targets from decoys, yet the false positives that
matter are *real-but-absent* proteins, whose RT and mass consistency resemble
present proteins and differ from decoys. The very features that could add depth are
exactly those a present and an absent-but-real protein share. The count-based
binomial is robust precisely because it uses only peptide multiplicity, where an
absent protein and a decoy are indistinguishable. We therefore *retain the binomial
scorer as default* and ship the LDA disabled @the2016percolator. This is a
transferable cautionary result: learned protein rescoring in MS1-only search trades
FDR integrity for apparent depth even under careful regularization, and entrapment
validation is essential to detect it. The RT fine-tune (§3.7) is adopted precisely
*because* it passes this same test where the LDA fails.

== Performance engineering

The end-to-end search was profiled on the benchmark matrix (20 cores). The original
implementation ran in *43 s*, of which *\~40 s* was protein inference. The cost was
the threshold-optimization grid (§2.6): \~99 (RT, mass) cells, each running a greedy
parsimony loop that re-scored every surviving protein on every round. Instrumentation
showed up to *17.7 million* exact-binomial evaluations per cell (incomplete-beta
continued fraction), \~1 billion across the grid.

Three changes removed this without altering output:
(i) the independent grid cells were run in *parallel* (rayon), with the winner
selected by the same deterministic rule;
(ii) within each cell the binomial was *memoized* on $(n_"theo", n_"obs")$,
collapsing \~18 M evaluations to a few thousand distinct ones, and the greedy loop
was restructured as a *best-first heap with lazy invalidation* — only the few
proteins that lose a peptide when a winner is chosen are re-scored and re-pushed,
turning the per-round full rescan into work proportional to peptides removed;
(iii) a release profile with link-time optimization and a single codegen unit.

The protein phase fell from *\~40 s to \~1 s* single-threaded (the algorithmic
gain alone, $approx 40 times$), and the full search from *43 s to 3.4 s*
($approx 13 times$). Output is *byte-identical* to the original across the whole
pipeline (identical protein lists and FDR), and remains deterministic across runs.

== Quantitation tracks a spike-in dilution; low-abundance differential sensitivity trails DirectMS1Quant

On the DirectMS1Quant UPS1-into-#emph[E. coli] dilution series, uno identifies *361
proteins* at 1% FDR (honest entrapment FDP $approx 1.8%$) and recovers the designed
structure of the experiment: spike-in (UPS1) protein abundance *tracks the dilution*
(median Spearman $rho = +0.71$ against log spike amount; 78% of UPS1 proteins with
$rho > 0.5$), while the constant #emph[E. coli] background and the absent human probe
stay flat (median $rho = +0.04$ and $+0.05$). Peptide-fold-change differential testing
recovers spike-in proteins with *high specificity* — zero #emph[E. coli] false
positives at the extreme (50 vs 0.1 fmol) contrast.

Differential *sensitivity* at the dilute end, however, trails the reference: uno flags
$approx 12$ changed UPS1 proteins at the extreme contrast versus 42 for DirectMS1Quant.
The gap is not in the differential test — the per-protein measurement count and
statistic match — but in the quantitative values feeding it, and it decomposes into two
separable causes. (i) #emph[Fill specificity]: koth's default box-sum integrates every
co-eluting hill in the extraction window, inflating dilute cells with off-pattern
interference; the averagine-projection estimator (Methods, #emph[Quantitation]) retains only the on-pattern
component and roughly *halves the false-positive differential calls* (e.g. 15 true / 17
false $arrow.r$ 12 true / 8 false under a common imputation; precision $47% arrow.r
60%$). (ii) #emph[Detection sensitivity]: tracing spike-in peptides confident at 50 fmol
down to 0.1 fmol, *63% remain detected*; of those lost, the signal is physically present
but sub-threshold — dilute peptides elute in only 1–2 scans, and koth's hill-builder
(`min_scans = 2`) discards their traces. Lowering the threshold (`min_scans = 1`)
recovers $approx 50%$ of them but multiplies the feature count $approx 6 times$ with
noise, and an end-to-end test confirms this is a net loss (differential precision
$27% arrow.r 9%$). Untargeted detection cannot separate a one-scan peptide from
one-scan noise; the resolved path is a *targeted, identification-seeded re-extraction* —
a matched-filter requant at the ($m/z$, RT) of peptides already identified in the
abundant samples, which confirms known targets rather than discovering new ones — scoped
in the validation roadmap. uno's quantitation is thus validated for identification-linked
abundance and dilution tracking, with trace-level differential-sensitivity parity
remaining open.

= Limitations

We report these candidly, as they delimit where MS1-only search and uno apply.

- *Inherent PSM-FDR floor.* Mass+RT cannot separate target/decoy at the
  single-peptide level (§3.2); identification depends on protein-level multiplicity.
  This is a property of MS1-only data, not of uno.
- *No ion mobility, no peptide-centric depth.* On Orbitrap data without ion
  mobility, peptide-centric extraction from raw hills does not add identifications:
  the (m/z, RT) coordinate is too degenerate without a third dimension. Ion-mobility
  data (timsTOF) is expected to lift this; uno predicts CCS via an analogous
  distilled model where IM is present.
- *Search-breadth expansion is counterproductive.* Adding missed cleavages or
  modifications inflates the theoretical count faster than the observed count and
  raises the decoy background, reducing identifications — consistent with DirectMS1's
  minimal digest.
- *TMT-labeled MS1 data.* uno currently yields $approx 0$ identifications on TMT MS1
  data; DeepLC RT for TMT-modified peptides correlates only $r approx 0.71$ with
  observed RT, too weak to seed calibration.
- *Entrapment resolution.* The single-probe entrapment bounds the empirical FDR only
  coarsely (§3.7); a size-matched probe and a multi-dataset panel are required to
  tighten absolute estimates.
- *Low-abundance differential sensitivity.* On a spike-in dilution series uno tracks
  abundance and calls changed proteins with high specificity, but at the most dilute
  contrasts it recovers fewer changed proteins than DirectMS1Quant. The cause is
  upstream feature detection — trace peptides yield 1–2 scan hills the feature finder
  discards — not the differential test; the scoped remedy is a targeted,
  identification-seeded matched-filter re-extraction rather than a global detection
  threshold, which trades specificity for recall.

= Validation roadmap

+ *Size-matched / multi-dataset entrapment* to tighten the absolute FDR estimate for
  both the binomial default and the RT fine-tune.
+ *Multi-dataset benchmark* ($>= 3$ datasets) spanning instruments/gradients/
  organisms; per dataset: uno vs ms1searchpy protein counts, FDR validity, variance
  ($n >= 10$ runs), and speed, including the fine-tune depth gain.
+ *Formal reproducibility experiment* ($n >= 10$ runs × both engines × each dataset).
+ *Quantitation.* Dilution tracking and differential specificity are validated here;
  the open items are trace-level differential #emph[sensitivity] via an
  identification-seeded matched-filter re-extraction, and CV/accuracy on a
  known-ratio benchmark.
+ Software-release hygiene for open-source expectations: a build that does not
  require a local sage checkout, README, example data, tests/CI, and an archived
  release with DOI.

= Conclusions

uno is a fast, deterministic, open-source Rust engine for MS1-only protein
identification that meets or exceeds the identification depth of ms1searchpy while
removing its run-to-run variability and its heavyweight Python/DeepLC dependency.
Its output is byte-identical across runs by construction; its distilled order-aware
RT model reaches DeepLC's MS1 depth in pure Rust; an on-data fine-tune of the RT
readout adds further depth (+51% here) while passing two-proteome entrapment; and
targeted performance engineering brings an end-to-end search to a few seconds with
unchanged results. The accompanying reverse-sequence benchmark provides a direct,
reusable measure of RT-model sequence-order sensitivity, and our entrapment analysis
shows that learned protein rescoring — while raising apparent counts — inflates FDR,
motivating uno's transparent binomial model and the entrapment-validated RT fine-tune
adopted in its place.

= Data and Software Availability

uno is open-source (MIT) at #link("https://github.com/")[\[URL\]]. Benchmark data:
PXD015669. RT distillation labels, trained models, and analysis scripts:
\[DOI/Zenodo — TODO\].

#bibliography("paper.bib")
