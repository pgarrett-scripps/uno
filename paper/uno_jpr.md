# uno: A Fast, Reproducible, Dependency-Light Rust Engine for MS1-Only Protein Identification

**Authors:** Patrick Garrett, et al.
**Target venue:** Journal of Proteome Research (Technical Note / Article)
**Draft status:** WORKING DRAFT — results in §3.1–3.6 are measured (incl. determinism fix and entrapment-FDR validation); §3.7 (speed) and items tagged `[PENDING]` require the remaining experiments in §5 (multi-dataset, n≥10 reproducibility, speed table) before submission.

---

## Abstract

MS1-only protein identification (the DirectMS1 approach) infers peptide and
protein identities from precursor mass and retention time (RT) alone, without
MS2 fragmentation, enabling very short LC gradients and high sample throughput.
The reference implementation, ms1searchpy, is written in Python and depends on
the deep-learning RT predictor DeepLC. We present **uno**, an open-source MS1-only
search engine written in Rust that consumes aligned feature matrices from the
koth feature finder. uno delivers three practical advantages over the existing
tooling at equal identification depth: (i) **reproducibility** — on identical
input, ms1searchpy returned 152–303 proteins across repeated runs (σ≈55) whereas
uno is deterministic (σ≈1); (ii) **speed** — the core search completes in seconds
of compiled Rust; and (iii) **dependency-light operation** — a distilled,
pure-Rust gradient-boosted RT model reproduces DeepLC's MS1 identification depth
(R²=0.980 vs DeepLC predictions) while running in-process ~30× faster and
removing the Python/DeepLC dependency entirely. We further introduce a
**reverse-sequence RT benchmark** that quantifies how well an RT model captures
sequence order — the property that distinguishes a peptide from its reversed
decoy — and show that uno's engineered features recover 72% of this signal that
composition-only models miss completely. Finally, a two-proteome **entrapment
validation** confirms uno's protein-level FDR is approximately honest (estimated
empirical FDR 1.7% at a reported 1%), and—using the same harness—demonstrates
that a learned (LDA) protein rescorer, despite raising apparent identifications,
inflates the empirical FDR to ~11%, underscoring why the deterministic
binomial scorer is retained. [PENDING: multi-dataset extension.]

**Keywords:** MS1-only proteomics; DirectMS1; retention time prediction;
reproducibility; Rust; open-source software; label-free quantification.

---

## 1. Introduction

Bottom-up proteomics conventionally identifies peptides from MS2 fragment
spectra. MS1-only identification (DirectMS1 [REF: Ivanov DirectMS1]) takes a
different route: it matches the accurate precursor mass of an LC–MS feature to a
theoretical peptide and confirms the assignment with predicted retention time,
omitting MS2 entirely. Because no time is spent acquiring fragment spectra, MS1-
only workflows support very short gradients (minutes) and correspondingly high
throughput, at the cost of an information-poor identification problem: mass and
RT alone are weak discriminators, and identity must be concentrated at the
protein level across multiple consistent peptides.

The established open-source tool for this approach, **ms1searchpy**
[REF: ms1searchpy], couples a precursor mass search to RT prediction by **DeepLC**
[REF: DeepLC] and a LightGBM-based rescoring stage, with protein-level inference
by a binomial model. It is effective but has two practical costs that matter at
scale. First, its identification counts are **not reproducible run to run**: the
DeepLC and LightGBM stages are stochastic and unseeded, so repeated searches of
the same features yield materially different protein lists (§3.1). Second, it
carries a **heavyweight Python/deep-learning dependency** (DeepLC, TensorFlow)
whose RT-prediction step dominates wall-clock time for a single run.

We present **uno**, an MS1-only identification and label-free quantitation engine
written in Rust. uno is one component of an all-Rust pipeline: the **koth**
feature finder performs alignment, deisotoping and XIC/LFQ; uno consumes koth's
aligned feature matrix and assigns peptide/protein identities, controls FDR, and
rolls intensities up into quant matrices. uno's design goals are determinism,
speed, and minimal external dependencies, while matching the identification
depth of ms1searchpy. The contributions of this work are:

1. A clean, open-source Rust reimplementation of DirectMS1-style MS1 search with
   **deterministic** output (§2).
2. A quantitative demonstration that **MS1-only identification counts are
   strongly implementation-dependent**: ms1searchpy spans 152–303 proteins on one
   fixed feature set, and uno sits at its median with ~50× lower variance (§3.1).
3. A **distilled, pure-Rust RT model** with engineered order-aware features that
   reproduces DeepLC's MS1 identification depth (R²=0.980) at ~30× lower cost and
   no Python dependency (§3.3).
4. A **reverse-sequence RT benchmark** that isolates an RT model's sequence-order
   sensitivity — the basis for decoy discrimination — and a measurement of how
   far feature-based models close the gap to a sequence CNN (§3.4).

---

## 2. Experimental Procedures (Methods)

### 2.1 Input and pipeline position

uno consumes the aligned feature matrix produced by koth: per-feature
`massCalib`, `mz`, `charge`, `rtApex`, ion mobility (`im`, where available),
`combined_score`, and per-sample intensities. RT is in koth's reference-run space,
so a single calibration suffices rather than per-file calibration.

### 2.2 Peptide database and decoys

In-silico tryptic digestion is performed with sage-core [REF: Sage] (enzyme
rules, missed cleavages, length/mass bounds, fixed/variable modifications), and
masses are recomputed in f64. Candidate lookup is a binary search over a
mass-sorted peptide vector. Decoys are generated by sequence shuffling with an
optional mass shift; we use shuffled (random-mass) decoys because they provide a
clean protein-level FDR null, whereas same-mass reversed decoys are
RT-indistinguishable "twins" of their targets under composition-based RT models
(see §3.4 and [REF: PLAN/decoy]).

### 2.3 Feature–peptide matching

Each feature is matched to all peptides within a ppm tolerance, enumerating
isotope-error offsets (0, ±1 Da/charge). The resulting peptide–feature matches
(PFMs) are reduced to the top-K best-mass candidates per feature to bound the
downstream cost without discarding the RT-correct assignment.

### 2.4 Retention-time model (this work, §3.3)

uno predicts RT with a gradient-boosted-tree model (the pure-Rust `gbdt` crate)
**distilled from DeepLC**: DeepLC predictions over a large, diverse tryptic digest
serve as training labels. Sequences are encoded by a 99-dimensional **order-aware
feature set**: amino-acid composition (20), N- and C-terminal residue identities
(first three, last three; 6), **position-binned composition** (three sequence
bins × 20 = 60), and 13 scalar descriptors (length, ln mass, Kyte–Doolittle /
Krokhin-SSRCalc / Eisenberg hydrophobicity sums, an N-terminus-weighted SSRCalc
term, proline/aromatic/basic counts, and a hydropathy moment). The model is
distilled offline once; at search time it predicts RT for all candidate peptides
in-process.

Because the frozen predictor lives in DeepLC's abstract RT space, uno calibrates
it to the run's observed gradient using a RANSAC-seeded bootstrap on confident
identifications [REF: PLAN/RANSAC].

### 2.5 Scoring and FDR

PFMs are scored by a deviation discriminant combining recalibrated mass error
(in σ units, via sigma-clipped re-estimation of the mass error distribution) and
calibrated RT deviation. Because mass and RT alone do not separate target from
decoy at the single-peptide level in this regime (§3.2), identity is
concentrated at the protein level: a binomial survival score
`-log10 P(X≥k | n, p)` evaluates whether a protein's k observed peptides exceed
the random-match background rate p estimated from decoys, followed by greedy
parsimony, target/decoy pair competition, and a target/decoy FDR walk. Mass/RT
acceptance thresholds are optimized on a grid to maximize target proteins at the
chosen FDR (mirroring ms1searchpy's `worker`).

All stages are implemented to be **deterministic**: there is no use of
process-randomized iteration order or unseeded RNGs in any result-affecting path,
so repeated runs on identical input produce byte-identical output (§3.1). An
optional learned (cross-validated LDA) protein rescorer is provided but **disabled
by default**, because entrapment validation shows it inflates the empirical FDR
(§3.6); the binomial model is the validated default.

### 2.6 Quantitation

Identified peptides/proteins are mapped onto koth's per-sample intensities to
emit peptide- and protein-level quant matrices; uno does not re-quantify (koth
owns alignment/XIC/LFQ).

### 2.7 Benchmark data and comparison

Benchmarking uses the HeLa dataset PXD015669 (run `QEHF1_09744_JB`, Q Exactive
HF), with features called by koth. ms1searchpy (v[REF], `-ml 1 -deeplc 1 -ts 2`,
1% FDR) and uno were run on the **identical** koth feature set. A reference
protein set of 268 identifications (ms1searchpy on biosaur2 features) and an
orthogonal DDA-derived protein list are used to assess overlap. All software was
run on [HARDWARE — to fill in].

---

## 3. Results and Discussion

### 3.1 MS1-only identification counts are strongly implementation-dependent

Running ms1searchpy on a single fixed koth feature set six times yielded
**152, 153, 155, 188, 203, and 303 proteins** at a nominal 1% FDR (σ≈55; the
logged random-match background p ranged 0.016–0.037). The high outlier (303) is
consistent with the stochastic LightGBM rescoring occasionally overfitting the
target/decoy boundary. uno, run on the same features through its deterministic
(pre-calibrated RT) path, returned **175–179 proteins (σ≈1)** — sitting at the
**median** of ms1searchpy's distribution with ~50× lower variance (Figure 1).

This reframes the comparison: rather than uno being "behind" any single
ms1searchpy run, the two engines are at equal identification depth, and uno
removes a previously-unreported source of run-to-run variability in MS1-only
search. [PENDING: extend to n≥10 runs and ≥3 datasets; Figure 1 = paired
violin/strip plot of per-engine count distributions.]

uno's determinism is by design and was verified by construction. Two sources of
run-to-run nondeterminism in early development were traced to Rust's
process-randomized hash-map iteration order: (i) the RANSAC RT-calibration seed
was assembled from a hash map and sampled positionally, and (ii) the protein
greedy-parsimony and FDR steps broke score ties by iteration order. Both were
made order-independent (sort the calibration seed by feature index; resolve ties
on accession). After these fixes, three independent runs on identical input are
**byte-identical** — including on the fast in-process RT path (§3.3), where on
the upgraded koth feature set uno deterministically reports 218 proteins / 169
reference-overlap (previously varying 184–217 across runs). Determinism is thus a
guaranteed property of the implementation, not a statistical observation.

### 3.2 Identity is concentrated at the protein level, not the peptide level

We confirmed directly that in this regime mass+RT cannot separate target from
decoy peptides: among matches with near-zero RT deviation (<0.01 min), target and
decoy peptides occur at ≈50:50. Most target peptides in a proteome are not
actually present in a given sample, so a random target match is no more likely
than a decoy match. Discrimination therefore arises only from **multiplicity** —
real, present proteins accumulate several RT-consistent peptides while decoy
proteins remain singletons — which the binomial protein model exploits. This
property is intrinsic to MS1-only search and bounds what any single-peptide
discriminant (including ML rescoring) can achieve. [Supports the limitations in
§4.]

### 3.3 A distilled pure-Rust RT model reaches DeepLC's identification depth

The RT predictor is the heaviest dependency in MS1-only search. We asked whether
a fast, dependency-free model can replace DeepLC without losing depth.

A naïve composition-feature gbdt (the prior 69-d Klammer set) trailed: R²=0.972
vs DeepLC predictions, and 141 proteins / 119 reference-overlap on the benchmark
(vs the DeepLC library's 176 / 152). Replacing the encoding with the 99-d
**order-aware** feature set (§2.4) and increasing tree depth raised the distilled
model to **R²=0.980 (RMSE 6.62), matching the LightGBM ceiling**, and restored
identification depth to **~180 proteins / ~152 reference-overlap — at DeepLC
parity** (Table 1). Critically, this is achieved **in-process in pure Rust**: RT
prediction for ~5×10⁵ candidate peptides takes ~4–5 s within the search, versus
~130 s for a separate DeepLC Python preprocessing step — a ~30× reduction for the
RT stage, with no Python dependency.

We found the **boosting backend was not the limiter** (the Rust gbdt crate
reaches R²=0.972 vs LightGBM's 0.975 on identical features); the **feature
encoding** was the lever. A still-richer 215-d feature set gave only marginal
improvement (R²=0.9824), indicating the 99-d set is near the practical ceiling for
hand-crafted features.

**Table 1.** RT model vs identification depth on the HeLa benchmark (1% FDR).

| RT model | R² vs DeepLC | Proteins / ref-overlap | RT-stage cost | Dependency |
|---|---|---|---|---|
| Composition gbdt (69-d) | 0.972 | 141 / 119 | ~4 s in-process | none (Rust) |
| **Order-aware gbdt (99-d, this work)** | **0.980** | **~180 / ~152** | **~4 s in-process** | **none (Rust)** |
| DeepLC library | (teacher) | 176 / 152 | ~130 s preprocessing | Python/TensorFlow |

[Note: protein counts on this single file are noise-dominated by the on-data
calibration step; reference-overlap (~152) is the robust metric and is reported
as the primary comparison. PENDING: multi-dataset confirmation.]

### 3.4 A reverse-sequence benchmark quantifies sequence-order capture

A peptide and its reverse have identical amino-acid composition but different
true retention time; an order-aware predictor must assign them different RTs,
while a composition-only model is mathematically forced to predict the same value
for both. This forward−reverse RT difference is precisely the signal that lets a
search engine distinguish a real peptide from a reversed decoy.

We constructed a benchmark of 46,304 held-out tryptic peptides and their reverses
(non-tryptic, out-of-distribution), predicted DeepLC RT for both, and measured
how well each model reproduced the forward−reverse RT shift (DeepLC std = 12.25
native units; Table 2). A composition-only model captured **0%** of this signal
(it predicts identical RT for a sequence and its reverse, by construction). uno's
order-aware model captured **72%** (Pearson r=0.68 between model and DeepLC
forward−reverse differences) while predicting reverse-sequence RT nearly as well
as forward (R²_rev=0.96, R²_fwd=0.98). The residual ~28% is the genuine gap
between hand-crafted features and DeepLC's sequence convolution, and is not
closed by additional features or tree depth.

**Table 2.** Reverse-sequence RT benchmark (46,304 held-out peptides + reverses).

| Model | R²_fwd | R²_rev | order-Δ correlation | order magnitude captured |
|---|---|---|---|---|
| Composition-only | 0.96 | 0.95 | 0.00 | 0% |
| **uno order-aware (this work)** | **0.98** | **0.96** | **0.68** | **72%** |

This benchmark is, to our knowledge, a novel and direct way to report RT-model
sequence-order sensitivity, and is reusable for any RT predictor.

### 3.5 Two-proteome entrapment confirms the protein-level FDR is honest

Target/decoy FDR can be optimistic if the decoys do not faithfully model the
random-coincidence behaviour of real-but-absent peptides. We tested uno's
reported FDR directly with a two-proteome entrapment: the HeLa (human) sample was
searched against a database of human targets plus an entrapment proteome —
*Saccharomyces cerevisiae* (6,066 proteins), which is absent from the sample —
together with shuffled decoys for both. Any accepted yeast protein is a known
false positive. We estimate the empirical false-discovery proportion as
`n_yeast_accepted × (N_human_db / N_yeast_db) / n_human_accepted`.

At a reported 1% protein FDR, uno's binomial scorer accepted **201 human and 1
yeast** protein, an estimated **empirical FDR of 1.7%** — i.e. approximately
honest, the small excess being within Poisson noise of a single entrapment hit
(Table 3). This validates uno's shuffled-decoy protein FDR on this dataset.
[PENDING: repeat with a larger/size-matched entrapment proteome to tighten the
estimate, and across the multi-dataset panel of §5.]

### 3.6 Learned protein rescoring raises apparent counts but inflates FDR

We evaluated whether a learned protein-level rescorer — combining the binomial
score with log peptide count, coverage, and RT/mass consistency — could increase
depth, as ms1searchpy's PSM-level machine learning does. To avoid reintroducing
run-to-run variance, we used a deterministic, closed-form Fisher LDA with 3-fold
cross-validation by accession hash (no RNG). It raised apparent identifications
(218→242 proteins; 169→178 reference-overlap) and was byte-identical across runs.

However, the **same entrapment test exposed the gain as illusory**: the LDA
accepted **7 yeast (false) proteins versus the binomial's 1**, an estimated
empirical FDR of **11.5%** at a nominal 1% (Table 3).

To determine whether this was avoidable overfitting, we hardened the rescorer with
textbook semi-supervised safeguards — high-confidence iterated positives
(q ≤ 0.01, refreshed; rather than a loose set), target/decoy pairs co-assigned to
the same cross-validation fold, and stronger covariance shrinkage. This **removed
the inflated apparent count (242 → 223 proteins) but did not restore FDR control**
(still 11.4% empirical; Table 3), and it accepted *fewer* reference-validated
proteins than the binomial. The inflation is therefore **not classical
overfitting** but an **objective mismatch**: a learned scorer is trained to
separate targets from *shuffled* decoys, yet the false positives that matter are
*real-but-absent* proteins (e.g. the yeast entrapment), whose retention-time and
mass consistency resemble present proteins and differ from shuffled decoys. The
very features that could add depth (RT/mass/score consistency) are exactly those
that real-but-absent proteins share with present ones, so optimizing
target-vs-decoy separation promotes absent real proteins. The count-based binomial
is robust precisely because it uses only peptide multiplicity, where an absent
protein and a shuffled decoy are indistinguishable. We therefore **retain the
binomial scorer as the default** and ship the LDA disabled. This is a transferable
cautionary result: learned protein rescoring in MS1-only search trades FDR
integrity for apparent depth even under careful regularization, and entrapment
validation is essential to detect it.

**Table 3.** Two-proteome entrapment (HeLa human sample vs human + yeast database).

| Scorer | Human accepted | Yeast (false) accepted | Estimated empirical FDR | Reported |
|---|---|---|---|---|
| **Binomial (default)** | 201 | 1 | **1.7%** | 1% |
| LDA rescorer (naïve) | 208 | 7 | 11.5% | 1% |
| LDA rescorer (hardened CV) | 209 | 7 | 11.4% | 1% |

### 3.7 Speed [PENDING — finalize on standardized hardware]

Indicative figures on [HARDWARE]: uno core search (DB build + matching +
discriminant scoring) ~3–5 s; RT prediction in-process ~4–5 s; ms1searchpy total
~40–45 s per run including its internal DeepLC stage. [PENDING: clean,
contention-free wall-clock table; cold vs warm (RT-library reuse); peak memory.]

---

## 4. Limitations

We report these candidly, as they delimit where MS1-only search and uno apply.

- **Inherent PSM-FDR floor.** Mass+RT cannot separate target/decoy at the
  single-peptide level (§3.2); identification depends on protein-level
  multiplicity. This is a property of MS1-only data, not of uno.
- **No ion mobility, no peptide-centric depth.** On Orbitrap data without ion
  mobility, peptide-centric extraction from raw hills does not add identifications
  (target/decoy enrichment ≈1.06× on this benchmark): the (m/z, RT) coordinate is
  too degenerate without a third dimension. Ion-mobility data (timsTOF) is
  expected to lift this; uno predicts CCS via an analogous distilled model where
  IM is present, but we do not yet have an IM dataset with ground truth.
- **Search-breadth expansion is counterproductive.** Adding missed cleavages or
  modifications inflates the theoretical count faster than the observed count and
  raises the decoy background, reducing identifications — consistent with
  DirectMS1's use of a minimal digest.
- **TMT-labeled MS1 data.** uno currently yields ~0 identifications on TMT MS1
  data; DeepLC RT for TMT-modified peptides correlates only r≈0.71 with observed
  RT, too weak to seed calibration.
- **The fast-RT identification gain is noise-limited on a single file.** Depth-8
  vs depth-5 improved clean RT metrics (R², order capture) but not protein count,
  which is saturated/noisy on one dataset; multi-dataset evaluation is required to
  resolve real differences.

---

## 5. Validation roadmap (pre-submission checklist) — TODO

These are required before submission and are not yet complete:

1. **[DONE — extend to multi-dataset] Entrapment / two-proteome FDR validation.**
   Implemented (§3.5): HeLa human sample searched against human + *S. cerevisiae*
   entrapment; binomial FDR estimated honest (1.7% at reported 1%), and the LDA
   rescorer shown to inflate it (11.5%). Remaining: repeat on a larger/size-matched
   entrapment proteome and across the §5.2 dataset panel to tighten the estimate.
2. **Multi-dataset benchmark (≥3 datasets)** spanning instruments/gradients/
   organisms; per dataset: uno vs ms1searchpy protein counts, FDR validity,
   variance (n≥10 runs each), and speed.
3. **Formal reproducibility experiment** (n≥10 runs × both engines × each
   dataset; mean±SD) → headline figure.
4. **Clean speed/memory benchmark** on standardized hardware (cold vs warm).
5. **Quantitation validation** if the paper claims quant (known-ratio sample,
   CV/accuracy) — otherwise scope to ID-only and state so.
6. Software release hygiene for JPR's open-source expectations: build that does
   not require a local sage checkout (resolve the path dependency), README,
   example data, tests/CI, archived release with DOI.

---

## 6. Conclusions

uno is a fast, deterministic, open-source Rust engine for MS1-only protein
identification that matches the identification depth of ms1searchpy while removing
its run-to-run variability and its heavyweight Python/DeepLC dependency. Its
output is byte-identical across runs by construction, its distilled order-aware RT
model reaches DeepLC's MS1 identification depth in pure Rust at ~30× lower cost,
and two-proteome entrapment confirms its protein-level FDR is approximately
honest. The accompanying reverse-sequence benchmark provides a direct, reusable
measure of RT-model sequence-order sensitivity, and our entrapment analysis shows
that learned protein rescoring — while raising apparent counts — inflates FDR,
motivating uno's transparent binomial model. [PENDING: final depth/variance/speed
claims contingent on the multi-dataset panel of §5.]

---

## Data and Software Availability

uno is open-source (MIT) at [URL]. Benchmark data: PXD015669. RT distillation
labels, trained models, and analysis scripts: [DOI/Zenodo — TODO].

## References

See `paper.bib`. Key references to include: DirectMS1 (Ivanov et al.);
ms1searchpy; DeepLC (Bouwmeester et al.); biosaur2; Sage (Lazear); Krokhin
SSRCalc; Kyte–Doolittle; Eisenberg hydrophobicity; the Rust language / `gbdt`
crate; mokapot/Percolator (semi-supervised rescoring); pyteomics.
