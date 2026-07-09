# RT model pretraining plan — make the in-house CNN generalize like DeepLC

## The problem (what the data forced us to conclude)

uno's URTC CNN (`models/rt_cnn_*.urtc`) is trained **only** on ProteomeTools/Prosit
synthetic iRT. It wins ~2× on that holdout but **loses ~2× in real search**:

| RT source (honest entrapment on run10, yeast-covered) | real IDs | run10 PWL-MAE |
|---|---|---|
| DeepLC generic (zero-shot)        | **827** | **0.76 min** |
| CNN ch64 + frozen-embed transfer  |   607   | ~2.0 min |
| CNN ch64 generic                  |   546   | 2.65 min |

The gap is **RT prediction accuracy on a real LC gradient**, not bootstrap mechanics
(final `rt_sigma` is a core-only artifact — DeepLC's core holds 89% of confident
peptides, the CNN's only 48%; the missing half lands 2–5 min off and never seeds).

Frozen-embedding run-specific transfer plateaus at ~2 min (the pretrained features cap
it; cleaning labels didn't help). From-scratch run-specific training is *worse* (3.77).
**The bottleneck is pretraining: the CNN never saw diverse real LC gradients.** DeepLC
did, which is why it generalizes zero-shot. So: pretrain the CNN on a diverse, real-LC,
common-scale corpus. Full evidence: memory `uno-rt-search-validation`.

## The fix: pretrain on the Chronologer harmonized corpus

DeepLC generalizes because of training breadth. The directly-usable equivalent is the
**Chronologer harmonized DB**: 2.25M unique peptides / 2.6M observations from **11
community datasets**, all mapped into one physicochemical scale — **Hydrophobic Index
(HI = %ACN at C18 elution)**. The common scale is the crux: it lets many gradients be
mixed into one training target, and uno's existing per-run PWL/RANSAC calibration maps
HI → observed minutes per run (exactly as it does for DeepLC/iRT today).

- Source: `searlelab/chronologer` (Apache-2.0), `data/Chronologer_DB_220308.gz`.
- Encoding: EncyclopeDIA `PeptideModSeq` (e.g. `ALSVLGC[+57.021464]GHTSSTK`), 6–50 AA,
  17 mod types, HI label column.
- Paper: Wilburn et al., *Deep learning from harmonized peptide libraries…* (biorxiv
  2023.05.30.542978 / PMC10312522).

## Steps

1. **Get + convert data** (done-tooling: `python/chronologer_to_labels.py`)
   ```
   git clone https://github.com/searlelab/chronologer
   python crates/uno-rt-train/python/chronologer_to_labels.py \
       chronologer/data/Chronologer_DB_220308.gz data/rt/chronologer/ --peek   # verify cols
   python crates/uno-rt-train/python/chronologer_to_labels.py \
       chronologer/data/Chronologer_DB_220308.gz data/rt/chronologer/
   ```
   Produces `data/rt/chronologer/{train,val,holdout}.tsv` (peptide-disjoint, HI labels).
   v1 keeps unmodified + carbamidomethyl-C only (URTC is AA-only). Expect a large
   fraction retained (most peptides are unmodified or Cam-C).

2. **Train the CNN on HI** (GPU venv; existing `python/train_rt.py`)
   ```
   python crates/uno-rt-train/python/train_rt.py --arch cnn --ch 64 --blocks 1,2,4 \
     --train data/rt/chronologer/train.tsv --val data/rt/chronologer/val.tsv \
     --out models/rt_cnn_chrono_ch64.urtc --epochs 80 --bs 4096 --patience 12
   ```
   Start with ch64 (capacity helped: ch32→ch64 seeds 302→418). Consider ch96/ch128 —
   2.25M diverse peptides can feed more capacity than 760k synthetic ones could.
   Optionally **mix** Prosit iRT in: convert Prosit `irt`→HI by affine-aligning shared
   peptides, then concatenate. Or **warm-start** from the existing Prosit model with
   `train_rt.py --init models/rt_cnn_ch64.urtc` (full-network fine-tune) to keep its
   sequence priors while adapting to real-LC breadth. Keep masking-after-every-conv
   (Rust inference depends on it). Verify Rust↔torch parity with `urtc_torch_predict.py`.

3. **VALIDATE IN SEARCH — not on holdout** (done-tooling: `scripts/validate_rt_in_search.sh`)
   ```
   cargo build --release -p uno-rt-train -p uno
   crates/uno-rt-train/scripts/validate_rt_in_search.sh models/rt_cnn_chrono_ch64.urtc
   ```
   Reports honest run10 IDs (yeast-covered entrapment) + run10 PWL-MAE. **This is the
   only metric that counts.** Holdout MAE already lied once. Target: PWL-MAE → ~0.76,
   honest IDs → ~827, within-1-min → ~89%.

4. **Iterate** on capacity / data mix / epochs **using the in-search number**. If it
   stalls well above 0.76 with full capacity + data, the architecture (not data) is the
   limit — revisit (deeper/attention, or accept DeepLC for production depth).

## Success / stop criteria
- **Ship** if a Chronologer-pretrained CNN reaches within ~5% of DeepLC's honest IDs
  (~785+) — that makes a pure-Rust, ~40× faster, no-Python RT path viable.
- **Stop / keep DeepLC** if it plateaus below ~700 honest after exhausting capacity and
  the Prosit+Chronologer data mix.

## Notes / gotchas
- **Entrapment libraries must cover yeast** or FDR is fake (the old `deeplc_pp.tsv`
  baseline of "775 @ 0%" was this bug). The harness predicts the full candidate set.
- Per-run calibration already handles absolute-RT scale, so HI vs iRT vs minutes as the
  training target only matters for *cross-dataset mixing* — pick HI and convert all in.
- Mods: extending URTC to Chronologer's 17 mod types (needs a featurizer change) is a
  bigger, separate win — do AA-only first, prove the generalization thesis, then add mods.
- Benchmark inputs currently live in `/tmp/cmp` (run10 matrix, entrap fasta, candidate
  set, `rtcalib.tsv`). Relocate them into the repo (e.g. `data/bench/run10/`) so the
  validation harness is reproducible across machines.
