# MS1 quant: uno output + downstream stats

uno is an ID + quant tool. For a multi-sample (koth-aligned) matrix it emits a
protein-group x sample ABUNDANCE matrix and stops there:

    uno --matrix intensity_matrix.tsv --matrix-min-runs N --fasta db.fasta \
        --rt-model <model> --profile-reassign --out-abundance abundance.tsv

That path does the quant-correctness work internally (design-free):
feature competition (1 peak -> 1 protein), profile-guided reassignment (route a
feature to the protein whose across-sample profile it tracks), and protein grouping
(indistinguishable proteins -> one group; unique + razor peptides). No design,
no differential testing, no FDR — those belong downstream.

## Downstream stats (this is NOT uno)
`de_validate.py` is an example downstream analysis (MSstats/limma-style) that takes
uno's abundance matrix + an experimental design and calls what changed. It also
scores against the IonStar species ground truth, so it doubles as our validation
harness. Design is inferred from run-name substrings; no design file needed here.

    python de_validate.py abundance.tsv --dilution A:1,B:2,C:3,D:4,E:5   # order only
    python de_validate.py abundance.tsv --contrast A,E

## Validation status (IonStar PXD003881, full real-MBR matrix, all 5 levels)
Downstream dilution DE on uno's abundance matrix: ecoli 54/100 @5% FDR (92%
precision), 58% @10%, 64% @20%. Human (true-negative) stays the false-positive
class. Recall is higher than the old in-uno peptide-level test because the
protein-level trend over 20 runs has more power, on uno's clean abundances.

`make_consensus.py` is a stand-in for koth_align (builds a matrix from per-run
features) kept for when raw data isn't available.

## DirectMS1Quant dilution benchmark (PXD033350: UPS1-into-E.coli)

Head-to-head vs DirectMS1Quant on its own spike-in series (8 levels 0.1-50 fmol x
3 reps, 24 runs). Pipeline: raw -> koth features -> koth_align (**`[lfq]
averagine_projection = true`**) -> `uno --matrix --profile-reassign` -> the two
downstream scripts here. DB = 48 UPS1 (spike) + E. coli (flat background) + human
(known-absent entrapment).

    # koth_align config: [lfq] averagine_projection = true   <-- the standard quant path
    uno --matrix intensity_matrix.tsv --matrix-min-runs 6 --fasta ups_ecoli.fasta \
        --rt-model <model> --profile-reassign --out-abundance ab.tsv --out-lfq lfq.tsv
    python dilution_track.py ab.tsv                 # spike tracks? background flat?
    python peptide_fc_de.py  lfq.tsv 50 0_1         # DirectMS1Quant-style contrast (indicative)

**Estimator choice — use averagine.** Under the current koth feature finder (after
the downward-chain-walk removal) the raw box-sum **regressed the honest protein-ID
count by 21%** (306 -> 243) while the averagine projection **held it** (305 -> 326).
Averagine projects each cell onto its mass-derived isotope pattern and rejects the
co-eluting interference the box-sum absorbs; it is mass-based and ID-FDR-free.
Dilution tracking itself is stable across both estimators, but averagine is the
standard path. (Box-sum copies kept as `results/*_boxsum.tsv` for the record.)

**Result (new koth + averagine, re-measured 2026-07-28; ID + tracking: strong).**
uno IDs **326** proteins @1% (honest human-entrapment FDP **0.34%**). Spike tracks
the dilution, background/absent stay flat (n = proteins quantified in >=8 of 24 runs):

| species | role       | n   | median rho (abundance vs log fmol) | frac(rho>0.5) |
|---------|------------|-----|------------------------------------|---------------|
| UPS1    | spike      | 18  | **+0.71**                          | 89%           |
| E. coli | flat bg    | 303 | +0.02                              | 4%            |
| human   | absent     | 5   | +0.25                              | 0%            |

This dilution trend is the **primary** quantitative result — it is robust to the
match-between-runs over-fill that makes the differential contrast fragile.

**Differential calling (indicative only, low-abundance sensitivity trails).**
DirectMS1Quant reports 42 changed UPS proteins at the extreme 50-vs-0.1 contrast;
uno's peptide-fold-change test recovers roughly half. We do NOT report a precision
figure: the floor-imputation DE is imputation-sensitive (e.g. 50-vs-0.1 gives
22 UPS / 23 E. coli / 60 human at BH-1%+2-fold — the FP swings with the feature set),
so it is treated as indicative. The real, robust signal for the gap is **detection**:

- Tracing UPS1 precursors confident at 50 fmol down to 0.1 fmol, the aligned-matrix
  cell stays **filled for ~82%** but **independent detection (koth per-cell q<=0.01)
  falls to ~1%**. The spread is match-between-runs over-fill; the lower trace is what
  bounds differential sensitivity. Of those lost, the signal is present but the
  peptide elutes in 1-2 scans and koth's hill-builder (`min_scans=2`) drops it.
- A blanket `min_scans=1` does NOT fix it (see `results/minscans_ab/README.md`):
  features x6.4, mass-anchor background p 0.15 -> 0.51, uno accepts **0 proteins @1%**
  (the search collapses). Untargeted detection cannot separate a one-scan peptide
  from one-scan noise. The clean fix is a targeted, ID-seeded matched-filter
  re-extraction (confirm known targets, do not discover) — scoped as future work.

Bottom line: uno identifies and tracks abundance as well as DirectMS1Quant, with
better determinism; trace-level differential sensitivity is the one open gap, with a
scoped remedy. Use the averagine path.

### 2026-07-28 re-measurement note (the earlier numbers did not reproduce)

Rebuilt uno (`cargo build --release`, commit 7c80199) and re-ran the benchmark on
both koth_align matrices. Commands (from the repo root):

    # target-only DB: 48 UPS1 + 4431 E. coli + 20193 human (UPS entries count as human)
    awk '/^>/{p=(!/^>DECOY_/)} p' \
      /mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/part1/sprot_ecoli_human_with_ups_replaced_shuffled.fasta \
      > benchmark/quant/results/ups_ecoli_human.fasta

    # M = align (default box-sum) and align_avg (averagine projection)
    ./target/release/uno \
      --matrix /mnt/labdata/pgarrett/PXD033350_DirectMS1Quant/uno_proc/$M/intensity_matrix.tsv \
      --matrix-min-runs 6 --fasta benchmark/quant/results/ups_ecoli_human.fasta \
      --rt-model models/rt_cnn_chrono_global_v3_ch32.urtc --profile-reassign \
      --out-abundance ab.tsv --out-lfq lfq.tsv --out-proteins proteins.tsv

    python benchmark/quant/dilution_track.py benchmark/quant/results/ab.tsv
    python benchmark/quant/peptide_fc_de.py  benchmark/quant/results/lfq.tsv 50 0_1
    python paper/scripts/figures/fig5_quant.py
    python paper/scripts/figures/figS5_quant_detail.py

Saved outputs: `benchmark/quant/results/` (`{ab,lfq,proteins}{,_avg}.tsv`, the uno
logs, `analysis_stdout.txt`, `fig5_stdout.txt`, `figS5_stdout.txt`, and
`dilution_track_median_fix.txt`).

What changed and why:

- **361 proteins @1% -> 306** (align) / 305 (align_avg). The 361 does not
  reproduce on either matrix.
- **honest entrapment FDP ~1.8% -> 0.07%** (1 human hit of 306 accepted). The
  earlier figure was ~25x too pessimistic.
- **median rho +0.75 -> +0.71, frac 81% -> 88%.** Two causes. First,
  `dilution_track.py` printed `sorted(r)[len(r)//2]`, the upper-middle order
  statistic, not a median; on the n=16 UPS set that reads +0.74 where the true
  median is +0.71. Fixed to `statistics.median` in this commit. Second, the human
  row (+0.12, 0%) came from the **averagine** matrix while the UPS and E. coli
  rows came from the box-sum matrix, so the old table mixed two matrices. On the
  box-sum matrix human is n=1 at +0.09.
- **"averagine roughly halves false positives, precision 47% -> 60%" is wrong.**
  Measured: it raises true positives 8 -> 14 and false positives 5 -> 8, precision
  62% -> 64%. It is a sensitivity lever, not a precision lever, and precision was
  already ~62% with the default box-sum.
- **the single "63% stay detected" does not reproduce under either definition.**
  It sits between the two traces above (83% matrix-filled, 6% independently
  detected). Report both; the honest pair is more informative than the single
  number, since it isolates match-between-runs over-fill as the cause.
