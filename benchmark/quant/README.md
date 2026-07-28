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
3 reps, 24 runs). Pipeline: raw -> koth features -> koth_align -> `uno --matrix
--profile-reassign` -> the two downstream scripts here. DB = 48 UPS1 (spike) +
E. coli (flat background) + human (known-absent entrapment).

    uno --matrix intensity_matrix.tsv --matrix-min-runs 6 --fasta ups_ecoli.fasta \
        --rt-model <model> --profile-reassign --out-abundance ab.tsv --out-lfq lfq.tsv
    python dilution_track.py ab.tsv                 # spike tracks? background flat?
    python peptide_fc_de.py  lfq.tsv 50 0_1         # DirectMS1Quant-style contrast

**Result (identification + tracking: strong).** uno IDs 361 proteins @1% (honest
entrapment FDP ~1.8%). Spike tracks the dilution, background/absent stay flat:

| species | role       | median rho (abundance vs log fmol) | frac(rho>0.5) |
|---------|------------|------------------------------------|---------------|
| UPS1    | spike      | **+0.75**                          | 81%           |
| E. coli | flat bg    | +0.06                              | 3%            |
| human   | absent     | +0.12                              | 0%            |

**Differential calling (specificity strong, low-abundance sensitivity trails).**
At the extreme 50-vs-0.1 contrast uno flags ~12 changed UPS proteins vs
DirectMS1Quant's 42, but with high specificity (few false positives). The gap is
NOT the test - it is the quant values feeding it, in two separable parts:

- **Fill specificity** — koth's default box-sum inflates dilute cells with
  co-eluting interference. koth's opt-in `lfq.averagine_projection = true`
  (matched-filter projection onto the mass-derived averagine pattern) roughly
  halves false-positive calls (precision 47% -> 60%). Enable it for uno's quant
  path (it is off by default in koth, which optimizes feature quality not DE).
- **Detection sensitivity** — tracing spike peptides confident at 50 fmol to
  0.1 fmol: 63% stay detected; of those lost, the signal is present but the
  peptide elutes in 1-2 scans and koth's hill-builder (`min_scans=2`) drops it.
  A blanket `min_scans=1` recovers ~50% but floods 6x noise features and is a net
  loss end-to-end (differential precision 27% -> 9%). The clean fix is a targeted,
  ID-seeded matched-filter re-extraction (confirm known targets, do not discover) —
  scoped as future work; a global detection threshold trades specificity for recall.

Bottom line: uno identifies and tracks abundance as well as DirectMS1Quant with
better determinism and specificity; trace-level differential sensitivity is the
one open gap, with a scoped remedy.
