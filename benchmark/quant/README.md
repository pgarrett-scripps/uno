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
