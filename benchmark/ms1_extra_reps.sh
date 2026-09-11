#!/usr/bin/env bash
# Extra ms1searchpy replicates, so its run-to-run spread is measured with the
# same n as uno's seed sweep (benchmark/seed_sweep.sh, 10 seeds). run.sh did 3
# replicates; this adds reps 4..10 on the named datasets. Runs are SERIAL: two
# concurrent ms1searchpy processes need ~15 GB and have frozen this box before.
set -euo pipefail
cd "$(dirname "$0")"
MS1="${MS1:-../../FeatureFinder/.venv/bin/ms1searchpy}"
FASTA=data/entrap_human_yeast.fasta
DATASETS="${DATASETS:-hela5min_09746 faims_LUM2_05433_VG}"
REPS="${REPS:-4 5 6 7 8 9 10}"
for DS in $DATASETS; do
  for i in $REPS; do
    OUT="results/ms1_${DS}_rep$i"
    [ -s "$OUT/${DS}.features_proteins_full.tsv" ] && { echo "skip $DS rep$i"; continue; }
    mkdir -p "$OUT"
    nice /usr/bin/time -v "$MS1" "data/${DS}.features.tsv" -d "$FASTA" -deeplc 1 -ad 1 -o "$OUT" \
      > "$OUT/ms1.log" 2>"$OUT/ms1.time"
    echo "done $DS rep$i"
  done
done
