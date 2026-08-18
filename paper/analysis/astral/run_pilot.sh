#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/patrick-garrett/Repos/uno
DATA=${ASTRAL_DATA:-/home/patrick-garrett/Data/PXD066701_Astral}
STEM=CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng
RAW="$DATA/raw/$STEM.raw"
MZML="$DATA/mzML/$STEM.mzML"
FEATURES="$DATA/features/$STEM.features.tsv"
LOGS="$DATA/logs"
RESULTS="$DATA/results"
FASTA="$ROOT/benchmark/data/entrap_human_yeast.fasta"
UNO="$ROOT/target/release/uno"
RT="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
MS1=/home/patrick-garrett/Repos/FeatureFinder/.venv/bin/ms1searchpy

mkdir -p "$DATA/mzML" "$DATA/features" "$LOGS" "$RESULTS"

test "$(stat -c %s "$RAW")" = 1698374785

sha256sum "$RAW" "$UNO" "$RT" "$FASTA" > "$LOGS/pilot_sha256.txt"
git -C "$ROOT" rev-parse HEAD > "$LOGS/uno_head.txt"
git -C /home/patrick-garrett/Repos/FeatureFinder rev-parse HEAD \
  > "$LOGS/ms1searchpy_head.txt"
git -C "$ROOT" diff --binary > "$LOGS/uno_worktree.patch"
git -C /home/patrick-garrett/Repos/FeatureFinder diff --binary \
  > "$LOGS/ms1searchpy_worktree.patch"

if [[ ! -s "$MZML" ]]; then
  docker run --rm -v "$DATA:/data" \
    quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
    ThermoRawFileParser -i "/data/raw/$STEM.raw" -o /data/mzML -f 2 \
    > "$LOGS/${STEM}.conversion.log" 2>&1
fi

if [[ ! -s "$FEATURES" ]]; then
  biosaur2 "$MZML" -o "$FEATURES" \
    > "$LOGS/${STEM}.biosaur2.log" 2>&1
fi

head -n 1 "$FEATURES" | tee "$LOGS/${STEM}.feature_header.txt"
wc -l "$FEATURES" | tee "$LOGS/${STEM}.feature_rows.txt"

/usr/bin/time -v "$UNO" \
  --features "$FEATURES" \
  --fasta "$FASTA" \
  --rt-model "$RT" \
  --config "$ROOT/benchmark/uno.toml" \
  --out-proteins "$RESULTS/uno_${STEM}.tsv" \
  > "$LOGS/${STEM}.uno.log" \
  2> "$LOGS/${STEM}.uno.time"

cp "$FEATURES" "$RESULTS/"
cp "$FASTA" "$RESULTS/"
(
  cd "$RESULTS"
  /usr/bin/time -v "$MS1" "$(basename "$FEATURES")" \
    -d "$(basename "$FASTA")" -deeplc 1 -ad 1 \
    > "$LOGS/${STEM}.ms1.log" \
    2> "$LOGS/${STEM}.ms1.time"
)

python3 "$ROOT/benchmark/entrapment.py" \
  uno "$RESULTS/uno_${STEM}.tsv" \
  ms1 "$RESULTS/${STEM}.features_proteins_full.tsv" \
  | tee "$RESULTS/${STEM}.entrapment.txt"
