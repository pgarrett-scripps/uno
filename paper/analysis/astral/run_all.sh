#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/patrick-garrett/Repos/uno
DATA=${ASTRAL_DATA:-/home/patrick-garrett/Data/PXD066701_Astral}
N_MS1=${N_MS1:-10}
FASTA="$ROOT/benchmark/data/entrap_human_yeast.fasta"
UNO="$ROOT/target/release/uno"
RT="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
MS1=/home/patrick-garrett/Repos/FeatureFinder/.venv/bin/ms1searchpy

STEMS=(
  CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng
  CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng
  CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng
)
SIZES=(1698374785 1876478721 1880709561)

mkdir -p "$DATA/mzML" "$DATA/features" "$DATA/logs" "$DATA/results/uno" \
  "$DATA/results/ms1"

for i in "${!STEMS[@]}"; do
  stem=${STEMS[$i]}
  raw="$DATA/raw/$stem.raw"
  mzml="$DATA/mzML/$stem.mzML"
  features="$DATA/features/$stem.features.tsv"

  test "$(stat -c %s "$raw")" = "${SIZES[$i]}"
  sha256sum "$raw" > "$DATA/logs/$stem.raw.sha256"

  if ! grep -q 'Processing completed 0 errors, 0 warnings' \
      "$DATA/logs/$stem.conversion.log" 2>/dev/null; then
    rm -f "$mzml"
    docker run --rm -v "$DATA:/data" \
      quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
      ThermoRawFileParser -i "/data/raw/$stem.raw" -o /data/mzML -f 2 \
      > "$DATA/logs/$stem.conversion.log" 2>&1
  fi

  if ! grep -q 'Feature detection is finished' \
      "$DATA/logs/$stem.biosaur2.log" 2>/dev/null; then
    rm -f "$features"
    biosaur2 "$mzml" -o "$features" \
      > "$DATA/logs/$stem.biosaur2.log" 2>&1
  fi

  head -n 1 "$features" > "$DATA/logs/$stem.feature_header.txt"
  wc -l "$features" > "$DATA/logs/$stem.feature_rows.txt"

  uno_out="$DATA/results/uno/$stem.tsv"
  if [[ ! -s "$uno_out" ]]; then
    /usr/bin/time -v "$UNO" \
      --features "$features" \
      --fasta "$FASTA" \
      --rt-model "$RT" \
      --config "$ROOT/benchmark/uno.toml" \
      --out-proteins "$uno_out" \
      > "$DATA/logs/$stem.uno.log" \
      2> "$DATA/logs/$stem.uno.time"
  fi
  python3 "$ROOT/benchmark/entrapment.py" uno "$uno_out" \
    > "$DATA/results/uno/$stem.entrapment.txt"

  for rep in $(seq -w 1 "$N_MS1"); do
    out="$DATA/results/ms1/$stem/rep$rep"
    mkdir -p "$out"
    full="$out/${stem}.features_proteins_full.tsv"
    if [[ ! -s "$full" ]]; then
      cp "$features" "$out/"
      cp "$FASTA" "$out/"
      (
        cd "$out"
        /usr/bin/time -v "$MS1" "${stem}.features.tsv" \
          -d "$(basename "$FASTA")" -deeplc 1 -ad 1 \
          > "$DATA/logs/$stem.ms1.rep$rep.log" \
          2> "$DATA/logs/$stem.ms1.rep$rep.time"
      )
    fi
    python3 "$ROOT/benchmark/entrapment.py" ms1 "$full" \
      > "$out/entrapment.txt"
  done
done
