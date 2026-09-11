#!/usr/bin/env bash
# Preregistered CW0710 input-series extension (see README.md). Strictly serial.
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
DATA=${DATA:-/home/patrick-garrett/Data/PXD066701_Astral}
LIST=$ROOT/paper/analysis/astral/input_series_files.tsv
UNO=$ROOT/target/release/uno
RT=$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc
FASTA=$ROOT/benchmark/data/entrap_human_yeast.fasta
MS1=${MS1:-$ROOT/../FeatureFinder/.venv/bin/ms1searchpy}
ENT=$ROOT/benchmark/entrapment.py
N_STOCK=3; N_PL=3
mkdir -p "$DATA"/{mzML,features,logs,results/uno,results/ms1,results/ms1_pl}
export DATA
tail -n +2 "$LIST" | while IFS=$'\t' read -r file bytes url; do
  stem=${file%.raw}; raw=$DATA/raw/$file; mzml=$DATA/mzML/$stem.mzML; feat=$DATA/features/$stem.features.tsv
  if [[ ! -f "$raw" || "$(stat -c %s "$raw")" != "$bytes" ]]; then echo "SKIP $stem: raw missing/incomplete" >> "$DATA/logs/input_series.progress"; continue; fi
  [[ -f "$DATA/logs/$stem.raw.sha256" ]] || sha256sum "$raw" > "$DATA/logs/$stem.raw.sha256"
  if ! grep -q 'Processing completed 0 errors, 0 warnings' "$DATA/logs/$stem.conversion.log" 2>/dev/null; then
    rm -f "$mzml"
    docker run --rm -v "$DATA:/data" quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
      ThermoRawFileParser -i "/data/raw/$file" -o /data/mzML -f 2 > "$DATA/logs/$stem.conversion.log" 2>&1
    grep -q 'Processing completed 0 errors, 0 warnings' "$DATA/logs/$stem.conversion.log" || { echo "FAIL convert $stem" >> "$DATA/logs/input_series.progress"; continue; }
  fi
  if ! grep -q 'Feature detection is finished' "$DATA/logs/$stem.biosaur2.log" 2>/dev/null; then
    rm -f "$feat"; nice biosaur2 "$mzml" -o "$feat" > "$DATA/logs/$stem.biosaur2.log" 2>&1
    grep -q 'Feature detection is finished' "$DATA/logs/$stem.biosaur2.log" || { echo "FAIL biosaur2 $stem" >> "$DATA/logs/input_series.progress"; continue; }
  fi
  head -n 1 "$feat" > "$DATA/logs/$stem.feature_header.txt"; wc -l "$feat" > "$DATA/logs/$stem.feature_rows.txt"
  uno_out=$DATA/results/uno/$stem.tsv
  if [[ ! -s "$uno_out" ]]; then
    nice /usr/bin/time -v "$UNO" --features "$feat" --fasta "$FASTA" --rt-model "$RT" --config "$ROOT/benchmark/uno.toml" \
      --out-proteins "$uno_out" > "$DATA/logs/$stem.uno.log" 2> "$DATA/logs/$stem.uno.time"
  fi
  python3 "$ENT" uno "$uno_out" > "$DATA/results/uno/$stem.entrapment.txt"
  run_ms1() { # outdir logprefix extra-args...
    local out=$1 lp=$2; shift 2; mkdir -p "$out"
    local full="$out/${stem}.features_proteins_full.tsv"
    if [[ ! -s "$full" ]]; then
      cp "$feat" "$FASTA" "$out/"
      ( cd "$out" && nice /usr/bin/time -v "$MS1" "${stem}.features.tsv" -d "$(basename "$FASTA")" -deeplc 1 -ad 1 "$@" \
          > "$DATA/logs/$lp.log" 2> "$DATA/logs/$lp.time" )
      rm -f "$out/${stem}.features.tsv" "$out/$(basename "$FASTA")" "$out"/*PFMs*.tsv
    fi
    python3 "$ENT" ms1 "$full" > "$out/entrapment.txt" 2>&1
  }
  for rep in $(seq -w 1 $N_STOCK); do run_ms1 "$DATA/results/ms1/$stem/rep$rep" "$stem.ms1.rep$rep"; done
  pl=$DATA/external/pl_diann_$stem.tsv
  if [[ -s "$pl" ]]; then
    for rep in $(seq -w 1 $N_PL); do run_ms1 "$DATA/results/ms1_pl/$stem/rep$rep" "$stem.ms1pl.rep$rep" -pl "$pl"; done
  else echo "NOPL $stem" >> "$DATA/logs/input_series.progress"; fi
  echo "DONE $stem" >> "$DATA/logs/input_series.progress"
done
echo "ALLDONE" >> "$DATA/logs/input_series.progress"
