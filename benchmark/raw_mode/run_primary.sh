#!/usr/bin/env bash
# benchmark/raw_mode/run_primary.sh
# Download primary-dataset .raw from PRIDE, convert with ThermoRawFileParser 1.4.5
# (docker, same as benchmark/README.md), then uno --raw (in-process koth features).
set -u
ROOT=/home/patrick-garrett/Repos/uno
P=/home/patrick-garrett/Data/primary_raw
UNO="$ROOT/target/release/uno"; FASTA="$ROOT/benchmark/data/entrap_human_yeast.fasta"
RT="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"; CFG="$ROOT/benchmark/uno.toml"
while read ds fname size url; do
  [ -z "$url" ] && { echo "no url for $ds"; continue; }
  echo "=== $ds $(date -u +%FT%TZ)"
  if [ ! -s "$P/raw/$fname" ] || [ "$(stat -c %s "$P/raw/$fname")" != "$size" ]; then
    curl -s -L --retry 5 -o "$P/raw/$fname" "$url" || { echo "download failed $ds"; continue; }
  fi
  echo "downloaded $(stat -c %s "$P/raw/$fname") of $size"
  base="${fname%.raw}"
  if [ ! -s "$P/mzML/$base.mzML" ]; then
    docker run --rm -v "$P/raw:/in" -v "$P/mzML:/out" quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
      ThermoRawFileParser -i "/in/$fname" -o /out -f 2 > "$P/logs/$ds.conversion.log" 2>&1
    echo "converted exit=$?"
  fi
  /usr/bin/time -v "$UNO" --raw "$P/mzML/$base.mzML" --fasta "$FASTA" --rt-model "$RT" --config "$CFG" \
    --out-proteins "$P/results/uno_raw_$ds.tsv" > "$P/logs/$ds.uno_raw.log" 2> "$P/logs/$ds.uno_raw.time"
  echo "uno exit=$? $(grep -E 'Elapsed|Maximum resident' "$P/logs/$ds.uno_raw.time" | tr '\n' ' ')"
done < "$ROOT/benchmark/raw_mode/urls.txt"
echo DONE

# Stage the protein tables and timing logs beside the biosaur2-arm outputs so
# summarize.py and the figure script find both arms in one place.
cp "$P"/results/uno_raw_*.tsv "$ROOT/benchmark/results/"
for t in "$P"/logs/*.uno_raw.time; do cp "$t" "$ROOT/benchmark/results/$(basename "$t" .uno_raw.time).uno_raw.time"; done
