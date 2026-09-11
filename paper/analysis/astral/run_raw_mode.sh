#!/usr/bin/env bash
# uno --raw (in-process koth feature finder) on the three Astral mzML files.
# Same FASTA, RT model, and uno config as the biosaur2-fed production runs
# (run_all.sh); the only difference is the feature finder. Serialized, timed.
set -u
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
D=/home/patrick-garrett/Data/PXD066701_Astral
UNO="$ROOT/target/release/uno"
FASTA="$ROOT/benchmark/data/entrap_human_yeast.fasta"
RT="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
CFG="$ROOT/benchmark/uno.toml"
OUT="$D/results/uno_raw"; mkdir -p "$OUT"
for f in "$D"/mzML/*.mzML; do
  n=$(basename "$f" .mzML)
  echo "=== $n $(date -u +%FT%TZ)"
  /usr/bin/time -v "$UNO" --raw "$f" --fasta "$FASTA" --rt-model "$RT" --config "$CFG" \
    --out-proteins "$OUT/$n.tsv" > "$D/logs/$n.uno_raw.log" 2> "$D/logs/$n.uno_raw.time"
  echo "exit=$? $(grep -E 'Elapsed|Maximum resident' "$D/logs/$n.uno_raw.time" | tr '\n' ' ')"
done
echo DONE
