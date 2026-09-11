#!/usr/bin/env bash
# Prespecified RT-model architecture analysis for the JPR Technical Note.
#
# The shipped sequence CNN and uno's legacy feature-engineered MLP are run from
# the same release binary, with identical feature tables, FASTA, digest, scoring
# configuration, and serialized timing. Existing outputs are reused unless
# FORCE=1, so an interrupted run can resume without mixing configurations.
set -euo pipefail

cd "$(dirname "$0")"
ROOT=".."
UNO="$ROOT/target/release/uno"
FASTA="data/entrap_human_yeast.fasta"
CONFIG="uno.toml"
OUT="results/rt_model_ablation"
TIME=/usr/bin/time
DATASETS=(
  hela5min_09746
  hela5min_09766
  hela5min_09734
  hela5min_09778
  hela5min_09950
  faims_LUM2_05433_VG
)

mkdir -p "$OUT"

parse_time() {
  awk '
    /User time/          {u=$NF}
    /Maximum resident/   {r=$NF/1024}
    /Elapsed .wall clock/{ n=split($NF,a,":");
                           w = (n==3)? a[1]*3600+a[2]*60+a[3] : a[1]*60+a[2] }
    END{ printf "%.1f %.1f %.0f\n", u, w, r }
  ' "$1"
}

score() {
  python3 - "$1" <<'PY'
import sys
sys.path.insert(0, ".")
import entrapment as e

rows = e.load(sys.argv[1])
out, _, _ = e.curve(rows, thresholds=(0.01, 0.05))
print(out[0.01], out[0.05])
PY
}

printf "dataset\tmodel\tdepth_1pct\tdepth_5pct\tuser_cpu_s\twall_s\tpeak_rss_mb\n" \
  > "$OUT/summary.tsv"
git rev-parse HEAD > "$OUT/git_head.txt"
sha256sum "$UNO" "$CONFIG" "$FASTA" \
  "$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc" \
  "$ROOT/models/rt_mlp.bin" > "$OUT/input_hashes.sha256"
git status --short > "$OUT/worktree_status.txt"

for ds in "${DATASETS[@]}"; do
  features="data/$ds.features.tsv"
  test -s "$features"
  for model in cnn_v3 legacy_mlp; do
    if [[ "$model" == cnn_v3 ]]; then
      rt_model="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
    else
      rt_model="$ROOT/models/rt_mlp.bin"
    fi
    result="$OUT/${ds}.${model}.tsv"
    log="$OUT/${ds}.${model}.log"
    timing="$OUT/${ds}.${model}.time"
    if [[ "${FORCE:-0}" == 1 || ! -s "$result" ]]; then
      "$TIME" -v "$UNO" \
        --features "$features" \
        --fasta "$FASTA" \
        --config "$CONFIG" \
        --rt-model "$rt_model" \
        --out-proteins "$result" \
        > "$log" 2> "$timing"
    fi
    read -r d1 d5 < <(score "$result")
    read -r user wall rss < <(parse_time "$timing")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$ds" "$model" "$d1" "$d5" "$user" "$wall" "$rss" \
      >> "$OUT/summary.tsv"
    printf "%s %-10s 1%%=%s 5%%=%s wall=%ss rss=%sMiB\n" \
      "$ds" "$model" "$d1" "$d5" "$wall" "$rss"
  done
done

