#!/usr/bin/env bash
# Apples-to-apples uno vs ms1searchpy benchmark, over multiple datasets.
#
# Held IDENTICAL per dataset: the biosaur2 single-run feature file, the human+yeast
# FASTA, digest params (uno.toml aligns uno to ms1's 7-30/0-MC). Each tool makes its
# own decoys; yeast is the entrapment false-positive class.
# Measured: entrapment depth @1%/@5%, user-CPU (load-immune), wall, peak RSS.
# ms1searchpy is non-deterministic -> N replicates; uno is deterministic -> 1.
# Runs are SERIALIZED so timings are not corrupted by concurrent load.
#
# Datasets = every data/*.features.tsv (run10 medium HeLa, deephela large 240min,
# 09744 small QEHF1). Override with DATASETS="run10 09744" to pick a subset.
set -u
cd "$(dirname "$0")"
ROOT=".."

FASTA="data/entrap_human_yeast.fasta"
RT_MODEL="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
UNO="$ROOT/target/release/uno"
# ms1searchpy console script from its installed venv (ms1searchpy + deeplc 3.1.13
# + pyteomics + lightgbm). Override with MS1=... if it lives elsewhere.
MS1="${MS1:-$ROOT/../FeatureFinder/.venv/bin/ms1searchpy}"
N_MS1="${N_MS1:-3}"                          # ms1 replicates (non-deterministic)
# Datasets: names of data/<name>.features.tsv. Default = all present.
if [ -z "${DATASETS:-}" ]; then
  DATASETS=""; for f in data/*.features.tsv; do DATASETS="$DATASETS $(basename "$f" .features.tsv)"; done
fi
TIME=/usr/bin/time

mkdir -p results
SUMMARY="results/summary.tsv"
printf "dataset\ttool\trep\tdepth_1pct\tdepth_5pct\tuser_cpu_s\twall_s\tpeak_rss_mb\n" > "$SUMMARY"

parse_time() {  # $1 = time-v logfile -> "user_s wall_s rss_mb"
  awk '
    /User time/          {u=$NF}
    /Maximum resident/   {r=$NF/1024}
    /Elapsed .wall clock/{ n=split($NF,a,":");
                           w = (n==3)? a[1]*3600+a[2]*60+a[3] : a[1]*60+a[2] }
    END{ printf "%.1f %.1f %.0f", u, w, r }
  ' "$1"
}
score1() {  # $1 = proteins tsv -> "d1 d5"
  python3 -c "
import sys; sys.path.insert(0,'.')
import entrapment as e
rows=e.load('$1'); out,_,_=e.curve(rows)
print(out[0.01], out[0.05])"
}

for DS in $DATASETS; do
  FEATURES="data/${DS}.features.tsv"
  [ -f "$FEATURES" ] || { echo "skip $DS (no $FEATURES)"; continue; }
  echo "############### DATASET: $DS ($(wc -l < "$FEATURES") features) ###############"

  echo "--- uno (deterministic, 1 run) ---"
  $TIME -v $UNO --features "$FEATURES" --fasta "$FASTA" \
    --rt-model "$RT_MODEL" --config uno.toml \
    --out-proteins "results/uno_${DS}.tsv" \
    > "results/uno_${DS}.log" 2>"results/uno_${DS}.time"
  read U_USER U_WALL U_RSS < <(parse_time "results/uno_${DS}.time")
  read U_D1 U_D5 < <(score1 "results/uno_${DS}.tsv")
  printf "%s\tuno\t1\t%s\t%s\t%s\t%s\t%s\n" "$DS" "$U_D1" "$U_D5" "$U_USER" "$U_WALL" "$U_RSS" >> "$SUMMARY"
  echo "uno[$DS]: @1%=$U_D1 @5%=$U_D5  user=${U_USER}s wall=${U_WALL}s rss=${U_RSS}MB"

  echo "--- ms1searchpy (-deeplc 1 -ad 1, ${N_MS1} reps) ---"
  for i in $(seq 1 "$N_MS1"); do
    OUT="results/ms1_${DS}_rep$i"; mkdir -p "$OUT"
    $TIME -v $MS1 "$FEATURES" -d "$FASTA" -deeplc 1 -ad 1 -o "$OUT" \
      > "$OUT/ms1.log" 2>"$OUT/ms1.time"
    PF="$OUT/${DS}.features_proteins_full.tsv"
    read M_USER M_WALL M_RSS < <(parse_time "$OUT/ms1.time")
    if [ -f "$PF" ]; then read M_D1 M_D5 < <(score1 "$PF"); else M_D1=NA; M_D5=NA; fi
    printf "%s\tms1\t%s\t%s\t%s\t%s\t%s\t%s\n" "$DS" "$i" "$M_D1" "$M_D5" "$M_USER" "$M_WALL" "$M_RSS" >> "$SUMMARY"
    echo "ms1[$DS] rep$i: @1%=$M_D1 @5%=$M_D5  user=${M_USER}s wall=${M_WALL}s rss=${M_RSS}MB"
  done
done

echo; echo "================= SUMMARY (results/summary.tsv) ================="
column -t "$SUMMARY"
