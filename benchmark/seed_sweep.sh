#!/usr/bin/env bash
# uno seed-to-seed spread on the six primary datasets.
#
# uno is byte-deterministic for a fixed seed; that is a property of seeded code,
# not of the algorithm. This measures what ms1searchpy's three replicates
# measure -- the spread of the depth over independent random draws -- by running
# uno once per seed (`--seed N`, N = 1..10; seed 0 = the shipped default already
# in results/uno_<DS>.tsv). Depths are scored with entrapment.py exactly as
# run.sh scores everything else. Results: results/seed_sweep/summary.tsv.
#
#   SEEDS="1 2 3" JOBS=3 benchmark/seed_sweep.sh
set -euo pipefail
cd "$(dirname "$0")"
ROOT=..; FASTA=data/entrap_human_yeast.fasta
RT=$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc; UNO=$ROOT/target/release/uno
OUT=results/seed_sweep; mkdir -p "$OUT"
SEEDS="${SEEDS:-1 2 3 4 5 6 7 8 9 10}"
DATASETS="${DATASETS:-hela5min_09746 hela5min_09766 hela5min_09734 hela5min_09778 hela5min_09950 faims_LUM2_05433_VG}"
JOBS="${JOBS:-3}"
export OUT UNO FASTA RT
run_one() {
  local DS=$1 SEED=$2
  local P="$OUT/uno_${DS}_seed${SEED}"
  [ -s "$P.tsv" ] && return 0
  nice /usr/bin/time -v "$UNO" --features "data/$DS.features.tsv" --fasta "$FASTA" --rt-model "$RT" \
    --config uno.toml --seed "$SEED" --out-proteins "$P.tsv" > "$P.log" 2> "$P.time"
  echo "done $DS seed $SEED"
}
export -f run_one
for DS in $DATASETS; do for S in $SEEDS; do echo "$DS $S"; done; done \
  | xargs -P "$JOBS" -n 2 bash -c 'run_one "$0" "$1"'

python3 - "$OUT" "$DATASETS" "$SEEDS" <<'PY'
import sys, os
sys.path.insert(0, '.')
import entrapment as e
out, datasets, seeds = sys.argv[1], sys.argv[2].split(), sys.argv[3].split()
rows = ["dataset\tseed\tdepth_1pct\tdepth_5pct"]
for ds in datasets:
    for seed in ['0'] + seeds:
        p = f"results/uno_{ds}.tsv" if seed == '0' else f"{out}/uno_{ds}_seed{seed}.tsv"
        if not os.path.exists(p):
            continue
        d, _, _ = e.curve(e.load(p))
        rows.append(f"{ds}\t{seed}\t{d[0.01]}\t{d[0.05]}")
open(f"{out}/summary.tsv", "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
PY
