#!/usr/bin/env bash
# uno vs ms1searchpy on the DirectMS1 FAIMS benchmark (PXD023977): 3 replicate
# 200 ng HeLa Lumos+FAIMS runs (LUM2_05433/34/35). Apples-to-apples: the SAME
# biosaur2 features, the SAME target-only human+E.coli DB (both tools make their
# own decoys), honest E. coli-entrapment FDP.
#
# Data (download once, see benchmark/faims/README.md):
#   $DATA/LUM2_054{33,34,35}_VG.features.tsv           biosaur2 features
#   $DATA/sprot_human_and_ecoli_shuffled.fasta         targets + DECOY_ (stripped below)
set -u
cd "$(dirname "$0")"
# Absolute paths throughout: the ms1searchpy step below runs inside $OUT, so a
# path relative to this directory would resolve against the wrong cwd there.
ROOT="$(cd ../.. && pwd)"
DATA="${DATA:-/mnt/labdata/pgarrett/PXD023977_FAIMS}"
UNO="$ROOT/target/release/uno"
RT="$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc"
MS1="${MS1:-$ROOT/../FeatureFinder/.venv/bin/ms1searchpy}"
OUT="${OUT:-faims_out}"; mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"

# Target-only DB (strip the pre-made DECOY_ entries; both engines add their own).
FASTA="$OUT/hye_target.fasta"
awk '/^>/{p=(!/^>DECOY_/)} p' "$DATA/sprot_human_and_ecoli_shuffled.fasta" > "$FASTA"
NH=$(grep -c '_HUMAN' "$FASTA"); NE=$(grep -c '_ECOLI' "$FASTA")
echo "DB: $NH human + $NE ecoli (target-only)"
export N_TARGET_DB=$NH N_ENTRAP_DB=$NE

for R in 05433 05434 05435; do
  F="$DATA/LUM2_${R}_VG.features.tsv"
  echo "### rep $R ###"
  # uno (deterministic)
  $UNO --features "$F" --fasta "$FASTA" --rt-model "$RT" --config "$ROOT/benchmark/uno.toml" \
       --out-proteins "$OUT/uno_${R}.tsv" > "$OUT/uno_${R}.log" 2>&1
  python3 score_entrap.py uno "$OUT/uno_${R}.tsv"
  # ms1searchpy (-deeplc 1 -ad 1). It writes next to the input + needs a writable
  # cwd, so copy the features into $OUT and run there.
  cp "$F" "$OUT/"; ( cd "$OUT" && $MS1 "LUM2_${R}_VG.features.tsv" -d "$(basename "$FASTA")" \
       -deeplc 1 -ad 1 > "ms1_${R}.log" 2>&1 )
  python3 score_entrap.py ms1 "$OUT/LUM2_${R}_VG.features_proteins_full.tsv"
done
echo "done -> $OUT/"
