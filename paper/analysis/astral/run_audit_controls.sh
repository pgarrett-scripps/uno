#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/patrick-garrett/Repos/uno
DATA=${ASTRAL_DATA:-/home/patrick-garrett/Data/PXD066701_Astral}
MS1=/home/patrick-garrett/Repos/FeatureFinder/.venv/bin/ms1searchpy
FASTA="$ROOT/benchmark/data/entrap_human_yeast.fasta"
STEM=CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng
SOURCE="$DATA/features/$STEM.features.tsv"
BASE="$DATA/audit_controls"

mkdir -p "$BASE"

# Stock ms1searchpy ignores biosaur2's isoerror and isoerror2 columns. Removing
# them exercises that behavior while leaving every standard input column intact.
STOCK="$BASE/$STEM.stock_columns.features.tsv"
"/home/patrick-garrett/Repos/FeatureFinder/.venv/bin/python3" -c 'import pandas as pd,sys; p=sys.argv[1]; o=sys.argv[2]; d=pd.read_csv(p,sep="\t"); d.drop(columns=[c for c in ("isoerror","isoerror2") if c in d],inplace=True); d.to_csv(o,sep="\t",index=False)' "$SOURCE" "$STOCK"

run_control() {
  label=$1
  feature=$2
  shift 2
  out="$BASE/$label"
  mkdir -p "$out"
  cp "$feature" "$out/"
  cp "$FASTA" "$out/"
  local_feature=$(basename "$feature")
  (
    cd "$out"
    /usr/bin/time -v "$MS1" "$local_feature" -d "$(basename "$FASTA")" \
      -ad 1 "$@" > "$label.log" 2> "$label.time"
  )
  full="$out/${local_feature%.tsv}_proteins_full.tsv"
  python3 "$ROOT/benchmark/entrapment.py" ms1 "$full" > "$out/entrapment.txt"
}

# One replicate is sufficient to diagnose whether either choice changes the
# order of magnitude; production stochasticity is already characterized by 10.
run_control stock_columns_deeplc "$STOCK" -deeplc 1
run_control stock_columns_default_rt "$STOCK" -deeplc 0
