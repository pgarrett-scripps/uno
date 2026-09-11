#!/usr/bin/env bash
# ms1searchpy RT-calibration controls on the three Astral files (see RT_CALIBRATION.md).
# Not part of the preregistered benchmark protocol; diagnostic controls only.
#   sweep   : shipped CLI options, rep1 only, 1 run each
#   pl      : -pl with the DIA-NN MS2 calibration list (needs build_diann_calib_lists.py), 3 files x 3 runs
#   robust  : local -robust_calib 1 patch (ms1searchpy_rt_calibration.patch), 3 files x 3 runs
# Runs inside a set were launched 3-way concurrent; wall times are NOT comparable to the
# serialized production timings.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
DATA=${DATA:-/home/patrick-garrett/Data/PXD066701_Astral}
OUT=${OUT:-$DATA/results/rt_calibration}
MS1=${MS1:-$ROOT/../FeatureFinder/.venv/bin/ms1searchpy}
ENT=$ROOT/benchmark/entrapment.py
FASTA=$DATA/results/entrap_human_yeast.fasta
STEMS=(CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng)
run() { # name stem args...
  local name=$1 stem=$2; shift 2; local out=$OUT/$name/$stem; mkdir -p "$out"
  cp "$DATA/features/$stem.features.tsv" "$FASTA" "$out/"
  ( cd "$out" && /usr/bin/time -v "$MS1" "$stem.features.tsv" -d entrap_human_yeast.fasta -deeplc 1 -ad 1 "$@" > stdout.log 2> stderr.log )
  python3 "$ENT" ms1 "$out/$stem.features_proteins_full.tsv" > "$out/entrapment.txt"
}
case ${1:-all} in
  sweep|all)
    s=${STEMS[0]}
    run sweep/v0_base "$s" -save_calib; run sweep/v1_es1 "$s" -es 1; run sweep/v2_ptol3 "$s" -ptol 3
    run sweep/v3_ci2 "$s" -ci 2 -csc 2; run sweep/v4_es1_ptol3 "$s" -es 1 -ptol 3; run sweep/v5_es1_ptol3_ci2 "$s" -es 1 -ptol 3 -ci 2 -csc 2 ;;&
  pl|all)
    for r in 1 2 3; do for s in "${STEMS[@]}"; do run pl/rep$r "$s" -pl "$DATA/external/pl_diann_$s.tsv"; done; done ;;&
  robust|all)
    for r in 1 2 3; do for s in "${STEMS[@]}"; do run robust/rep$r "$s" -robust_calib 1; done; done ;;
esac
