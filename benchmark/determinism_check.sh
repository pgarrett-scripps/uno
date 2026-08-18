#!/usr/bin/env bash
# Two-run byte-diff determinism check.
#
# The manuscript claims uno is byte-deterministic on identical input. The unit
# tests in uno-core cover the specific reducers that were found to be order-
# dependent; this covers the claim as stated, end to end, on a real feature file.
#
# It is a script rather than a #[test] because the inputs it needs -- a benchmark
# feature file, a FASTA, and an RT model -- are gitignored or large, so a cargo
# test could only skip silently, which is the failure mode three existing tests in
# this repo already have.
#
# RUN THIS BEFORE THE RE-RUN. A failure here invalidates paper.typ's determinism
# claim, not just this script.
#
#   benchmark/determinism_check.sh [features.tsv] [fasta]
set -euo pipefail
cd "$(dirname "$0")/.."

# Defaults mirror benchmark/run.sh, so this checks the configuration the paper
# reports rather than some other one.
FEATURES="${1:-benchmark/data/09744.features.tsv}"
FASTA="${2:-benchmark/data/entrap_human_yeast.fasta}"
MODEL="models/rt_cnn_chrono_global_v3_ch32.urtc"
CONFIG="benchmark/uno.toml"
OUT="${TMPDIR:-/tmp}/uno-determinism-$$"

for f in "$FEATURES" "$FASTA" "$MODEL" "$CONFIG"; do
  if [ -z "$f" ] || [ ! -f "$f" ]; then
    echo "determinism_check: missing input '$f'" >&2
    echo "  usage: $0 [features.tsv] [fasta]" >&2
    exit 2
  fi
done

mkdir -p "$OUT"
trap 'rm -rf "$OUT"' EXIT

echo "building release binary..."
cargo build --release --quiet

run() {   # $1 = run label
  ./target/release/uno \
    --features "$FEATURES" \
    --fasta "$FASTA" \
    --rt-model "$MODEL" \
    --config "$CONFIG" \
    --out-proteins "$OUT/$1.proteins.tsv" \
    >"$OUT/$1.log" 2>&1
}

echo "run 1..."; run a
echo "run 2..."; run b

rc=0
if cmp -s "$OUT/a.proteins.tsv" "$OUT/b.proteins.tsv"; then
  echo "  OK    proteins: byte-identical ($(wc -l < "$OUT/a.proteins.tsv") lines)"
else
  echo "  FAIL  proteins: outputs differ"
  diff "$OUT/a.proteins.tsv" "$OUT/b.proteins.tsv" | head -20
  rc=1
fi

if [ $rc -eq 0 ]; then
  echo "determinism_check: PASS (two runs on $FEATURES are byte-identical)"
else
  echo "determinism_check: FAIL -- paper.typ's determinism claim does not hold" >&2
fi
exit $rc
