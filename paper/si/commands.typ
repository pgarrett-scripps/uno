// ADAPTED (not verbatim) from benchmark/README.md, benchmark/run.sh,
// benchmark/faims/run_faims.sh and benchmark/quant/README.md. The scripts loop
// over datasets and replicates and set their paths internally; this listing
// unrolls one iteration of each and defines every variable it uses, so the
// commands run as printed. Paths are relative to the repository root.
//
// Two deliberate differences from the scripts. run.sh scores by importing
// `entrapment` inline via `python3 -c`; the equivalent CLI form is shown here.
// run_faims.sh copies each feature file into $OUT and runs ms1searchpy from
// there, which is reproduced below, because ms1searchpy writes next to its input
// and needs a writable working directory.
```sh
# ---------------------------------------------------------------- variables
ROOT=$PWD                                   # the repository root
DS=hela5min_09746                           # one of the six (see Table S2)
OUT=$ROOT/benchmark/results/ms1_${DS}       # ms1searchpy writes here
RT=$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc
MS1=/path/to/ms1searchpy                    # the patched 2.8.12, see Section S1
mkdir -p "$OUT"

# 0. inputs: raw file -> mzML -> single-run biosaur2 features
#    (benchmark/README.md; run once per dataset, output named per Table S2)
docker run --rm -v $PWD:/data \
  quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
  ThermoRawFileParser -i /data/RUN.raw -o /data -f 2
biosaur2 RUN.mzML -o $ROOT/benchmark/data/${DS}.features.tsv

# 1. identification benchmark, both engines, identical input (benchmark/run.sh)
/usr/bin/time -v $ROOT/target/release/uno \
  --features $ROOT/benchmark/data/${DS}.features.tsv \
  --fasta    $ROOT/benchmark/data/entrap_human_yeast.fasta \
  --rt-model $RT \
  --config   $ROOT/benchmark/uno.toml \
  --out-proteins $ROOT/benchmark/results/uno_${DS}.tsv

# ms1searchpy writes beside its input, so it is run from $OUT on a copy.
cp $ROOT/benchmark/data/${DS}.features.tsv "$OUT/"
( cd "$OUT" && /usr/bin/time -v $MS1 ${DS}.features.tsv \
    -d $ROOT/benchmark/data/entrap_human_yeast.fasta -deeplc 1 -ad 1 )

# 2. entrapment scoring, identical scorer for both engines
python3 $ROOT/benchmark/entrapment.py \
  uno $ROOT/benchmark/results/uno_${DS}.tsv \
  ms1 $OUT/${DS}.features_proteins_full.tsv

# 3. FAIMS benchmark (benchmark/faims/run_faims.sh), E. coli entrapment
DATA=/path/to/PXD023977                     # the deposit, downloaded
FOUT=$ROOT/benchmark/faims/faims_out; mkdir -p "$FOUT"
R=05433                                     # replicate: 05433, 05434 or 05435

# Target-only DB: strip the deposit's pre-made DECOY_ entries, since both
# engines generate their own. The two counts set the entrapment scale factor.
awk '/^>/{p=(!/^>DECOY_/)} p' $DATA/sprot_human_and_ecoli_shuffled.fasta \
  > $FOUT/hye_target.fasta
export N_TARGET_DB=$(grep -c '_HUMAN' $FOUT/hye_target.fasta)
export N_ENTRAP_DB=$(grep -c '_ECOLI' $FOUT/hye_target.fasta)

$ROOT/target/release/uno --features $DATA/LUM2_${R}_VG.features.tsv \
  --fasta $FOUT/hye_target.fasta --rt-model $RT \
  --config $ROOT/benchmark/uno.toml --out-proteins $FOUT/uno_${R}.tsv
cp $DATA/LUM2_${R}_VG.features.tsv "$FOUT/"
( cd "$FOUT" && $MS1 LUM2_${R}_VG.features.tsv -d hye_target.fasta \
    -deeplc 1 -ad 1 )

cd $ROOT/benchmark/faims
python3 score_entrap.py uno $FOUT/uno_${R}.tsv
python3 score_entrap.py ms1 $FOUT/LUM2_${R}_VG.features_proteins_full.tsv

# 4. quantitation (benchmark/quant/README.md), PXD033350, 24 runs
QDATA=/path/to/PXD033350                    # the deposit, downloaded
QRES=$ROOT/benchmark/quant/results; mkdir -p "$QRES"

# 4a. per-run koth features, then cross-run alignment. averagine_projection is
#     NOT koth's default and must be set; the box-sum arm it replaces regresses
#     the entrapment-scored protein count (Section S7).
printf '[lfq]\naveragine_projection = true\n' > $QRES/align.toml
koth_align $QDATA/koth_runs -o $QDATA/align_avg -c $QRES/align.toml

# 4b. target-only DB: 48 UPS1 + E. coli + human (the known-absent probe)
awk '/^>/{p=(!/^>DECOY_/)} p' \
  $QDATA/sprot_ecoli_human_with_ups_replaced_shuffled.fasta \
  > $QRES/ups_ecoli_human.fasta

# 4c. --matrix-min-runs 6 is REQUIRED, not a tuning choice: without the floor the
#     mass-anchor background rises and the search returns nothing (Section 3.9).
$ROOT/target/release/uno \
  --matrix $QDATA/align_avg/intensity_matrix.tsv \
  --matrix-min-runs 6 --fasta $QRES/ups_ecoli_human.fasta \
  --rt-model $RT --profile-reassign \
  --out-abundance $QRES/ab.tsv --out-lfq $QRES/lfq.tsv \
  --out-proteins $QRES/proteins.tsv

python3 $ROOT/benchmark/quant/dilution_track.py $QRES/ab.tsv
python3 $ROOT/benchmark/quant/peptide_fc_de.py  $QRES/lfq.tsv 50 0_1
```
