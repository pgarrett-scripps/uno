#!/usr/bin/env bash
# Re-run ONLY uno (both arms) on the six primary datasets with the current
# target/release/uno; ms1searchpy rows in results/summary.tsv are kept verbatim.
# 2026-09-04: needed because the July uno outputs no longer match the binary.
set -u
cd "$(dirname "$0")"
ROOT=..; FASTA=data/entrap_human_yeast.fasta; RT=$ROOT/models/rt_cnn_chrono_global_v3_ch32.urtc; UNO=$ROOT/target/release/uno
P=/home/patrick-garrett/Data/primary_raw
parse_time() { awk '/User time/{u=$NF} /Maximum resident/{r=$NF/1024} /Elapsed .wall clock/{n=split($NF,a,":"); w=(n==3)?a[1]*3600+a[2]*60+a[3]:a[1]*60+a[2]} END{printf "%.1f %.1f %.0f",u,w,r}' "$1"; }
score1() { python3 -c "
import sys; sys.path.insert(0,'.')
import entrapment as e
rows=e.load('$1'); out,_,_=e.curve(rows); print(out[0.01], out[0.05])"; }
grep -v -P "^\S+\tuno\t" results/summary.tsv > results/summary.tsv.tmp
for DS in hela5min_09746 hela5min_09766 hela5min_09734 hela5min_09778 hela5min_09950 faims_LUM2_05433_VG; do
  nice /usr/bin/time -v $UNO --features data/$DS.features.tsv --fasta $FASTA --rt-model $RT --config uno.toml --out-proteins results/uno_$DS.tsv > results/uno_$DS.log 2> results/uno_$DS.time
  read U_USER U_WALL U_RSS < <(parse_time results/uno_$DS.time); read D1 D5 < <(score1 results/uno_$DS.tsv)
  printf "%s\tuno\t1\t%s\t%s\t%s\t%s\t%s\n" "$DS" "$D1" "$D5" "$U_USER" "$U_WALL" "$U_RSS" >> results/summary.tsv.tmp
  echo "uno[$DS] @1%=$D1 @5%=$D5 wall=${U_WALL}s"
  base=$(awk -v d=$DS '$1==d{print $2}' raw_mode/urls.txt | sed 's/\.raw$//')
  nice /usr/bin/time -v $UNO --raw $P/mzML/$base.mzML --fasta $FASTA --rt-model $RT --config uno.toml --out-proteins results/uno_raw_$DS.tsv > $P/logs/$DS.uno_raw.log 2> results/$DS.uno_raw.time
  cp results/uno_raw_$DS.tsv $P/results/ 2>/dev/null; cp results/$DS.uno_raw.time $P/logs/$DS.uno_raw.time 2>/dev/null
  echo "raw[$DS] @1%=$(score1 results/uno_raw_$DS.tsv | cut -d' ' -f1)"
done
# restore original row order: dataset blocks in the order of the old file, uno row first
python3 - <<'PY'
rows=[l.rstrip("\n") for l in open("results/summary.tsv.tmp")]
hdr=rows[0]; body=rows[1:]
order=[]; 
for l in open("results/prev_july/summary.tsv"):
    d=l.split("\t")[0]
    if d!="dataset" and d not in order: order.append(d)
out=[hdr]
for d in order:
    out+= [r for r in body if r.startswith(d+"\tuno\t")] + [r for r in body if r.startswith(d+"\tms1\t")]
open("results/summary.tsv","w").write("\n".join(out)+"\n")
PY
rm results/summary.tsv.tmp; echo RERUN_DONE
