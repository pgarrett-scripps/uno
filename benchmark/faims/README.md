# FAIMS benchmark: uno vs ms1searchpy (DirectMS1 PXD023977)

DirectMS1/ms1searchpy's own FAIMS "2000-protein" benchmark: 3 replicate 200 ng
HeLa Lumos+FAIMS runs. Apples-to-apples — SAME biosaur2 features, SAME target-only
human+E.coli DB (both engines make their own decoys), honest E. coli-entrapment FDP.

## Data (download once)
ProteomeXchange PXD023977 (PRIDE):
  https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD023977/
Grab: LUM2_054{33,34,35}_VG.features.tsv, sprot_human_and_ecoli_shuffled.fasta.
Kept at /mnt/labdata/pgarrett/PXD023977_FAIMS (override with DATA=...).

## Run
    ./run_faims.sh                 # DATA=... MS1=... OUT=... to override

## Result (2026-07-27, honest E. coli-entrapment, mean of 3 reps)
| metric | uno  | ms1searchpy |
|--------|------|-------------|
| @1% FDR| 1591 | 1432  (+11% uno) |
| @5% FDR| 2060 | 1883  (+9%  uno) |

uno matches/beats the reference method on its own FAIMS data and is deterministic
(ms1searchpy swings rep-to-rep: 1659/1136/1503 @1%). uno ~2060 @5% = DirectMS1's
headline "2000". Complements the Orbitrap (no-FAIMS) result in ../README.md.
