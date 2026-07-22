# Canonical benchmark results (DirectMS1 / ms1searchpy paper datasets)

5-min HeLa from **PXD015669** (DirectMS1, Anal Chem 2020) + 1 FAIMS run from **PXD023977**
(ms1searchpy/DirectMS1Quant "2000 proteins"). Same biosaur2 features fed to both engines;
human+yeast entrapment; ms1 with **DeepLC** (`-deeplc 1 -ad 1`), 3 reps; uno deterministic, 1.

| dataset | uno @1% | ms1 @1% (med, range) | uno @5% | ms1 @5% med | uno wall | ms1 wall | uno RSS | ms1 RSS |
|---|---|---|---|---|---|---|---|---|
| hela5min_09746 | **862** | 837 (823–853) | **1209** | 1122 | 33s | 195s | 3.3G | 6.9G |
| hela5min_09766 | **875** | 804 (789–913) | **1100** | 989 | 33s | 232s | 3.3G | 7.4G |
| hela5min_09734 | **1035** | 880 (770–951) | **1245** | 1073 | 32s | 203s | 3.2G | 7.0G |
| hela5min_09778 | **1222** | 1173 (1147–1220) | **1653** | 1510 | 40s | 239s | 3.9G | 7.5G |
| hela5min_09950 | **1440** | 1279 (1169–1469) | **1740** | 1740 | 37s | 245s | 3.4G | 7.4G |
| faims_LUM2_05433_VG | **2072** | 1483 (1481–2081) | **2916** | 2319 | 45s | 234s | 4.1G | 7.7G |

**uno wins @1% on all six datasets**; wins or ties @5%. ~5–6× faster wall, ~2× less RAM,
byte-deterministic vs ms1's ±13–40% run-to-run variance. On the reference tools' own data.

FAIMS note: uno reads the run as one (ignores the −50/−65/−80 V CV fractions) and still hits
2072 @1%, matching the DirectMS1Quant headline — CV is orthogonal upside uno currently discards.
