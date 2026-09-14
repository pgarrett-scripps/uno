# timsTOF IM result — 13 September 2026

**Neither IM setting produced a consistent identification gain.** The primary
three-sigma filter improved 4 of 9 matched searches, worsened 4, and tied 1.
The stricter two-sigma filter improved 6 of 9, but its losses offset the gains;
neither met the predefined consistency rule.

Median present-species protein accessions across three seeds, at **1% empirical
entrapment FDP**, are shown below. These are retrospective accession counts,
not UNO's internally accepted q-value counts or pooled protein unions.

| Acquisition | IM off | IM 3σ (primary) | IM 2σ (sensitivity) |
|---|---:|---:|---:|
| HYE A2 | 263 | 256 | 210 |
| HYE B2 | 69 | 0 | 139 |
| HYE C2 | 381 | 217 | 395 |

The old bundled model expects 70 features, while UNO supplies 100; its feature
range therefore excludes the supplied charge. A replacement Sequoia model,
trained on the public [Meier CCS data](https://proteomicsml.org/datasets/ionmobility/Meier_TIMS.html),
reduces median relative IM error for doubly charged peptides in the separate
A1 calibration run from 10.4% to 1.2%. Better prediction did not translate into
stable protein gains: the primary filter rejected only about 5–6% of candidate
matches, and several searches still returned empty protein lists. The stricter
filter also turned one 413-accession baseline into an empty result.

All 27 paired searches used fresh KOTH 0.5.0 features, the existing Bruker
quality defaults, the same HYE-plus-absent-Arabidopsis database, and seeds 0–2.
An initial setup accidentally omitted those quality defaults; all completed
and interrupted unfiltered attempts are preserved separately, including the
one nonempty IM result. The setup was corrected to frozen main's existing
thresholds before any standard-workflow results were obtained; the model,
IM widths, seeds, database, and outcome rule were unchanged.

An independent calculation matched all **81 threshold evaluations** (1%, 3%,
5%). IM-off matches frozen main byte-for-byte, all 76 Rust tests pass, and
58,213 predictions are unchanged by the final model-width guard. The search
binary's exact source is commit `1223f92`; subsequent changes add the Sequoia
width guard, setup correction, and reporting without replacing that binary.

Keep IM experimental on `codex/im-timstof`; the published main branch and
`paper-2026-09-13` remain unchanged. The [machine-readable results](../benchmark/im_timstof_20260913.json)
contain every seed, internal-q counts, entrapment results, and provenance;
[the protocol](im-timstof-plan.md) records calibration and decision details.
