# External Orbitrap Astral validation

Pilot data live outside the repository at
`/home/patrick-garrett/Data/PXD066701_Astral`.

## Preregistered dataset

- Deposit: PXD066701, Ryan Kelly laboratory, Brigham Young University.
- Sample: 10 ng Pierce HeLa digest.
- Instrument: Orbitrap Astral with FAIMS Pro at one fixed compensation voltage
  (-45 V).
- LC: dual-column short-cycle method (11.4 min total cycle, near-100% duty
  cycle), as described in the deposited project metadata.
- Selected complete matched set:
  - `CW0710_OldFAIMS_10ng__RA4_rep1_ch1_DIA10ng.raw`
  - `CW0710_OldFAIMS_10ng__RA4_rep4_ch2_DIA10ng.raw`
  - `CW0710_OldFAIMS_10ng__RA4_rep6_ch2_DIA10ng.raw`

Selection was made from file names and acquisition metadata before viewing any
uno or ms1searchpy result. Files will not be replaced based on search outcome.

## Pilot gate

The first file is processed through ThermoRawFileParser 1.4.5 and biosaur2
0.2.28. Production searches proceed only if conversion finishes, the feature
table contains the columns required by both engines, and both engines complete
RT calibration and protein inference. Failure is reported as an instrument or
feature-finder compatibility result rather than silently selecting another run.

## Comparison

Both engines receive the identical biosaur2 feature table, human-plus-yeast
target FASTA, and digest. uno uses the paper's frozen RT model. ms1searchpy uses
DeepLC and the same local quality-column patch as the main benchmark. All timed
runs are serialized and captured with `/usr/bin/time -v`. The common
`benchmark/entrapment.py` scorer reports human depth at 1% and 5% entrapment q.

The working tree is currently dirty, so this pilot records both the repository
HEAD and binary hashes. No Astral values enter the manuscript until the unified
benchmark is rerun from a committed, pinned build.

The completed production run and focused comparator controls are assessed in
`AUDIT.md`. `run_audit_controls.sh` tests stock-equivalent feature columns with
both DeepLC and ms1searchpy's default RT model; these are diagnostic controls,
not additional benchmark replicates.

## Data provenance

Official project and file API responses are preserved as
`project_metadata.json` and `file_metadata.json` in the external data root.
PRIDE reports no source checksum, so completed downloads are verified against
the deposited byte size and assigned a local SHA-256 checksum.
