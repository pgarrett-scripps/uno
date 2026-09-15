# IM and LFQ integration

## Included work

This integration combines the completed IM working-copy snapshot (`f3d09b0`), the paired timsTOF study (`291e793`), and the LFQ consolidation (`17ff014`) on top of main's publication-freeze commit (`2fff764`). The original study commits remain in the merged history, including the 31 recovered LFQ research files documented in `lfq-consolidation-20260914.md`.

The completed bundled CNN, soft IM rescoring, charge-aware prediction cache, CLI, and Python training workflow take precedence over the earlier optional `UNO_IM_GATE` experiment. The experimental gate was deliberately not restored: the completed implementation uses IM as a classifier feature after RT shortlisting and explicitly avoids IM cutoffs. The earlier model-width validation, legacy tree loading, Sequoia training CLI, study scripts, and public training artifacts are retained. The old LFQ checkout's remaining training difference was an obsolete CLI help comment; its five deferred IM file versions are covered by the earlier study history and local recovery snapshot.

The manuscript keeps the IM model/search sections and the KOTH 0.6.0 LFQ findings. Its joint-search paragraph combines the updated LFQ results with the exploratory, dataset-dependent interpretation and the warning that pooled identification does not establish detection in every sample. Shared statistics and IM/LFQ tables were regenerated from their recorded inputs. The entrapment assumption was joined to its explanatory paragraph without changing its words.

## Dependency and benchmark boundary

UNO uses a local path dependency on the sibling KOTH checkout. That checkout currently reports `koth_ff` 0.7.0, so Cargo.lock was reconciled from 0.5.0 to 0.7.0 to make the current local build reproducible with `--locked`; no registry package versions were changed. This records the existing sibling source version and does not rerun, replace, or relabel historical benchmark evidence. The LFQ manuscript continues to describe its frozen 0.6.0 study. The broader KOTH release rerun remains a separate task described in `koth-next-release-rerun-plan.md`.

## Validation

- Rust: 86 workspace tests passed (8 CLI, 77 core, 1 LightGBM binding), with no failures; remaining crate and documentation test targets also passed. The vendored LightGBM directory was supplied in LD_LIBRARY_PATH for the binding test runner.
- Python: 8 RT/IM workflow tests passed with the newly built UNO binary, including CLI training, evaluation, artifact validation, lineage checks, and Python/Rust prediction parity. All 8 matrix tests also passed.
- Manuscript: `just paper` and Word export passed. `just verify` passed formatting, extractor checks, prose rules, declarations for 180 statistics and 57 assets, output staleness, and PDF overflow checks. Existing advisory prose/declaration warnings remain visible.
- Word export: 53 native equations, 17 tables, and 11 images were exported; the document was rendered with the bundled LibreOffice renderer for layout inspection.
- These are integration and regression checks. The historical full benchmark campaigns were not rerun.

## Recovery

Before integration, all 43 changed/untracked main files were snapshotted and their hashes verified. The earlier LFQ snapshot retains all 139 changed/untracked entries from that checkout. Local recovery material is stored under the ignored `experiments_archive/` directory; it is not uploaded with the source repository. Obsolete worktree contents are preserved there before their registrations and merged branches are removed.
