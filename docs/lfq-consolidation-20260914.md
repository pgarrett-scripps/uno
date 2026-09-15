# LFQ consolidation while IM development continues

Integration follow-up: see [IM and LFQ integration](integration-20260914.md) for the completed combination, deferred-IM decisions, and combined validation. The record below describes the earlier isolated LFQ stage.

## Scope and provenance

The `codex/lfq-consolidation` branch starts at committed `main` (`2fff76492352a974baffb8fb8e11d58049576f6e`) and incorporates both existing KOTH 0.6.0 LFQ commits, `db3eab7` and `c314b19`. It also recovers 31 previously untracked LFQ research scripts, result summaries, provenance records, and notes from `experiments/lfq-matrix-fixes/uno`, whose HEAD was `500d66a7214b42a17becef8e3ac1c58e5f3c68c0`.

Recovery preserves the original bytes and executable modes. Historical result files retain their original scope, paths, and claims; this consolidation did not rerun benchmarks or promote diagnostic experiments into new manuscript results. Scripts that launch studies or write result directories were syntax-checked without executing them.

## Source inventory

| Disposition | Entries |
| --- | ---: |
| Already identical to committed main | 99 |
| Unique LFQ files recovered here | 31 |
| IM differences preserved for later integration | 5 |
| Local data/results/cache symlinks left in the original checkout | 4 |

The original LFQ checkout remains intact. Before recovery, all 139 changed or untracked entries were archived with SHA-256 checks for regular files, symlink targets, the source HEAD, a binary tracked-change patch, and a worktree inventory. Archive contents and source stability were verified before copying. Ignored datasets and build directories were neither copied nor removed.

Local recovery archive (outside all experiment worktrees): `experiments_archive/lfq-consolidation-20260915T004344Z`. Its `manifest.json` records every disposition and checksum; `dirty-files.tar.gz` holds the file snapshot and `tracked.patch` records the tracked changes relative to the original LFQ HEAD. This ignored archive is machine-local and is not included in Git commits.

## Deferred IM differences

The following LFQ-checkout files differ from committed main because they contain IM model, training, or gating changes. Their complete original versions remain in the recovery archive and original checkout; none were copied into this LFQ branch.

- `Cargo.lock`
- `crates/uno-core/src/im.rs`
- `crates/uno-core/src/score.rs`
- `crates/uno-rt-train/Cargo.toml`
- `crates/uno-rt-train/src/main.rs`

All tracked files under `crates/`, plus `Cargo.toml` and `Cargo.lock`, remain identical to the starting main commit. The active main checkout and `codex/im-timstof` checkout were not modified by this consolidation.

## Validation

- Matrix suite: all 8 tests passed using `python -m unittest discover -s benchmark -p test_matrix_study.py -v` with the existing root Python environment.
- Recovered files: all 31 match the snapshot checksums; 16 Python scripts compile and 8 JSON records parse successfully. Seven Markdown notes were preserved byte for byte.
- `just paper`: passed in the isolated worktree, using the existing paper Python environment with environment synchronization and network access disabled.
- `just verify`: formatting, extractor tests, prose checks (0 errors, 24 warnings), stats declarations (164 values), asset declarations (54 assets), and PDF overflow checks passed. The overall command exits nonzero solely because the fresh worktree has no `paper.docx`; Word export remains a final integration check.
- External benchmark inputs absent from the isolated worktree are reported as unverified by the declaration checks. No expensive benchmark runs or deep numerical re-derivations were performed.
- The 31 recovered files are additive and introduce no Rust changes, so this consolidation did not rerun Rust tests. Rust and IM tests are required after the eventual combined integration.

The built LFQ manuscript has 4,857 main-text words and 2,950 SI words (7,807 total); Flesch–Kincaid grades are 12.1 and 11.6 respectively. These measurements describe this branch's existing KOTH manuscript, not the ongoing IM manuscript.

## Final integration boundary

Keep this branch isolated until IM development finishes. At that point, reconcile shared manuscript sources, statistics, and asset metadata with the completed IM changes, review the five deferred IM file differences, regenerate combined outputs, and run the relevant Rust/Python and manuscript checks before updating main. Worktree removal, branch deletion, and remote updates are deferred; no existing branch or checkout was removed.

## Recovered inventory

| File | SHA-256 |
| --- | --- |
| `benchmark/check_joint_bundle.py` | `17d32beb06c4f0c79d48b12c85585687292d60ff42cee5327af1b184085a45db` |
| `benchmark/finalize_joint_gap.py` | `8d645556c80d76fe31505a50a5aadb1fee4d813b48f854bd705562b1667f0265` |
| `benchmark/finalize_joint_matrix.py` | `19a03d32ba8a8310777b76c70540058b3b77f0b4fb1e8a76b1f38f4e15921330` |
| `benchmark/joint_cell_quality_pilot.py` | `081395d27e47613970574fb3dbcbb28d607e635d9a9ab2f223fad02fe9fd8de2` |
| `benchmark/joint_cell_quality_summary.json` | `f24eaa845722d8da7095b99f3c5dd354a19a39d2eaad4209f64aaf60941aa8b3` |
| `benchmark/joint_evidence_summary.json` | `d4b6733fe6edba2b2dd84337ba77cac581ace4c4fc6030e78d266ab94c9964e4` |
| `benchmark/joint_evidence_validation.py` | `35be9958198bcdfa049e3928e7ecf479d14ca1b293ea2b1864707785d0eae4f1` |
| `benchmark/joint_gap_alternatives.py` | `6a4f20abe607c764eee398f2a3510635a5ee1c3e53d1e465588e7020e7bb8b70` |
| `benchmark/joint_gap_astral_controls.py` | `9e5c8abd3c89a4a54bb829ab5e6995a1e5aae2f47a904057453d7ad56ff29ce8` |
| `benchmark/joint_gap_audit.py` | `45bc7b71763c736d457603bc49ef97e0e1349ffb88325920e49bde3e1c36e4d6` |
| `benchmark/joint_gap_pilot.py` | `68d7f4b49386cc2ca0524d1904e7d21d201f076a6bee7944c2563d61e6311fa1` |
| `benchmark/joint_gap_retention.py` | `d6787d583ba413164c5be1fd664d24ccad8aff8d7eacbda5cc488928e5b887bf` |
| `benchmark/joint_gap_stage_audit.py` | `9818dcbb4daa9765b046454dda2972fa6dc029f6a9e5cb6077c096c998a3dafd` |
| `benchmark/joint_gap_summary.json` | `1d0a9b1ef9472d55e9b31fad4e5d29d82f28f10466116b3005ec72cb92f85446` |
| `benchmark/joint_gap_validation.py` | `1413fb70c8d86b27c79d89172591c4ec9c83bb6d3a8a6d6aa4dd3f6b61c113a0` |
| `benchmark/joint_matrix_pilot.py` | `924181b39243ef76a3309ad3161044956bdcf4737b2344d621217a2bca2ab60d` |
| `benchmark/joint_matrix_summary.json` | `f17dc19fda60bf21f09e5720f2c629f5f13898b1ebc40f89b164322fa34a25d4` |
| `benchmark/lfq_matrix_fixes.py` | `338ed3a1719fe309031ff0719184fd297c97eb20b2e1a8071d3021f4f8ea3b81` |
| `benchmark/lfq_matrix_fixes_pilot.json` | `556b17f32914566631cb492a7fe03c6e32d53c0f8dd38bd4c18b70c69a53f79f` |
| `benchmark/matrix_release_evidence_links.json` | `83e97e1bcc31cc90a8f1f0c4cb1712b470d1b92c6f40e374105363b08bbabee0` |
| `benchmark/matrix_release_manuscript_baseline.json` | `b193491dd770601aa556311bd8d9d31fe671b4513930ae51b5cd396fd3b2bfea` |
| `benchmark/matrix_release_manuscript_delivery.json` | `e9a989128e6f45555b4449f20171c1bb630e01d1b680e82091e55231506c4155` |
| `benchmark/matrix_seed_ensemble.py` | `aec60a7739b9cc6028f4f94d9dbe6ef861bc48b4d5cbc00d6a32ecca3ec2d498` |
| `benchmark/matrix_seed_ensemble_check.py` | `e4c90692556b87511e79d439b13e391951369393c47aa64425bed0244af5d65d` |
| `docs/joint-cell-quality-results.md` | `e89574efeb76d5dffb3652554a811c28e515b65c88ff5a25dbe92cf0ec694c00` |
| `docs/joint-evidence-investigation.md` | `d6b1f96e5800fe964d8e71f4753a8afecb7f07a76bd079aa51cd195aef379d74` |
| `docs/joint-gap-results.md` | `62415e06bc0c05d6bda29e0121d559ebc595fae9f5aef7ca122ac38322a984c3` |
| `docs/joint-matrix-search-results.md` | `88551447f0467c495a3d18876365fa164ff7b64eb315a6bfcf83ce02ed8cb3df` |
| `docs/lfq-matrix-fix-plan.md` | `66a2c478a1652c3dd460932e0349074b7a719b408efb7048a65c8ebac442476b` |
| `docs/lfq-matrix-fix-results.md` | `f95db251905df549032717335d70274a1640d1f8c94dc87384f5f12c5fd03854` |
| `docs/matrix-seed-ensemble-20260913.md` | `0c6e117c170aece4155a81912837c3e60369a32b0b5d9da9dede3cd05edcdcb9` |
