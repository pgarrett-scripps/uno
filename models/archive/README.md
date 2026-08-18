# Retired RT models

Trained artifacts that nothing in the tree loads any more. They are kept rather
than deleted because each cost training data and hours to produce, cannot be
regenerated from what is in this repository, and the Zenodo deposit planned for
release (`paper/SUBMISSION_TODO.md` B1) is specified to carry the trained models.

Moved here 2026-07-31. None is referenced by any `.rs`, `.toml`, `.sh`, `.py`,
`.typ` or `.md` outside this directory.

| file | what it was |
|---|---|
| `rt_deeplc.gbdt` | the first pure-Rust RT model, a gbdt distilled from DeepLC over the 99-d order-aware feature set. Superseded by the CNN. |
| `rt_posmlp_irt.bin` | positional-encoding MLP trained on iRT-scale labels. |
| `rt_cnn_irt.urtc` | CNN trained on iRT-scale labels, before the move to Chronologer real-LC data. |
| `rt_cnn_ed.urtc` | CNN trained on the earlier development split. |
| `rt_cnn_ft.urtc` | CNN with the on-data readout fine-tuning that was purged in `a91c6a9`. The paper subsection describing it was deleted (decision D7). |

The shipped default is `../rt_cnn_chrono_global_v3_ch32.urtc`; `../rt_cnn_chrono_ch16.urtc`
is the speed dial. See `../rt_cnn_README.md`.
