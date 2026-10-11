# Experimental IM sequence CNN

`im_cnn_public_20260914.uimc` uses the same `RTCnn` network as retention-time
training, with 12 globals instead of 11: charge/4 is the added input. The
embedding, residual dilated convolutions and head are trained on measured
public IM-related labels, not on outputs distilled from an external predictor.
It is embedded in the search binary and enabled by default after RT shortlisting.
Use `--no-im` to disable it or `--im-model` to override its weights.

- Training: 123,565 observations / 100,000 sequences from the earlier public
  Meier TIMS CCS corpus preparation; fixed carbamidomethyl C, lengths 7–30,
  charges 2–4. CCS was converted to proportional inverse mobility.
- Validation: 12,284 observations / 10,000 disjoint sequences. Twenty epochs,
  seed 0, emb32/ch32, dilations 1/2/4, kernel 5, head64; best epoch 19.
- Independent public test: 45,929 observations / 37,101 disjoint sequences;
  MAE 0.01987 on the proportional IM scale, RMSE 0.02865, R² 0.96381.
- The previous tree model, trained on a larger public training split, has
  slightly lower error on that same test: MAE 0.01911, R² 0.96865. The sequence
  architecture is useful for consistent training/adaptation, but is not a
  demonstrated overall prediction improvement over that tree.
- A 1,000-observation order diagnostic increased MAE from 0.01897 to 0.02429
  when internal residue order was reversed while composition and termini were
  held constant. Rust/PyTorch maximum prediction discrepancy: 2.51e-7.

Search calibrates proportional IM per charge and uses the relative residual
only as a downstream classifier feature. Python is required only for training;
the artifact loads with pure Rust and includes architecture, output
normalization, label units, training settings and peptide lineage hashes.

See [training and usage](../docs/im-models.md). Prediction accuracy and
search identification depth are separate outcomes; the experimental search
comparison retains all three runs and three seeds for each arm.

A 27-search comparison found six CNN gains, two losses, and one tie, with higher medians in all three files; the predefined consistency rule was not met.
