//! Phase 7: quant rollup. koth already did alignment + XIC + LFQ, so uno only
//! assigns identities onto its per-sample intensities and rolls them up to
//! peptides and proteins. No requantification.

use crate::database::PeptideDb;
use crate::input::FeatureMatrix;
use crate::score::Identification;
use std::collections::HashMap;

/// A quant row: an entity (peptide or protein) and its per-sample intensities.
pub struct QuantRow {
    pub key: String,
    pub decoy: bool,
    pub q_value: f64,
    pub intensities: Vec<f64>,
}

/// Roll identified features up to peptides: one row per (winning) feature ID,
/// carrying that feature's koth intensities. Only target IDs at or below
/// `max_q` are emitted.
pub fn peptide_quant(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    ids: &[Identification],
    max_q: f64,
) -> Vec<QuantRow> {
    ids.iter()
        .filter(|id| !id.decoy && id.q_value <= max_q)
        .map(|id| QuantRow {
            key: db.peptides[id.peptide_idx as usize].sequence.to_string(),
            decoy: id.decoy,
            q_value: id.q_value,
            intensities: matrix.features[id.feature_idx as usize].intensities.clone(),
        })
        .collect()
}

/// Roll identified features up to proteins: per-sample intensity is the sum of
/// the protein's identified peptide-features. `protein_q` maps accession ->
/// protein q-value; only target proteins at or below `max_q` are emitted.
pub fn protein_quant(
    matrix: &FeatureMatrix,
    db: &PeptideDb,
    ids: &[Identification],
    protein_q: &HashMap<String, f64>,
    max_q: f64,
) -> Vec<QuantRow> {
    let n = matrix.n_samples();
    let mut acc: HashMap<&str, Vec<f64>> = HashMap::new();
    for id in ids {
        if id.decoy {
            continue;
        }
        let f = &matrix.features[id.feature_idx as usize];
        for prot in &db.peptides[id.peptide_idx as usize].proteins {
            let v = acc.entry(prot.as_ref()).or_insert_with(|| vec![0.0; n]);
            for (s, &x) in f.intensities.iter().enumerate() {
                v[s] += x;
            }
        }
    }
    acc.into_iter()
        .filter_map(|(acc_name, intensities)| {
            let q = *protein_q.get(acc_name)?;
            if q <= max_q {
                Some(QuantRow {
                    key: acc_name.to_string(),
                    decoy: false,
                    q_value: q,
                    intensities,
                })
            } else {
                None
            }
        })
        .collect()
}
