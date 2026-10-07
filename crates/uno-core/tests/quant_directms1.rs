use std::fs;
use uno_core::quant::directms1::{
    analyze,
    input::{self, ReferenceOptions},
    uno_adapter, Config, Dataset, InputPeptide,
};
#[test]
fn malformed_matrix_returns_errors_without_panics() {
    let cfg = Config::default();
    let data = Dataset {
        samples: vec!["A".into()],
        n_s1: 3,
        rows: vec![],
    };
    assert!(analyze(&data, &cfg).is_err());
    let mut data = Dataset {
        samples: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        n_s1: 2,
        rows: vec![],
    };
    assert!(analyze(&data, &cfg)
        .unwrap_err()
        .to_string()
        .contains("six eligible"));
    data.rows = vec![InputPeptide {
        key: "x".into(),
        sequence: "PEPTIDE".into(),
        proteins: vec!["P".into()],
        values: vec![Some(-1.); 4],
    }];
    assert!(analyze(&data, &cfg)
        .unwrap_err()
        .to_string()
        .contains("positive"));
}
#[test]
fn reference_refuses_reusing_a_run() {
    let tmp = tempfile::tempdir().unwrap();
    let p = tmp.path().join("a_proteins_full.tsv");
    fs::write(&p, "dbname\nP\n").unwrap();
    let e = input::reference(
        &[p.clone(), p.clone()],
        &[p.clone(), p],
        &Config::default(),
        &ReferenceOptions::default(),
    )
    .unwrap_err();
    assert!(e.to_string().contains("multiple samples"));
}
#[test]
fn uno_evidence_retains_decoys_and_drops_zero_signal() {
    let tmp = tempfile::tempdir().unwrap();
    let headers="schema_version\trun\tfeature_id\tpeptide_id\tsequence\tdecoy\tproteins\tcharge\tion_mobility\tobserved_mass\trt_minutes\tintensity\tidentity_evidence\n";
    let mut manifest = "sample\tcondition\tevidence\tproteins\n".to_string();
    for i in 0..4 {
        fs::write(
            tmp.path().join(format!("p{i}.tsv")),
            "accession\tdecoy\tq_value\nP\t0\t0.001\nP\t1\t0.8\n",
        )
        .unwrap();
        let rows=format!("{headers}1\tr{i}\t0\t0\tPEPTIDE\t0\tP\t2\t\t900\t12\t100\tuno_ms1_rt_mass_assignment\n1\tr{i}\t1\t1\tEDITPEP\t1\tP\t2\t\t900\t12\t50\tuno_ms1_rt_mass_assignment\n1\tr{i}\t2\t2\tZEROROW\t0\tP\t2\t\t900\t12\t0\tuno_ms1_rt_mass_assignment\n");
        fs::write(tmp.path().join(format!("e{i}.tsv")), rows).unwrap();
        manifest += &format!("r{i}\tS{}\te{i}.tsv\tp{i}.tsv\n", if i < 2 { 1 } else { 2 });
    }
    let path = tmp.path().join("samples.tsv");
    fs::write(&path, manifest).unwrap();
    let (data, hashes) = uno_adapter::load(&path, &Config::default(), 0.01).unwrap();
    assert_eq!(data.rows.len(), 2);
    assert_eq!(data.n_s1, 2);
    assert_eq!(hashes.len(), 9);
    assert!(data.rows.iter().any(|r| r.proteins == vec!["DECOY_P"]));
    assert!(data
        .rows
        .iter()
        .all(|r| r.values.iter().all(Option::is_some)));
}
#[test]
fn tied_intensity_panel_is_explicitly_unsupported() {
    let rows = (0..12)
        .map(|i| InputPeptide {
            key: i.to_string(),
            sequence: format!("PEPTIDE{i}"),
            proteins: vec![if i < 6 { "P" } else { "DECOY_P" }.into()],
            values: vec![Some(100.); 4],
        })
        .collect();
    let data = Dataset {
        samples: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        n_s1: 2,
        rows,
    };
    assert!(analyze(&data, &Config::default())
        .unwrap_err()
        .to_string()
        .contains("quintile"));
}

#[test]
fn exclusive_feature_gate_is_symmetric_and_keeps_alternative_observations() {
    let tmp = tempfile::tempdir().unwrap();
    let header = "schema_version\trun\tfeature_id\tpeptide_id\tsequence\tdecoy\tproteins\tcharge\tion_mobility\tobserved_mass\trt_minutes\tintensity\tidentity_evidence\n";
    let mut design = "sample\tcondition\tevidence\tproteins\n".to_string();
    for i in 0..4 {
        fs::write(
            tmp.path().join(format!("p{i}.tsv")),
            "accession\tdecoy\tq_value\nP\t0\t0.001\n",
        )
        .unwrap();
        // The strongest feature is shared by target and decoy. Both must lose
        // that signal; the target has an independent lower-intensity feature.
        let mut text = format!("{header}1\tr\t0\t0\tPEPTIDE\t0\tP\t2\t\t900\t12\t100\tuno_ms1_rt_mass_assignment\n1\tr\t0\t1\tEDITPEP\t1\tP\t2\t\t900\t12\t100\tuno_ms1_rt_mass_assignment\n1\tr\t1\t0\tPEPTIDE\t0\tP\t2\t\t900\t12\t50\tuno_ms1_rt_mass_assignment\n");
        if i < 2 {
            text.push_str(
                "1\tr\t2\t2\tONESIDE\t0\tP\t2\t\t900\t12\t30\tuno_ms1_rt_mass_assignment\n",
            );
        }
        fs::write(tmp.path().join(format!("e{i}.tsv")), text).unwrap();
        design.push_str(&format!(
            "r{i}\tS{}\te{i}.tsv\tp{i}.tsv\n",
            if i < 2 { 1 } else { 2 }
        ));
    }
    let path = tmp.path().join("samples.tsv");
    fs::write(&path, design).unwrap();
    let (baseline, _) = uno_adapter::load(&path, &Config::default(), 0.01).unwrap();
    assert_eq!(baseline.rows.len(), 3);
    let policy = uno_adapter::AdmissionPolicy {
        exclusive_features: true,
        min_per_condition: 2,
        ..Default::default()
    };
    let (data, _, audit) =
        uno_adapter::load_with_policy(&path, &Config::default(), 0.01, &policy).unwrap();
    assert_eq!(data.rows.len(), 1);
    assert_eq!(data.rows[0].sequence, "PEPTIDE");
    assert_eq!(data.rows[0].values, vec![Some(50.); 4]);
    assert_eq!(audit.precursors_before_condition_filter, 2);
    assert_eq!(audit.precursors_after_condition_filter, 1);
    assert!(audit
        .runs
        .iter()
        .all(|r| r.ambiguous_features == 1 && r.excluded_ambiguous_assignments == 2));
    let impossible = uno_adapter::AdmissionPolicy {
        min_per_condition: 3,
        ..policy
    };
    assert!(uno_adapter::load_with_policy(&path, &Config::default(), 0.01, &impossible).is_err());
}

#[test]
fn native_missing_mobility_uses_pandas_nan_key_without_merging_zero() {
    let tmp = tempfile::tempdir().unwrap();
    let mut paths = Vec::new();
    for i in 0..4 {
        let base = tmp.path().join(format!("run{i}"));
        let base = base.to_str().unwrap();
        for suffix in ["_proteins.tsv", "_proteins_full.tsv"] {
            fs::write(format!("{base}{suffix}"), "dbname\nP\n").unwrap();
        }
        fs::write(
            format!("{base}_PFMs_ML.tsv"),
            "seqs\tqpreds\tplen\nPEPTIDE\t0\t7\n",
        )
        .unwrap();
        fs::write(format!("{base}_PFMs.tsv"),
            "sequence\tproteins\tcharge\tion_mobility\tIntensity\tRT\nPEPTIDE\tP\t2\t\t100\t1\nPEPTIDE\tP\t2\tNaN\t200\t1\nPEPTIDE\tP\t2\t0.0\t50\t1\n").unwrap();
        paths.push(std::path::PathBuf::from(format!(
            "{base}_proteins_full.tsv"
        )));
    }
    let (data, _) = input::reference(
        &paths[..2],
        &paths[2..],
        &Config::default(),
        &ReferenceOptions::default(),
    )
    .unwrap();
    assert_eq!(data.rows.len(), 2);
    let absent = data.rows.iter().find(|r| r.key == "PEPTIDE2nan").unwrap();
    assert_eq!(absent.values, vec![Some(200.); 4]);
    let zero = data.rows.iter().find(|r| r.key == "PEPTIDE20.0").unwrap();
    assert_eq!(zero.values, vec![Some(50.); 4]);
}

/// Four runs of one fixture: feature 0 is claimed by PEPTIDE (protein P) and a
/// decoy EDITPEP (decoy P); PEPTIDE has a second, closer feature 1; LOWCONF
/// is a low-scoring, out-of-ensemble assignment on feature 2.
fn option_fixture(tmp: &std::path::Path, with_scores: bool) -> std::path::PathBuf {
    let ev = "schema_version\trun\tfeature_id\tpeptide_id\tsequence\tdecoy\tproteins\tcharge\tion_mobility\tobserved_mass\trt_minutes\tintensity\tidentity_evidence\n\
1\tr\t0\t0\tPEPTIDE\t0\tP\t2\t\t900.000000\t12\t100\tuno_ms1_rt_mass_assignment\n\
1\tr\t0\t1\tEDITPEP\t1\tP\t2\t\t900.000000\t12\t100\tuno_ms1_rt_mass_assignment\n\
1\tr\t1\t0\tPEPTIDE\t0\tP\t2\t\t900.000000\t12\t50\tuno_ms1_rt_mass_assignment\n\
1\tr\t2\t2\tLOWCONF\t0\tP\t2\t\t900.000000\t12\t40\tuno_ms1_rt_mass_assignment\n";
    let pr = "accession\tdecoy\tbinomial_score\tq_value\nP\t0\t10.000\t0.001\nP\t1\t5.000\t0.5\n";
    let sc = "peptide\tdecoy\tproteins\tfeature_idx\tml\tdeviation\trank_fraction\tin_ensemble\tml_selected\tgroup_anchor_run\tgroup_anchor_observation\n\
PEPTIDE\t0\tP\t0\t5\t20\t0.25\t1\t1\t\t\n\
PEPTIDE\t0\tP\t1\t4\t5\t0.5\t1\t1\t\t\n\
EDITPEP\t1\tP\t0\t1\t30\t0.75\t0\t1\t\t\n\
LOWCONF\t0\tP\t2\t0.5\t1\t1\t0\t1\t\t\n";
    let ft = "massCalib\tintensityApex\tintensitySum\n900.0\t1e2\t1e3\n900.0\t5e1\t5e2\n900.0\t4e1\t4e2\n";
    fs::write(tmp.join("e.tsv"), ev).unwrap();
    fs::write(tmp.join("f.tsv"), ft).unwrap();
    fs::write(tmp.join("s.tsv"), sc).unwrap();
    let mut design = String::from("sample\tcondition\tevidence\tproteins\tfeatures");
    design.push_str(if with_scores {
        "\tpeptide_scores\n"
    } else {
        "\n"
    });
    for i in 0..4 {
        // Each sample needs its own evidence file.
        fs::copy(tmp.join("e.tsv"), tmp.join(format!("e{i}.tsv"))).unwrap();
        fs::write(tmp.join(format!("p{i}.tsv")), pr).unwrap();
        design.push_str(&format!(
            "r{i}\tS{}\te{i}.tsv\tp{i}.tsv\tf.tsv{}\n",
            1 + i / 2,
            if with_scores { "\ts.tsv" } else { "" }
        ));
    }
    let path = tmp.join("samples.tsv");
    fs::write(&path, design).unwrap();
    path
}

#[test]
fn quant_options_each_change_only_their_own_step() {
    use uno_adapter::{AdmissionPolicy, IntensityKind, PeptideConfidence, Reduction};
    let tmp = tempfile::tempdir().unwrap();
    let path = option_fixture(tmp.path(), true);
    let run = |policy: AdmissionPolicy| {
        let (data, _, audit) =
            uno_adapter::load_with_policy(&path, &Config::default(), 0.01, &policy).unwrap();
        let values: std::collections::BTreeMap<String, f64> = data
            .rows
            .iter()
            .map(|r| (r.key.clone(), r.values[0].unwrap()))
            .collect();
        (values, audit.runs[0].clone())
    };
    let (reference, audit) = run(AdmissionPolicy::default());
    assert_eq!(
        reference.into_iter().collect::<Vec<_>>(),
        vec![
            ("EDITPEP|2||1".into(), 100.),
            ("LOWCONF|2||0".into(), 40.),
            ("PEPTIDE|2||0".into(), 100.)
        ]
    );
    assert_eq!(audit.precursors_with_multiple_observations, 1);
    // The reference ignores extra manifest columns entirely.
    let bare = tempfile::tempdir().unwrap();
    let plain = option_fixture(bare.path(), false);
    assert_eq!(
        uno_adapter::load(&plain, &Config::default(), 0.01)
            .unwrap()
            .0
            .rows
            .len(),
        3
    );

    let (v, audit) = run(AdmissionPolicy {
        feature_competition: true,
        ..Default::default()
    });
    assert!(!v.contains_key("EDITPEP|2||1") && v["PEPTIDE|2||0"] == 100.);
    assert_eq!(audit.removed_by_competition, 1);

    let (v, audit) = run(AdmissionPolicy {
        peptide_confidence: PeptideConfidence::Ensemble,
        ..Default::default()
    });
    assert_eq!(v.keys().collect::<Vec<_>>(), vec!["PEPTIDE|2||0"]);
    assert_eq!(audit.removed_by_confidence, 2);

    // q: PEPTIDE 0, 0; EDITPEP and LOWCONF both 1/3 after the monotone pass.
    let (v, _) = run(AdmissionPolicy {
        peptide_confidence: PeptideConfidence::Qvalue,
        peptide_qvalue: 0.34,
        ..Default::default()
    });
    assert_eq!(v.len(), 3);
    let (v, _) = run(AdmissionPolicy {
        peptide_confidence: PeptideConfidence::Qvalue,
        peptide_qvalue: 0.33,
        ..Default::default()
    });
    assert_eq!(v.keys().collect::<Vec<_>>(), vec!["PEPTIDE|2||0"]);

    let (v, _) = run(AdmissionPolicy {
        intensity: IntensityKind::Integrated,
        ..Default::default()
    });
    assert_eq!((v["PEPTIDE|2||0"], v["LOWCONF|2||0"]), (1000., 400.));

    let (v, _) = run(AdmissionPolicy {
        reduction: Reduction::MinDeviation,
        ..Default::default()
    });
    assert_eq!(v["PEPTIDE|2||0"], 50.);
}

#[test]
fn quant_options_refuse_missing_inputs() {
    use uno_adapter::{AdmissionPolicy, PeptideConfidence};
    let tmp = tempfile::tempdir().unwrap();
    let path = option_fixture(tmp.path(), false);
    let policy = AdmissionPolicy {
        peptide_confidence: PeptideConfidence::Ensemble,
        ..Default::default()
    };
    assert!(uno_adapter::load_with_policy(&path, &Config::default(), 0.01, &policy).is_err());
}
