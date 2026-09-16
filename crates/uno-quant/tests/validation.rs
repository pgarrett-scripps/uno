use std::fs;
use uno_quant::{
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
