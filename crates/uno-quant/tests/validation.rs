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
