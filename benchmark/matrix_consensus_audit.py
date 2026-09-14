#!/usr/bin/env python3
"""Reproduce consensus failure mechanisms against the actual KOTH source.

This is a diagnostic, not a corrected algorithm. It compiles the consensus
implementation with identity-alignment fixtures, removing only serde/logging
dependencies. Assertions describe the current problematic behavior.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


STUBS = r'''
#![allow(dead_code, unused_variables)]
macro_rules! info { ($($arg:tt)*) => {}; }
mod log { pub(crate) use info; }
struct ConsensusConfig {
    min_member_combined_score: f64, min_group_size: usize,
    min_seed_combined_score: f64,
}
struct LfqConfig {
    consensus: ConsensusConfig, mz_ppm: f64,
    rt_window_pct: f64, im_tolerance: f64,
}
struct Feature { charge: u8, rt: f64 }
impl Feature {
    fn rt_apex(&self) -> f64 { self.rt }
    fn im_apex(&self) -> f64 { 0.0 }
}
struct ScoredFeature {
    feature: Feature, mz: f64, combined_score: f64,
    theoretical_pattern: Vec<f64>,
}
impl ScoredFeature { fn monoisotopic_mz(&self) -> f64 { self.mz } }
struct RunInput { name: String, features: Vec<ScoredFeature> }
impl RunInput { fn rt_range(&self) -> (f64, f64) { (0.0, 100.0) } }
struct Identity;
impl Identity {
    fn correct_mz(&self, x:f64, _:f64) -> f64 { x }
    fn warp_rt(&self, x:f64) -> f64 { x }
    fn correct_im(&self, x:f64, _:f64) -> f64 { x }
}
struct AlignmentResult {
    reference_idx: usize, alignments: std::collections::HashMap<String, Identity>,
}
fn config(min_group_size:usize, min_seed_combined_score:f64) -> LfqConfig {
    LfqConfig { consensus: ConsensusConfig {
        min_member_combined_score:0.5, min_group_size, min_seed_combined_score,
    }, mz_ppm:10.0, rt_window_pct:0.01, im_tolerance:0.05 }
}
fn feature(mass:f64, rt:f64, score:f64) -> ScoredFeature {
    ScoredFeature { feature:Feature{charge:2, rt},
        mz:mass/2.0+1.007_276_466_621, combined_score:score,
        theoretical_pattern:vec![0.6,0.3,0.1] }
}
fn interleaved(score:f64) -> (Vec<RunInput>, AlignmentResult) {
    let runs=vec![
        RunInput{name:"a".into(),features:vec![feature(1000.0,10.0,score),
                                                   feature(1000.001,40.0,0.9)]},
        RunInput{name:"b".into(),features:vec![feature(1000.002,10.0,0.9)]},
    ];
    let al=AlignmentResult{reference_idx:0,
        alignments:std::collections::HashMap::from([("b".into(),Identity)])};
    (runs,al)
}
fn cf(mass:f64, rt:f64, score:f64, slots:Vec<Option<u32>>) -> ConsensusFeature {
    ConsensusFeature { ref_mz:mass/2.0+1.007_276_466_621, ref_rt:rt, ref_im:0.0,
        charge:2, theoretical_pattern:vec![0.6,0.3,0.1],neutral_mass:mass,
        seed_combined_score:score,seed_run_idx:0,
        n_contributing_runs:slots.iter().filter(|x|x.is_some()).count(),
        per_run_feature:slots }
}
#[test]
fn seed_filter_before_merge_erases_valid_support() {
    let (runs,al)=interleaved(0.7);
    let actual=build_consensus(&runs,&al,&config(1,0.75));
    assert_eq!(actual.iter().find(|c|c.ref_rt==10.0).unwrap().n_contributing_runs,1);
    let delayed=build_consensus(&runs,&al,&config(1,0.0));
    let survivor=delayed.iter().find(|c|c.ref_rt==10.0).unwrap();
    assert_eq!(survivor.n_contributing_runs,2);
    assert!(survivor.seed_combined_score>=0.75);
}
#[test]
fn default_group_floor_before_merge_erases_replicated_feature() {
    let (runs,al)=interleaved(0.9);
    assert!(build_consensus(&runs,&al,&config(2,0.75)).is_empty());
    let delayed=build_consensus(&runs,&al,&config(1,0.75));
    assert_eq!(delayed.iter().filter(|c|c.n_contributing_runs>=2).count(),1);
}
#[test]
fn transitive_merge_exceeds_rt_and_mass_tolerances() {
    let groups=vec![cf(1000.0,10.0,0.99,vec![Some(0),None,None]),
        cf(1000.015,11.5,0.9,vec![None,Some(0),None]),
        cf(1000.030,13.0,0.9,vec![None,None,Some(0)])];
    let merged=merge_consensus(groups,&config(1,0.75),1.0);
    assert_eq!(merged.len(),1);
    assert_eq!(merged[0].n_contributing_runs,3);
    // Endpoint is 30 ppm and 3 minutes from the retained seed;
    // configured pairwise tolerances are 20 ppm and 2 minutes.
    assert!((1000.030-merged[0].neutral_mass)/merged[0].neutral_mass*1e6>20.0);
    assert!((13.0-merged[0].ref_rt)>2.0);
}
#[test]
fn owner_seed_score_does_not_select_best_feature_in_each_run() {
    let projected=vec![
        ProjectedFeature{run_idx:0,feature_idx:0,ref_mz:501.0,ref_rt:10.0,
            ref_im:0.0,charge:2,neutral_mass:1000.0,combined_score:0.99,
            theoretical_pattern:vec![]},
        ProjectedFeature{run_idx:1,feature_idx:10,ref_mz:501.0,ref_rt:10.0,
            ref_im:0.0,charge:2,neutral_mass:1000.0,combined_score:0.51,
            theoretical_pattern:vec![]},
        ProjectedFeature{run_idx:1,feature_idx:20,ref_mz:501.0,ref_rt:10.0,
            ref_im:0.0,charge:2,neutral_mass:1000.0,combined_score:0.90,
            theoretical_pattern:vec![]}];
    let groups=vec![emit_group(&projected,&[0,1],2),emit_group(&projected,&[2],2)];
    let merged=merge_consensus(groups,&config(1,0.75),1.0);
    assert_eq!(merged.len(),1);
    assert_eq!(merged[0].per_run_feature[1],Some(10)); // 0.51 wins over 0.90
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    original = args.source.read_text()
    # Keep production types and functions verbatim; stub only external types,
    # serde derivations and logging. Include build_consensus, not a reimplementation.
    body = original[original.index('/// One row in the multi-run consensus'):]
    body = body.split('#[cfg(test)]')[0]
    assert 'pub fn build_consensus(' in body and 'fn merge_consensus(' in body
    source = STUBS + '\nuse std::collections::HashMap;\nuse std::time::Instant;\n' + body
    with tempfile.TemporaryDirectory(prefix='koth-consensus-audit-') as tmp:
        rs = Path(tmp) / 'audit.rs'
        binary = Path(tmp) / 'audit'
        rs.write_text(source)
        subprocess.run(['rustc', '--edition=2021', '--test', str(rs), '-o', str(binary)],
                       check=True, capture_output=True, text=True)
        result = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
    report = {
        'source': str(args.source.resolve()),
        'source_sha256': hashlib.sha256(original.encode()).hexdigest(),
        'harness_sha256': hashlib.sha256(source.encode()).hexdigest(),
        'scope': 'Synthetic mechanism reproductions; not measured identification recovery.',
        'test_output': result.stdout,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
