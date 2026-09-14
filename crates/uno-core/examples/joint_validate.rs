//! Verify the final quality-filter integration leaves default pilot input exact.
use anyhow::{ensure, Result};
use uno_core::{joint::restrict_observations, long_input::Bundle, MatchConfig};

fn main() -> Result<()> {
    let path = std::env::args().nth(1).expect("usage: joint_validate BUNDLE");
    let mut matrix = Bundle::open(std::path::Path::new(&path))?.load_joint(1)?;
    let original = matrix.clone();
    restrict_observations(&mut matrix, &MatchConfig::default());
    ensure!(matrix.features == original.features, "default quality filtering changed quantitative rows");
    let a = matrix.joint_evidence.as_ref().unwrap();
    let b = original.joint_evidence.as_ref().unwrap();
    ensure!(a.len() == b.len(), "changed group count");
    let mut observations = 0;
    for (a, b) in a.iter().zip(b) {
        ensure!(a.len() == b.len(), "changed group membership");
        for (a, b) in a.iter().zip(b) {
            ensure!(a.feature == b.feature && a.run_id == b.run_id && a.observation_id == b.observation_id
                && a.bootstrap_eligible == b.bootstrap_eligible, "changed original scoring evidence");
            observations += 1;
        }
    }
    println!("{}", serde_json::json!({"default_scoring_input_exact":true,"groups":a.len(),"observations":observations}));
    Ok(())
}
