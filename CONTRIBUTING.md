# Contributing to uno

Thanks for your interest in `uno`. It is a Rust MS1-only search engine: it
searches LC-MS peptide features, scores candidates with retention-time and
ion-mobility prediction, and estimates FDR with target-decoy competition.
Contributions of all kinds are welcome, from bug reports and documentation
fixes to new features.

## Reporting bugs and requesting features

Please use [GitHub issues](https://github.com/pgarrett-scripps/uno/issues).

- For bugs, include what you ran (command and config), what you expected, what
  happened, and enough detail to reproduce (OS, Rust version, input type:
  feature table, mzML, Bruker `.d` or Thermo `.raw`).
- For feature requests, describe the use case and the outcome you want.

Search existing issues first to avoid duplicates.

## Development setup

```sh
git clone https://github.com/pgarrett-scripps/uno.git
cd uno
cargo build --locked
cargo test --workspace --locked
```

uno requires Rust **1.97+** and builds on Linux x86_64. `rust-toolchain.toml`
pins the development compiler.

## Code style

Before submitting, run:

```sh
just check
```

This runs `cargo fmt --all -- --check`, `cargo clippy --workspace --all-targets
--locked -- -D warnings` and `cargo test --workspace --locked`. CI enforces all
three, and Clippy treats warnings as errors.

## Adding tests

- Put unit tests inline in the module they cover, inside a `#[cfg(test)]`
  module.
- Small test inputs live in `tests/data` and `tests/fixtures`.

New behavior should come with tests. Bug fixes should ideally include a test
that fails before the fix and passes after. Changes to scoring or FDR should
also be checked on a real acquisition, since unit tests cannot show a change
in identification counts or false discovery proportion.

## Pull request process

1. Branch from `main`.
2. Keep commits focused, one logical change per commit where practical.
3. Make sure `just check` passes locally.
4. Open the PR against `main` and make sure CI is green.
5. Reference any related issue in the PR description.

## License and Code of Conduct

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE), and that you will follow the
[Code of Conduct](CODE_OF_CONDUCT.md).
