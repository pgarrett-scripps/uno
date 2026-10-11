# uno release helpers. `just --list` shows every recipe.

default:
    @just --list

# Write Cargo.toml's workspace version and today's date into CITATION.cff (run at release time)
cite-sync date=`date +%F`:
    #!/usr/bin/env bash
    set -euo pipefail
    v=$(sed -n '/^\[workspace.package\]/,/^\[/s/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
    sed -i -E "s/^version: .*/version: \"$v\"/; s/^date-released: .*/date-released: \"{{date}}\"/" CITATION.cff
    echo "CITATION.cff: version $v, date-released {{date}}"

# Fail if CITATION.cff's version differs from Cargo.toml's (CI runs the same check)
cite-check:
    #!/usr/bin/env bash
    set -euo pipefail
    crate=$(sed -n '/^\[workspace.package\]/,/^\[/s/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
    cff=$(sed -n 's/^version: "\{0,1\}\([^"]*\)"\{0,1\}$/\1/p' CITATION.cff)
    [ "$crate" = "$cff" ] || { echo "CITATION.cff version $cff != Cargo.toml $crate: run just cite-sync"; exit 1; }
    echo "CITATION.cff matches Cargo.toml ($crate)"

# Format, lint and test the workspace (fmt and clippy are not yet clean on main)
check:
    cargo fmt --all -- --check
    cargo clippy --workspace --all-targets --locked -- -D warnings
    cargo test --workspace --locked

# List what each publishable crate would ship
package-list:
    for p in uno-ms-lightgbm-sys uno-ms-core uno-ms; do echo "== $p"; cargo package --list --allow-dirty -p $p; done
