# Database Residue Screen

A Python framework for reproducible screening of deposited nucleotide records
against declared library-preparation reference panels.

**Status: implemented and tested components; final scientific application pending.**
The package provides traceback-validated matching, pinned dinucleotide controls,
resumable BLAST execution with a production database gate, synthetic software QA,
original-panel selection, eligibility-aware matched-parent aggregation, dated
taxonomy reconciliation, conditional control-precision reports and result-table
validation. It retains historical audit and flank-check commands separately.
These implementations do not establish biological specificity, complete null
precision, or a coherent research database. Final freeze is deferred until the
application and its evidence are complete; working runs remain versioned.

## Installation

Python 3.11 or later is required. Use a project-specific environment:

```bash
python -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m unittest discover -s tests -v
```

BLAST+ is optional for contract checks and required for database inspection,
searches and integration tests. uShuffle is pinned to an exact upstream commit
because the Python-compatible 1.1.3 build was not available from PyPI. Installing
it requires a C compiler and build dependencies; its source and binary provenance
are recorded with every generated control panel.
Its executable version is recorded explicitly; an installed version is not
assumed to match a previous analysis. The scientific Python environment and the
Word-document editing environment should remain separate.

## Quick Start

```bash
residue-screen example .local/example
residue-screen preflight .local/example/study.json --output .local/example/preflight
residue-screen validate-output .local/example/study.json .local/example/preflight
residue-screen freeze .local/example/study.json --output .local/example/design.lock.json
```

These commands generate synthetic fixtures and a **synthetic G0-only lock**.
They do not authenticate research constructs or demonstrate detector performance.
The lock explicitly states that production search readiness and scientific
result validation are false. Output locations must be new; commands never
overwrite an existing study, report bundle, or lock.

Use `python -m database_residue_screen` instead of `residue-screen` when needed.

## Existing Analysis Archives

```bash
residue-screen import-legacy /path/to/core_nt_database_search --output .local/study
residue-screen preflight .local/study/study.json --output .local/study/preflight
residue-screen validate-output .local/study/study.json .local/study/preflight
```

The legacy adapter is explicit and optional. The core contracts do not require
the original directory layout. The adapter imports physical queries, source
references and mechanically recoverable components into a new study. It checks
local source agreement, mapped-decoy invariants, and the structure of 17 active
historical BLAST tables. Unauthenticated queries remain exploratory. Derived
standalone components are reference queries, not newly authenticated molecules.
Original FASTAs, BLAST files, historical results and manuscripts are untouched.

`preflight` returns exit code **2** while a design is blocked and still writes
the full report when `--output` is supplied. A well-structured report can validly
say that a study is not ready; `validate-output` verifies the report, not the
scientific acceptability of the study.

For the original-reference publication scope, use `profile: "reference_screen"`
in the study JSON and `select-originals` or `prepare-reference-inputs` to record
the selected inputs. Original `primary` and `exploratory` rows both belong to
that reference population; added `reference` fragments are excluded by default.
This profile does not require proprietary kit/junction authentication or new
component comparisons. Source identity, available construction-basis evidence,
database acceptance and honest interpretation remain required. The default
legacy `authenticated_residues` profile retains its stricter design gates.

## Database Evidence

```bash
residue-screen inspect-database /path/to/core_nt/core_nt \
  --study-root .local/study --output .local/study/evidence/database
```

Point the study's `inputs.database` to the resulting descriptor relative to the
study root. The inspector records runtime metadata, supplied JSON metadata,
volume presence and hashes for small files. Large files are unverified by default.
`--hash-all` requests a potentially expensive full content-hash pass; it is not
an internal BLAST database-integrity test. Metadata disagreements and incomplete
hash coverage block a verified identity. Do not overwrite conflicting evidence
to make the check pass.

`search` now defaults to **production**, requiring `--acceptance` from
`accept-database`. Acceptance requires full content hashes, consistent volumes,
`blastdbcheck -full`, exact accession/length probes and explicit BLAST taxonomy
support files. It is not triggered automatically by a completed download.
`--mode working` is an explicit software-development mode; its outputs cannot
pass production figure-input export. See [docs/reanalysis.md](docs/reanalysis.md).

## Repository Layout

```text
src/database_residue_screen/  reusable contracts, validation, adapters and CLI
src/database_residue_screen/schemas/  versioned JSON Schema contracts
tests/                       synthetic positive and failure-path tests
docs/                        contracts, design decisions and publication boundary
examples/                    guidance for generated synthetic studies
.github/workflows/           CI for package installation and tests
.local/                      ignored study inputs, reports and private evidence
```

See [docs/contracts.md](docs/contracts.md) for inputs and outputs,
[docs/design.md](docs/design.md) for pre-freeze decisions and current limits, and
[docs/publication.md](docs/publication.md) for what belongs in GitHub.

## Scientific Safeguards

- Same-HSP identity and coverage are checked together; subject-level `qcovs`
  cannot establish full-query support.
- Physical records, condition-level canonical units and global sequence unions
  are separate quantities with explicit parent relationships.
- The reference-screen profile accepts the original supplier-reference and
  industry-knowledge panel without a physical-authentication claim. Missing
  proprietary details are not invented and do not silently exclude originals.
- Physical component authentication is outside that profile. Existing component
  validation remains available for studies that explicitly retain that scope.
- Unknown topology, record type, taxonomy and project metadata remain explicit.
- Search failure is not a successful zero-hit query. The runner requires the
  declared per-query blocks and complete BLAST footer, validates every HSP and
  commits only completed immutable batch attempts. Failed attempts remain archived.
- Production searches require an accepted immutable database identity and the
  selected executable, not merely equal volume dates or an accessible prefix.
- Linear-genomic results require explicit DNA, genomic class and linear topology,
  with known RNA/circular/engineered/targeted exclusions and unresolved categories.
- Taxonomy is reconciled against pinned files at exact ranks; unknown categories
  stay in accession-level denominators and source disagreements remain explicit.
- Benchmark fixtures do not define publication query-length cutoffs. Short-query
  applicability and fixed-panel retained outcomes are declared separately.
- Draft approval, missing database evidence, unknown control version and absent
  taxonomy still block the relevant final gates. Final freeze follows results.
- Old result files and figure annotations are never accepted as new validated
  framework output simply because they can be parsed.

## License and Data

The framework code is MIT licensed; see [LICENSE](LICENSE). Third-party packages
retain their own licenses. This license does not grant redistribution rights for
manufacturer sequence documents, existing research panels, manuscripts, database
volumes or historical analysis data. Those files are intentionally not included.

## Working Reanalysis

See [docs/reanalysis.md](docs/reanalysis.md) for commands, output contracts and
the distinction between synthetic evaluation, retained-data descriptions,
completed exploratory searches, and final scientific acceptance.