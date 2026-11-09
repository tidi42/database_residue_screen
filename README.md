# Database Residue Screen

A Python framework for reproducible screening of deposited nucleotide records
against one declared library-preparation reference panel.

**Version 1.0.0.** Release 1.0.0 is the version used for the Paper 3 analysis of NCBI core_nt (`studies/paper3`).
The package provides traceback-validated strict and tolerant matching, pinned
dinucleotide-shuffled controls, resumable BLAST execution with a production
database gate and a-priori length-aware testability, eligibility-aware
matched-parent aggregation, exact-version record metadata and dated taxonomy
reconciliation, positional and flank analyses, an exhaustive exact scan, a
panel-derived planted-sequence benchmark, conditional control-precision reports
and result-table validation. Every condition of the panel is analysed equally and
enters its group union and the global union. These implementations do not
establish biological specificity, complete null precision, or a coherent research
database.

## Installation

Python 3.11 or later is required; the test suite runs on Python 3.11, 3.12,
3.13, 3.14, and 3.15. Use a project-specific environment:

```bash
python -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m unittest discover -s tests -v
```

BLAST+ is optional for contract checks and required for database inspection,
searches and integration tests. `pdftotext` (poppler-utils) is required to build
panels from PDF sources and by the test suite. uShuffle is pinned to an exact
upstream commit because the Python-compatible 1.1.3 build was not available from PyPI. Installing
it requires a C compiler, Python development headers (e.g. `python3-dev`) and build dependencies; its source and binary provenance
are recorded with every generated control panel.
Its executable version is recorded explicitly; an installed version is not
assumed to match a previous analysis. The scientific Python environment and the
Word-document editing environment should remain separate.

## Quick Start

```bash
residue-screen example .local/example
residue-screen preflight .local/example/study.json --output .local/example/preflight
residue-screen validate-output .local/example/study.json .local/example/preflight
```

These commands generate synthetic fixtures and a checksummed preflight report.
They do not authenticate research constructs or demonstrate detector performance.
Output locations must be new; commands never overwrite an existing study or
report.

Use `python -m database_residue_screen` instead of `residue-screen` when needed.

Errors return exit code **1**; success returns **0**. `preflight` returns exit
code **2** while a design is blocked (`design_ready` is false) and still writes
the full report when `--output` is supplied. A well-structured report can validly
say that a study is not ready; `validate-output` verifies the report, not the
scientific acceptability of the study.

For the reference-screen publication scope, use `profile: "reference_screen"`
in the study JSON. Every query of the panel belongs to the reference population;
`prepare-panel` binds every metadata row to its canonical parent and matched
control draws. This profile does not require proprietary kit/junction
authentication or new component comparisons. Source identity, available
construction-basis evidence, database acceptance and honest interpretation remain
required. The default `authenticated_residues` profile additionally requires
complete component decompositions and an independent validation truth set.

## Running a Complete Study

```bash
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3
residue-screen study-status /data/runs/paper3 --verify --verify-code
```

`run-study` executes the `screen` pipeline: pinned panel sources, panel and
controls, database acceptance, BLAST search, analysis, metadata, taxonomy,
positional analysis, record cases, flank census, exact scan, target-limit arms,
benchmark, inference, validation, publication metrics, reports and figures. Every
stage records its request, code identity, tools, inputs and outputs, so a
completed run can be verified and reused stage by stage. See
[docs/study-workflow.md](docs/study-workflow.md) for the stages and the
configuration, [docs/study-engine.md](docs/study-engine.md) for the engine and
[docs/reports.md](docs/reports.md) for the publication tables.

## Study Panels from Pinned Sources

```bash
residue-screen fetch-sources studies/paper3/sources.json --directory .local/paper3-sources
residue-screen build-panel studies/paper3/panel_definition.tsv \
  --sources .local/paper3-sources --output .local/paper3-panel
```

A study panel is versioned as a definition (IDs, rules, source locators, segment
coordinates and SHA-256 values, no sequences) plus a source manifest (URLs and
SHA-256 values). `fetch-sources` downloads the vendor and catalogue documents;
`build-panel` derives every query from them and fails on any missing or changed
source or any hash disagreement. The resulting `queries.fasta`/`queries.tsv`
are the query inputs of `controls` and `prepare-panel`. Vendor documents and
the derived sequences stay outside the repository. See
[docs/contracts.md](docs/contracts.md#panel-build) for the rules.

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

`search` defaults to **production**, requiring `--acceptance` from
`accept-database`. Acceptance requires full content hashes, consistent volumes,
`blastdbcheck -full`, exact accession/length probes and explicit BLAST taxonomy
support files. It is not triggered automatically by a completed download.
`--mode working` is an explicit software-development mode; its outputs cannot
pass production figure-input export. See [docs/analysis.md](docs/analysis.md).

## Repository Layout

```text
src/database_residue_screen/  reusable contracts, analysis stages, validation and CLI
src/database_residue_screen/schemas/  versioned JSON Schema contracts
studies/paper3/              Paper 3 study: configuration, panel definition, source manifest (hashes, no sequences), figure scripts
tests/                       synthetic positive and failure-path tests
docs/                        contracts, analysis commands, design decisions and publication boundary
examples/                    guidance for generated synthetic studies
.github/workflows/           CI for package installation and tests
.local/                      ignored study inputs, reports and private evidence
```

See [docs/contracts.md](docs/contracts.md) for inputs and outputs,
[docs/analysis.md](docs/analysis.md) for commands and analysis rules,
[docs/design.md](docs/design.md) for design decisions and current limits, and
[docs/publication.md](docs/publication.md) for what belongs in GitHub.

## Scientific Safeguards

- Same-HSP identity and coverage are checked together; subject-level `qcovs`
  cannot establish full-query support.
- Physical records, condition-level canonical units and global sequence unions
  are separate quantities with explicit parent relationships. Every condition
  enters its group and the global union; condition-level estimates are kept.
- The reference-screen profile accepts a supplier-reference and industry-knowledge
  panel without a physical-authentication claim. Missing proprietary details are
  not invented and do not silently exclude queries.
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
- Draft approval, missing database evidence, unknown control version and an
  undeclared taxonomy snapshot block the design gate (`design_ready`).
- Result directories are validated by contract, hashes and recomputation; a file
  that merely parses is never accepted as validated framework output.

## License and Data

The framework code is MIT licensed; see [LICENSE](LICENSE). Third-party packages
retain their own licenses. This license does not grant redistribution rights for
manufacturer sequence documents, research panels, manuscripts, database volumes
or research results. Those files are intentionally not included.

## Analysis Commands

See [docs/analysis.md](docs/analysis.md) for commands, output contracts and the
distinction between synthetic evaluation, completed searches and scientific
acceptance.

## Length-Aware Search and Testability

The default `search` keeps the production E-value fixed (E=10 unless declared
otherwise). Runtime XML calibration computes retention and a-priori testability
per length, rule and actual database size; it is recorded in the search request,
and `analyse-search` refuses a search without it. Primary denominators contain
only `testable` matched parents. `censored_untestable` and `chance_limited` have
positive=-1 and separate attrition, never zero-filled denominators.

The declared uniform-iid chance bound uses both strands and all qualifying spans
(tolerant: a conservative union bound). Defaults: <=0.01 testable; (0.01,0.05]
chance-limited; >0.05 or incomplete score retention untestable. Configure
`--testable-chance-limit` and `--informativeness-limit` before analysis.
`--threshold-mode adaptive` is an optional diagnostic producing chance-saturated
output; it is OFF by default and never overrides testability exclusions.
`--search-space` permits fixed, VecScreen-like E-value scaling across releases;
actual database chance burden still changes and is always recomputed.
See `length_sensitivity.tsv` for the primary testable cohort, the full panel
inventory and the baseline exact-score-calibrated length panel, each with
explicitly labelled ungated fractions (testability gate removed);
`length_testability.tsv` for per-length classes and their reasons; and
`common_rule_contrasts.tsv` for strict/tolerant comparisons on parents eligible for
both rules. `precision-spec` takes the draw ladder (for example 5 10 15) and the
retention ladder (for example 500 2000) explicitly; `assess-precision` evaluates
those draw prefixes and also writes `conditional_contrasts.tsv` (grouped Hoeffding
bounds with Holm adjustment).

Use `--study study.json` to validate/hash declared search parameters and check
analysis position grids. The implemented failure policy is
`abort_run_resume_retries_failed_batches`: resume retries failed/incomplete
batches, and analysis requires a complete run. The schema value
`retry_then_block_matched_parent` is reported as unimplemented by preflight and
rejected when binding a search.

See [length-aware search details](docs/length-aware-search.md) and
[construct validity limits](docs/construct-validity.md). Vendor documents,
derived sequences, database files and run outputs are not part of this source tree.
A pinned source-install environment (Python 3.15.2) is supplied in `requirements.lock`.
