# Input and Output Contracts

Study contract version: **1.0.0**. Study-relative paths must resolve inside the
study directory; absolute paths and symlink escapes are rejected there. Explicit
search/inspection/acceptance commands take an external BLAST prefix and record
its identity. New analysis bundles use `reference_screen_v1`, independently of
historical preflight and retained-output schemas.

## Inputs

| Input | Required structure and meaning |
| --- | --- |
| Study JSON | Exact schema version; study/scope identity; approval status and record; paths; analysis, search, control, validation, metadata, figure and companion policies. Unknown configuration keys are rejected. |
| Query FASTA | Unique stable physical IDs, nonempty ACGT sequences. Separate standalone barcode, adapter/primer and composite representations. |
| Query metadata TSV | `query_id`, `group`, `condition`, `construct_class`, `family_id`, `role`, `source_id`, `provenance_status`, `architecture_status`. IDs must exactly match the FASTA population. |
| Source TSV | `source_id`, `reference`, `version`, `evidence_file`, `evidence_sha256`, `verification`, `redistribution`. A source extract agreeing with a query is not automatically verified kit applicability. |
| Component TSV | `composite_id`, `component_id`, `component_kind`, `start`, `end`, `strand`. One-based inclusive query coordinates. Standalone sequence or its reverse complement must exactly occupy the stated interval. Primary composites need complete nonoverlapping barcode/adapter decomposition. |
| Database JSON | Name/type, date/version, declared identity status, evidence path/hash, sequence/base counts. BLAST identity evidence must reconcile runtime and supplied metadata and state complete hash coverage. The software does not independently certify user-supplied verification declarations. |
| Subject metadata TSV | `raw_id`, `accession_version`, `length`, `topology`, `record_type`, `taxids`, `assembly_id`, `project_id`, `metadata_source`, `metadata_date`. Unknown values are not silently converted to biological facts. |
| Validation subject FASTA | Stable subject IDs; IUPAC ambiguity is allowed, unlike query FASTA. Exact matching treats ambiguity literally, not as a wildcard that guarantees a match. |
| Validation case TSV | `case_id`, `query_id`, `subject_id`, `truth`, `split`, `family_id`, `subject_source`, `error_class`. Independent exact-sequence checks reject contradictory strict truth labels. Tuning requires separate development/evaluation families and sources; sequence aliases are also checked. |

Actual contracts are shipped under
[src/database_residue_screen/schemas](../src/database_residue_screen/schemas).
TSV header order is flexible on input but field names must match exactly, with
no missing/extra/duplicate columns or malformed row widths. Numeric conversion
rejects non-finite values. Empty metadata tables retain their headers and remain
explicitly incomplete rather than fabricating records.

The optional study `profile` is `authenticated_residues` (legacy default) or
`reference_screen`. The latter selects original primary/exploratory rows without
requiring proprietary architecture authentication, complete component
decomposition or a new biological truth set. Valid optional truth fixtures still
cannot contain contradictions. Source-file/evidence integrity and the unrelated
database, null, taxonomy and approval checks remain required.

New extended record tables retain the original subject columns and add
`molecule`, `origin`, `targeted`, `organelle` and `context`; their schema is in
[eligibility.py](../src/database_residue_screen/eligibility.py). Molecule and
record class are distinct. Exact accession, length and source fields are retained.
Legacy subject tables are accepted with explicit unknowns for missing fields,
not inferred DNA/topology. Annotation keeps original BLAST and record taxids
separate and records mismatch/unresolved states.

## Preflight Outputs

```text
preflight/
  manifest.json
  preflight.json
  query_inventory.tsv
  component_coverage.tsv
```

The manifest records the study hash, validator version/code checksum, and exact artifact roles, paths, formats,
row counts and SHA-256 values. The report contains input hashes, blocker/warning
codes, inventory counts and detailed component coverage. Query inventory extends
the metadata with `canonical_id` (full SHA-256 of the lexicographically smaller
sequence/reverse complement) and length. Component coverage distinguishes exact
sequence decomposition from verified physical architecture.

`validate-output` requires every declared artifact, rejects undeclared files,
checks hashes and row contracts, and recomputes the preflight against current
inputs. Stale source evidence or an altered study invalidates the bundle. A
bundle describing a blocked study is still structurally valid.

## Freeze Lock

The `freeze` command writes a new JSON lock only after all G0 contract checks
pass. It records scope, study and evidence hashes and explicitly sets
`production_search_ready` and `scientific_results_validated` to false.

G0 checks declaration/file completeness and tested mechanical relationships. It
does not authenticate a manufacturer claim, replace an author decision, approve
a clinical/scientific study, measure search sensitivity, or certify inferential
results. A generated synthetic example is never a research design approval.

## Search Output Contract

The search runner is implemented. The raw BLAST
fields are the legacy 18 fields plus `qseq`, `sseq`, `btop`, `sstrand`, and
`qcovhsp`, as exported by [matching.py](../src/database_residue_screen/matching.py).
All scoring/masking/database/query settings must be captured with the run, and
real and control targets must use the same scientific settings.

The older standalone `search_status` row schema defines `query_id`, `role`,
`parent_id`, `query_sha256`, `search_id`, `batch_id`, `status`, `exit_code`,
`retained_rows`, `parameters_sha256`, `database_sha256`, and `error`.
Successful zero requires exit code zero, no error and zero retained rows.
Failure and incomplete execution are distinct states. A control draw retains
its own identity and canonical parent even when its sequence repeats.

`validate-table search_status` checks individual rows only. It does **not**
certify the completeness of a run, query membership, successful execution,
retrieval stability or control-parent matching. Full search-bundle validation
is performed by the runner; non-preflight stages are rejected by
`validate-output` rather than falsely certified.

The runner's actual sequence-level `status.tsv` uses `STATUS_SCHEMA` from
[search.py](../src/database_residue_screen/search.py), with per-query HSP/target
counts and completion states. Its immutable batch/attempt files are validated
by the runner and `analyse-search`; they are not the older standalone logical
status schema. Production mode requires an accepted database artifact; working
mode is explicit. Failed/incomplete attempt files cannot supply zeros.

## Downstream Release Structure

Implemented `reference_screen_v1` bundles contain hashed analysis policy,
sequence-level events, logical query/parent/draw/seed mappings, search status,
real-event projections, eligibility, matched outcomes/contrasts and
population-specific density/position/overlap tables. Taxonomy adds exact-rank
accession assignments and composition tables. Precision reports are separate
and reference the source analysis and prespecified pilot hashes.

`validate-analysis` recomputes logical weights, eligibility and derived tables;
`--require-production` rejects working outputs and mismatched database provenance.
`figure-inputs` exports one production-origin source release without historical
fallback. These commands certify declared artifact consistency and arithmetic,
not exhaustive retrieval, biological specificity, scientific approval or a
completed final freeze. See [reanalysis.md](reanalysis.md) for the full workflow.