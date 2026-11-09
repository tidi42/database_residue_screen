# Input and Output Contracts

Study contract version: **1.0.0**. Study-relative paths must resolve inside the
study directory; absolute paths and symlink escapes are rejected there. Explicit
search/inspection/acceptance commands take an external BLAST prefix and record
its identity. Completed analyses use the result contract `reference_screen_v2`
with `annotated_v1` event tables (file-format versions).

## Inputs

| Input | Required structure and meaning |
| --- | --- |
| Study JSON | Exact schema version; study/scope identity; approval status and record; paths; analysis, search, control, validation, metadata, figure and companion policies. Unknown configuration keys are rejected. |
| Query FASTA | Unique stable physical IDs, nonempty ACGT sequences. Separate standalone barcode, adapter/primer and composite representations. |
| Query metadata TSV | `query_id`, `group`, `condition`, `construct_class`, `family_id`, `source_id`, `provenance_status`, `architecture_status`. IDs must exactly match the FASTA population. Every query has a nonempty family and source, provenance `verified` or `synthetic`, and composites a `verified` or `synthetic` architecture. |
| Source TSV | `source_id`, `reference`, `version`, `evidence_file`, `evidence_sha256`, `verification`, `redistribution`. A source extract agreeing with a query is not automatically verified kit applicability. |
| Component TSV | `composite_id`, `component_id`, `component_kind`, `start`, `end`, `strand`. One-based inclusive query coordinates. Standalone sequence or its reverse complement must exactly occupy the stated interval. In the `authenticated_residues` profile every composite needs a complete nonoverlapping barcode/adapter decomposition. |
| Database JSON | Name/type, date/version, declared identity status, evidence path/hash, sequence/base counts. BLAST identity evidence must reconcile runtime and supplied metadata and state complete hash coverage. The software does not independently certify user-supplied verification declarations. |
| Subject metadata TSV | `raw_id`, `accession_version`, `length`, `topology`, `record_type`, `taxids`, `assembly_id`, `project_id`, `metadata_source`, `metadata_date`. Unknown values are not silently converted to biological facts. |
| Validation subject FASTA | Stable subject IDs; IUPAC ambiguity is allowed, unlike query FASTA. Exact matching treats ambiguity literally, not as a wildcard that guarantees a match. |
| Validation case TSV | `case_id`, `query_id`, `subject_id`, `truth`, `split`, `family_id`, `subject_source`, `error_class`. Independent exact-sequence checks reject contradictory strict truth labels. Tuning requires separate development/evaluation families and sources; sequence aliases are also checked. |
| Panel source manifest JSON | `kind: panel_sources`; per source `name` (local file name), `url`, `retrieved`, `sha256`, `bytes`, `purpose`, `redistribution`. URLs and hashes only; no document content. |
| Panel definition TSV | `query_id`, `group`, `condition`, `family`, `construct_class`, `length`, `sequence_sha256`, `kit_or_set`, `orientation`, `segments`, `transformation`, `rule`, `element`, `sources`, `verification`, `index_name`, `index_side`. No sequence column: `rule` and `element` name how and from which source element each query is derived; `sources` assigns manifest names to the rule's roles (`role=name;...`, optional `reference=` evidence). |

Actual contracts are shipped under
[src/database_residue_screen/schemas](../src/database_residue_screen/schemas).
TSV header order is flexible on input but field names must match exactly, with
no missing/extra/duplicate columns or malformed row widths. Numeric conversion
rejects non-finite values. Empty metadata tables retain their headers and remain
explicitly incomplete rather than fabricating records.

The optional study `profile` is `authenticated_residues` (default) or
`reference_screen`. In both, every query row with a sequence belongs to the one
study panel, and its group must be listed in `analysis.groups`; every query
counts in the panel results. `reference_screen` does not require
proprietary architecture authentication, complete component decomposition or a
new biological truth set. Valid optional truth fixtures still cannot contain
contradictions. Source-file/evidence integrity and the unrelated database, null,
taxonomy and approval checks remain required.

Extended record tables retain the base subject columns and add
`molecule`, `origin`, `targeted`, `organelle` and `context`; their schema is in
[eligibility.py](../src/database_residue_screen/eligibility.py). Molecule and
record class are distinct. Exact accession, length and source fields are retained.
Base subject tables (the study input contract) are accepted with explicit unknowns
for the added fields, not inferred DNA/topology. Annotation keeps BLAST and record taxids
separate and records mismatch/unresolved states.

## Panel Build

`fetch-sources MANIFEST --directory DIR` downloads every manifest URL and keeps a
file only if its SHA-256 equals the pinned value. An existing file with that hash
is accepted; a different existing file is never replaced and a changed download
is not stored. `build-panel DEFINITION --sources DIR --output OUT` reads the
`sources.json` beside the definition, verifies every referenced source, extracts
PDF text with `pdftotext -layout` (poppler) and derives each query by its rule:

| Rule | Derivation |
| --- | --- |
| `ont_native_barcode` | ONT native table, forward column; reverse column and Dorado `SQK-NBD114-96` entry must agree |
| `ont_native_context_k40` | ONT forward native flank statement around the barcode, within the Dorado kit flanks |
| `ont_native_context_k46` | Union of the ONT forward and reverse-complemented reverse flank statements; equals the Dorado kit flanks |
| `ont_rapid_context_90` | ONT Rapid flank statement (= Dorado `SQK-RBK114-96` flanks) around the barcode at kit position RBnn |
| `porechop_adapter` | Named Porechop definition; `verified` only if it occurs in the current ONT document |
| `illumina_index` | Bases in adapter of an index table; sample-sheet columns confirm the column layout |
| `illumina_index_adapter` | Index embedded in the document's `front[i7]rear` or `front[i5]rear` adapter formula |
| `illumina_formula_segment` | Constant front or rear of that formula |
| `documented_oligo` | Oligo after a locator (`anchor > anchor #n`), with a declared `U_to_T` or `A_tail_base` transformation |
| `mgi_filter_adapter` | Default `【adapterN】` value on one manual page (`page:parameter`) |
| `pacbio_barcode` | Official barcode file, else the declared mirror (`verified_secondary_source`); both must agree |
| `pacbio_barcode_hairpin` | Barcode + catalogue hairpin; every documented barcoded adapter must read barcode-hairpin-reverse complement |
| `catalogue_record` | FASTA catalogue record by accession (e.g. UniVec) |

Derived length, segments, verification status and SHA-256 must all equal the
definition; otherwise nothing is written. Outputs are `queries.fasta` (one line
per sequence), `queries.tsv` (the query metadata contract without condition
roles; `family_id` is the canonical ID of the carried barcode/index or of the
oligo itself) and `panel_build.json` (source and text-extract hashes, counts,
canonical parents and per-query cross-source notes, without sequences).
`queries.fasta` and `queries.tsv` are the query inputs of `controls` and
`prepare-panel`.

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
inputs. Stale source evidence or an altered study invalidates the preflight
outputs. Outputs describing a blocked study are still structurally valid.

## Design Readiness

`design_ready` is true when the preflight has no blocker; `preflight` then exits
with 0, otherwise with 2. Readiness checks declaration/file completeness and
tested mechanical relationships. It does not authenticate a manufacturer claim,
replace an author decision, approve a clinical/scientific study, measure search
sensitivity, or certify inferential results. A generated synthetic example is
never a research design approval.

## Search Output Contract

The search runner is implemented. The raw BLAST
fields are 18 standard tabular fields plus `qseq`, `sseq`, `btop`, `sstrand`, and
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
retrieval stability or control-parent matching. Full search validation is
performed by the runner and `analyse-search`; non-preflight stages are rejected
by `validate-output` rather than falsely certified.

The runner's actual sequence-level `status.tsv` uses `STATUS_SCHEMA` from
[search.py](../src/database_residue_screen/search.py), with per-query HSP/target
counts and completion states. Its immutable batch/attempt files are validated
by the runner and `analyse-search`; they are not the older standalone logical
status schema. Production mode requires an accepted database artifact; working
mode is explicit. Failed/incomplete attempt files cannot supply zeros.

## Downstream Release Structure

Completed analyses (`reference_screen_v2`) contain hashed analysis policy,
sequence-level events, logical query/parent/draw/seed mappings, search status,
real-event projections, eligibility, matched outcomes/contrasts (every condition
enters its group and the global union), length sensitivity, common-rule contrasts,
per-length testability and population-specific density/position/overlap tables.
Taxonomy adds exact-rank accession assignments and composition tables. Precision
reports are separate and reference the source analysis and the precision
specification hashes.

`validate-analysis` recomputes logical weights, eligibility and derived tables;
`--require-production` rejects working outputs and mismatched database provenance.
`figure-inputs` exports one production-origin source release. These commands
certify declared artifact consistency and arithmetic, not exhaustive retrieval,
biological specificity or scientific approval. See [analysis.md](analysis.md) for
the commands and analysis rules.