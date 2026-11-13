# Publication metrics, reports and figures

The last four stages of the `screen` pipeline ([study-workflow.md](study-workflow.md))
turn the final analysis and the stage outputs into the tables and figures of a
publication. They search nothing and annotate nothing; every input is read from
a completed upstream attempt and checked against its recorded hashes. Values are
written at full precision; rounding happens only in text and figures. Every
output directory has a `manifest.json` with the SHA-256 of each file and of each
input. Code: [reporting/](../src/database_residue_screen/reporting/).

## Counting rules

- **Unit.** A canonical parent: one reference sequence (both orientations) with
  its matched control draws. Physical queries that share a sequence share a parent.
- **Endpoint-native cohort.** For one endpoint (strict or tolerant), record
  population (all records or linear genomic) and measure (detection or breadth):
  the parents whose reference and every control draw are a priori testable,
  completed and resolved.
- **Common cohort.** Parents in the native cohorts of both endpoints; used for
  strict-versus-tolerant comparisons.
- **Excluded states.** `censored_untestable`, `chance_limited`,
  `partially_score_censored`, incomplete and unresolved evidence count neither as
  0 nor as 1; they appear only as attrition.
- **Unions.** Every panel condition enters its platform group and the global union.

## `metrics` stage

| File | Content |
|---|---|
| `metrics.json` | Counts of the analysed panel, chance-limit pairs with their testable lengths, exact-score panel, zero-search accounting, global detection. |
| `cohorts.tsv` | Native and common cohorts per scope (global, group, display class, condition), endpoint, population and measure: eligible parents, reference positives, control rates. |
| `table1_source.tsv` | Detection rows of every scope for both cohorts (source of Table 1). |
| `display_classes.tsv` | Display classes of Table 1 (`[metrics] display_classes`; default: built-in Paper 3 classes) with members, parents and lengths. Labels are unique within a group. |
| `cohort_membership.tsv` | Cohort membership of every parent. |
| `attrition.tsv` | Per group, endpoint, population and measure: planned parents and the count at every exclusion step. |
| `testability.tsv` | A priori class per query length and endpoint under every chance-limit pair, with the panel at that length. |
| `chance_limits.tsv` | The same evidence re-gated with the baseline pair of the analysis policy and the `[metrics] chance_limit_variants`. |
| `sensitivity_panels.tsv` | Ungated panel (no testability gate) and exact-score-retainable length panel. |
| `conversion.tsv` | Strict-to-tolerant conversion in the common detection cohort against each parent's own controls. |
| `breadth_lower_bound.tsv` | Mean breadth of resolved parents and a lower bound that keeps parents with unresolved breadth at their resolved accessions. |
| `zero_searches.tsv`, `zero_search_units.tsv` | Per query length and search role: searches with hits, capped searches and zero results, separated into zeros forced by the score threshold (the exact score cannot reach E = 10) and ordinary zeros; the same states per counting unit. |

## `reports` stage

One directory per report (`reports.json` lists them):

| Directory | Content |
|---|---|
| `components/` | Component inventory of composite queries, within-set barcode distances and barcode-segment readout (segments from the panel definition). |
| `cap_review/` | Every production search that reached the target limit, joined to its parents with its review state. |
| `target_limit/` | Per parent, endpoint and population: completion, retained and unresolved accessions, detection and breadth in each target-limit arm and in the production search; summary of agreement. |
| `eligibility/` | Accession-level eligibility flow, taxonomy context and population of the record fields used by the eligibility rules. |
| `positional/` | Terminal-window events under all positional variants (checked against the positional stage) and record type by taxonomy. |
| `nesting/` | Per parent: strict events against distinct loci, and the conditions that share a locus (nested queries are co-positive). |
| `benchmark_misses/` | Tolerant truth-positive pairs not accepted in the `[reports] benchmark_arm`, at every database size, with the reason (no seed, seed without HSP, HSP not qualifying) and denominators. |
| `census/` | Flank-census class counts per family and window status. |
| `exact_scan_validation/` | Exact scan against the strict BLAST endpoint: concordance by query length and discordant classes. |
| `inference/` | Conditional reference-versus-control bounds per group and measure (grouped Hoeffding bounds with Holm adjustment) and control precision at the `[reports] draw_prefixes`. |
| `record_cases/` | Sequence-level checks of the selected record cases (window rebuilt from the retrieved record and compared base by base). |

## `figure_sources` and `figures` stages

`figure_sources` writes every table read by the R figure scripts (contract
`publication_figure_sources_v1`) and a manifest with the counts quoted in the
figures. `figures` runs `build_figures.sh` from `[figures] scripts`
([studies/paper3/figures](../studies/paper3/figures/)) with the configured
`Rscript`. The build verifies the figure sources, renders every figure as PDF and
PNG with its CSV data, and writes `figures_manifest.json`. PNG and CSV files are
byte-deterministic; PDF files embed a creation time, so the manifest records a
PDF content digest with the date blanked. Required R packages: ggplot2,
patchwork, dplyr, tidyr, readr, jsonlite, digest, ragg and scales.
