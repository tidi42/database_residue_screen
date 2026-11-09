# Study workflow: the `screen` pipeline

One configuration file runs a complete reference-panel screen of one BLAST
nucleotide database, from the pinned panel sources to the publication figures:

```bash
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3
residue-screen study-status /data/runs/paper3 --verify --verify-code
```

The engine (records, attempts, reuse, production and working modes) is described
in [study-engine.md](study-engine.md). This page lists the stages of the `screen`
pipeline ([pipeline/screen.py](../src/database_residue_screen/pipeline/screen.py)
and [pipeline/reporting.py](../src/database_residue_screen/pipeline/reporting.py)),
their configuration sections and their outputs. Every stage writes only below
its attempt's `output/` directory; the paths below are relative to it.

## Stages

| Stage | Depends on | Configuration | Purpose and main outputs |
|---|---|---|---|
| `sources` | | `[sources]` | Hash-checked snapshot of the pinned panel sources (`sources/files/`, `sources/sources.json`, `sources/collection.json`). Missing files are downloaded into `[sources] directory`. A changed file is accepted only with `allow_source_drift = true` and only if every derived query still matches its definition SHA-256. |
| `panel` | sources | `[panel]`, `[subset]` | Derives every query of the panel definition (`build/`), selects the study subset (`queries.fasta`, `queries.tsv`, `segments.tsv`, `selection.json`) and derives the run-time materials of later stages from the built panel and its sources: census flank elements, barcode catalogue and benchmark templates (`census_elements.json`, `census_catalogue.tsv`, `benchmark_templates.json`). `panel.json` records counts, hashes and accepted source drift. |
| `controls` | panel | `[controls]` | Dinucleotide-preserving shuffles (uShuffle, k = 2) of every canonical parent; the seed of draw *d* of parent *p* is derived from `sha256("ushuffle-k2-v1:{seed}:{p}:{d}")`, so each parent's draws depend only on its own sequence (`controls/`). |
| `design` | panel, controls | `[controls]`, `[analysis_policy]` | Physical and logical queries, distinct search sequences and the hashed analysis policy (`search_panel/`). |
| `database` | | `[database]` | Inventory with full file hashes, acceptance (`blastdbcheck -full`, accession and length probes, BLAST taxonomy support, download evidence) and volume audit (`inventory/`, `acceptance/`, `integrity/`, `probes.json`). |
| `search` | design, database | `[search]` | Resumable batch search of every distinct sequence with the per-length effective search space (`search/`). With `resume = true`, completed batches of an unfinished attempt with an identical request are reused. |
| `analysis` | design, search | | Strict and tolerant events, matched-parent outcomes and testability (`analysis/`). |
| `metadata` | analysis | `[ncbi]`, `[metadata]` | Exact-version NCBI ESummary records of every retained accession, with retrieval times (`metadata/`). |
| `events` | analysis, metadata | | Eligibility and record populations added to the analysis (`analysis/`). |
| `taxonomy` | events | `[taxonomy]` | Pinned taxdump snapshot and exact-rank taxonomy of every event (`snapshot/`, `analysis/`). This analysis is the final analysis read by all later stages. |
| `insdc` | taxonomy, metadata | `[ncbi]`, `[insdc]` | Dated INSDC GenBank headers (keywords, primers, source features) of the retained accessions (`insdc/`). |
| `positional` | taxonomy, metadata, insdc | `[positional]` | Terminal-window analysis over the declared windows and minimum record lengths (`positional_sensitivity.tsv`, `positional.json`). |
| `record_cases` | taxonomy, metadata, search, database, panel | `[record_cases]` | Deterministic selection of one record per declared category, with local flank retrieval and sequence verification (`record_cases/`). |
| `census` | taxonomy, design, panel, database | `[census]` | Flank classes of the retained strict windows of the declared families (`census/`). |
| `exact_scan` | design, database | `[exact_scan]` | Exhaustive exact occurrences, on both strands, of every search sequence of at least `minimum_length` nt in every database volume (`scan/`). |
| `scan_validation` | exact_scan, search | `[scan_validation]`, `[exact_scan]` | Concordance of the exact scan with the strict BLAST endpoint (`scan_validation/`). |
| `cap_review` | search, exact_scan | `[exact_scan]` | For every search that reached the target limit: weakest retained score against the minimum qualifying scores, with cap-free strict accession counts from the exact scan (`cap_review.tsv`, `cap_review.json`). |
| `target_limit` | design, database, search | `[target_limit]`, `[search]` | Every reference sequence searched once per configured target limit with otherwise identical settings; per-parent comparison of consecutive arms and of the production-limit arm with the production search (`query_sets/`, `search_<limit>/`, `comparison/`). |
| `benchmark` | design, panel, database | `[benchmark]` | One hashed planted-sequence fixture set with precomputed truth, calibrated and searched at every configured database size (`fixtures/`, `dimensions/`, `calibration/<size>/`, `runs/<size>/`, `benchmark.json`). See [length-aware-search.md](length-aware-search.md). |
| `inference` | design, taxonomy | `[inference]` | Precision specification bound to the search panel and conditional control-precision assessment at the declared draw prefixes (`specification/`, `precision/`). |
| `validation` | taxonomy | `[validation]` | Hashes, policy, logical weights, eligibility, derived tables and taxonomy of the final analysis, recomputed (`validation.json`); `require_production` requires production search evidence. |
| `metrics` | design, search, taxonomy | `[metrics]` | Central cohorts and denominators, chance-limit and testability sensitivity panels. See [reports.md](reports.md). |
| `reports` | metrics and the stages above | `[reports]` | Publication report tables. See [reports.md](reports.md). |
| `figure_sources` | metrics, reports and run outputs | | Every table read by the figure scripts. |
| `figures` | figure_sources | `[figures]` | Figures rendered by the R scripts in `[figures] scripts` (tool `Rscript`). |

Every acquisition and analysis stage from `panel` onward declares `[subset]`, so
attempts of a dry run on a subset are never reused for the complete study.

## Configuration

[studies/paper3/study.toml](../studies/paper3/study.toml) is the complete
configuration of the Paper 3 study. Machine-specific paths and the NCBI contact
address come from environment variables (listed at the top of the file); an
optional `NCBI_API_KEY` is recorded as `<redacted>`. The study configuration
holds rules and hashes only: queries, flank elements, barcode catalogues and
benchmark templates are derived at run time from the pinned sources.

[studies/paper3/study_dryrun.toml](../studies/paper3/study_dryrun.toml) is the
same configuration on one query of every panel family, two exact-scan volumes and
one benchmark shard per fixture panel at one database size. A test checks that
the two files differ only in these keys.

Tools (`[tools]`): `blastn`, `blastdbcmd`, `blastdbcheck`, `makeblastdb` and
`blast_formatter` from one BLAST+ installation, a C compiler `cc` (the benchmark
compiles the NCBI length-adjustment routine from the pinned source archive),
`pdftotext` and `Rscript`. Their identities are recorded with every stage that
uses them.

## Typical invocations

```bash
# dry run in working mode (uncommitted changes allowed and recorded)
residue-screen run-study studies/paper3/study_dryrun.toml --output runs/dryrun --working

# production run from a clean checkout of the release
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3

# stop after the searches, or rebuild only the figures
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3 --until search
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3 --only figures

# check every record, hash and the stage code against the checkout
residue-screen study-status /data/runs/paper3 --verify --verify-code
```

A run directory inside the repository must be git-ignored (for example `runs/`).

## Time and resources

On 60 threads the search stage of the Paper 3 panel (1,384 canonical parents,
15 control draws each) against core_nt 2026-09-12 takes about two days; the
target-limit arms, the exact scan and the benchmark add several hours each. The
database stage hashes every database file (about 1 TB for core_nt).
