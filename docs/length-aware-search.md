# Length-aware search amendment (v19)

This method amendment requires reanalysis of retained outputs with an explicit a-priori calibration.
New searches are needed only for unavailable evidence or explicitly diagnostic arms. A completed
fixed-threshold search that predates runtime calibration is classified by `calibrate-search` (see
"Historical fixed-threshold searches" below) without rereading or repeating any search. It is not evidence that an earlier
run used these settings. Manuscript numbers cannot be regenerated without the
archived panel, controls, database and output bundles.

## Execution and database scope

Supply any local nucleotide BLAST database prefix (NCBI download or `makeblastdb`)
and its `inspect-database` evidence. No core_nt dimensions are used by `search`.
`blastdbcmd -info` supplies bases and sequence count; a fresh runtime inspection
must agree with the inventory. Production searches retain the existing database
acceptance gate. Protein databases and remote searches are outside this runner.
Physical volumes must be resolvable in the prefix directory; aliases spanning
external directories must be materialised as a local database before inspection.

```bash
residue-screen inspect-database /data/my_nucleotide_db --output run/db-evidence
residue-screen search panel/search_queries.fasta /data/my_nucleotide_db \
  --evidence run/db-evidence/database_identity.json --mode working \
  --study study.json --output run/search --threads 8
residue-screen analyse-search panel run/search --output run/analysis
```

## Statistical contract

Scoring is fixed: blastn-short, +1 match, -3 mismatch, gap open 5, extension 2,
DUST off, soft masking off, both strands. XML probes against a generated tiny
database obtain K, lambda and H from the actual executable and verify each local
effective search space and emitted E-value. The Python length adjustment follows
NCBI's public-domain `BLAST_ComputeLengthAdjustment` with alpha/lambda=1/H,
beta=0, including the 20-iteration bound and integer verification. Source:
[NCBI blast_stat.c](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_stat.c).
The older source-archive adapter remains pinned to the exact BLAST+ 2.17.0 source
SHA256; new tests cover direct/archive XML and the production-scale reference
values. No source download/compiler is required for ordinary searches.

For each length L and actual database (n bases, N sequences), with length
adjustment ell, space=(L-ell)*max(1,n-N*ell), and E(S)=K*space*exp(-lambda*S).
An explicit `--search-space` overrides that space and is labelled in settings;
it can be declared as a fixed, VecScreen-like search space for release-invariant
E-value thresholds. It cannot make actual chance expectations release-invariant.
There is no extra factor of two in the nucleotide context length.

Strict exact score is L. Tolerant ungapped spans range from ceil(0.9L) to L;
with span a, at most floor(a/20) substitutions are allowed. The minimum qualifying
score is min[a-4*floor(a/20)] over those integer spans. For 19/20/21/22 nt these
scores are 18/16/16/16, respectively. Full exact E-values alone would miss
partial score censoring of tolerant matching.

Fixed mode is the default: no per-length threshold raising. Each length is
calibrated against the actual database release, size and fixed scoring contract.
Word size is min(7,L), with BLAST's floor of four. Isolated exact witnesses must
be recovered. No score threshold guarantees every heuristic HSP representation.
The legacy boolean names `strict_testable` and `tolerant_fully_testable` are
retained as score-retention compatibility aliases; use `*_testability` for
primary inference. Explicit `*_score_retainable` fields remove this ambiguity.

## Protocol v2: chance-aware testability (20 October 2026)

Raising thresholds corrects score censoring but turns short-query zeros into
chance-saturated positives for references AND controls. This is a refinement of
the author's proposal: calibration must establish informativeness as well as
retention. The three states are defined before observing outcomes, per length,
rule, actual database size and scoring:

* `censored_untestable`: the minimum qualifying score cannot pass production E,
  no exact witness/seed is available, OR the model chance bound exceeds the
  declared informativeness ceiling (default 0.05).
* `chance_limited`: retainable and chance bound in (0.01,0.05].
* `testable`: retainable and chance bound <=0.01.

Both limits are configurable; require 0 < testable limit <= informativeness
ceiling. The narrow intermediate band makes the user's hard ceiling explicit;
it differs from the audit's suggested 0.05-to-1 band. Equal limits give a binary
classification. These are conservative protocol choices, not fitted results.
Archive their sensitivity and freeze them before recomputing publication tables.

Strict chance bound is 2*n/4^L. Tolerant bound sums, for each span a>=ceil(.9L),
2*n*(L-a+1)*4^-a*sum_{k=0..floor(a/20)} choose(a,k)*3^k.
This double-counts overlapping spans and is a conservative union bound UNDER
uniform iid assumptions, not a bound for arbitrary genome composition, biological
false positives, or a replacement for controls. `E(S)` counts chance HSPs at a
score; it is recorded separately from predicate-qualifying chance burden.

References and controls inherit the same class. Both excluded classes have
positive=-1 even with retained matches; retained_accessions remains descriptive.
They are excluded from detection AND breadth primary denominators and reported
separately in contrasts/precision. Raw evidence is preserved. Partial tolerant
score censoring at 19–25 nt (including ONT 24 nt) is recorded separately in the
calibration, although its primary class is `censored_untestable`.

At the study database size the default model classifies strict 23 nt and tolerant
29 nt as chance-limited; strict 24 nt and tolerant 32 nt are testable. These are
protocol-v2 model boundaries, not newly measured study outcomes. Recompute each
release/size/scoring: the approximate exact boundary rises 1.7 nt per tenfold
increase in database size. Current runner scoring remains the fixed validated
blastn-short contract; other scoring requires recalibration and validation.

Optional `--threshold-mode adaptive` uses max(base E,1.1*E(minimum tolerant score))
with optional `--maximum-evalue`. This DIAGNOSTIC, OFF by default, produces
chance-saturated output at short lengths and never promotes testability.
`scripts/sanity_short_queries.py` also defaults to fixed E=10; opt into its raised
arm with `--diagnostic-raised-thresholds`.

`length_sensitivity.tsv` includes primary_testable_parents, full_panel_inventory
and baseline_exact_calibrated_length_panel (baseline E=10). Cohort inventory and
actual eligible denominator are separate: full-panel and calibrated-length
inventory do not count unavailable outcomes as zeros. For the paper retain the
1,288 and 836 inventory sensitivities, and separately reanalyse historical v1
fractions with explicit labels. Do not substitute inventory for primary N.

Proven singleton dinucleotide null spaces are flagged by bounded exact enumeration
in `parents.tsv` and analysis policies. Repeated sampled draws are not a proof.
The default policy retains those parents and sampling weights, reports their count,
and does not imply a broadly adequate biological null.

## Study binding and failure policy

`--study` validates task, base word size, base E, masking, strand, both target limits,
threshold mode, optional E ceiling, both chance limits and fixed search space against the request, and hashes the study.
`analyse-search` checks its position grids. Without `--study`, the request is an
explicit unbound working specification and makes no claim to execute study.json.
The study is a G0 design declaration, not a final scientific freeze. Control
budgets/seed/algorithm, pilot draw and retention ladders, stability tolerance,
max_mc_standard_error and validation database_size_policy remain declarative until
bound by their separate control/pilot/benchmark artifacts; they are not silently
claimed to drive those stages. Equality of all other stage parameters is not
certified by search binding. Unequal real/control caps are rejected for a joint run.

Implemented failure policy: abort the run, preserve failed/incomplete attempts,
retry unfinished batches on identical-request resume, and require all queries
completed before analysis. There is no selective matched-parent failure blocking.
Legacy `retry_then_block_matched_parent` remains readable in old study files but
is rejected by executable binding. Exit codes: 0 success, 2 a completed preflight
with blocked design, 1 operational/validation error (including refused freeze).

## Current benchmark and subset entry points

```bash
residue-screen calibrate-benchmark ncbi-blast-2.17.0+-src.tar.gz /path/to/bin \
  run/db-evidence/database_identity.json panel/search_queries.fasta --output calibration
residue-screen current-benchmark prepare benchmark --panel panel --calibration calibration
residue-screen current-benchmark run benchmark --workers 4
python -m database_residue_screen.subset_search --help
```

Calibration lengths, original-parent count, group names and database dimensions
come from supplied evidence/panels. Subset execution requires identical per-length
scientific settings on overlapping lengths; assembled requests merge the required
calibrations. Batch and final completion files use pending-file rename.
Uncalibrated historical search requests are rejected by `analyse-search` unless a
hash-bound migration from `calibrate-search` is supplied; they never produce silent zeros.
Historical v1/v2 table files have fewer columns than this amendment; preserve
original bundles and regenerate into a new directory, never overwrite them.

## v20 review additions

**Reasons for each class.** `censored_untestable` covers two different situations:
the weakest qualifying alignment cannot pass the production threshold (score
censoring, or no runtime exact witness), and a retainable alignment whose chance
bound exceeds the informativeness ceiling (chance saturation, e.g. a 19-nt exact
query at core_nt size). Every calibration row therefore records
`strict_testability_reason` and `tolerant_testability_reason`
(`score_not_retainable`, `no_exact_witness`,
`chance_bound_above_informativeness_limit`, joined by `+`; or
`chance_bound_between_limits` / `chance_bound_within_testable_limit`). Policy
validation recomputes classes, reasons, retention flags, chance bounds and score
E-values from each row's own inputs, so hand-edited calibrations are rejected.
`length_testability.tsv` lists every length and rule with these inputs.

**Historical fixed-threshold searches.** Protocol-v2 classes depend only on query
lengths, database dimensions, the fixed scoring contract and the production E-value.
For a completed search without runtime calibration that used the validated
blastn-short contract:

```bash
residue-screen calibrate-search run/search --queries panel/search_queries.fasta --output run/migration
residue-screen analyse-search panel run/search --output run/analysis \
  --calibration-migration run/migration/calibration_migration.json
```

`calibrate-search` refuses searches that already carry a calibration, a different
scoring contract, changed database evidence or a different BLAST+ version, and it
reads no HSP. The migration is bound by hash to the request, completion, queries and
database evidence; `validate-analysis` re-verifies it. The historical request and
outputs are never edited.

**Sensitivity outputs.** `length_sensitivity.tsv` adds
`historical_zero_filled_positive/_parents/_rate`: the protocol-v1 fraction obtained
by counting completed untestable searches as zeros, for the full inventory and the
exact-score-calibrated length inventory (historically 1,288 and 836 parents). They
are labelled sensitivities, never primary estimates. `common_rule_contrasts.tsv`
gives strict/tolerant detection on parents eligible for both rules, conversion of
strict-negative references and controls, and the conversion expected from each
strict-negative reference's own controls (parent clustering retained), for the
protocol-v2 cohort and, as a diagnostic only, under historical zero-filled semantics.

**Post hoc conditional contrasts.** `assess-precision` also writes
`conditional_contrasts.tsv` for the final precision stage: per group, strict
all-record detection and breadth, D = reference mean - control generator mean, the
stage-protected Hoeffding interval and tail bound min(1, 2K exp(-2D^2/(B^2 sum w^2))),
Holm adjustment over all group-by-measure slots (non-estimable slots retained as one),
and observed versus generator-expected reference totals. Applied to the printed
counts of the historical table it reproduces the published detection bounds. These
are exploratory Monte Carlo bounds for fixed sequences, not exchangeability,
contamination or generalisation tests.

**Primary event views (v20 consolidation).** With a detectability calibration,
`*_density.tsv`, `*_positions.tsv` and `*_overlap.tsv` apply the strict a-priori
applicability gate before aggregation; `taxonomic_composition.tsv` applies the
native strict/tolerant gate separately. `events.tsv` and `real_events.tsv` retain
all supplied evidence for explicitly labelled historical sensitivities. Do not
present an unfiltered inventory as a primary testability-restricted view.
If any strict-negative reference has no strict-negative own control draws,
`expected_reference_conversions_under_own_controls` is `not_estimated`, with the
missing-reference count reported; a partial sum is not an expected total.
