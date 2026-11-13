# Length-aware search and testability

Every search records an explicit a-priori calibration of each query length
against the actual database and the fixed scoring contract before any batch
runs. `analyse-search` uses that runtime calibration and refuses a search
without it. Classes depend only on query length, database dimensions, scoring
and the declared limits, never on observed reference or control output.

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

Scoring is fixed: blastn-short, +1 match, -3 mismatch, gap open 5, gap extend 2,
DUST off, soft masking off, both strands (`SCORING_CONTRACT` in
[length_aware.py](../src/database_residue_screen/length_aware.py); `search_settings`
adds only the production E-value, target limit and optional search space). XML
probes against a generated tiny database obtain K, lambda and H from the actual
executable and verify each local effective search space and emitted E-value. The
Python length adjustment follows NCBI's public-domain `BLAST_ComputeLengthAdjustment`
with alpha/lambda=1/H, beta=0, including the 20-iteration bound and integer
verification. Source:
[NCBI blast_stat.c](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_stat.c).
The source-archive adapter used by `calibrate-benchmark` is pinned to the exact
BLAST+ 2.17.0 source SHA256; tests cover direct/archive XML statistics and the
production-scale reference values. No source download/compiler is required for
ordinary searches.

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
The boolean fields `strict_testable` and `tolerant_fully_testable` record score
retention with a recovered exact witness; use `*_testability` for primary
inference. The equal `*_score_retainable` fields state this meaning explicitly.

## Chance-aware testability

Raising thresholds corrects score censoring but turns short-query zeros into
chance-saturated positives for references AND controls, so calibration must
establish informativeness as well as retention. The three states are defined
before observing outcomes, per length, rule, actual database size and scoring:

* `censored_untestable`: the minimum qualifying score cannot pass production E,
  no exact witness/seed is available, OR the model chance bound exceeds the
  declared informativeness ceiling (default 0.05).
* `chance_limited`: retainable and chance bound in (0.01,0.05].
* `testable`: retainable and chance bound <=0.01.

Both limits are configurable; require 0 < testable limit <= informativeness
ceiling. Equal limits give a binary classification. These are conservative
declared choices, not fitted results; declare them in the study settings before
the run.

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

At the size of core_nt (about 1.08e12 bases) the default model classifies strict
23 nt and tolerant 29 nt as chance-limited; strict 24 nt and tolerant 32 nt are
testable. These are model boundaries, not measured outcomes. Recompute each
release/size/scoring: the approximate exact boundary rises 1.7 nt per tenfold
increase in database size. Current runner scoring remains the fixed validated
blastn-short contract; other scoring requires recalibration and validation.

Optional `--threshold-mode adaptive` uses max(base E,1.1*E(minimum tolerant score))
with optional `--maximum-evalue`. This DIAGNOSTIC, OFF by default, produces
chance-saturated output at short lengths and never promotes testability.
`scripts/sanity_short_queries.py` also defaults to fixed E=10; opt into its raised
arm with `--diagnostic-raised-thresholds`.

`length_sensitivity.tsv` includes primary_testable_parents, full_panel_inventory
and baseline_exact_calibrated_length_panel (lengths whose exact score is retained
at the baseline E, default the production E, with a recovered exact witness).
Cohort inventory and actual eligible denominator are separate: full-panel and
calibrated-length inventory do not count unavailable outcomes as zeros. The
`ungated_*` columns repeat each cohort without the testability gate (completed
untestable searches count as zeros; incomplete searches never do) as an explicitly
labelled sensitivity. Do not substitute inventory for primary N.

Proven singleton dinucleotide null spaces are flagged by bounded exact enumeration
in `parents.tsv` and analysis policies. Repeated sampled draws are not a proof.
The default policy retains those parents and sampling weights, reports their count,
and does not imply a broadly adequate biological null.

## Study binding and failure policy

`--study` validates task, base word size, base E, masking, strand, both target limits,
threshold mode, optional E ceiling, both chance limits and fixed search space against the request, and hashes the study.
`analyse-search` checks its position grids. Without `--study`, the request is an
explicit unbound working specification and makes no claim to execute study.json.
The study is a design declaration, not scientific acceptance. Control
budgets/seed/algorithm, precision draw and retention ladders, stability tolerance,
max_mc_standard_error and validation database_size_policy remain declarative until
bound by their separate control/precision/benchmark artifacts; they are not silently
claimed to drive those stages. Equality of all other stage parameters is not
certified by search binding. Unequal real/control caps are rejected for a joint run.

Implemented failure policy: abort the run, preserve failed/incomplete attempts,
retry unfinished batches on identical-request resume, and require all queries
completed before analysis. There is no selective matched-parent failure blocking.
The schema value `retry_then_block_matched_parent` is reported by preflight as
unimplemented and rejected by executable binding. Exit codes: 0 success, 2 a
completed preflight with blocked design, 1 operational/validation error.

## Benchmark and subset entry points

```bash
residue-screen calibrate-benchmark ncbi-blast-2.17.0+-src.tar.gz /path/to/bin \
  run/db-evidence/database_identity.json panel/search_queries.fasta --output calibration
residue-screen benchmark prepare fixtures --panel panel --extra-lengths 19 20 21 22 \
  --templates panel/benchmark_templates.json
residue-screen benchmark run fixtures --calibration calibration --output runs/core_nt --workers 4
python -m database_residue_screen.subset_search query-sets panel --output sets
python -m database_residue_screen.subset_search assemble panel --sources run/a run/b --output run/assembled
```

The benchmark is generated once from the real parents of one panel: panel A
plants every parent exactly and reverse-complemented (internal and terminal) with
a no-seed negative; panel B holds boundary templates at the panel lengths with all
boundary variants; panel C plants the tolerant-boundary and combined
clip/substitution variants of every parent; panel D repeats the boundary templates
at the declared `--extra-lengths`, together with the declared templates (named panel
queries extended by repeating a declared tail; the study panel stage writes them).
Panels A/B and C/D use separate fixed seeds, so extra lengths never change the
other fixtures. `--shard-limit N` keeps the first N shards of each panel for a dry run. Fixtures and the independent
sequence truth of every query-subject pair are hashed before any search; `run`
verifies those hashes and the pinned binaries and searches every shard in a local,
a production-calibrated and a permissive diagnostic arm. The calibration must cover
the panel and extra lengths; database dimensions come from the calibration, so one
fixture set runs at several database sizes (the study `benchmark` stage calibrates
and runs each configured size).

Calibration lengths, real-parent count, group names and database dimensions
come from supplied evidence/panels. Subset assembly requires identical per-length
scientific settings on overlapping lengths; assembled requests merge the required
calibrations. Batch and final completion files use pending-file rename. Search
requests without runtime calibration are rejected by `analyse-search`; they never
produce silent zeros. Regenerate results into a new directory; never overwrite them.

## Testability reasons and sensitivity outputs

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

**Sensitivity outputs.** `length_sensitivity.tsv` adds
`ungated_positive/_parents/_rate`: the fraction obtained without the testability
gate, counting completed untestable searches as zeros, for the full panel
inventory, the exact-score-calibrated length panel and the testable cohort. They
are labelled sensitivities, never primary estimates. `common_rule_contrasts.tsv`
gives strict/tolerant detection on parents eligible for both rules, conversion of
strict-negative references and controls, and the conversion expected from each
strict-negative reference's own controls (parent clustering retained), for the
a-priori testable cohort (`primary_testable`) and, as a diagnostic only, ungated
(`ungated_diagnostic`).

**Conditional contrasts.** `assess-precision` also writes
`conditional_contrasts.tsv` for the largest declared draw prefix: per group, strict
all-record detection and breadth, D = reference mean - control generator mean, the
stage-protected relative-entropy Hoeffding interval over n seed streams and the tail bound
min(1, 2K exp(-n kl(control/B || reference/B))),
Holm adjustment over all group-by-measure slots (non-estimable slots retained as one),
and observed versus generator-expected reference totals. These are Monte Carlo
bounds for fixed sequences, not exchangeability, contamination or generalisation
tests.

**Event views.** With a detectability calibration,
`*_density.tsv`, `*_positions.tsv` and `*_overlap.tsv` apply the strict a-priori
applicability gate before aggregation; `taxonomic_composition.tsv` applies the
native strict/tolerant gate separately. Every panel condition contributes to these
group and global views. `events.tsv` and `real_events.tsv` retain all supplied
evidence for explicitly labelled sensitivities. Do not present an unfiltered
inventory as a primary testability-restricted view.
If any strict-negative reference has no strict-negative own control draws,
`expected_reference_conversions_under_own_controls` is `not_estimated`, with the
missing-reference count reported; a partial sum is not an expected total.
