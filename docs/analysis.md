# Analysis Commands and Rules

The reference-screen profile compares deposited nucleotide records to one
declared supplier-reference/industry-knowledge panel. It does not authenticate
proprietary kit junctions, add component comparisons, evaluate raw reads or
prove biological contamination. Every condition of the panel is analysed
equally: it keeps its own condition-level estimates and enters its group union
and the global union. No preparation command starts a database search.

## Prepare Inputs Without Database Access

```bash
residue-screen controls queries.fasta --draws 15 --seed 42 --output controls
residue-screen prepare-panel queries.fasta queries.tsv controls --draws 15 --output panel
residue-screen precision-spec panel --draw-ladder 5 10 15 --retention-ladder 500 2000 --output precision-spec
```

`prepare-panel` binds every row of the query metadata, whose IDs must exactly
match the FASTA, to its canonical parent and control draws. Identical and
reverse-complement sequences are searched once; every physical alias, condition
and control draw keeps its own logical identity and seed. Control QC must belong
to the exact source FASTA. `--sample-parents N` selects a result-independent
stratified subset of canonical parents for small runs (round-robin over groups;
smallest canonical hash per construct, length and composition stratum).
`precision-spec` binds the panel hash to the explicitly declared draw and retention
ladders. All destinations are new and immutable.

## Analysis Policy

Every new panel contains a hashed `analysis_policy.json`. `--policy` accepts an
explicit alternative, validated by [policy.py](../src/database_residue_screen/policy.py).
It separates matching, query applicability, record population and interpretation.

- Strict: one full-query, exact, ungapped HSP. Secondary: >=95% identity and
  >=90% same-HSP span, gap-free. Subject-level `qcovs` cannot establish this.
- Every panel query remains in the search manifest and outcome tables.
  Below the configured 7-base word size, applicability is `below_blast_word_size`.
  Longer references default to `operational_retained_search`, not a sensitivity
  guarantee. Explicit higher length exclusions must be declared in the policy.
- Benchmark results do not filter publication outcomes. Planted-sequence
  benchmark metrics are software and detectability QA, not biological accuracy or
  exhaustive core_nt recall. Reverse-complement search is orientation handling.
- `all_records` includes identified records regardless of type or topology.
  `linear_genomic` requires explicit DNA, genomic class and linear topology;
  excludes known RNA/transcript, circular, synthetic/vector/engineered and
  targeted/amplicon records; unknown required fields remain unresolved.
- `biomol=genomic` alone is not DNA evidence. Unknown engineered/targeted origin
  is not proof of a natural record; residual uncertainty is retained as a limit.
- Position defaults are a <=500-bp window on records >=10 kb, with all
  100/500/1,000-bp and 1/10/50-kb sensitivities. Stored record ends are not
  necessarily chromosome ends. N2 is a geometric expectation, not trimming failure.

## Database Acceptance and Search

Download completion is not extraction or database acceptance. Safely extract one
verified archive release into a clean directory, record archive-to-extracted-file
provenance, and do not mix volumes or overwrite an existing installation. That
setup is a separate operation, not part of input preparation.

```bash
residue-screen inspect-database /path/to/new/core_nt --hash-all --output database-evidence
residue-screen accept-database /path/to/new/core_nt \
  --evidence database-evidence/database_identity.json --probes accession-probes.json \
  --taxonomy-support /path/to/blast-taxdb --download-evidence /path/to/download/evidence \
  --output accepted-database
residue-screen search panel/search_queries.fasta /path/to/new/core_nt \
  --evidence database-evidence/database_identity.json \
  --acceptance accepted-database/accepted.json --targets 2000 --threads 8 \
  --batch-size 16 --timeout-seconds 3600 --output search
```

The probe JSON has `version: "1.0.0"` and a `probes` list. Each entry requires
`accession_version`, `expected_length`, `allow_absent` and an evidence `basis`.
An allowed absent accession is recorded, not treated as a successful lookup;
at least one exact identity/length probe must pass. Wrong IDs always fail.

Acceptance verifies complete file hashes, metadata, volume counts/builds,
`blastdbcheck -full`, probes and executable identities. It archives failures and
publishes `accepted.json` only on success. For core_nt it also requires the
pinned downloader's completed checksum evidence. Local extracted hashes identify
the checked installation; they do not replace archive/extraction provenance.

BLAST's taxonomy support (`taxdb.bti` and `taxdb.btd`) is separate from the dated
taxdump used for research lineage summaries. Supply it explicitly; the installed
checker can fail metadata tests without it even when sequence checks pass.
`--blastn`, `--blastdbcmd` and `--blastdbcheck` select the intended tools.
Full acceptance can read the entire database and take substantial time.

`search` defaults to production and fails without matching acceptance. Only
software fixtures/development should use `--mode working`; those outputs cannot
be exported as production figure inputs. Requests capture scientific settings,
code, executable and database identity. Accepted files must remain immutable;
startup/batch guards check membership, size and nanosecond mtime rather than
rehashing hundreds of GB per batch. These are accidental-mutation guards, not
protection against an actor deliberately changing data and restoring metadata.

Search resumes only identical requests and checksums completed attempts. Runtime
allowance may change, but scientific settings, batching, inputs and runner/parser
identity may not. Failed/incomplete attempts never contribute zeros. Each batch
records elapsed time, output bytes and the process-lifetime child RSS high-water
mark (KiB on Linux), not a per-query memory estimate. Full resource and cap
adequacy still require representative real-database runs.

Dedicated searches of a declared sequence set, for example all references of
the panel searched in two arms that differ only in the target limit, start from
`python -m database_residue_screen.subset_search query-sets panel --output sets`,
which writes `references.fasta` and `controls.fasta` with hashes. Search each arm
with the same FASTA, `--batch-size`, `--threads` and search space; identical
inputs give identical batches. `subset_search assemble panel --sources ...`
builds an analysable search over a panel from completed searches with identical
scientific settings, as verified per-query projections without running BLAST.
`cap-review` compares the weakest retained subject of every cap-saturated query
with the minimum qualifying scores.

## Annotate and Recompute

```bash
residue-screen analyse-search panel search --output analysis
residue-screen fetch-metadata analysis --output record-metadata
residue-screen annotate-events analysis record-metadata --output annotated
residue-screen pin-taxonomy /path/to/ncbi-taxdump --source-url https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ \
  --retrieved-utc 2026-09-15T12:00:00Z --output taxonomy-snapshot
residue-screen annotate-taxonomy annotated taxonomy-snapshot --output reconciled
residue-screen assess-precision reconciled precision-spec/precision_specification.json --output precision
residue-screen validate-analysis reconciled --require-production
residue-screen figure-inputs reconciled --output figure-inputs
```

The taxonomy timestamp above is an example, not an actual acquisition record.
Use the real retrieval timestamp; `--release-date` is optional and must not be
inferred from the download date. No taxonomy download is started by these commands.

`fetch-metadata` archives exact-version ESummary responses. Missing metadata,
wrong versions and length conflicts remain explicit. DNA/RNA and record class
are separate. The extended subject table adds `molecule`, `origin`, `targeted`,
`organelle` and `context`; base subject tables (the study input contract) remain
readable but cannot acquire DNA evidence by inference. Additional explicit record
annotations need archived source evidence; titles do not silently assign topology
or origin.

`annotate-events` rebuilds outcomes, contrasts, eligibility and derived tables.
It does not copy broad contrasts and rename them genomic. Sequence-level HSPs
stay in `events.tsv`; `logical_queries.tsv` preserves every physical query and
every control draw including its seed. `event_memberships.tsv` maps searched
sequences to their real parent, group, condition and role without multiplying
control draw weights. `real_events.tsv` projects real events once per parent and group.

`outcomes.tsv` has separate match endpoint and record population, applicability,
completion, resolution, retained-accession count and target-limit state.
`positive=-1` means unresolved/noncompleted, never zero. Known positives with
unresolved additional breadth stay positive but are not fully resolved.
`contrasts.tsv` uses equally weighted canonical parents, separate condition/group
and global unions, and fully resolved matched subsets with explicit attrition.
Unresolved controls are not selectively removed from a parent's denominator.

The real-event source tables are `{population}_density.tsv`,
`{population}_positions.tsv` and `{population}_overlap.tsv`. Units are retained
HSPs/normalised accessions, not genomes or independent biological events.
Taxonomy adds `accession_taxonomy.tsv` and `taxonomic_composition.tsv`.

## Taxonomy and Precision Rules

Taxonomy snapshots archive hashes of `nodes.dmp`, `names.dmp`, `merged.dmp` and
`delnodes.dmp`. Merge chains are resolved, cycles/invalid parents rejected,
deleted IDs retained as unresolved, and ranks are exact. The default reports
domain and superkingdom separately without silently mapping one rank onto the
other. Multiple assignments need unanimous rank agreement. BLAST and later
exact-version record taxids remain distinct sources; conflicts are not overwritten.
Viruses/synthetic/vector lineages and unresolved categories remain explicit.
Percentages include unresolved buckets in the named accession denominator;
the global union is deduplicated across groups and HSPs.

The precision specification has no default ladders: `precision-spec` requires the
draw ladder (nested draw prefixes per parent, for example 5/10/15) and the
retention ladder (target limits, for example 500/2,000). The draw budget defaults
to the largest prefix. The declared targets are a 0.01 detection-rate half-width
and a breadth half-width of max(1 accession, 5% of the control mean); they are
declared criteria, not achieved precision or biological uncertainty.

`assess-precision` evaluates the declared draw prefixes of a completed analysis
bound to that specification. It reports conditional bounds with seed streams as
the independent units: the relative-entropy (Chernoff) form of the Hoeffding bound
(Hoeffding 1963, Theorem 1) on equally sized seed streams, and the additive
Hoeffding bound with weights grouped by seed otherwise. The relative-entropy form
uses the same assumptions (bounded outcomes, independent streams) and is much
tighter for means near zero. The report also lists non-zero zero-success bounds, per-parent diversity,
unsearched prefixes and attrition. Coverage assumes independence between distinct
seed streams and is simultaneous across the declared prefixes within a stratum
only. For the largest prefix it also writes `conditional_contrasts.tsv`: per group,
strict all-record detection and breadth reference-minus-generator differences
with these intervals and tail bounds, Holm-adjusted over all
group-by-measure slots.
Data-dependent unresolved-parent subsets get diagnostics, not certified coverage.
Target-capped breadth remains a retained endpoint even when Monte Carlo precision
is good. Null precision does not certify retrieval or biological specificity.

## Validation

`validate-analysis` checks result membership/hashes, search provenance,
policy, recomputed logical outcomes/weights and derived tables. With taxonomy it
also recomputes assignments and denominators against the pinned snapshot.
`--require-production` additionally checks production provenance. `figure-inputs`
uses that gate and exports one exact table release, not significance stars or
finished figures. Arithmetic validation is not final scientific acceptance.

Only completed analyses with the current result contract are accepted, and
`analyse-search` refuses searches without their runtime length calibration.
Rebuild panels and results instead of retrofitting hashes, and never merge
quarantined searches. The preflight `validate-output` command is separate from
result validation.

## Database Failure Policy

If volume builds differ or accession lookup returns an unrelated record,
quarantine searches made against that directory. Do not edit metadata to hide
the disagreement or assume only a harmless JSON file is stale. Preserve the
directory, obtain a complete coherent release in a separate location, verify
hashes and representative retrieval, then repeat the real/control searches and
dependent validation. Do not merge a replacement snapshot with old results.