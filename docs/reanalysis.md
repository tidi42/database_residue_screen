# Working Reanalysis

The reference-screen profile compares deposited nucleotide records to the
author's original supplier-reference/industry-knowledge panel. It does not
authenticate proprietary kit junctions, add component comparisons, evaluate
raw reads or prove biological contamination. Final freeze follows completed
results and review. No preparation command starts a database search.

## Prepare Inputs Without Database Access

For generic inputs:

```bash
residue-screen select-originals queries.fasta queries.tsv --output originals
residue-screen controls originals/queries.fasta --draws 100 --seed 42 --output controls
residue-screen prepare-panel originals/queries.fasta originals/queries.tsv controls \
  --draws 100 --output full-panel
residue-screen prepare-panel originals/queries.fasta originals/queries.tsv controls \
  --draws 100 --pilot-parents 16 --output pilot-panel
residue-screen pilot-spec pilot-panel --output pilot-spec
```

`select-originals` selects roles `primary` and `exploratory`, excludes added
`reference` rows and records IDs, exact sequence hashes, canonical keys and
source-file hashes. `prepare-panel` also excludes `reference` rows by default;
`--selection all` is an explicit change of population, not the publication default.
Control QC must belong to the exact source FASTA. Never substitute the hash of
an expanded panel for an original-only panel.

For the historical research layout, `prepare-reference-inputs` additionally
compares the selected imported sequences and IDs to the four original construct
FASTAs before generating controls:

```bash
residue-screen prepare-reference-inputs /path/to/research-root /path/to/imported/inputs \
  --draws 100 --seed 42 --output .local/prepared-reference-inputs
```

All destinations are new and immutable. The output includes originals, controls,
full/pilot panels, scope/source evidence and a fixed-stage pilot specification.
Existing controls, source panels and old analyses remain unchanged. No new
biological source authentication or final research approval is asserted.

## Analysis Policy

Every new panel contains a hashed `analysis_policy.json`. `--policy` accepts an
explicit alternative, validated by [policy.py](../src/database_residue_screen/policy.py).
It separates matching, query applicability, record population and interpretation.

- Strict: one full-query, exact, ungapped HSP. Secondary: >=95% identity and
  >=90% same-HSP span, gap-free. Subject-level `qcovs` cannot establish this.
- Every original query remains in the search manifest and outcome tables.
  Below the configured 7-base word size, applicability is `below_blast_word_size`.
  Longer references default to `operational_retained_search`, not a sensitivity
  guarantee. Explicit higher length exclusions must be declared in the policy.
- The benchmark's 24/32-bp fixture thresholds do not filter publication outcomes.
  Synthetic benchmark metrics remain software QA, not biological accuracy or a
  core_nt retrieval benchmark. Reverse-complement search is orientation handling.
- `all_records` includes identified records regardless of type or topology.
  `linear_genomic` requires explicit DNA, genomic class and linear topology;
  excludes known RNA/transcript, circular, synthetic/vector/engineered and
  targeted/amplicon records; unknown required fields remain unresolved.
- `biomol=genomic` alone is not DNA evidence. Unknown engineered/targeted origin
  is not proof of a natural record; residual uncertainty is retained as a limit.
- Position defaults are a <=500-bp window on records >=10 kb, with all
  100/500/1,000-bp and 1/10/50-kb sensitivities. Stored record ends are not
  necessarily chromosome ends. N2 is a geometric expectation, not trimming failure.

## Database Acceptance and Search (Later)

Download completion is not extraction or database acceptance. Preserve the old
database. Safely extract one verified archive release into a clean directory,
record archive-to-extracted-file provenance, and do not mix volumes. That setup
is a separately authorised operation, not part of input preparation.

```bash
residue-screen inspect-database /path/to/new/core_nt --hash-all --output database-evidence
residue-screen accept-database /path/to/new/core_nt \
  --evidence database-evidence/database_identity.json --probes accession-probes.json \
  --taxonomy-support /path/to/blast-taxdb --download-evidence /path/to/download/evidence \
  --output accepted-database
residue-screen search pilot-panel/search_queries.fasta /path/to/new/core_nt \
  --evidence database-evidence/database_identity.json \
  --acceptance accepted-database/accepted.json --targets 500 --threads 8 \
  --batch-size 16 --timeout-seconds 3600 --output pilot-search
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
adequacy still require representative real-database pilots.

## Annotate and Recompute (Later)

```bash
residue-screen analyse-search pilot-panel pilot-search --output analysis
residue-screen fetch-metadata analysis --output record-metadata
residue-screen annotate-events analysis record-metadata --output annotated
residue-screen pin-taxonomy /path/to/ncbi-taxdump --source-url https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ \
  --retrieved-utc 2026-09-15T12:00:00Z --output taxonomy-snapshot
residue-screen annotate-taxonomy annotated taxonomy-snapshot --output reconciled
residue-screen assess-precision reconciled pilot-spec/pilot.json --output precision
residue-screen validate-analysis reconciled --require-production
residue-screen figure-inputs reconciled --output figure-inputs
```

The taxonomy timestamp above is an example, not an actual acquisition record.
Use the real retrieval timestamp; `--release-date` is optional and must not be
inferred from the download date. No taxonomy download is started by these commands.

`fetch-metadata` archives exact-version ESummary responses. Missing metadata,
wrong versions and length conflicts remain explicit. DNA/RNA and record class
are separate. The extended subject table adds `molecule`, `origin`, `targeted`,
`organelle` and `context`; legacy tables remain readable but cannot acquire DNA
evidence by inference. Additional explicit record annotations need archived
source evidence; titles do not silently assign topology or origin.

`annotate-events` rebuilds outcomes, contrasts, eligibility and derived tables.
It does not copy broad contrasts and rename them genomic. Sequence-level HSPs
stay in `events.tsv`; `logical_queries.tsv` preserves physical originals and
every control draw including its seed. `event_memberships.tsv` maps searched
sequences to original parent/group/condition and role without multiplying control
draw weights. `real_events.tsv` projects real events once per parent and group.

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

The default pilot specification declares fixed 10/50/100 stages, maximum 100
draws per parent, retention limits 500/2,000/10,000, a 0.01 detection-rate
half-width and a breadth half-width target of max(1 accession, 5% of the control
mean). These are declared implementation defaults for review before the pilot,
not achieved precision or biological uncertainty. No optional stopping is enabled.

`assess-precision` evaluates nested draw prefixes from a completed panel bound
to that specification. It reports conditional Hoeffding bounds with weights
grouped by shared seed, non-zero zero-success bounds, per-parent diversity,
unsearched stages and attrition. Coverage assumes independence between distinct
seed streams and is simultaneous across declared stages within a stratum only.
Data-dependent unresolved-parent subsets get diagnostics, not certified coverage.
Target-capped breadth remains a retained endpoint even when Monte Carlo precision
is good. Null precision does not certify retrieval or biological specificity.

## Validation and Historical Compatibility

`validate-analysis` checks result membership/hashes, original search provenance,
policy, recomputed logical outcomes/weights and derived tables. With taxonomy it
also recomputes assignments and denominators against the pinned snapshot.
`--require-production` additionally checks production provenance. `figure-inputs`
uses that gate and exports one exact table release, not significance stars or
finished figures. Arithmetic validation is not final scientific acceptance.

Old panels/results lacking the new policy, seed or result contracts are preserved
as historical artifacts. Rebuild new panels and results; do not retrofit hashes
or merge old quarantined searches. `reanalyse-retained` and `context-evidence`
remain retrospective audit tools, not required new component/biological studies.
The preflight `validate-output` command is separate from result validation.

## Database Failure Policy

If volume builds differ or accession lookup returns an unrelated record,
quarantine searches made against that directory. Do not edit metadata to hide
the disagreement or assume only a harmless JSON file is stale. Preserve the
directory, obtain a complete coherent release in a separate location, verify
hashes and representative retrieval, then repeat the real/control searches and
dependent validation. Do not merge a replacement snapshot with old results.