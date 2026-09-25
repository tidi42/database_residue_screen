# Pre-freeze Decisions and Boundaries

## Reference-Screen Publication Profile

`profile: "reference_screen"` adopts the original author reference panels built
from available supplier references and industry knowledge. It does not require
physical kit/junction authentication, proprietary disclosure, newly added
component queries or raw-read validation. Both original `primary` and
`exploratory` rows form its reference population; added `reference` fragments
are excluded. Source-file hashes and available construction-basis evidence remain
required, without inventing verified supplier information.

The implementation and current commands are described in
[reanalysis.md](reanalysis.md). They include original-input reconciliation,
valid controls, versioned analysis policy, record eligibility, dated taxonomy,
matched outcomes and precision diagnostics. Production searches require accepted
database evidence. Final scientific acceptance and freeze remain after results.
The stricter authenticated-residue profile remains the default for older study
contracts, so its historical requirements are not silently weakened.

## Implemented Decisions

- Strict full-query ungapped HSP support is the proposed primary endpoint.
  Identity >=95% and same-HSP query span >=90% with zero gaps is the named
  secondary endpoint. These thresholds are explicit defaults, not optimised
  detector-performance claims.
- Canonical-parent units are counted within conditions, with a separate global
  canonical union and equal parent weights. Physical aliases remain traceable.
- Finite-panel descriptive inference is the only supported design mode in this
  version. Platform-wide hypothesis tests and automatic significance stars are
  not enabled. Broader inference needs its own reviewed contract and tests.
- Pinned whole-query dinucleotide shuffling is implemented with invariant QC,
  logical draw/seed preservation and fixed-stage precision reporting. Adequacy
  requires actual comparable searches; generated controls alone do not establish it.
- New component comparisons are outside the reference-screen profile. Existing
  component tools are separate historical/optional analyses; whole-query shuffling
  is not a barcode-specific conditional null.
- Broad record matching is the primary reference-screen outcome. Linear-genomic
  DNA positions are a narrower endpoint with explicit molecule, topology, known
  exclusions and unresolved states, not an automatic replacement population.
- Automatic event collapse is disabled until validated. Descriptive accession
  overlap does not imply sequencing-platform origin, independent projects or an
  independently defined population-level overlap test.
- Metadata rules are explicit: dated taxonomy, ambiguous multiple taxids,
  disabled title fallback, and missing project/assembly linkage reported.
- Synthetic tests validate software behaviour, not biological sensitivity or
  specificity. They do not impose implicit length exclusions on the publication
  panel or require a new wet-laboratory/raw-read validation study.
- Proposed sensitivity ladders and numerical tolerances are draft design
  settings, not evidence of achieved stability. Author approval must confirm
  the population, scope and precision criteria before a research freeze.

## What Software Cannot Resolve Alone

When physical kit/junction authenticity is claimed, it requires source evidence
and scientific review. This claim is outside the reference-screen publication
scope; existing verification labels are retained as provenance, not eligibility
gates for that profile. The schema cannot establish physical truth from a label.
Likewise, an approval record proves that a file was provided, not that its author
has actually endorsed a scientific claim. Do not set verified/approved merely to
make the gate green.

Database accessibility and header counts do not certify file identity. Conflicting
JSON metadata, missing content hashes or unrecorded updates must be reconciled
against actual files/provenance before freezing. Record annotation and a few
successfully retrieved flanks do not substitute for a complete subject metadata
or biological adjudication dataset.

## Execution Order

1. Reconcile and hash original references; declare scope, counting units,
  applicability and metadata/taxonomy/precision policies.
2. Prepare matched controls and a result-independent representative pilot panel.
  Validate the software on local fixtures without claiming biological accuracy.
3. Install and accept a coherent database separately, then run retention/resource
  pilots and execute comparable reference/control searches under a versioned
  supported specification. Do not merge incompatible results.
4. Reconcile completion, annotate records/taxonomy, assess conditional null
  precision and recompute declared endpoint populations with explicit missingness.
5. Validate one result release, restore supported figures and documents, complete
  remaining author requirements, then perform final freeze and publication.

The pre-freeze code in this repository implements contracts and tests that
support this order. It does not claim completion of the later search, matched
null, contextual evaluation or statistical stages.

## Figure Policy

The intended result release can restore the earlier manuscript's scientific
panel designs. Values and uncertainty must be recomputed; unsupported contrasts
must remain unestimated. Workflow and audit graphics can move to supplementary
or revision material. Original images may be retained as explicitly historical
material, never as evidence that reanalysis is complete.

No manuscript, old estimate, biological validation claim, or figure annotation
is bundled with the public framework source. Public deposition and author/journal
declarations remain distinct requirements from software testing.