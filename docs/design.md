# Design Decisions and Boundaries

## Reference-Screen Publication Profile

`profile: "reference_screen"` adopts one reference panel built from available
supplier references and industry knowledge. It does not require physical
kit/junction authentication, proprietary disclosure, newly added component
queries or raw-read validation. Every query of the panel belongs to its reference
population and every condition is analysed equally. Source-file hashes and
available construction-basis evidence remain required, without inventing
verified supplier information.

The implementation and commands are described in [analysis.md](analysis.md).
They include valid controls, versioned analysis policy, record eligibility, dated
taxonomy, matched outcomes and precision diagnostics. Production searches require
accepted database evidence. Scientific acceptance follows the results. The
stricter authenticated-residue profile remains the default study profile; it
additionally requires complete component decompositions and an independent
validation truth set.

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
  logical draw/seed preservation and precision reporting at explicitly declared
  draw prefixes. Adequacy requires actual comparable searches; generated controls
  alone do not establish it.
- New component comparisons are outside the reference-screen profile. Component
  support tools are separate optional analyses; whole-query shuffling is not a
  barcode-specific conditional null.
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
- Sensitivity ladders and numerical tolerances are declared design settings,
  not evidence of achieved stability. Author approval must confirm the population,
  scope and precision criteria before the production run.

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
against actual files/provenance before production. Record annotation and a few
successfully retrieved flanks do not substitute for a complete subject metadata
or biological adjudication dataset.

## Execution Order

1. Hash the reference panel; declare scope, counting units, applicability and
  metadata/taxonomy/precision policies.
2. Prepare matched controls and, if needed, a result-independent representative
  sample of parents. Validate the software on local fixtures without claiming
  biological accuracy.
3. Install and accept a coherent database separately, then execute comparable
  reference/control searches under one versioned specification. Do not merge
  incompatible results.
4. Reconcile completion, annotate records/taxonomy, assess conditional null
  precision and recompute declared endpoint populations with explicit missingness.
5. Validate one result release, produce figures and documents from it, complete
  remaining author requirements, then publish.

The code in this repository implements contracts and tests that support this
order. It does not claim completion of the search, matched null, contextual
evaluation or statistical stages of any study.

## Figure Policy

Figures are produced from one validated result release. Values and uncertainty
come from that release; unsupported contrasts remain unestimated. Workflow and
audit graphics belong in supplementary material.

No manuscript, estimate, biological validation claim, or figure annotation is
included in the public framework source. Public deposition and author/journal
declarations remain distinct requirements from software testing.