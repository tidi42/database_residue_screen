# GitHub Publication Boundary

Target: https://github.com/tidi42/database_residue_screen

The repository is deliberately separate from the surrounding research archive.
Its `origin` must point to the target above. Do not add or commit unrelated files
from the parent repository.

## Include

- Installable Python package and versioned schemas.
- Focused tests, synthetic fixtures/generator and runnable CLI examples.
- Input/output contract, scientific-scope and reproducibility documentation.
- MIT license for the new framework code and dependency license attribution.
- CI that runs without private inputs, database volumes, credentials or research
  result files.

## Keep Local

- Manuscripts, supplements, backups, reviews with unpublished author information,
  and historical scientific figures.
- Full FASTA/query panels until their provenance and redistribution conditions
  have been verified; generated control panels and BLAST output.
- Database volumes, taxonomic archives, local execution logs, absolute machine
  paths, private metadata, credentials and environment directories.
- Local study/preflight/evidence bundles under ignored `.local/` or `runs/`.

The absence of research inputs from GitHub is deliberate. They can later be
published as a separately reviewed, versioned research-data release with legal
and scientific provenance. The MIT code license does not license third-party
sequences or manufacturer documents.

## Code-Only Publication

On 25 September 2026, the author explicitly authorized committing and pushing
the framework code separately from the research archive. This supersedes the
earlier code-publication hold, not the study's scientific acceptance gates.
The development release includes reusable source, schemas, synthetic tests,
packaging, CI and generic usage documentation. Study data, reports, reviews,
figures, manuscript files and local environments remain excluded.

Publishing code does not freeze or approve any scientific application, establish
biological validity, or authorize publishing third-party reference panels.
Research-data and scientific-result releases require separate approval and
their applicable provenance, verification and freeze checks.

## Release Checks

1. Run unit tests and the documented synthetic CLI workflow from an installed
   wheel, not only through a workspace `PYTHONPATH`.
2. Keep the code-release checks separate from scientific acceptance. Do not add
  private study inputs, output bundles or freeze evidence to the commit.
3. Inspect the exact proposed Git index and built package file list for data,
   secrets, unexpected large files, local paths and unrelated content.
4. Commit only the reviewed framework files and record the tested version.
5. Push normally to the authorised remote; do not force-push or rewrite history.

None of these steps authorises changing the underlying scientific results to
obtain a passing freeze status.