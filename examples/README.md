# Synthetic Examples

Generate a new self-contained example with:

```bash
residue-screen example .local/example
```

The generated input contains a barcode, adapter, composite, complete component
coordinates, synthetic provenance, subject records, exact/reverse/negative truth
cases and explicit draft-analysis defaults. The example is preapproved only as
a software fixture. It is not a manufacturer sequence panel, an empirical BLAST
benchmark, a research approval or biological ground truth for contamination.

The generator and CLI workflow are tested in CI. Outputs remain in ignored local
storage so the repository does not accumulate study data or duplicated reports.