# Paper 3: library-preparation sequences in NCBI core_nt

This directory holds everything that defines the Paper 3 study except the data:

| File | Content |
| --- | --- |
| `study.toml` | Complete study configuration of the `screen` pipeline (see [docs/study-workflow.md](../../docs/study-workflow.md)). |
| `study_dryrun.toml` | The same configuration on one query of every family, two exact-scan volumes and one benchmark shard per fixture panel. |
| `panel_definition.tsv` | The 1,456 physical queries: rule, source locators, segments and the SHA-256 of each sequence. |
| `sources.json` | The 15 pinned public documents (URL, retrieval date 2026-10-20, SHA-256) from which every query is derived. |
| `figures/` | R scripts of the publication figures. |

No file contains a sequence. Run the study with

```bash
residue-screen run-study studies/paper3/study.toml --output /data/runs/paper3
```

after setting the environment variables listed at the top of `study.toml`. The
panel alone can be built with

```bash
residue-screen fetch-sources studies/paper3/sources.json --directory .local/paper3-sources
residue-screen build-panel studies/paper3/panel_definition.tsv \
  --sources .local/paper3-sources --output .local/paper3-panel
```

## Panel

| Group | Families (queries) |
| --- | --- |
| ONT (392) | `ont_native_barcode_24` 96, `ont_native_barcode_context_k40` 96, `ont_native_barcode_context_k46` 96, `ont_rapid_barcode_context_90` 96, `ont_ligation_adapter` 2, `ont_rapid_adapter` 2, `ont_pcr_adapter` 3, `ont_cdna_adapter` 1 |
| Illumina (868) | `illumina_i7_bare` 218, `illumina_i7_embedded` 218, `illumina_i5_bare` 210, `illumina_i5_embedded` 210, `illumina_nextera` 3, `illumina_smallrna` 2, `illumina_truseq` 7 |
| BGI/DNBSEQ (2) | `mgi_dnbseq_filter_adapter` 2 |
| PacBio (194) | `pacbio_barcode_16` 96, `pacbio_barcode_plus_hairpin_61` 96, `pacbio_smrtbell_adapter` 2 |

The panel has 1,384 canonical parents (identical or reverse-complement sequences
counted once).

- Rapid barcode contexts use the Dorado `SQK-RBK114-96` (`RBK_1_96`) barcodes;
  positions 26, 39, 40, 48, 54 and 60 differ from the generic ONT barcode table,
  which `panel_build.json` records per query.
- `verified_secondary_source` (14): four ONT oligos of earlier kits defined only
  by Porechop, and ten PacBio barcodes (BC1005, BC1014, BC1025, BC1038, BC1041,
  BC1042, BC1051, BC1086, BC1089, BC1090) that are absent from the official
  Sequel 96 v2 file and taken from a third-party mirror of the Sequel 384 v1 set.
- `search_representation` (96): PacBio barcode followed by the blunt hairpin,
  the first 61 nt of the documented 16 + 45 + 16 barcoded adapter.
- `C2_primer_forward` belongs to condition `pacbio_smrtbell_adapter` with
  construct class `primer`.

Vendor pages change. The `sources` stage uses matching local copies and downloads
missing files; a file whose SHA-256 differs from its pin is refused unless
`allow_source_drift = true`, and then accepted only if every query derived from
it still matches its definition SHA-256 (recorded in `panel.json`). Published
Illumina sequences require the notice "Oligonucleotide sequences © 2025
Illumina, Inc. All rights reserved."
