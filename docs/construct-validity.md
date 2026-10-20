# Construct validity boundary

`examples/` and `example.py` contain synthetic software fixtures, not vendor
panels. The actual four submitted FASTAs and source CSVs are absent. The legacy
adapter imports sequences unchanged and annotates identifiable components;
`recover_components` is not a vendor construct synthesiser. Do not “correct” an
unseen historical sequence or infer an authentic molecule from concatenation.

The v19 audit alongside the patch checks the manuscript-described classes against
Illumina 1000000002694 v22, ONT Chemistry Technical Document
CHTD_500_v1_revAV_24July2026, PacBio manufacturer documentation and MGI adapter
documents. The described Illumina i7/i5 formulas match the manufacturer's Nextera
and Illumina UD formulas; each actual index/version still needs checking. Bases
in adapter, sample-sheet bases and index-read orientation are separate fields.

ONT Native NB and Rapid RB identifiers do not imply an interchangeable insert in
identical flanks. Official NB01 forward is the reverse complement of RB01; that
relationship does not hold uniformly for every numbered barcode. A 23-nt motif
plus a barcode is a designed search representation until its kit/junction is
verified. Published Rapid flanks plus a 24-nt tag sum to 90 nt, not the described
111 nt; additional sequence needs its own source and orientation evidence.

PacBio historical 16-nt barcode sets and 45-nt blunt hairpins must not be assumed
to form the actual 61-nt commercial adapter for every BC1001–BC1096 identifier.
The historical manufacturer shared protocol explicitly lists a 77-base oligo
for a 16-base barcoded adapter (16 + 45 + 16); its spreadsheet supplies 24
such examples. A 61-base barcode+hairpin query is therefore a one-sided
representation, not the complete oligo in that design. Modern SMRTbell barcoded adapter plate 3.0 uses 10-bp barcodes and 5-bp padding.
This is evidence against treating different kit generations as interchangeable,
not proof that the submitted designed query is sequence-wrong. The C2 oligo is a
sequencing primer and must remain separate from ligated adapters.

Required per-query release: physical ID; submitted 5′–3′ sequence and SHA256;
canonical sequence/orientation; component coordinates; kit and revision; exact
vendor locator/archived source hash; index bases-in-adapter; platform-specific
sample-sheet convention; verified/inferred/uncheckable verdict. A barcode or
index name alone cannot authenticate sequence identity (also the Paper 2 concern).
