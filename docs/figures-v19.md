# Figure reconstruction specification, v19

No manuscript plot source, editable figure project or publication source tables
were supplied in this repository. This file is a reconstruction specification;
it does not claim that the DOCX-embedded figures have been regenerated.

| Figure | Required change and validation |
|---|---|
| 1C | Place a-priori testability before completion: score retention plus actual-database chance bound. Show `censored_untestable` and `chance_limited` separately with no zero/one inference; only testable parents enter the original completion/resolution branches. Partial tolerant 19–25-nt censoring is a reason flag, not a new primary denominator. |
| 2A / Table 1 | Use testable primary parents; show excluded-class inventory separately as NE. Strict/tolerant common comparisons require both rules testable. Add ‡ for partial tolerant 19–25-nt censoring including ONT 24 nt. Supplement: full panel (1,288) and baseline calibrated lengths (836), separate inventory from actual denominator. Regenerate all affected numbers from retained outputs; no raised-threshold production run. |
| 3A | Regenerate from taxonomic composition with 29278=`vector`, 32630=`synthetic`, residual 81077=`artificial_other`. Reconcile broad union and eligibility separately; don't merely rename every former vector as artificial. The historical total 1,837 and union 3,555 can only be checked against archived events. |
| 4A | Keep historical E=10/diagnostic curves labelled historical. Show exact and minimum-tolerant E curves, fixed production threshold, chance bounds and declared testability classes; identify the database dimensions, scoring and cap. |
| 4B–D | Species labels *Streptococcus suis*, *Mus musculus* in italic glyphs inside the actual image. Caption/footer: “searched (canonical) query orientation”. Stored-record coordinates remain 1-based forward-record coordinates. If the shown alignment is reverse to the record, reverse-complement the whole displayed subject including flanks and reverse its coordinate ordering; never reverse bases alone. Compare rendered Q/S and flanks against verified case packets. A circular vector's stored origin is not a physical end. |
| S4 | Print `78,996 nt`, `213,866 nt`, `32,882 nt` etc. using comma thousands separators consistently. Italicise the two species in the image, preserve accession versions. Apply the same orientation rule as Fig. 4. |
| S1–S3 | Preserve explicit NE/unknown categories and state each denominator, parent/control weights, cap and transformation. No inferred missing publication statistics. |
| Graphical abstract | Replace headline fractions with testable-primary results (pending reanalysis); add testability first and explain that terminal events can cluster by accession/project and occur in designed amplicon/clone ends. Establish icon attribution/licence from the actual source artwork. |

Author data needed: original figure scripts/projects; SD1 event, alias and case
packets; actual palette and icon origins; used 5/10/15 and cap specifications.
Existing context/provenance functions are not the missing full publication driver.
Post-hoc Hoeffding/Holm and cap-stability manuscript outputs cannot be reconstructed
by substituting the generic precision module. No new significance labels are added.

v20 source tables: Figure 4A and Figure S1A guides come from `length_testability.tsv`
(class, reason, minimum qualifying score, its E-value and chance bound per length and
rule); Figure 2A paired strict/tolerant marks and Table S2 common cohorts come from
`common_rule_contrasts.tsv`; Table S3 comes from `conditional_contrasts.tsv`; Table S9
and the historical 1,288/836 fractions come from `length_sensitivity.tsv`
(`historical_zero_filled_*`, labelled sensitivity only).
