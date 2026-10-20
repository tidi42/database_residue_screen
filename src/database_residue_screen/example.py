import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .contracts import schema, write_table
from .design import digest
from .sequences import reverse_complement


def draft_study(study_id: str, scope: str = "research") -> dict:
    return dict(
        schema_version="1.0.0", study_id=study_id, scope=scope, decision_status="draft", approval_file="",
        inputs=dict(queries="inputs/queries.fasta", query_metadata="inputs/queries.tsv", components="inputs/components.tsv",
                    sources="inputs/sources.tsv", database="inputs/database.json", subjects="inputs/subjects.tsv",
                    validation_cases="inputs/validation_cases.tsv", validation_subjects="inputs/validation_subjects.fasta"),
        analysis=dict(primary_groups=["ONT", "Illumina", "PacBio"], primary_endpoint="strict_full_query_ungapped",
                      secondary_endpoint="identity95_span90_ungapped", condition_unit="canonical_parent_within_condition",
                      pooled_unit="global_canonical_union", parent_weight="equal", inference="descriptive_fixed_panel",
                      multiplicity="no_hypothesis_tests", nesting="disabled_until_validated", overlap="descriptive_only",
                      position_record_type="genomic", position_topology="linear",
                      unknown_record_policy="exclude_primary_position_report_missing", position_windows=[100, 500, 1000],
                      minimum_record_lengths=[1000, 10000, 50000], low_complexity_policy="report_separately_no_posthoc_filter"),
        search=dict(task="blastn-short", word_size=7, evalue=10.0, dust=False, strand="both", real_target_limit=500,
                    control_target_limit=500, retention_ladder=[500, 2000, 10000], relative_stability_tolerance=0.01,
                    traceback="qseq_sseq_btop", completion_policy="per_query_manifest_atomic_batches",
                    failure_policy="abort_run_resume_retries_failed_batches"),
        controls=dict(algorithm="ushuffle", algorithm_version="", k=2, seed=42, sampling="with_replacement",
                      parent_identity="retain_and_report", generation_unit="global_canonical_parent",
                      draw_ladder=[10, 50, 100], max_mc_standard_error=0.01,
                      component_null="not_claimed_descriptive_segment_comparison"),
        validation=dict(unit="query_subject_pair", independent_exact_reference="both_strand_exhaustive", tuning=False,
                        strict_recovery_minimum=1, tolerant_recovery_minimum=0.95, precision_minimum=0.95,
                        required_error_classes=["exact", "reverse", "partial", "substitution", "indel", "no_match", "low_complexity"],
                        database_size_policy="production_effective_search_space"),
        metadata=dict(taxonomy_snapshot="", multi_taxid="explicit_ambiguous_not_duplicated", title_fallback="disabled",
                      project_missingness="report_not_infer_independence", topology_unknown="exclude_primary_position",
                      case_selection="stratified_illustrative_no_population_rate"),
        figures=dict(restoration="v3_design_corrected_data", historical_results="never_use_as_current",
                     unsupported_contrasts="not_estimated", audit_panels="supplement_or_revision_appendix",
                     source_tables="same_versioned_event_release"),
        companion=dict(reuse="undecided", reconciliation_file="", citation_status="pending"),
    )


def save_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def save_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        SeqIO.write([SeqRecord(Seq(sequence), id=name, description="") for name, sequence in records.items()], handle, "fasta")


def create_example(destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=False)
    inputs = destination / "inputs"
    inputs.mkdir()
    barcode = "ACGTCCGTAACC"
    adapter = "GATTCAGGCTAGTTCA"
    sequences = {"tag": barcode, "adapter": adapter, "composite": barcode + adapter}
    save_fasta(inputs / "queries.fasta", sequences)
    queries = [dict(query_id=name, group="Example", condition=name, construct_class="barcode" if name == "tag" else name,
                    family_id="synthetic-family", role="primary", source_id="synthetic-source", provenance_status="synthetic",
                    architecture_status="synthetic" if name == "composite" else "not_applicable") for name in sequences]
    write_table(inputs / "queries.tsv", queries, schema("queries"))
    components = [dict(composite_id="composite", component_id="tag", component_kind="barcode", start=1, end=len(barcode), strand="+"),
                  dict(composite_id="composite", component_id="adapter", component_kind="adapter", start=len(barcode) + 1,
                       end=len(barcode + adapter), strand="+")]
    write_table(inputs / "components.tsv", components, schema("components"))
    write_table(inputs / "sources.tsv", [dict(source_id="synthetic-source", reference="Generated synthetic fixture; no library provenance claim",
                                             version="fixture-v1", evidence_file="inputs/queries.fasta",
                                             evidence_sha256=digest(inputs / "queries.fasta"), verification="synthetic", redistribution="synthetic")], schema("sources"))
    subjects = {"EXAMPLE_POS.1": "GGG" + barcode + adapter + "TTT", "EXAMPLE_REV.1": "GGG" + reverse_complement(barcode) + "TTT",
                "EXAMPLE_NEG.1": "T" * 40}
    save_fasta(inputs / "validation_subjects.fasta", subjects)
    subject_rows = [dict(raw_id=name, accession_version=name, length=len(sequence), topology="linear", record_type="synthetic",
                         taxids="", assembly_id="", project_id="", metadata_source="generated fixture", metadata_date="2000-01-01")
                    for name, sequence in subjects.items()]
    write_table(inputs / "subjects.tsv", subject_rows, schema("subjects"))
    cases = [dict(case_id="exact", query_id="tag", subject_id="EXAMPLE_POS.1", truth="strict_positive", split="evaluation",
                  family_id="synthetic-family", subject_source="fixture-positive", error_class="exact"),
             dict(case_id="reverse", query_id="tag", subject_id="EXAMPLE_REV.1", truth="strict_positive", split="evaluation",
                  family_id="synthetic-family", subject_source="fixture-reverse", error_class="reverse"),
             dict(case_id="negative", query_id="tag", subject_id="EXAMPLE_NEG.1", truth="negative", split="evaluation",
                  family_id="synthetic-family", subject_source="fixture-negative", error_class="no_match")]
    write_table(inputs / "validation_cases.tsv", cases, schema("validation_cases"))
    save_json(inputs / "database.json", dict(name="synthetic-example", kind="synthetic_fasta", build_date="2000-01-01",
                                             blast_version="not_executed", identity_status="synthetic",
                                             evidence_file="inputs/validation_subjects.fasta", evidence_sha256=digest(inputs / "validation_subjects.fasta"),
                                             sequence_count=len(subjects), base_count=sum(map(len, subjects.values()))))
    with (destination / "synthetic_decisions.txt").open("x", encoding="utf-8") as handle:
        handle.write("Synthetic configuration for software tests only. Not research approval or a production-ready search design.\n")
    study = draft_study("synthetic-example", "synthetic_example")
    study.update(decision_status="approved", approval_file="synthetic_decisions.txt")
    study["analysis"]["primary_groups"] = ["Example"]
    study["controls"]["algorithm_version"] = "not_executed_in_synthetic_example"
    study["metadata"]["taxonomy_snapshot"] = "synthetic-no-taxonomy"
    study["validation"]["required_error_classes"] = ["exact", "reverse", "no_match"]
    study["companion"] = dict(reuse="none", reconciliation_file="", citation_status="not_applicable")
    save_json(destination / "study.json", study)
    return destination / "study.json"
