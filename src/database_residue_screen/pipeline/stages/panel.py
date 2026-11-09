"""Stage `panel`: build the panel from the source snapshot, select the study subset and derive the vendor-derived
materials of later stages (census flank elements, barcode catalogue, benchmark templates)."""
from pathlib import Path

from ...constructs import Library, build_panel, load_manifest, read_definition, write_fasta
from ...contracts import read_json, read_table, schema, write_table
from ...design import digest
from ...example import save_json
from ...materials import CATALOGUE_SCHEMA, SEGMENT_SCHEMA, derive_catalogue, derive_elements, derive_templates, segment_rows, select_subset
from ...sequences import read_fasta


def panel(context):
    """Derive every query, verify it against its definition hash, select the study subset and derive run-time materials."""
    section, subset = context.config["panel"], context.config["subset"]
    snapshot = context.upstream["sources"].output / "sources"
    pdftotext = context.tool("pdftotext")
    definition = Path(section["definition"])
    build = build_panel(definition, snapshot / "files", context.output / "build", snapshot / "sources.json", pdftotext)
    entries = {entry["query_id"]: entry for entry in read_definition(definition)}
    sequences = read_fasta(context.output / "build" / "queries.fasta")
    selected, selection = select_subset(sequences, subset["parents"], subset["query_ids"])
    write_fasta(context.output / "queries.fasta", {name: sequences[name] for name in selected})
    metadata = {row["query_id"]: row for row in read_table(context.output / "build" / "queries.tsv", schema("queries"))}
    write_table(context.output / "queries.tsv", [metadata[name] for name in selected], schema("queries"))
    write_table(context.output / "segments.tsv", segment_rows(entries, selected), SEGMENT_SCHEMA)
    save_json(context.output / "selection.json", selection)
    elements = derive_elements(section["census_elements"], entries, sequences)
    save_json(context.output / "census_elements.json", dict(kind="census_flank_elements", elements=elements))
    library = Library(snapshot / "files", load_manifest(snapshot / "sources.json"), pdftotext)
    catalogue, catalogue_summary = derive_catalogue(section["census_catalogue"], section["census_catalogue_checks"], library, entries, sequences)
    write_table(context.output / "census_catalogue.tsv", catalogue, CATALOGUE_SCHEMA)
    templates = derive_templates(section["benchmark_templates"], sequences)
    save_json(context.output / "benchmark_templates.json", dict(kind="benchmark_templates", templates=templates))
    drift = read_json(snapshot / "collection.json")["drift"]
    accepted = [dict(name=row["name"], pinned_sha256=row["pinned_sha256"], sha256=row["sha256"], origin=row["origin"],
                     status="accepted_every_derived_sequence_matches_its_definition_sha256") for row in drift]
    for row in accepted:
        context.log.warning("drifted source %s accepted: every derived sequence matches its definition SHA-256", row["name"])
    files = ("queries.fasta", "queries.tsv", "segments.tsv", "selection.json", "census_elements.json", "census_catalogue.tsv",
             "benchmark_templates.json")
    report = dict(kind="study_panel", definition=str(definition), definition_sha256=digest(definition),
                  build={key: build[key] for key in ("queries", "distinct_sequences", "canonical_parents", "by_group", "by_family",
                                                     "by_verification", "files")},
                  selection={key: value for key, value in selection.items() if key not in ("query_ids", "parent_ids")},
                  source_drift=accepted,
                  census_elements={item["name"]: dict(length=item["length"], sha256=item["sha256"], source_query=item["source_query"],
                                                      queries_checked=item["queries_checked"]) for item in elements},
                  census_catalogue=catalogue_summary,
                  benchmark_templates=[{key: item[key] for key in ("name", "query", "length", "repeat_from", "sha256")} for item in templates],
                  materials_scope="Derived from the complete built panel and its pinned sources; the subset selects only the study queries",
                  files={name: digest(context.output / name) for name in files})
    save_json(context.output / "panel.json", report)
    return dict(panel_queries=build["queries"], panel_parents=build["canonical_parents"], study_queries=selection["selected_queries"],
                study_parents=selection["selected_parents"], selection=selection["mode"], census_elements=len(elements),
                catalogue_entries=len(catalogue), benchmark_templates=len(templates), source_drift=len(accepted))
