import csv
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

from .contracts import ContractError, schema, write_table
from .design import digest
from .example import draft_study, save_fasta, save_json
from .matching import BLAST_FIELDS, parse_legacy_hsp, strict_hsp, tolerant_hsp
from .sequences import canonical_id, read_fasta, reverse_complement


PANEL_FILES = ("constructs.fasta", "constructs_ont_extended.fasta", "constructs_illumina_extended.fasta", "constructs_pacbio_extended.fasta")
DECOY_FILES = ("decoys.fasta", "decoys_ont_extended.fasta", "decoys_illumina_extended.fasta", "decoys_pacbio_extended.fasta")
REAL_BLAST_FILES = ("ont_barcode_blast_db.tsv", "ont_adapter_barcode_blast_db.tsv", "ont_rapid_barcode_blast_db.tsv",
                    "ont_ligation_adapter_blast_db.tsv", "ont_rapid_adapter_blast_db.tsv", "ont_pcr_adapter_blast_db.tsv",
                    "ont_cdna_adapter_blast_db.tsv", "illumina_truseq_blast_db.tsv", "illumina_nextera_blast_db.tsv",
                    "illumina_smallrna_blast_db.tsv", "illumina_bgi_blast_db.tsv", "illumina_extended_constructs_blast_db.tsv",
                    "pacbio_extended_constructs_blast_db.tsv")
DECOY_BLAST_FILES = ("decoys_blast_db.tsv", "decoys_ont_extended_blast_db.tsv", "decoys_illumina_extended_blast_db.tsv", "decoys_pacbio_extended_blast_db.tsv")
SOURCE_FILES = ("barcodes_only.fasta", "barcodes_with_ont_adapter.fasta", "illumina_ud_index_setA.csv", "illumina_ud_index_setB.csv",
                "Nextera_Index 1_i7_Adapters.csv", "Nextera_Index 2_i5_Adapters.csv", "pacbio_barcodes_repo.csv")
COMPOSITE_CONDITIONS = {"ont_adapter_barcode", "ont_rapid_barcode", "illumina_i7_index_adapter", "illumina_i5_index_adapter", "pacbio_barcode_adapter"}
BARCODE_CONDITIONS = {"ont_barcode", "illumina_i7_index", "illumina_i5_index", "pacbio_barcode"}


def group_of(condition: str) -> str:
    if condition == "illumina_bgi":
        return "BGI/DNBSEQ"
    try:
        return {"ont": "ONT", "illumina": "Illumina", "pacbio": "PacBio"}[condition.split("_")[0]]
    except KeyError as error:
        raise ContractError(f"Unknown legacy condition: {condition}") from error


def query_classification(identifier: str, condition: str) -> tuple[str, str]:
    if identifier == "C2_primer_forward":
        return "primer", "pacbio_c2_primer"
    kind = "composite" if condition in COMPOSITE_CONDITIONS else "barcode" if condition in BARCODE_CONDITIONS else "adapter"
    return kind, condition


def counterpart_id(query_id: str, condition: str) -> str | None:
    if condition == "ont_adapter_barcode" and query_id.endswith("_adapter_plus_barcode"):
        return query_id.removesuffix("_adapter_plus_barcode")
    if condition == "ont_rapid_barcode" and query_id.startswith("RBK"):
        return "NB" + query_id[3:]
    if condition.startswith("illumina_"):
        return query_id.replace("_index_adapter_forward", "_index_forward")
    if condition == "pacbio_barcode_adapter":
        return query_id.replace("_barcode_adapter_forward", "_barcode_forward")
    return None


def recover_components(sequences: dict[str, str], metadata: list[dict]) -> tuple[list[dict], list[str], list[dict]]:
    query_map = {row["query_id"]: row for row in metadata}
    aliases = defaultdict(list)
    for identifier, sequence in sequences.items():
        aliases[canonical_id(sequence)].append(identifier)
    components, added, unresolved = [], [], []
    for query in list(metadata):
        if query["construct_class"] != "composite":
            continue
        query_id = query["query_id"]
        barcode_id = counterpart_id(query_id, query["condition"])
        if barcode_id not in sequences:
            unresolved.append(dict(query_id=query_id, reason="bare_barcode_counterpart_missing"))
            continue
        composite, barcode = sequences[query_id], sequences[barcode_id]
        positions = [(position, strand) for strand, oriented in (("+", barcode), ("-", reverse_complement(barcode)))
                     for position in range(len(composite) - len(oriented) + 1) if composite[position:position + len(oriented)] == oriented]
        if barcode == reverse_complement(barcode):
            positions = [(position, strand) for position, strand in positions if strand == "+"]
        if len(positions) != 1:
            unresolved.append(dict(query_id=query_id, reason="barcode_footprint_ambiguous_or_absent"))
            continue
        position, strand = positions[0]
        query["family_id"] = query_map[barcode_id]["family_id"]
        components.append(dict(composite_id=query_id, component_id=barcode_id, component_kind="barcode",
                               start=position + 1, end=position + len(barcode), strand=strand))
        for start, end in ((0, position), (position + len(barcode), len(composite))):
            if start == end:
                continue
            sequence = composite[start:end]
            key = canonical_id(sequence)
            matches = sorted(identifier for identifier in aliases[key] if query_map[identifier]["construct_class"] in {"adapter", "primer"})
            if matches:
                identifier = matches[0]
                orientation = "+" if sequences[identifier] == sequence else "-"
                kind = query_map[identifier]["construct_class"]
            else:
                identifier = "derived_component_" + key.removeprefix("seq_")
                orientation, kind = "+", "adapter"
                sequences[identifier] = sequence
                record = dict(query_id=identifier, group=query["group"], condition="derived_standalone_component", construct_class=kind,
                              family_id=key, role="reference", source_id=query["source_id"], provenance_status="unverified",
                              architecture_status="not_applicable")
                metadata.append(record)
                query_map[identifier] = record
                aliases[key].append(identifier)
                added.append(identifier)
            components.append(dict(composite_id=query_id, component_id=identifier, component_kind=kind, start=start + 1, end=end, strand=orientation))
    return components, added, unresolved


def source_crosscheck(root: Path, sequences: dict[str, str]) -> dict:
    checked, mismatches, missing = 0, [], []

    def compare(identifier, sequence):
        nonlocal checked
        checked += 1
        if identifier not in sequences:
            missing.append(identifier)
        elif sequences[identifier] != sequence.strip().upper():
            mismatches.append(identifier)

    for filename in SOURCE_FILES:
        path = root / "panels" / filename
        if not path.is_file():
            continue
        if filename.endswith(".fasta"):
            for identifier, sequence in read_fasta(path).items():
                if identifier in sequences:
                    compare(identifier, sequence)
            continue
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if filename.startswith("illumina_ud"):
                    name = re.sub(r"[^A-Za-z0-9]", "", row["index_name"])
                    for index in ("i7", "i5"):
                        compare(name + "_" + index + "_index_forward", row[index + "_bases_in_adapter"])
                elif filename.startswith("Nextera_Index"):
                    index = "i7" if "1_i7" in filename else "i5"
                    name = re.sub(r"[^A-Za-z0-9]", "", row[index + "_index_name"])
                    compare(name + "_" + index + "_index_forward", row["bases_in_adapter"])
                else:
                    compare(row["barcode_name"] + "_barcode_forward", row["barcode_sequence"])
    return dict(comparisons=checked, sequence_mismatches=mismatches, missing_expected_queries=missing,
                meaning="Agreement with local source extracts only; not independent manufacturer or kit authentication")


def audit_decoys(root: Path, sequences: dict[str, str]) -> dict:
    mapped, excluded, failures, repeated, identical = 0, 0, 0, 0, 0
    parents, seen_sequences, seen_ids = Counter(), defaultdict(set), set()
    invalid = []
    for filename in DECOY_FILES:
        path = root / "panels" / filename
        if not path.is_file():
            continue
        for record in SeqIO.parse(path, "fasta"):
            match = re.fullmatch(r"(.+)_decoy[0-9]+", record.id)
            if match is None:
                invalid.append(record.id)
                continue
            parent = match.group(1)
            if parent not in sequences:
                excluded += 1
                continue
            if record.id in seen_ids:
                invalid.append(record.id)
            seen_ids.add(record.id)
            sequence, original = str(record.seq).upper(), sequences[parent]
            if len(sequence) != len(original) or Counter(sequence) != Counter(original):
                invalid.append(record.id)
            actual = Counter(sequence[index:index + 2] for index in range(len(sequence) - 1))
            expected = Counter(original[index:index + 2] for index in range(len(original) - 1))
            failures += actual != expected
            repeated += sequence in seen_sequences[parent]
            identical += sequence == original
            seen_sequences[parent].add(sequence)
            parents[parent] += 1
            mapped += 1
    return dict(mapped_draws=mapped, historical_parent_draws_excluded=excluded, dinucleotide_failures=failures,
                invalid_length_bases_or_ids=invalid, draws_per_parent=dict(Counter(parents.values())),
                parents_without_draws=sorted(set(sequences) - set(parents)), parent_identical_draws=identical,
                repeated_within_parent=repeated, eligible_for_new_null=False)


def audit_blast(path: Path, sequences: dict[str, str], real: bool) -> dict:
    widths, rows, malformed, selected = Counter(), 0, [], 0
    old_ids, strict_ids, tolerant_ids = set(), set(), set()
    old_rows, strict_rows, tolerant_rows = 0, 0, 0
    with path.open(encoding="utf-8", newline="") as handle:
        for number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            rows += 1
            widths[len(fields)] += 1
            try:
                row = parse_legacy_hsp(fields)
                if real and row["qseqid"] in sequences:
                    if row["qlen"] != len(sequences[row["qseqid"]]):
                        raise ContractError("HSP qlen differs from current FASTA")
                    selected += 1
                    old = row["pident"] == 100 and row["qcovs"] == 100 and row["mismatch"] == 0
                    strict, tolerant = strict_hsp(row), tolerant_hsp(row)
                    old_rows += old
                    strict_rows += strict
                    tolerant_rows += tolerant
                    if old:
                        old_ids.add(row["qseqid"])
                    if strict:
                        strict_ids.add(row["qseqid"])
                    if tolerant:
                        tolerant_ids.add(row["qseqid"])
            except (ValueError, KeyError) as error:
                malformed.append(dict(line=number, error=str(error)))
    return dict(file=path.name, rows=rows, column_counts=dict(widths), malformed_rows=len(malformed), malformed_examples=malformed[:10],
                selected_current_query_rows=selected if real else None,
                legacy_positive_query_ids=len(old_ids), strict_positive_query_ids=len(strict_ids), tolerant_positive_query_ids=len(tolerant_ids),
                legacy_perfect_rows=old_rows, strict_rows=strict_rows, tolerant_rows=tolerant_rows,
                missing_new_traceback_fields=list(BLAST_FIELDS[18:]), completion="not_assessed",
                usable_as_validated_framework_search=False)


def import_legacy(root: Path, destination: Path) -> tuple[Path, dict]:
    root = root.resolve()
    required = [("panels", name) for name in PANEL_FILES + DECOY_FILES + SOURCE_FILES]
    required += [("blast", name) for name in REAL_BLAST_FILES + DECOY_BLAST_FILES]
    files = []
    for folder, name in required:
        path = root / folder / name
        files.append(dict(path=f"{folder}/{name}", present=path.is_file(), bytes=path.stat().st_size if path.is_file() else None,
                          sha256=digest(path) if path.is_file() else None))
    missing_panels = [name for name in PANEL_FILES if not (root / "panels" / name).is_file()]
    if missing_panels:
        raise ContractError("Cannot import without current physical panels: " + ", ".join(missing_panels))
    destination.mkdir(parents=True, exist_ok=False)
    inputs, evidence = destination / "inputs", destination / "evidence"
    inputs.mkdir()
    evidence.mkdir()
    sequences, metadata, sources = {}, [], []
    for filename in PANEL_FILES:
        path = root / "panels" / filename
        panel = read_fasta(path)
        descriptions = {record.id: record.description for record in SeqIO.parse(path, "fasta")}
        source_id = "legacy_" + path.stem
        shutil.copyfile(path, evidence / filename)
        sources.append(dict(source_id=source_id, reference=f"Local legacy panel {filename}; source authenticity not verified",
                            version="unrecorded", evidence_file="evidence/" + filename, evidence_sha256=digest(path),
                            verification="unverified", redistribution="unknown"))
        for identifier, sequence in panel.items():
            if identifier in sequences:
                raise ContractError("Duplicate physical ID across panels: " + identifier)
            match = re.search(r"(?:^|\s)condition=([^\s]+)", descriptions[identifier])
            if match is None:
                raise ContractError("Missing condition metadata: " + identifier)
            condition = match.group(1)
            kind, condition = query_classification(identifier, condition)
            metadata.append(dict(query_id=identifier, group=group_of(condition), condition=condition, construct_class=kind,
                                 family_id=canonical_id(sequence), role="exploratory", source_id=source_id,
                                 provenance_status="unverified", architecture_status="unverified" if kind == "composite" else "not_applicable"))
            sequences[identifier] = sequence
    original_sequences = dict(sequences)
    source_check = source_crosscheck(root, original_sequences)
    decoys = audit_decoys(root, original_sequences)
    components, added, unresolved = recover_components(sequences, metadata)
    save_fasta(inputs / "queries.fasta", sequences)
    write_table(inputs / "queries.tsv", metadata, schema("queries"))
    write_table(inputs / "components.tsv", components, schema("components"))
    write_table(inputs / "sources.tsv", sources, schema("sources"))
    write_table(inputs / "subjects.tsv", [], schema("subjects"))
    write_table(inputs / "validation_cases.tsv", [], schema("validation_cases"))
    save_fasta(inputs / "validation_subjects.fasta", {})
    save_json(inputs / "database.json", dict(name="core_nt", kind="blast", build_date="", blast_version="",
                                             identity_status="unverified", evidence_file="", evidence_sha256="", sequence_count=0, base_count=0))
    study = draft_study("legacy-core-nt")
    save_json(destination / "study.json", study)
    outputs = [audit_blast(root / "blast" / filename, original_sequences, filename in REAL_BLAST_FILES)
               for filename in REAL_BLAST_FILES + DECOY_BLAST_FILES if (root / "blast" / filename).is_file()]
    group_counts = []
    original_metadata = [row for row in metadata if row["query_id"] in original_sequences]
    for group in ("ONT", "Illumina", "BGI/DNBSEQ", "PacBio"):
        selected = [row for row in original_metadata if row["group"] == group]
        group_counts.append(dict(group=group, physical_queries=len(selected),
                                 condition_canonical_entries=len({(row["condition"], canonical_id(original_sequences[row["query_id"]])) for row in selected}),
                                 global_canonical_sequences=len({canonical_id(original_sequences[row["query_id"]]) for row in selected})))
    failures = []
    for path in sorted((root / "logs").glob("*failed*.txt")):
        failures.append(dict(file=str(path.relative_to(root)), nonempty_lines=sum(bool(line.strip()) for line in path.read_text().splitlines()),
                             sha256=digest(path), meaning="Historical failure list; not a complete per-query execution manifest"))
    for record in files:
        if record["present"] and digest(root / record["path"]) != record["sha256"]:
            raise ContractError("Source changed during import: " + record["path"])
    report = dict(schema_version="1.0.0", imported_profile="legacy-core-nt", source_root=str(root), source_files=files,
                  missing_files=[record["path"] for record in files if not record["present"]],
                  original_physical_queries=len(original_sequences), group_inventory=group_counts, local_source_crosscheck=source_check,
                  condition_remapping=[dict(query_id="C2_primer_forward", previous_condition="pacbio_smrtbell_adapter", condition="pacbio_c2_primer")],
                  derived_standalone_queries=added, component_rows=len(components), unresolved_composites=unresolved,
                  component_interpretation="Exact decomposition of recorded sequences only; unverified physical architecture and provenance",
                  legacy_decoy_audit=decoys, legacy_blast_audit=outputs, failure_lists=failures,
                  excluded_blast_files=sorted(path.name for path in (root / "blast").glob("*.tsv") if path.name not in REAL_BLAST_FILES + DECOY_BLAST_FILES),
                  source_files_unchanged=True, research_design_ready=False,
                  unresolved=["authenticated primary queries and physical junctions", "source redistribution permissions",
                              "database identity evidence", "scientific decision approval", "independent evaluation truth data",
                              "taxonomy snapshot and record metadata", "shuffle version and resource/precision pilot",
                              "Paper 2 input reconciliation", "modern tracebacks and per-query search completion"])
    save_json(destination / "legacy_inventory.json", report)
    return destination / "study.json", report
