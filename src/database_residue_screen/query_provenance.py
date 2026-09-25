import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from .contracts import ContractError, read_json, read_table, schema
from .design import digest
from .legacy import PANEL_FILES, SOURCE_FILES, counterpart_id
from .sequences import canonical_id, read_fasta, reverse_complement


def sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def component_map(query: dict, sequence: str, components: list[dict], sequences: dict[str, str]) -> tuple[list, list]:
    candidates = []
    if query["construct_class"] != "composite":
        role = "barcode" if query["construct_class"] == "barcode" else "scaffold" if query["construct_class"] in {"adapter", "primer"} else "unknown"
        pieces = [dict(start=1, end=len(sequence), component_id=query["query_id"], strand="+", role=role,
                       assignment_status="recorded_query_class" if role != "unknown" else "unknown_query_class")]
    else:
        barcode_id = counterpart_id(query["query_id"], query["condition"])
        barcode = sequences.get(barcode_id)
        if barcode:
            for strand, oriented in (("+", barcode), ("-", reverse_complement(barcode))):
                if strand == "-" and oriented == barcode:
                    continue
                candidates.extend(dict(query_id=query["query_id"], component_id=barcode_id, start=position + 1,
                                       end=position + len(barcode), strand=strand)
                                  for position in range(len(sequence) - len(barcode) + 1)
                                  if sequence[position:position + len(barcode)] == oriented)
        if len(candidates) != 1 or not components:
            pieces = [dict(start=1, end=len(sequence), component_id="unknown", strand="unknown", role="unknown",
                           assignment_status="barcode_boundary_ambiguous" if len(candidates) > 1 else "component_evidence_missing")]
        else:
            pieces = []
            candidate = candidates[0]
            recorded_barcodes = [row for row in components if row["component_kind"] == "barcode"]
            if len(recorded_barcodes) != 1 or any(recorded_barcodes[0][key] != candidate[key] for key in ("component_id", "start", "end", "strand")):
                raise ContractError("Archived barcode coordinates disagree with the unique local sequence footprint")
            for row in sorted(components, key=lambda item: (item["start"], item["end"], item["component_id"])):
                component_sequence = sequences.get(row["component_id"])
                if component_sequence is None or row["strand"] not in {"+", "-"}:
                    raise ContractError("Archived component lacks its sequence or orientation")
                expected = component_sequence if row["strand"] == "+" else reverse_complement(component_sequence)
                if not 1 <= row["start"] <= row["end"] <= len(sequence) or sequence[row["start"] - 1:row["end"]] != expected:
                    raise ContractError("Archived component sequence does not match the submitted query interval")
                pieces.append(dict(start=row["start"], end=row["end"], component_id=row["component_id"], strand=row["strand"],
                                   role="barcode" if row["component_kind"] == "barcode" else "scaffold",
                                   assignment_status="unique_local_barcode_footprint" if row["component_kind"] == "barcode" else "recorded_non_barcode_construct_segment"))
            position, filled = 1, []
            for piece in pieces:
                if piece["start"] < position:
                    raise ContractError("Archived components overlap; no double-counted segment assignment is allowed")
                if piece["start"] > position:
                    filled.append(dict(start=position, end=piece["start"] - 1, component_id="unknown", strand="unknown", role="unknown", assignment_status="unassigned_query_bases"))
                filled.append(piece)
                position = piece["end"] + 1
            if position <= len(sequence):
                filled.append(dict(start=position, end=len(sequence), component_id="unknown", strand="unknown", role="unknown", assignment_status="unassigned_query_bases"))
            pieces = filled
    result = []
    for index, piece in enumerate(pieces, 1):
        fragment = sequence[piece["start"] - 1:piece["end"]]
        result.append(dict(query_id=query["query_id"], segment_id=f"{query['query_id']}:segment_{index:02d}", **piece,
                           component_sequence=fragment, component_sequence_sha256=sequence_digest(fragment),
                           physical_architecture_authenticated=False))
    if sum(row["end"] - row["start"] + 1 for row in result) != len(sequence):
        raise ContractError("Segment assignments do not cover every submitted query base")
    return result, candidates


def local_extracts(root: Path) -> tuple[dict, dict]:
    by_query, files = defaultdict(list), {}
    for name in SOURCE_FILES:
        path = root / "panels" / name
        if not path.is_file():
            continue
        files[str(path.resolve())] = digest(path)
        if path.suffix == ".fasta":
            for index, record in enumerate(SeqIO.parse(path, "fasta"), 1):
                by_query[record.id].append(dict(file=str(path.relative_to(root)), record=index, column="FASTA_sequence", sequence=str(record.seq).upper()))
            continue
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for index, row in enumerate(csv.DictReader(handle), 2):
                if name.startswith("illumina_ud"):
                    identifier = re.sub(r"[^A-Za-z0-9]", "", row["index_name"])
                    for tag in ("i7", "i5"):
                        by_query[identifier + "_" + tag + "_index_forward"].append(dict(file=str(path.relative_to(root)), record=index,
                            column=tag + "_bases_in_adapter", sequence=row[tag + "_bases_in_adapter"].strip().upper()))
                elif name.startswith("Nextera_Index"):
                    tag = "i7" if "1_i7" in name else "i5"
                    identifier = re.sub(r"[^A-Za-z0-9]", "", row[tag + "_index_name"])
                    by_query[identifier + "_" + tag + "_index_forward"].append(dict(file=str(path.relative_to(root)), record=index,
                        column="bases_in_adapter", sequence=row["bases_in_adapter"].strip().upper()))
                else:
                    by_query[row["barcode_name"] + "_barcode_forward"].append(dict(file=str(path.relative_to(root)), record=index,
                        column="barcode_sequence", sequence=row["barcode_sequence"].strip().upper()))
    return by_query, files


def declared_source(query: dict) -> dict:
    source = query["source_id"]
    if source == "legacy_constructs_ont_extended":
        return dict(declared_reference="Porechop adapters.py; barcode counterparts from the local ONT panel",
                    declared_url="https://github.com/rrwick/Porechop", declared_version="unrecorded",
                    construction_evidence="scripts/build_ont_extended_panel.py", primary_source_status="declared_unverified")
    if source == "legacy_constructs_illumina_extended":
        return dict(declared_reference="Illumina UD Indexes adapter page, document 1000000002694; local physical-index CSVs",
                    declared_url="unrecorded", declared_version="v22_as_declared_in_builder",
                    construction_evidence="scripts/build_illumina_extended_panel.py", primary_source_status="declared_unverified")
    if source == "legacy_constructs_pacbio_extended":
        reference = "PacificBiosciences/barcodes; local Sequel 96-plex CSV"
        if query["construct_class"] == "composite":
            reference += "; NCBI UniVec NGB00972.1; author-recorded barcode + adapter concatenation"
        elif query["query_id"] in {"SMRTbell_blunt_adapter_forward", "C2_primer_forward"}:
            reference = "NCBI UniVec " + ("NGB00973.1" if query["query_id"] == "C2_primer_forward" else "NGB00972.1")
        return dict(declared_reference=reference, declared_url="https://github.com/PacificBiosciences/barcodes" if query["construct_class"] in {"barcode", "composite"} else "unrecorded",
                    declared_version="UniVec_accession_version_declared; barcode_repository_commit_unrecorded",
                    construction_evidence="scripts/build_pacbio_extended_panel.py", primary_source_status="declared_unverified")
    return dict(declared_reference="Inherited local reference panel; supplier-source locator not individually archived",
                declared_url="unrecorded", declared_version="unrecorded", construction_evidence="panels/constructs.fasta",
                primary_source_status="not_individually_recorded")


def provenance_rows(root: Path, panel: Path, study_inputs: Path, reference_basis: Path) -> tuple[list, list, list, dict]:
    basis = read_json(reference_basis)
    inputs = {str(reference_basis.resolve()): digest(reference_basis)}
    for name, expected in basis["source_panel_sha256"].items():
        if digest(Path(name)) != expected:
            raise ContractError("Original panel changed since the reference-input preparation")
        inputs[name] = expected
    if digest(Path(basis["imported_source_table"])) != basis["imported_source_table_sha256"]:
        raise ContractError("Archived source declarations changed")
    inputs[basis["imported_source_table"]] = basis["imported_source_table_sha256"]
    physical = read_table(panel / "physical_queries.tsv", schema("queries"))
    sequences = read_fasta(panel / "physical_queries.fasta")
    archived_sequences = read_fasta(study_inputs / "queries.fasta")
    components = read_table(study_inputs / "components.tsv", schema("components"))
    grouped = defaultdict(list)
    for row in components:
        grouped[row["composite_id"]].append(row)
    for name in ("queries.fasta", "components.tsv"):
        inputs[str((study_inputs / name).resolve())] = digest(study_inputs / name)
    originals = {}
    for name in PANEL_FILES:
        path = root / "panels" / name
        for index, record in enumerate(SeqIO.parse(path, "fasta"), 1):
            if record.id in originals:
                raise ContractError("Duplicate physical ID in original source panels")
            originals[record.id] = dict(source_file=str(path.relative_to(root)), source_file_sha256=digest(path),
                                        source_record=index, source_header=record.description, sequence=str(record.seq).upper())
    if set(originals) != set(sequences) or {row["query_id"] for row in physical} != set(sequences) or len(physical) != len(sequences):
        raise ContractError("Physical-query register does not exactly match original source records")
    extracts, extract_files = local_extracts(root)
    inputs.update(extract_files)
    rows, segments, alternatives = [], [], []
    for query in physical:
        identifier, sequence = query["query_id"], sequences[query["query_id"]]
        original = originals[identifier]
        if original["sequence"] != sequence or archived_sequences.get(identifier) != sequence:
            raise ContractError("Physical query sequence differs from its original or archived local record")
        local = extracts.get(identifier, [])
        counterpart = counterpart_id(identifier, query["condition"]) if query["construct_class"] == "composite" else None
        inherited = extracts.get(counterpart, []) if counterpart else []
        for source in local:
            if source["sequence"] != sequence:
                raise ContractError("Local sequence extract disagrees with submitted physical query")
        for source in inherited:
            if source["sequence"] != sequences.get(counterpart):
                raise ContractError("Composite's barcode source disagrees with the unchanged counterpart")
        mapped, candidates = component_map(query, sequence, grouped[identifier], archived_sequences)
        sources = local if local else inherited
        declared = declared_source(query)
        builder = root / declared["construction_evidence"]
        inputs[str(builder.resolve())] = digest(builder)
        rows.append(dict(query_id=identifier, parent_id=canonical_id(sequence), sequence_id=canonical_id(sequence),
                         alias_orientation="+" if sequence <= reverse_complement(sequence) else "-",
                         **{name: query[name] for name in ("group", "condition", "construct_class", "source_id")},
                         **original, length=len(sequence), sequence_sha256=sequence_digest(sequence),
                         derivation="author_assembled_composite" if query["construct_class"] == "composite" else "local_extract_transcription" if local else "inherited_local_reference",
                         local_evidence_status="exact_full_query_agreement" if local else "exact_barcode_counterpart_agreement" if inherited else "original_panel_identity_only",
                         local_extracts=";".join(f"{item['file']}:{item['record']}:{item['column']}" for item in sources),
                         **declared, source_retrieval_date="unrecorded", physical_architecture_authenticated=False,
                         segment_status="explicit_unknown" if any(row["role"] == "unknown" for row in mapped) else "locally_resolved_not_physically_authenticated"))
        for row in mapped:
            row.update(parent_id=canonical_id(sequence), alias_orientation=rows[-1]["alias_orientation"],
                       evidence_file=str((study_inputs / "components.tsv").resolve()) if query["construct_class"] == "composite" else original["source_file"],
                       evidence_sha256=inputs[str((study_inputs / "components.tsv").resolve())] if query["construct_class"] == "composite" else original["source_file_sha256"])
        segments.extend(mapped)
        alternatives.extend(candidates)
    return rows, segments, alternatives, inputs