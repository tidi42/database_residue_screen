"""Record-case table: sequence-level checks of the record cases selected by the record_cases stage.

For every selected case the displayed window is rebuilt from the retrieved forward record window and checked against
the searched query (searched canonical orientation, HSP strand) and, when the case packet is present, against the
retrieved FASTA. Flanks are annotated with the longest stretch adjacent to the match that also occurs in another panel
reference (both orientations; sequences from the run's design outputs only). The strict testability of the query length
is attached. Exact sequence verification supports triage, not biological adjudication.
"""
import re
from collections import Counter
from pathlib import Path

from ..contracts import ContractError
from ..sequences import canonical_id, read_fasta, reverse_complement
from .common import read_tsv, write_json, write_tsv

REQUIRED = ("case_id", "category", "accession_version", "query_id", "parent_id", "group", "qseq", "sseq", "sstart", "send", "lower",
            "upper", "record_length", "distance_to_end", "display_left_flank", "display_right_flank", "display_orientation",
            "display_segment_roles", "forward_sequence", "flank_start", "topology", "molecule", "record_type", "origin", "context",
            "organism", "title", "taxids", "linear_eligibility")
MINIMUM_ADJACENT = 8


def binomial(name: str) -> bool:
    return re.fullmatch(r"[A-Z][a-z]+ [a-z]+", name or "") is not None


def adjacent_panel_match(flank: str, from_end: bool, references: dict[str, str]) -> tuple[int, str]:
    """Longest stretch of the flank next to the match (its 3' end if ``from_end``) occurring in a panel reference."""
    for size in range(len(flank), MINIMUM_ADJACENT - 1, -1):
        piece = flank[len(flank) - size:] if from_end else flank[:size]
        owners = sorted(name for name, sequence in references.items() if piece in sequence or piece in reverse_complement(sequence))
        if owners:
            return size, owners[0]
    return 0, "none"


def packet_sequence(record_cases: Path, case: dict) -> str | None:
    evidence = case.get("evidence_file", "")
    if not evidence:
        return None
    packet = (record_cases / evidence).parent / "retrieved.fasta"
    if not packet.is_file():
        return None
    return "".join(line.strip() for line in packet.read_text(encoding="utf-8").splitlines() if not line.startswith(">")).upper()


def write_record_cases(destination: Path, *, record_cases: Path, physical_fasta: Path, search_fasta: Path, policy: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    cases = read_tsv(record_cases / "record_cases.tsv") if (record_cases / "record_cases.tsv").is_file() else []
    cases = [case for case in cases if case.get("accession_version")]
    if cases:
        missing = [name for name in REQUIRED if name not in cases[0]]
        if missing:
            raise ContractError(f"record_cases.tsv lacks columns: {', '.join(missing)}")
    physical, searched = read_fasta(physical_fasta), read_fasta(search_fasta)
    calibration = {row["query_length"]: row for row in policy.get("detectability", [])}
    rows = []
    for case in cases:
        query = physical.get(case["query_id"])
        if query is None:
            raise ContractError(f"Record case {case['case_id']}: query {case['query_id']} is not a panel physical query")
        canonical_query = searched.get(case["parent_id"])
        plus = int(case["sstart"]) <= int(case["send"])
        forward = case["forward_sequence"].upper()
        offset = int(case["lower"]) - int(case["flank_start"])
        window = forward[offset:offset + int(case["upper"]) - int(case["lower"]) + 1]
        shown = case["display_left_flank"] + case["sseq"] + case["display_right_flank"]
        displayed = forward if case["display_orientation"] == "+" else reverse_complement(forward)
        packet = packet_sequence(record_cases, case)
        others = {name: sequence for name, sequence in physical.items() if canonical_id(sequence) != case["parent_id"]}
        left_size, left_owner = adjacent_panel_match(case["display_left_flank"], True, others)
        right_size, right_owner = adjacent_panel_match(case["display_right_flank"], False, others)
        strict = calibration.get(len(case["qseq"]), {})
        rows.append(dict(
            case_id=case["case_id"], category=case["category"], accession_version=case["accession_version"], organism=case["organism"],
            organism_italic=binomial(case["organism"]), query_id=case["query_id"], parent_id=case["parent_id"],
            group=case.get("group", ""), condition=case.get("condition", ""), query_length=len(query),
            query_parent_matches=canonical_id(query) == case["parent_id"], searched_sequence_found=canonical_query is not None,
            hsp_full_length_exact=canonical_query is not None and case["qseq"] == case["sseq"] == canonical_query,
            hsp_strand="plus" if plus else "minus",
            record_window_equals_hsp=window == (case["sseq"] if plus else reverse_complement(case["sseq"])),
            displayed_window_equals_record=shown in displayed, packet_checked=packet is not None,
            packet_identical=None if packet is None else packet == forward, display_orientation=case["display_orientation"],
            topology=case["topology"], molecule=case["molecule"], record_type=case["record_type"], record_length=int(case["record_length"]),
            distance_to_end=int(case["distance_to_end"]), left_flank_nt=len(case["display_left_flank"]),
            right_flank_nt=len(case["display_right_flank"]), left_flank_clipped=case.get("left_flank_clipped", ""),
            right_flank_clipped=case.get("right_flank_clipped", ""), left_adjacent_panel_nt=left_size, left_adjacent_panel_query=left_owner,
            right_adjacent_panel_nt=right_size, right_adjacent_panel_query=right_owner,
            strict_testability=strict.get("strict_testability", "not_calibrated"),
            strict_expected_chance_bound=strict.get("strict_expected_chance_bound")))
    fields = ("case_id", "category", "accession_version", "organism", "organism_italic", "query_id", "parent_id", "group", "condition",
              "query_length", "query_parent_matches", "searched_sequence_found", "hsp_full_length_exact", "hsp_strand",
              "record_window_equals_hsp", "displayed_window_equals_record", "packet_checked", "packet_identical", "display_orientation",
              "topology", "molecule", "record_type", "record_length", "distance_to_end", "left_flank_nt", "right_flank_nt",
              "left_flank_clipped", "right_flank_clipped", "left_adjacent_panel_nt", "left_adjacent_panel_query", "right_adjacent_panel_nt",
              "right_adjacent_panel_query", "strict_testability", "strict_expected_chance_bound")
    failed = [row["case_id"] for row in rows if not (row["query_parent_matches"] and row["searched_sequence_found"] and row["record_window_equals_hsp"]
                                                     and row["displayed_window_equals_record"] and row["packet_identical"] is not False)]
    if failed:
        raise ContractError(f"Record cases do not reproduce at sequence level: {', '.join(failed)}")
    write_tsv(destination / "record_case_table.tsv", rows, fields, boolean=("True", "False"))
    status_path = record_cases / "case_category_status.tsv"
    categories = read_tsv(status_path) if status_path.is_file() else []
    write_tsv(destination / "case_category_status.tsv", categories,
              list(categories[0]) if categories else ("category", "selection_subcategory", "eligible_candidates", "attempted_accessions",
                                                      "selected_candidate", "status"))
    summary = dict(cases=len(rows), categories=dict(sorted(Counter(row["category"] for row in rows).items())),
                   packets_checked=sum(row["packet_checked"] for row in rows),
                   unavailable_categories=sorted(row["category"] for row in categories if row.get("status", "selected") != "selected"),
                   strict_testability=dict(sorted(Counter(row["strict_testability"] for row in rows).items())),
                   statement="Exact sequence verification supports triage, not biological adjudication or a validation rate")
    write_json(destination / "summary.json", summary)
    return summary
