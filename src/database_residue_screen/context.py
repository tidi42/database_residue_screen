"""Record-context helpers: verified flank retrieval and component support of an accepted traceback HSP."""
from .contracts import ContractError
from .sequences import accession_version


def validate_flank_record(record, expected_accession: str, start: int, end: int) -> str:
    returned = accession_version(record.id.split(":", 1)[0])
    if returned != expected_accession:
        raise ContractError(f"Accession lookup mismatch: requested {expected_accession}, returned {returned}")
    if len(record.seq) != end - start + 1:
        raise ContractError("Retrieved flank has an unexpected length")
    if record.id != f"{expected_accession}:{start}-{end}":
        raise ContractError("Retrieved FASTA header does not confirm the requested range")
    return str(record.seq).upper()


def traceback_segment_coverage(hsp: dict, component: dict, alias_orientation: str = "+") -> dict:
    length = hsp["qlen"]
    start, end = component["start"], component["end"]
    if alias_orientation not in {"+", "-"} or not 1 <= start <= end <= length:
        raise ContractError("Invalid physical-query segment or alias orientation")
    aligned_query, aligned_subject = hsp["qseq"], hsp["sseq"]
    if (len(aligned_query) != len(aligned_subject) or "-" in aligned_query + aligned_subject
            or len(aligned_query) != abs(hsp["qend"] - hsp["qstart"]) + 1
            or not 1 <= min(hsp["qstart"], hsp["qend"]) <= max(hsp["qstart"], hsp["qend"]) <= length):
        raise ContractError("Component annotation requires a valid single ungapped accepted HSP")
    step = 1 if hsp["qend"] >= hsp["qstart"] else -1
    covered, differences = [], []
    for offset, (query_base, subject_base) in enumerate(zip(aligned_query, aligned_subject)):
        canonical_position = hsp["qstart"] + step * offset
        position = canonical_position if alias_orientation == "+" else length - canonical_position + 1
        if start <= position <= end:
            covered.append(position)
            if query_base != subject_base:
                differences.append(position)
    size, matching = len(covered), len(covered) - len(differences)
    complete = size == end - start + 1
    state = "not_covered" if not size else ("complete" if complete else "partial") + ("_exact" if not differences else "_mismatched")
    return dict(covered_bases=size, matching_bases=matching, mismatches=len(differences),
                component_length=end - start + 1, coverage_fraction=size / (end - start + 1),
                complete_segment=int(complete), support_state=state,
                covered_start=min(covered, default=0), covered_end=max(covered, default=0),
                mismatch_positions=";".join(map(str, sorted(differences))))
