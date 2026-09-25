import io
import os
import re
import subprocess
from pathlib import Path

from Bio import SeqIO

from .contracts import ContractError
from .context import validate_flank_record
from .design import digest
from .example import save_json
from .sequences import reverse_complement


def selection_specification() -> dict:
    return dict(version="record_triage_v1", matching_rule="strict_single_ungapped_full_query",
                category_order=["terminal", "synthetic_or_vector", "internal_or_ambiguous"],
                terminal_window=500, terminal_minimum_length=10000, flank_bases=50,
                maximum_attempts_per_category=3, maximum_selected_per_category=1,
                sort_order=["accession_version", "parent_id", "lower", "upper", "canonical_strand"],
                no_repeated_selected_accession=True, internal_preferred_to_ambiguous=True,
                synthetic_rule="Archived title explicitly declares a cloning/expression vector or synthetic construct; taxon/context alone is insufficient",
                incomplete_or_unidentified="retained_in_candidate_register_not_eligible_for_verified_panel",
                selection_basis="Purposive deterministic triage illustrations; not representative biological validation",
                new_blast_searches=0, external_retrieval_authorized=False, final_freeze=False)


def category_flags(event: dict, metadata: dict | None, linear_decision: tuple[str, str], completed: bool = True) -> dict:
    exact = (not event["accession_version"].startswith("unresolved:") and event.get("identity_status") == "exact_version"
             and metadata is not None and metadata.get("accessionversion") == event["accession_version"]
             and int(metadata.get("slen", -1)) == event["slen"] and completed)
    title = metadata.get("title", "") if metadata else ""
    declaration = re.match(r"^(?:cloning vector\b|expression vector\b|synthetic construct\b|synthetic DNA\b)", title, flags=re.IGNORECASE)
    genomic = exact and linear_decision[0] == "included"
    return dict(evidence_complete=bool(exact), terminal=bool(genomic and event["slen"] >= 10000 and event["distance_to_end"] <= 500),
                synthetic_or_vector=bool(exact and declaration), internal=bool(genomic and event["distance_to_end"] > 500),
                ambiguous=bool(exact and linear_decision[0] == "unresolved"),
                synthetic_support="archived_title_explicit_context_declaration" if declaration else "none",
                exclusion_reason="" if exact else "incomplete_or_missing_exact_version_length_metadata")


def choose_cases(candidates: list[dict], retrieve, specification: dict) -> tuple[list, list, list]:
    selected, attempts, categories, used = [], [], [], set()
    order = lambda row: tuple(row[name] for name in specification["sort_order"])
    for category in specification["category_order"]:
        if category == "internal_or_ambiguous":
            field = "internal" if any(row["internal"] and row["accession_version"] not in used for row in candidates) else "ambiguous"
        else:
            field = category
        eligible = sorted((row for row in candidates if row[field]), key=order)
        attempted_accessions, chosen = set(), None
        for rank, row in enumerate(eligible, 1):
            if row["accession_version"] in used or row["accession_version"] in attempted_accessions:
                attempts.append(dict(category=category, candidate_id=row["candidate_id"], rank=rank, attempted=False,
                                     status="skipped_duplicate_accession", reason="selected_earlier_or_attempted_within_category"))
                continue
            if len(attempted_accessions) >= specification["maximum_attempts_per_category"]:
                break
            attempted_accessions.add(row["accession_version"])
            result = retrieve(row, category, len(attempted_accessions))
            attempts.append(dict(category=category, candidate_id=row["candidate_id"], rank=rank, attempted=True,
                                 status=result["status"], reason=result.get("error", "")))
            if result["status"] == "sequence_verified_not_biologically_adjudicated":
                chosen = dict(row, category=category, selection_subcategory=field, **result)
                selected.append(chosen)
                used.add(row["accession_version"])
                break
        categories.append(dict(category=category, selection_subcategory=field, eligible_candidates=len(eligible),
                               attempted_accessions=len(attempted_accessions), selected_candidate=chosen["candidate_id"] if chosen else "",
                               status="selected" if chosen else "unavailable_no_verified_candidate"))
    return selected, attempts, categories


def verified_flanks(record, accession: str, hsp: dict, flank_bases: int = 50) -> dict:
    lower, upper = min(hsp["sstart"], hsp["send"]), max(hsp["sstart"], hsp["send"])
    start, end = max(1, lower - flank_bases), min(hsp["slen"], upper + flank_bases)
    sequence = validate_flank_record(record, accession, start, end)
    aligned = sequence[lower - start:upper - start + 1]
    oriented = aligned if hsp["sstart"] <= hsp["send"] else reverse_complement(aligned)
    if oriented != hsp["sseq"]:
        raise ContractError("Retrieved HSP interval disagrees with the accepted aligned subject")
    left, right = sequence[:lower - start], sequence[upper - start + 1:]
    minus = hsp["sstart"] > hsp["send"]
    return dict(flank_start=start, flank_end=end, forward_sequence=sequence,
                display_left_flank=reverse_complement(right) if minus else left,
                display_right_flank=reverse_complement(left) if minus else right,
                left_flank_clipped=start == 1 and lower - 1 < flank_bases,
                right_flank_clipped=end == hsp["slen"] and hsp["slen"] - upper < flank_bases,
                display_orientation="-" if minus else "+", hsp_subject_sequence=oriented,
                lower=lower, upper=upper, sequence_verified=True)


def retrieve_case(candidate: dict, hsp: dict, metadata: dict, prefix: Path, binary: Path,
                  binary_sha256: str, destination: Path, flank_bases: int = 50) -> dict:
    if digest(binary) != binary_sha256:
        raise ContractError("Pinned blastdbcmd executable changed")
    destination.mkdir(parents=True, exist_ok=False)
    lower, upper = min(hsp["sstart"], hsp["send"]), max(hsp["sstart"], hsp["send"])
    start, end = max(1, lower - flank_bases), min(hsp["slen"], upper + flank_bases)
    command = [str(binary), "-db", str(prefix), "-entry", candidate["accession_version"], "-target_only",
               "-strand", "plus", "-range", f"{start}-{end}", "-outfmt", "%f"]
    result = dict(status="not_checked", error="", command=command, candidate_id=candidate["candidate_id"],
                  hsp=hsp, metadata=metadata, database_prefix=str(prefix), executable_sha256=binary_sha256)
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=60,
                                 env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"))
        (destination / "retrieved.fasta").write_text(process.stdout, encoding="utf-8")
        (destination / "stderr.txt").write_text(process.stderr, encoding="utf-8")
        result.update(exit_code=process.returncode)
        if process.returncode:
            raise ContractError("Local exact-version retrieval failed: " + process.stderr.strip())
        record = SeqIO.read(io.StringIO(process.stdout), "fasta")
        result.update(verified_flanks(record, candidate["accession_version"], hsp, flank_bases))
        result["status"] = "sequence_verified_not_biologically_adjudicated"
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        result.update(status="retrieval_or_verification_failed", error=str(error))
    result["files"] = {path.name: digest(path) for path in destination.iterdir() if path.is_file()}
    save_json(destination / "case_evidence.json", result)
    return result