import json
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .contracts import ContractError, read_json, read_table, schema, write_table
from .design import digest
from .eligibility import annotate_record, read_records, record_schema
from .example import save_json
from .panels import LOGICAL_SCHEMA
from .policy import load_policy
from .reanalysis import ANNOTATED_EVENT_SCHEMA, derive_event_tables, read_events, write_completed_analysis
from .search import STATUS_SCHEMA


SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def summary_record(record: dict, requested: set[str], date: str) -> dict | None:
    accession = record.get("accessionversion", "")
    if accession not in requested:
        return None
    topology = record.get("topology", "unknown")
    topology = topology if topology in {"linear", "circular"} else "unknown"
    biomol = record.get("biomol", "").lower()
    record_type = "genomic" if biomol == "genomic" else "transcript" if biomol in {"mrna", "rrna", "trna", "ncrna", "transcribed rna"} else "unknown"
    moltype = record.get("moltype", "").lower()
    molecule = "DNA" if moltype in {"dna", "genomic dna"} else "RNA" if moltype in {"rna", "genomic rna", "mrna", "rrna", "trna", "ncrna"} else "unknown"
    if record_type == "transcript":
        molecule = "conflicting" if molecule == "DNA" else "RNA"
    origin = {32630: "synthetic", 29278: "vector", 81077: "artificial_other"}.get(record.get("taxid"), "unknown")
    if origin in {"synthetic", "vector"}:
        record_type = origin
    return dict(raw_id=accession, accession_version=accession, length=int(record["slen"]), topology=topology,
                record_type=record_type, taxids=str(record.get("taxid", "")), assembly_id="", project_id=str(record.get("projectid", "")),
                molecule=molecule, origin=origin,
                targeted="unknown", organelle="unknown", context="unknown",
                metadata_source="NCBI nuccore ESummary accession-version exact match", metadata_date=date)



def summary_with_retry(request, attempts=3):
    """Bounded retries; each failed attempt is retained as evidence."""
    failures = []
    email = os.environ.get("NCBI_EMAIL")
    if email:
        request = Request(request.full_url + "&" + urlencode(dict(email=email)), headers=dict(request.header_items()))
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=60) as response:
                payload = json.load(response)
            if not isinstance(payload.get("result"), dict) or payload.get("error"):
                raise ValueError("ESummary did not return a valid result")
            return payload, failures
        except (OSError, ValueError) as error:
            failures.append(dict(attempt=attempt, reason="retrieval_attempt_failed", error=str(error),
                                 retrieved_utc=datetime.now(timezone.utc).isoformat()))
            if attempt == attempts:
                # Caller archives the terminal failure separately; attach full audit.
                failure = OSError("ESummary retries exhausted")
                failure.attempts = failures
                raise failure from error
            time.sleep(0.4 * 2**(attempt - 1))


def verify_metadata_archive(directory: Path, depth: int = 0) -> dict:
    """Verify an archive and, for a chained cache, every upstream link by hash (bounded depth)."""
    evidence = read_json(directory / "metadata.json")
    if (evidence.get("record_contract") != "explicit_molecule_v1" or
            digest(directory / "subjects.tsv") != evidence["subject_table_sha256"]):
        raise ContractError("Cached metadata contract or subject table changed")
    for batch in evidence["batches"]:
        if Path(batch["file"]).name != batch["file"] or digest(directory / batch["file"]) != batch["sha256"]:
            raise ContractError("Cached metadata response changed")
        if batch.get("reused_from") and digest(Path(batch["reused_from"])) != batch["sha256"]:
            raise ContractError("Upstream cached metadata response changed: " + batch["reused_from"])
    upstream = evidence.get("reused_metadata")
    if upstream is not None:
        if depth >= 5:
            raise ContractError("Metadata cache chain is too deep")
        source = Path(upstream["directory"])
        if digest(source / "metadata.json") != upstream["manifest_sha256"]:
            raise ContractError("Upstream metadata manifest changed")
        verify_metadata_archive(source, depth + 1)
    return evidence


def fetch_metadata(analysis: Path, destination: Path, batch_size: int = 100, reuse_metadata: Path | None = None) -> dict:
    if not 1 <= batch_size <= 100:
        raise ContractError("Metadata batches must contain 1 to 100 accessions")
    source = read_json(analysis / "analysis.json")
    if digest(analysis / "events.tsv") != source["files"]["events.tsv"]:
        raise ContractError("Source events changed since analysis")
    events = read_events(analysis, source)
    expected_lengths = {}
    for row in events:
        accession = row["accession_version"]
        if accession.startswith("unresolved:"):
            continue
        if accession in expected_lengths and expected_lengths[accession] != row["slen"]:
            raise ContractError("One exact accession has inconsistent BLAST sequence lengths")
        expected_lengths[accession] = row["slen"]
    requested = sorted(expected_lengths)
    destination.mkdir(parents=True, exist_ok=False)
    today = datetime.now(timezone.utc).date().isoformat()
    rows, problems, batches = {}, [], []
    reused = None
    if reuse_metadata is not None:
        evidence = verify_metadata_archive(reuse_metadata)
        for index, batch in enumerate(evidence["batches"], 1):
            target = destination / f"cached_esummary_{index:04d}.json"
            shutil.copyfile(reuse_metadata / batch["file"], target)
            batches.append(dict(batch, file=target.name, reused_from=str((reuse_metadata / batch["file"]).resolve())))
        for row in read_records(reuse_metadata / "subjects.tsv"):
            accession = row["accession_version"]
            if accession in expected_lengths and row["length"] == expected_lengths[accession]:
                if accession in rows:
                    raise ContractError("Duplicate exact accession in cached metadata")
                rows[accession] = row
        shutil.copyfile(reuse_metadata / "metadata.json", destination / "source_metadata.json")
        reused = dict(directory=str(reuse_metadata.resolve()), manifest_sha256=digest(reuse_metadata / "metadata.json"),
                      matching_records=len(rows), source_manifest="source_metadata.json")
    missing = sorted(set(requested) - rows.keys())
    previous = 0.0
    for offset in range(0, len(missing), batch_size):
        identifiers = set(missing[offset:offset + batch_size])
        url = SUMMARY_URL + "?" + urlencode(dict(db="nuccore", id=",".join(sorted(identifiers)), retmode="json", tool="database_residue_screen"))
        delay = max(0.0, 0.4 - (time.monotonic() - previous))
        if delay:
            time.sleep(delay)
        previous = time.monotonic()
        batch_path = destination / f"esummary_{offset // batch_size + 1:04d}.json"
        try:
            request = Request(url, headers={"User-Agent": "database-residue-screen/0.1 research-metadata"})
            payload, attempts = summary_with_retry(request)
            problems.extend(dict(requested_batch=sorted(identifiers), **attempt) for attempt in attempts)
            save_json(batch_path, payload)
            batch_info = dict(file=batch_path.name, sha256=digest(batch_path), requested=sorted(identifiers), url=url,
                              retrieved_utc=datetime.now(timezone.utc).isoformat())
            batches.append(batch_info)
            result = payload.get("result", {})
            if not isinstance(result, dict):
                raise ContractError("Unexpected ESummary result structure")
            for identifier in result.get("uids", []):
                record = result[identifier]
                observed = summary_record(record, identifiers, today)
                if observed is None:
                    returned = record.get("accessionversion", "")
                    changed_versions = [accession for accession in identifiers if returned.rsplit(".", 1)[0] == accession.rsplit(".", 1)[0] and returned != accession]
                    if changed_versions:
                        problems.extend(dict(accession=accession, returned=returned, reason="version_mismatch") for accession in changed_versions)
                    else:
                        problems.append(dict(requested_batch=sorted(identifiers), returned=returned, reason="version_changed_or_unrequested"))
                elif observed["length"] != expected_lengths[observed["accession_version"]]:
                    problems.append(dict(accession=observed["accession_version"], reason="length_conflict", blast=expected_lengths[observed["accession_version"]], metadata=observed["length"]))
                elif observed["accession_version"] in rows:
                    raise ContractError("Duplicate exact accession returned by metadata service")
                else:
                    rows[observed["accession_version"]] = observed
        except (OSError, ValueError, KeyError) as error:
            problems.extend(dict(requested_batch=sorted(identifiers), **attempt) for attempt in getattr(error, "attempts", []))
            problems.append(dict(requested_batch=sorted(identifiers), reason="retrieval_failed", error=str(error)))
        print(f"Metadata: {min(offset + batch_size, len(missing))}/{len(missing)} uncached requested, {len(rows)} exact-version records matched", flush=True)
    write_table(destination / "subjects.tsv", list(rows.values()), record_schema())
    report = dict(source_events_sha256=digest(analysis / "events.tsv"), requested_accessions=len(requested),
                  matched_exact_version_records=len(rows), missing_accessions=sorted(set(requested) - set(rows)), problems=problems,
                  batches=batches, subject_table_sha256=digest(destination / "subjects.tsv"),
                  record_contract="explicit_molecule_v1", classification="Explicit ESummary moltype/biomol/topology; separate origin and unknown context; no title fallback",
                  limitation="Genomic biomol and linear topology do not establish organismal provenance, exclude every synthetic record, or authenticate a residue. No dated taxonomy lineage reconciliation performed.")
    if reused is not None:
        report["reused_metadata"] = reused
    save_json(destination / "metadata.json", report)
    return report


def annotate_events(analysis: Path, metadata: Path, destination: Path) -> dict:
    source, evidence = read_json(analysis / "analysis.json"), read_json(metadata / "metadata.json")
    for name, checksum in source["files"].items():
        if Path(name).name != name or digest(analysis / name) != checksum:
            raise ContractError("Source analysis artifact changed")
    if source.get("source_search") and (Path(source["source_search"]) / "quarantine.json").exists():
        raise ContractError("Source search is quarantined")
    if (digest(analysis / "events.tsv") != evidence["source_events_sha256"]
            or digest(metadata / "subjects.tsv") != evidence["subject_table_sha256"]):
        raise ContractError("Metadata does not belong to the unchanged event release")
    if evidence.get("reused_metadata") is not None:
        reused = evidence["reused_metadata"]
        if Path(reused["source_manifest"]).name != reused["source_manifest"] or digest(metadata / reused["source_manifest"]) != reused["manifest_sha256"]:
            raise ContractError("Reused metadata provenance changed")
    for batch in evidence["batches"]:
        if Path(batch["file"]).name != batch["file"] or digest(metadata / batch["file"]) != batch["sha256"]:
            raise ContractError("Metadata evidence changed")
    subject_rows = read_records(metadata / "subjects.tsv")
    subjects = {row["accession_version"]: row for row in subject_rows}
    if len(subjects) != len(subject_rows):
        raise ContractError("Duplicate accession in record metadata")
    problems = {row["accession"]: row["reason"] for row in evidence.get("problems", []) if "accession" in row}
    events = [annotate_record(row, subjects.get(row["accession_version"]), problems.get(row["accession_version"], ""))
              for row in read_events(analysis, source)]
    provenance = dict(source, metadata_evidence_sha256=digest(metadata / "metadata.json"), metadata_limitations=evidence["limitation"],
                      source_analysis_sha256=digest(analysis / "analysis.json"))
    for name in ("taxonomy_snapshot", "taxonomy_ranks", "taxonomy_source_policy"):
        provenance.pop(name, None)
    provenance["taxonomy_status"] = "requires_reconciliation_after_metadata_annotation"
    if source.get("result_contract") in {"reference_screen_v1", "reference_screen_v2"}:
        logical = read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)
        statuses = read_table(analysis / "search_status.tsv", STATUS_SCHEMA)
        policy = load_policy(analysis / "analysis_policy.json")
        if policy != source["analysis_policy"]:
            raise ContractError("Analysis policy differs from its report")
        return write_completed_analysis(events, logical, statuses, policy, destination, provenance)
    destination.mkdir(parents=True, exist_ok=False)
    for name, checksum in source["files"].items():
        if name not in {"events.tsv", "density.tsv", "positions.tsv", "overlap.tsv"}:
            shutil.copyfile(analysis / name, destination / name)
    write_table(destination / "events.tsv", events, ANNOTATED_EVENT_SCHEMA)
    derived = derive_event_tables(events, sorted({row["group"] for row in events}), destination)
    report = dict(provenance, event_contract="annotated_v1",
                  derived=derived, files={path.name: digest(path) for path in destination.glob("*.tsv")})
    report["stage"] = "annotated_retained_output"
    save_json(destination / "analysis.json", report)
    return report