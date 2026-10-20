"""Flank census of retained strict reference hits (protocol v2, descriptive).

For declared query families every retained strict HSP window is retrieved from
the local database with fixed flanks, checked to contain the physical query
exactly at the HSP coordinates, oriented to the family's declared orientation
and classified by exact flank rules (secondary: at most one mismatch per flank
element). Unresolved identifiers and failed retrievals are listed, never
dropped. Counts are windows and accessions; the census assigns context, not
biological origin, and adds no search.
"""
import json
import os
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .contracts import ContractError, read_json, read_table, table_schema, write_table
from .design import digest
from .panels import LOGICAL_SCHEMA
from .reanalysis import ANNOTATED_EVENT_SCHEMA
from .sequences import read_fasta, reverse_complement

FLANK = 50
K46_FRONT, K40_FRONT, NATIVE_REAR = "ATTGCTAAGGTTAA", "AAGGTTAA", "CAGCACCT"
RAPID_FRONT, RAPID_REAR_START = "GCTTGGGTGTTTAACC", "GTTTTCGCATTTATCG"
NEXTERA_A, NEXTERA_B = "TCGTCGGCAGCGTC", "GTCTCGTGGGCTCGG"
MU_PHAGE_TAXID = 10677
FAMILIES = {
    "ont_native_barcode_24": dict(conditions={"ont_barcode"}, orientation="query"),
    "ont_rapid_rear_50": dict(physical_ids={"RAPID_adapter_forward"}, orientation="query"),
    "pacbio_smrtbell_45": dict(physical_ids={"SMRTbell_blunt_adapter_forward"}, orientation="query"),
    "illumina_tn5_19": dict(physical_ids={"ILLUMINA_NEXTERA_core_motif"}, orientation="reverse_complement",
                            note="classified in the mosaic-end orientation AGATGTGTATAAGAGACAG"),
}
WINDOW_SCHEMA = table_schema(("family", "physical_id", "accession_version", "strand", "status", "census_class", "census_class_hamming1",
                              "upstream", "downstream", "taxids", "testability"),
                             ("lower", "upper", "record_length", "query_length", "upstream_length", "downstream_length"))
SUMMARY_SCHEMA = table_schema(("family", "rule", "census_class"), ("windows", "accessions", "family_windows", "family_accessions"))
CATALOGUE_SCHEMA = table_schema(("set", "name", "sequence"))


def matches(observed: str, expected: str, mismatches: int) -> bool:
    return len(observed) == len(expected) and sum(a != b for a, b in zip(observed, expected)) <= mismatches


def in_catalogue(observed: str, catalogue: set[str], mismatches: int) -> bool:
    if mismatches == 0:
        return observed in catalogue
    return any(matches(observed, item, mismatches) for item in catalogue if len(item) == len(observed))


def classify(family: str, up: str, down: str, catalogues: dict, mu: bool, mismatches: int = 0) -> str:
    def tail(n):
        return up[-n:] if len(up) >= n else None

    def head(n):
        return down[:n] if len(down) >= n else None

    def eq(observed, expected):
        return observed is not None and matches(observed, expected, mismatches)

    if family == "ont_native_barcode_24":
        if eq(tail(14), K46_FRONT) and eq(head(8), NATIVE_REAR):
            return "K46"
        if eq(tail(8), K40_FRONT) and eq(head(8), NATIVE_REAR):
            return "K40_only"
        if eq(tail(16), RAPID_FRONT) and eq(head(16), RAPID_REAR_START):
            return "rapid_context"
        if eq(tail(8), K40_FRONT):
            return "native_front_only"
        if eq(head(8), NATIVE_REAR):
            return "native_rear_only"
        return "insufficient_flank" if tail(16) is None or head(16) is None else "none"
    if family == "ont_rapid_rear_50":
        barcodes = catalogues.get("ont_barcodes", set())
        barcode = tail(24)
        listed = barcode is not None and (in_catalogue(barcode, barcodes, mismatches) or in_catalogue(reverse_complement(barcode), barcodes, mismatches))
        if listed and len(up) >= 40 and matches(up[-40:-24], RAPID_FRONT, mismatches):
            return "rapid_front_and_barcode"
        if listed:
            return "catalogued_barcode_only"
        if mu:
            return "phage_mu_record"
        return "insufficient_flank" if tail(40) is None else "none"
    if family == "pacbio_smrtbell_45":
        barcodes = catalogues.get("pacbio_barcodes", set())
        for observed in (tail(16), head(16)):
            if observed is not None and (in_catalogue(observed, barcodes, mismatches) or in_catalogue(reverse_complement(observed), barcodes, mismatches)):
                return "barcode_adjacent"
        return "insufficient_flank" if tail(16) is None or head(16) is None else "none"
    if family == "illumina_tn5_19":
        if eq(tail(14), NEXTERA_A):
            return "nextera_A"
        if eq(tail(15), NEXTERA_B):
            return "nextera_B"
        return "insufficient_flank" if tail(15) is None else "none"
    raise ContractError("Unknown census family: " + family)


def mu_lineage(taxids: str, taxonomy) -> bool:
    if taxonomy is None:
        return False
    for raw in taxids.replace(",", ";").split(";"):
        if raw.strip().isdigit():
            resolved, _ = taxonomy.resolve(int(raw))
            if resolved is not None and MU_PHAGE_TAXID in taxonomy.lineage(resolved):
                return True
    return False


def window_requests(analysis: Path, panel: Path) -> list[dict]:
    report = read_json(analysis / "analysis.json")
    for name, checksum in report["files"].items():
        if digest(analysis / name) != checksum:
            raise ContractError("Analysis artifact changed: " + name)
    definition = read_json(panel / "panel.json")
    for name, checksum in definition["files"].items():
        if digest(panel / name) != checksum:
            raise ContractError("Panel artifact changed: " + name)
    physical = read_fasta(panel / "physical_queries.fasta")
    logical = read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)
    members = defaultdict(list)
    for row in logical:
        if row["role"] != "real":
            continue
        for family, rule in FAMILIES.items():
            if row["condition"] in rule.get("conditions", set()) or row["logical_id"] in rule.get("physical_ids", set()):
                members[row["sequence_id"]].append((family, row))
    policy = report["analysis_policy"]
    testability = {row["query_length"]: row["strict_testability"] for row in policy.get("detectability", [])}
    requests, seen = [], set()
    for event in read_table(analysis / "real_events.tsv", ANNOTATED_EVENT_SCHEMA):
        if event["strict"] != "yes":
            continue
        for family, row in members.get(event["query_id"], []):
            lower, upper = min(event["sstart"], event["send"]), max(event["sstart"], event["send"])
            canonical_plus = event["sstart"] <= event["send"]
            physical_plus = canonical_plus == (row["orientation"] == "+")
            key = family, row["logical_id"], event["accession_version"], lower, upper
            if key in seen:
                continue
            seen.add(key)
            sequence = physical[row["logical_id"]]
            requests.append(dict(family=family, physical_id=row["logical_id"], query=sequence, accession_version=event["accession_version"],
                                 lower=lower, upper=upper, record_length=event["slen"], strand="+" if physical_plus else "-",
                                 taxids=event.get("taxids", ""), testability=testability.get(len(sequence), "not_calibrated")))
    return requests


def retrieve(requests: list[dict], prefix: Path, blastdbcmd: str, scratch: Path, batch: int = 2000) -> None:
    resolvable = [item for item in requests if not item["accession_version"].startswith("unresolved:")]
    for item in requests:
        if item["accession_version"].startswith("unresolved:"):
            item.update(status="unresolved_identifier", window="")
    for offset in range(0, len(resolvable), batch):
        chunk = resolvable[offset:offset + batch]
        entries = scratch / f"entries_{offset:07d}.txt"
        entries.write_text("".join(f"{item['accession_version']} {max(1, item['lower'] - FLANK)}-"
                                   f"{min(item['record_length'], item['upper'] + FLANK)} plus\n" for item in chunk), encoding="ascii")
        process = subprocess.run([blastdbcmd, "-db", str(prefix), "-entry_batch", str(entries), "-target_only", "-outfmt", "%a %s"],
                                 capture_output=True, text=True, env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"))
        lines = [line for line in process.stdout.splitlines() if line.strip()]
        if process.returncode or len(lines) != len(chunk):
            for item in chunk:
                item.setdefault("status", "retrieval_failed")
                item.setdefault("window", "")
            raise ContractError(f"blastdbcmd batch retrieval returned {len(lines)} of {len(chunk)} windows: {process.stderr.strip()[:300]}")
        for item, line in zip(chunk, lines):
            accession, _, window = line.partition(" ")
            if accession.split(".")[0] != item["accession_version"].split(".")[0]:
                raise ContractError("Batch retrieval order or identity mismatch for " + item["accession_version"])
            item.update(window=window.upper(), status="retrieved")


def orient_and_classify(item: dict, catalogues: dict, taxonomy) -> dict:
    if item["status"] != "retrieved":
        return dict(item, census_class=item["status"], census_class_hamming1=item["status"], upstream="", downstream="")
    start = max(1, item["lower"] - FLANK)
    left = item["window"][:item["lower"] - start]
    core = item["window"][item["lower"] - start:item["upper"] - start + 1]
    right = item["window"][item["upper"] - start + 1:]
    if item["strand"] == "-":
        left, core, right = reverse_complement(right), reverse_complement(core), reverse_complement(left)
    if core != item["query"]:
        return dict(item, status="integrity_failed", census_class="integrity_failed", census_class_hamming1="integrity_failed",
                    upstream=left, downstream=right)
    if FAMILIES[item["family"]]["orientation"] == "reverse_complement":
        left, right = reverse_complement(right), reverse_complement(left)
    mu = mu_lineage(item["taxids"], taxonomy)
    return dict(item, status="verified", upstream=left, downstream=right,
                census_class=classify(item["family"], left, right, catalogues, mu, 0),
                census_class_hamming1=classify(item["family"], left, right, catalogues, mu, 1))


def run_census(analysis: Path, panel: Path, prefix: Path, blastdbcmd: str, catalogue: Path, destination: Path,
               taxonomy_directory: Path | None = None) -> dict:
    catalogues = defaultdict(set)
    for row in read_table(catalogue, CATALOGUE_SCHEMA):
        catalogues[row["set"]].add(row["sequence"].upper())
    taxonomy = None
    if taxonomy_directory is not None:
        from .taxonomy import load_taxonomy
        taxonomy = load_taxonomy(taxonomy_directory)[0]
    requests = window_requests(analysis, panel)
    destination.mkdir(parents=True, exist_ok=False)
    scratch = destination / "retrieval"
    scratch.mkdir()
    retrieve(requests, prefix, blastdbcmd, scratch)
    rows = [orient_and_classify(item, catalogues, taxonomy) for item in requests]
    windows = [dict(family=r["family"], physical_id=r["physical_id"], accession_version=r["accession_version"], strand=r["strand"],
                    status=r["status"], census_class=r["census_class"], census_class_hamming1=r["census_class_hamming1"],
                    upstream=r["upstream"], downstream=r["downstream"], taxids=r["taxids"], testability=r["testability"],
                    lower=r["lower"], upper=r["upper"], record_length=r["record_length"], query_length=len(r["query"]),
                    upstream_length=len(r["upstream"]), downstream_length=len(r["downstream"])) for r in rows]
    windows.sort(key=lambda r: (r["family"], r["physical_id"], r["accession_version"], r["lower"]))
    write_table(destination / "census_windows.tsv", windows, WINDOW_SCHEMA)
    summary = []
    for family in sorted({row["family"] for row in windows}):
        selected = [row for row in windows if row["family"] == family]
        family_accessions = len({row["accession_version"] for row in selected})
        for rule, field in (("exact", "census_class"), ("hamming1", "census_class_hamming1")):
            counts, accessions = Counter(), defaultdict(set)
            for row in selected:
                counts[row[field]] += 1
                accessions[row[field]].add(row["accession_version"])
            summary.extend(dict(family=family, rule=rule, census_class=name, windows=count, accessions=len(accessions[name]),
                                family_windows=len(selected), family_accessions=family_accessions) for name, count in sorted(counts.items()))
    write_table(destination / "census_summary.tsv", summary, SUMMARY_SCHEMA)
    report = dict(kind="flank_census_v1", created_utc=datetime.now(timezone.utc).isoformat(), analysis=str(analysis.resolve()),
                  analysis_sha256=digest(analysis / "analysis.json"), panel_sha256=digest(panel / "panel.json"),
                  database_prefix=str(prefix), blastdbcmd=blastdbcmd, blastdbcmd_sha256=digest(Path(blastdbcmd)),
                  catalogue_sha256=digest(catalogue), taxonomy=str(taxonomy_directory) if taxonomy_directory else None,
                  flank=FLANK, families={name: {k: sorted(v) if isinstance(v, set) else v for k, v in rule.items()} for name, rule in FAMILIES.items()},
                  windows=len(windows), status_counts=dict(Counter(row["status"] for row in windows)),
                  integrity_complete=all(row["status"] in {"verified", "unresolved_identifier"} for row in windows),
                  code_sha256=digest(Path(__file__)),
                  files={name: digest(destination / name) for name in ("census_windows.tsv", "census_summary.tsv")},
                  interpretation="Descriptive context of retained strict hits; classes are exact flank rules, not proof of origin")
    (destination / "census.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report
