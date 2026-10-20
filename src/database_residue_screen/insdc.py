"""Targeted/amplicon sensitivity from structured INSDC annotation (protocol v2, secondary).

The primary eligibility rule does not populate targeted/amplicon or engineered
status automatically. This module adds a declared secondary sensitivity: GenBank
headers are fetched per exact accession.version (rettype=gb, seq_start=1,
seq_stop=1, so the source feature and KEYWORDS are present without the full
sequence), archived with hashes, and a record is flagged targeted only from
structured fields: KEYWORDS TLS / Targeted Locus Study / AMPLICON / BARCODE, a
source /PCR_primers qualifier, or ESummary tech targeted/barcode. Free-text
titles are never used. The fetch is a new dated snapshot, not historical evidence.
"""
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .contracts import ContractError, read_json, table_schema, write_table
from .design import digest

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TARGETED_KEYWORDS = ("TLS", "TARGETED LOCUS STUDY", "AMPLICON", "BARCODE")
TARGETED_TECH = {"targeted", "barcode"}
ANNOTATION_SCHEMA = table_schema(("accession_version", "status", "keywords", "pcr_primers", "esummary_tech", "targeted_insdc", "evidence",
                                  "archive_file"))


def parse_genbank_headers(text: str) -> dict[str, dict]:
    """Parse KEYWORDS and source /PCR_primers per record version from GenBank flat text."""
    records = {}
    for block in re.split(r"^//\s*$", text, flags=re.M):
        version = re.search(r"^VERSION\s+(\S+)", block, re.M)
        if not version:
            continue
        keywords = re.search(r"^KEYWORDS\s+(.*?)(?=^\S)", block, re.M | re.S)
        keyword_text = " ".join(keywords.group(1).split()) if keywords else ""
        source = re.search(r"^\s{5}source\s.*?(?=^\s{5}\S|^ORIGIN|^CONTIG|\Z)", block, re.M | re.S)
        primers = bool(source and "/PCR_primers=" in source.group(0))
        records[version.group(1)] = dict(keywords=keyword_text.rstrip("."), pcr_primers=primers)
    return records


def targeted_status(record: dict | None, tech: str) -> tuple[str, str]:
    if record is None:
        return "unknown", "no_exact_version_header"
    tokens = {token.strip().upper() for token in re.split(r"[;.]", record["keywords"]) if token.strip()}
    hits = sorted(token for token in tokens if token in TARGETED_KEYWORDS)
    evidence = [f"KEYWORDS:{hit}" for hit in hits]
    if record["pcr_primers"]:
        evidence.append("source:/PCR_primers")
    if tech.lower() in TARGETED_TECH:
        evidence.append(f"esummary_tech:{tech}")
    return ("yes" if evidence else "no"), ";".join(evidence)


def fetch_batch(identifiers: list[str], attempts: int = 3) -> tuple[str, list[dict]]:
    parameters = dict(db="nuccore", id=",".join(identifiers), rettype="gb", retmode="text", seq_start=1, seq_stop=1)
    if os.environ.get("NCBI_EMAIL"):
        parameters["email"] = os.environ["NCBI_EMAIL"]
    request = Request(EFETCH + "?" + urlencode(parameters), headers={"User-Agent": "database-residue-screen"})
    failures = []
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=120) as response:
                return response.read().decode("utf-8", errors="replace"), failures
        except OSError as error:
            failures.append(dict(attempt=attempt, error=str(error), retrieved_utc=datetime.now(timezone.utc).isoformat()))
            if attempt == attempts:
                raise ContractError("EFetch retries exhausted for a batch") from error
            time.sleep(0.4 * 2 ** (attempt - 1))
    raise AssertionError("unreachable")


def fetch_insdc_annotation(accessions: list[str], technologies: dict[str, str], destination: Path, batch_size: int = 100,
                           fetch=fetch_batch, delay: float = 0.4) -> dict:
    identifiers = sorted({a for a in accessions if not a.startswith("unresolved:")})
    destination.mkdir(parents=True, exist_ok=False)
    rows, batches, problems = [], [], []
    for index, offset in enumerate(range(0, len(identifiers), batch_size)):
        chunk = identifiers[offset:offset + batch_size]
        name = f"efetch_gb_header_{index:04d}.txt"
        try:
            text, failures = fetch(chunk)
        except ContractError as error:
            problems.append(dict(batch=name, error=str(error)))
            rows.extend(dict(accession_version=a, status="retrieval_failed", keywords="", pcr_primers="unknown",
                             esummary_tech=technologies.get(a, ""), targeted_insdc="unknown", evidence="", archive_file="") for a in chunk)
            continue
        problems.extend(dict(batch=name, **failure) for failure in failures)
        (destination / name).write_text(text, encoding="utf-8")
        parsed = parse_genbank_headers(text)
        for accession in chunk:
            record = parsed.get(accession)
            status = "exact_version" if record else "version_not_returned"
            targeted, evidence = targeted_status(record, technologies.get(accession, ""))
            rows.append(dict(accession_version=accession, status=status, keywords=record["keywords"] if record else "",
                             pcr_primers=("yes" if record["pcr_primers"] else "no") if record else "unknown",
                             esummary_tech=technologies.get(accession, ""), targeted_insdc=targeted, evidence=evidence,
                             archive_file=name))
        batches.append(dict(file=name, sha256=digest(destination / name), accessions=len(chunk),
                            retrieved_utc=datetime.now(timezone.utc).isoformat()))
        time.sleep(delay)
    rows.sort(key=lambda row: row["accession_version"])
    write_table(destination / "insdc_annotation.tsv", rows, ANNOTATION_SCHEMA)
    report = dict(kind="insdc_targeted_annotation_v1", created_utc=datetime.now(timezone.utc).isoformat(), endpoint=EFETCH,
                  parameters=dict(db="nuccore", rettype="gb", retmode="text", seq_start=1, seq_stop=1), requested=len(identifiers),
                  status_counts={s: sum(r["status"] == s for r in rows) for s in sorted({r["status"] for r in rows})},
                  targeted_counts={s: sum(r["targeted_insdc"] == s for r in rows) for s in ("yes", "no", "unknown")},
                  rule=dict(keywords=list(TARGETED_KEYWORDS), source_qualifier="/PCR_primers", esummary_tech=sorted(TARGETED_TECH),
                            titles="never used"),
                  batches=batches, problems=problems, code_sha256=digest(Path(__file__)),
                  annotation_sha256=digest(destination / "insdc_annotation.tsv"),
                  scope="New dated snapshot for a declared secondary sensitivity; not historical study evidence")
    (destination / "insdc.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report


def read_annotation(directory: Path) -> dict[str, dict]:
    report = read_json(directory / "insdc.json")
    if digest(directory / "insdc_annotation.tsv") != report["annotation_sha256"]:
        raise ContractError("INSDC annotation changed since retrieval")
    from .contracts import read_table
    return {row["accession_version"]: row for row in read_table(directory / "insdc_annotation.tsv", ANNOTATION_SCHEMA)}


def technologies_from_metadata(metadata: Path) -> dict[str, str]:
    import glob
    technologies = {}
    for path in sorted(glob.glob(str(metadata / "*esummary_*.json"))):
        for value in json.loads(Path(path).read_text()).get("result", {}).values():
            if isinstance(value, dict) and value.get("accessionversion"):
                technologies[value["accessionversion"]] = str(value.get("tech") or "")
    return technologies
