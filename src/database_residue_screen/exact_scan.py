"""Exhaustive exact occurrence scan of a nucleotide BLAST database (cap- and heuristic-free).

Every record of every volume is decoded with blastdbcmd and matched against all
query sequences of at least a declared length on both strands with an
Aho-Corasick automaton that reports overlapping and nested occurrences. The scan
is an independent second implementation of the strict endpoint: it does not use
BLAST seeds, E-values or target caps. Ambiguity codes never match. Per-volume
record and base totals must sum to the database dimensions in the hashed
database evidence, otherwise the scan is not accepted.

The unit is the database record (OID), as in BLAST. blastdbcmd prints one line
per defline, so an OID with merged deflines (identical sequences of several
accessions) appears on consecutive lines; it is scanned and counted once and
reported under its first accession, which is the subject id BLAST reports. The
further accessions are listed as aliases.

Each volume is finished atomically: its occurrence table is written under a
temporary name and renamed, then a completion record with totals and the table
hash is renamed into place. A scan that was interrupted can be resumed into the
same directory only with an identical frozen request; completed volumes are
reused after hash verification and unfinished volumes are scanned again.
"""
import csv
import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .contracts import ContractError, read_json, table_schema, write_table
from .design import digest
from .sequences import read_fasta

OCCURRENCE_SCHEMA = table_schema(("query_id", "accession_version", "alias_accessions", "strand", "volume"),
                                 ("oid", "start", "end", "record_length", "query_length"))
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SCAN_FORMAT = "%o %a %s"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def scan_patterns(queries: dict[str, str], minimum_length: int) -> dict[str, list[tuple[str, str]]]:
    """Map each searched string to (query_id, strand) owners; palindromes are recorded once on '+'."""
    patterns: dict[str, list[tuple[str, str]]] = {}
    for name, sequence in sorted(queries.items()):
        if len(sequence) < minimum_length:
            continue
        if set(sequence) - set("ACGT"):
            raise ContractError("Exact scan patterns must be unambiguous ACGT: " + name)
        patterns.setdefault(sequence, []).append((name, "+"))
        reverse = reverse_complement(sequence)
        if reverse != sequence:
            patterns.setdefault(reverse, []).append((name, "-"))
    return patterns


def automaton(patterns: dict[str, list[tuple[str, str]]]):
    import ahocorasick

    machine = ahocorasick.Automaton()
    for text, owners in patterns.items():
        machine.add_word(text, (len(text), tuple(owners)))
    machine.make_automaton()
    return machine


def _record_rows(record: list, volume: str) -> list[dict]:
    oid, accession, sequence, aliases, hits = record
    alias_text = ";".join(aliases)
    return [dict(query_id=name, accession_version=accession, alias_accessions=alias_text, oid=oid, strand=strand, volume=volume,
                 start=end - length + 2, end=end + 1, record_length=len(sequence), query_length=length) for end, length, name, strand in hits]


def scan_lines(lines, machine, volume: str):
    """Yield occurrence rows and finally a (records, bases, deflines) total from 'oid accession sequence' lines (one per defline)."""
    records = bases = deflines = 0
    record = None
    for line in lines:
        deflines += 1
        oid_text, _, rest = line.rstrip("\n").partition(" ")
        accession, _, sequence = rest.partition(" ")
        if not oid_text.isdigit() or not accession or not sequence:
            raise ContractError("Malformed blastdbcmd record line in " + volume)
        oid, sequence = int(oid_text), sequence.upper()
        if record is not None and oid == record[0]:
            if sequence != record[2]:
                raise ContractError(f"Deflines of OID {oid} in {volume} carry different sequences")
            if accession != record[1] and accession not in record[3]:
                record[3].append(accession)
            continue
        if record is not None:
            if oid < record[0]:
                raise ContractError("blastdbcmd OIDs are not ascending in " + volume)
            yield from _record_rows(record, volume)
        hits = [(end, length, name, strand) for end, (length, owners) in machine.iter(sequence) for name, strand in owners]
        record = [oid, accession, sequence, [], hits]
        records += 1
        bases += len(sequence)
    if record is not None:
        yield from _record_rows(record, volume)
    yield records, bases, deflines


def scan_volume(job: tuple) -> dict:
    blastdbcmd, volume_prefix, patterns, destination = job
    machine = automaton(patterns)
    volume = Path(volume_prefix).name
    output = Path(destination) / f"{volume}.occurrences.tsv"
    process = subprocess.Popen([blastdbcmd, "-db", volume_prefix, "-entry", "all", "-outfmt", SCAN_FORMAT],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1 << 24,
                               env={"LC_ALL": "C", "BLAST_USAGE_REPORT": "false"})
    rows, totals = [], None
    for item in scan_lines(process.stdout, machine, volume):
        if isinstance(item, tuple):
            totals = item
        else:
            rows.append(item)
    stderr = process.stderr.read()
    if process.wait() != 0 or totals is None:
        raise ContractError(f"blastdbcmd failed for {volume}: {stderr.strip()}")
    rows.sort(key=lambda row: (row["query_id"], row["accession_version"], row["start"], row["strand"]))
    partial = output.with_name(output.name + ".partial")
    write_table(partial, rows, OCCURRENCE_SCHEMA)
    os.replace(partial, output)
    record = dict(volume=volume, records=totals[0], bases=totals[1], deflines=totals[2], occurrences=len(rows), file=output.name,
                  sha256=digest(output))
    pending = Path(destination) / f"{volume}.complete.json.partial"
    pending.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    os.replace(pending, Path(destination) / f"{volume}.complete.json")
    return record


def completed_volume(destination: Path, volume: str) -> dict | None:
    path = destination / f"{Path(volume).name}.complete.json"
    if not path.is_file():
        return None
    record = read_json(path)
    if record["volume"] != Path(volume).name or digest(destination / record["file"]) != record["sha256"]:
        raise ContractError("Completed exact-scan volume changed: " + record["volume"])
    return record


def run_exact_scan(queries_fasta: Path, database_prefix: Path, evidence: Path, destination: Path, blastdbcmd: str,
                   workers: int = 8, minimum_length: int = 23, resume: bool = False) -> dict:
    if workers < 1 or minimum_length < 1:
        raise ContractError("Exact scan needs positive workers and minimum length")
    runtime = read_json(evidence)["runtime_metadata"]
    queries = read_fasta(queries_fasta)
    patterns = scan_patterns(queries, minimum_length)
    if not patterns:
        raise ContractError("No query reaches the declared minimum scan length")
    request = dict(kind="exact_occurrence_scan_request_v1", queries_sha256=digest(queries_fasta), database_prefix=str(database_prefix),
                   evidence_sha256=digest(evidence), blastdbcmd=blastdbcmd, blastdbcmd_sha256=digest(Path(blastdbcmd)),
                   minimum_length=minimum_length, patterns=len(patterns), volume_names=list(runtime["volume_names"]),
                   code_sha256=digest(Path(__file__)))
    if destination.exists():
        if not resume:
            raise ContractError("Exact-scan output exists; use a new directory or resume an identical request")
        if not (destination / "scan_request.json").is_file() or read_json(destination / "scan_request.json") != request:
            raise ContractError("Exact-scan resume requires an identical frozen request")
        if (destination / "scan.json").exists():
            raise ContractError("Exact scan is already complete")
        for leftover in destination.glob("*.partial"):
            leftover.unlink()
    else:
        destination.mkdir(parents=True, exist_ok=False)
        (destination / "scan_request.json").write_text(json.dumps(request, indent=1) + "\n", encoding="utf-8")
    directory = database_prefix.parent
    reused = {name: completed_volume(destination, name) for name in runtime["volume_names"]}
    for name, record in reused.items():
        if record is None and (destination / f"{Path(name).name}.occurrences.tsv").exists():
            (destination / f"{Path(name).name}.occurrences.tsv").unlink()
    jobs = [(blastdbcmd, str(directory / name), patterns, str(destination)) for name in runtime["volume_names"] if reused[name] is None]
    started = datetime.now(timezone.utc).isoformat()
    with ProcessPoolExecutor(workers) as pool:
        scanned = list(pool.map(scan_volume, jobs))
    volumes = sorted(scanned + [record for record in reused.values() if record is not None], key=lambda row: row["volume"])
    records, bases = sum(row["records"] for row in volumes), sum(row["bases"] for row in volumes)
    accepted = records == runtime["sequence_count"] and bases == runtime["base_count"]
    report = dict(kind="exact_occurrence_scan_v1", started_utc=started, completed_utc=datetime.now(timezone.utc).isoformat(),
                  queries_sha256=digest(queries_fasta), database_prefix=str(database_prefix), evidence_sha256=digest(evidence),
                  blastdbcmd=blastdbcmd, blastdbcmd_sha256=digest(Path(blastdbcmd)), minimum_length=minimum_length,
                  scanned_queries=len({name for owners in patterns.values() for name, _ in owners}),
                  patterns=len(patterns), records=records, bases=bases, expected_records=runtime["sequence_count"],
                  expected_bases=runtime["base_count"], totals_match_database=accepted,
                  deflines=sum(row["deflines"] for row in volumes),
                  occurrences=sum(row["occurrences"] for row in volumes), volumes=volumes, code_sha256=digest(Path(__file__)),
                  request_sha256=digest(destination / "scan_request.json"), reused_volumes=sorted(record["volume"] for record in reused.values() if record),
                  interpretation="All exact occurrences on both strands; independent of BLAST seeds, E-values and target caps. "
                                 "Records are database OIDs reported under their first accession (the BLAST subject id); "
                                 "merged-defline accessions are listed as aliases")
    (destination / "scan.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    if not accepted:
        raise ContractError("Scanned records/bases differ from the database evidence; the scan is not accepted")
    return report


def read_occurrences(scan: Path) -> list[dict]:
    report = read_json(scan / "scan.json")
    rows = []
    for volume in report["volumes"]:
        path = scan / volume["file"]
        if digest(path) != volume["sha256"]:
            raise ContractError("Exact scan output changed: " + volume["file"])
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows.append(dict(row, **{key: int(row[key]) for key in ("oid", "start", "end", "record_length", "query_length")}))
    return rows
