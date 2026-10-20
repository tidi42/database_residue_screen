"""Exhaustive exact occurrence scan of a nucleotide BLAST database (cap- and heuristic-free).

Every record of every volume is decoded with blastdbcmd and matched against all
query sequences of at least a declared length on both strands with an
Aho-Corasick automaton that reports overlapping and nested occurrences. The scan
is an independent second implementation of the strict endpoint: it does not use
BLAST seeds, E-values or target caps. Ambiguity codes never match. Per-volume
record and base totals must sum to the database dimensions in the hashed
database evidence, otherwise the scan is not accepted.
"""
import csv
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .contracts import ContractError, read_json, table_schema, write_table
from .design import digest
from .sequences import read_fasta

OCCURRENCE_SCHEMA = table_schema(("query_id", "accession_version", "strand", "volume"), ("start", "end", "record_length", "query_length"))
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SCAN_FORMAT = "%a %s"


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


def scan_lines(lines, machine, volume: str):
    """Yield occurrence rows and finally a (records, bases) total from 'accession sequence' lines."""
    records = bases = 0
    for line in lines:
        accession, _, sequence = line.rstrip("\n").partition(" ")
        if not accession or not sequence:
            raise ContractError("Malformed blastdbcmd record line in " + volume)
        sequence = sequence.upper()
        records += 1
        bases += len(sequence)
        for end, (length, owners) in machine.iter(sequence):
            for name, strand in owners:
                yield dict(query_id=name, accession_version=accession, strand=strand, volume=volume, start=end - length + 2,
                           end=end + 1, record_length=len(sequence), query_length=length)
    yield records, bases


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
    write_table(output, rows, OCCURRENCE_SCHEMA)
    return dict(volume=volume, records=totals[0], bases=totals[1], occurrences=len(rows), file=output.name, sha256=digest(output))


def run_exact_scan(queries_fasta: Path, database_prefix: Path, evidence: Path, destination: Path, blastdbcmd: str,
                   workers: int = 8, minimum_length: int = 23) -> dict:
    if workers < 1 or minimum_length < 1:
        raise ContractError("Exact scan needs positive workers and minimum length")
    runtime = read_json(evidence)["runtime_metadata"]
    queries = read_fasta(queries_fasta)
    patterns = scan_patterns(queries, minimum_length)
    if not patterns:
        raise ContractError("No query reaches the declared minimum scan length")
    destination.mkdir(parents=True, exist_ok=False)
    directory = database_prefix.parent
    jobs = [(blastdbcmd, str(directory / name), patterns, str(destination)) for name in runtime["volume_names"]]
    started = datetime.now(timezone.utc).isoformat()
    with ProcessPoolExecutor(workers) as pool:
        volumes = sorted(pool.map(scan_volume, jobs), key=lambda row: row["volume"])
    records, bases = sum(row["records"] for row in volumes), sum(row["bases"] for row in volumes)
    accepted = records == runtime["sequence_count"] and bases == runtime["base_count"]
    report = dict(kind="exact_occurrence_scan_v1", started_utc=started, completed_utc=datetime.now(timezone.utc).isoformat(),
                  queries_sha256=digest(queries_fasta), database_prefix=str(database_prefix), evidence_sha256=digest(evidence),
                  blastdbcmd=blastdbcmd, blastdbcmd_sha256=digest(Path(blastdbcmd)), minimum_length=minimum_length,
                  scanned_queries=len({name for owners in patterns.values() for name, _ in owners}),
                  patterns=len(patterns), records=records, bases=bases, expected_records=runtime["sequence_count"],
                  expected_bases=runtime["base_count"], totals_match_database=accepted,
                  occurrences=sum(row["occurrences"] for row in volumes), volumes=volumes, code_sha256=digest(Path(__file__)),
                  interpretation="All exact occurrences on both strands; independent of BLAST seeds, E-values and target caps")
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
                rows.append(dict(row, **{key: int(row[key]) for key in ("start", "end", "record_length", "query_length")}))
    return rows
