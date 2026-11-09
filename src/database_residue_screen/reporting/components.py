"""Components report: component inventory, within-set barcode distances and barcode-segment readout.

Segments come from the panel definition ``segments`` column (``role:start-end[=element]``, 1-based, physical
orientation). Sequences come only from the run outputs (the design's physical queries), never from this package.

* component_inventory.tsv: every declared segment of every physical query;
* within_set_distances.tsv: pairwise Hamming distances within each barcode/index set (Illumina i7 and i5 separately,
  ONT, PacBio), per length, distinct sequences only;
* segment_readout.tsv: for every tolerant HSP of a composite query (more than one segment, one of them a barcode or
  index), the subject bases aligned to the barcode/index segment compared with the queried barcode and its whole set:
  queried_exact, queried_1_mismatch, other_catalogued_exact (scaffold-preserving swap), uncatalogued,
  partial_segment or gapped; split by strict / tolerant-only and the a-priori testability of the query length.

The readout is descriptive; it does not authenticate a physical junction or a library origin.
"""
import itertools
import re
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError
from ..matching import strict_hsp, tolerant_hsp
from ..sequences import read_fasta, reverse_complement
from .common import read_tsv, write_json, write_tsv
from .raw_search import search_hsps

BARCODE_ROLES = ("barcode", "index")
SEGMENT = re.compile(r"^([a-z_]+):(\d+)-(\d+)(?:=(.+))?$")


def parse_segments(text: str) -> list[dict]:
    result = []
    for item in text.split(";"):
        match = SEGMENT.match(item)
        if not match:
            raise ContractError(f"Invalid segment declaration: {item!r}")
        result.append(dict(role=match.group(1), start=int(match.group(2)), end=int(match.group(3)), element=match.group(4) or ""))
    return result


def catalogue_set(entry: dict) -> str:
    group = entry["group"].lower()
    if group.startswith("illumina"):
        return "illumina_" + (entry.get("index_side") or "index")
    return re.sub(r"[^a-z0-9]+", "_", group).strip("_")


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def readout(observed: str, queried: str, catalogue: set[str]) -> str:
    if observed == queried:
        return "queried_exact"
    others = {sequence for sequence in catalogue if sequence != queried and len(sequence) == len(observed)}
    if observed in others or reverse_complement(observed) in others:
        return "other_catalogued_exact"
    if len(observed) == len(queried) and hamming(observed, queried) == 1:
        return "queried_1_mismatch"
    return "uncatalogued"


def definition_entries(definition: Path, physical: dict[str, str]) -> list[dict]:
    entries = []
    for row in read_tsv(definition):
        if row["query_id"] not in physical:
            continue
        sequence = physical[row["query_id"]]
        segments = parse_segments(row["segments"])
        if len(sequence) != int(row["length"]) or any(segment["end"] > len(sequence) for segment in segments):
            raise ContractError(f"{row['query_id']}: definition length or segments disagree with the searched physical query")
        entries.append(dict(row, sequence=sequence, parsed=segments))
    missing = sorted(set(physical) - {entry["query_id"] for entry in entries})
    if missing:
        raise ContractError(f"Physical queries without a panel definition: {', '.join(missing[:10])}")
    return entries


def write_components(destination: Path, *, definition: Path, physical_fasta: Path, logical: list[dict], search: Path,
                     policy: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    physical = read_fasta(physical_fasta)
    entries = definition_entries(definition, physical)
    classes = {row["query_length"]: row for row in policy.get("detectability", [])}
    inventory, catalogues, composites = [], defaultdict(set), {}
    for entry in entries:
        for segment in entry["parsed"]:
            inventory.append(dict(query_id=entry["query_id"], group=entry["group"], condition=entry["condition"], family=entry["family"],
                                  construct_class=entry["construct_class"], query_length=len(entry["sequence"]), role=segment["role"],
                                  start=segment["start"], end=segment["end"], segment_length=segment["end"] - segment["start"] + 1,
                                  element=segment["element"]))
        barcodes = [segment for segment in entry["parsed"] if segment["role"] in BARCODE_ROLES]
        for segment in barcodes:
            catalogues[catalogue_set(entry)].add(entry["sequence"][segment["start"] - 1:segment["end"]])
        if len(barcodes) == 1 and len(entry["parsed"]) > 1:
            segment = barcodes[0]
            composites[entry["query_id"]] = dict(segment, barcode=entry["sequence"][segment["start"] - 1:segment["end"]],
                                                 catalogue=catalogue_set(entry), family=entry["family"], length=len(entry["sequence"]))
    write_tsv(destination / "component_inventory.tsv", inventory)
    distances = []
    for name, members in sorted(catalogues.items()):
        by_length = defaultdict(set)
        for sequence in members:
            by_length[len(sequence)].add(sequence)
        for length, sequences in sorted(by_length.items()):
            counts = Counter(hamming(a, b) for a, b in itertools.combinations(sorted(sequences), 2))
            distances.extend(dict(set=name, length=length, members=len(sequences), hamming=distance, pairs=count)
                             for distance, count in sorted(counts.items()))
    write_tsv(destination / "within_set_distances.tsv", distances, ("set", "length", "members", "hamming", "pairs"))
    owners = defaultdict(list)
    for row in logical:
        if row["role"] == "real" and row["logical_id"] in composites:
            owners[row["sequence_id"]].append(row)
    rows = []
    for _, hsp in search_hsps(search, set(owners)):
        if not tolerant_hsp(hsp):
            continue
        strict = strict_hsp(hsp)
        for owner in sorted(owners[hsp["qseqid"]], key=lambda row: row["logical_id"]):
            segment = composites[owner["logical_id"]]
            start, end = segment["start"], segment["end"]
            if owner["orientation"] == "-":
                start, end = segment["length"] - segment["end"] + 1, segment["length"] - segment["start"] + 1
            observed = ""
            if "-" in hsp["qseq"] or "-" in hsp["sseq"]:
                state = "gapped"
            elif min(hsp["qstart"], hsp["qend"]) > start or max(hsp["qstart"], hsp["qend"]) < end:
                state = "partial_segment"
            else:
                first = min(hsp["qstart"], hsp["qend"])
                observed = hsp["sseq"][start - first:end - first + 1].upper()
                if owner["orientation"] == "-":
                    observed = reverse_complement(observed)
                state = readout(observed, segment["barcode"], catalogues[segment["catalogue"]])
            calibration = classes.get(segment["length"], {})
            rows.append(dict(physical_id=owner["logical_id"], group=owner["group"], condition=owner["condition"], family=segment["family"],
                             accession=hsp["sseqid"], endpoint_state="strict" if strict else "tolerant_only",
                             strict_testability=calibration.get("strict_testability", "not_calibrated"),
                             tolerant_testability=calibration.get("tolerant_testability", "not_calibrated"),
                             queried_barcode=segment["barcode"], observed_segment=observed, readout=state))
    rows.sort(key=lambda row: (row["physical_id"], row["accession"], row["endpoint_state"], row["readout"], row["observed_segment"]))
    fields = ("physical_id", "group", "condition", "family", "accession", "endpoint_state", "strict_testability", "tolerant_testability",
              "queried_barcode", "observed_segment", "readout")
    write_tsv(destination / "segment_readout.tsv", rows, fields)
    counts = Counter((row["group"], row["condition"], row["endpoint_state"], row["strict_testability"], row["readout"]) for row in rows)
    write_tsv(destination / "readout_summary.tsv",
              [dict(group=key[0], condition=key[1], endpoint_state=key[2], strict_testability=key[3], readout=key[4], hsps=count)
               for key, count in sorted(counts.items())],
              ("group", "condition", "endpoint_state", "strict_testability", "readout", "hsps"))
    minimum = {}
    for row in distances:
        key = f"{row['set']}|{row['length']}"
        minimum[key] = min(minimum.get(key, row["hamming"]), row["hamming"])
    summary = dict(physical_queries=len(entries), declared_segments=len(inventory), composite_queries=len(composites),
                   catalogue_sets={name: len(members) for name, members in sorted(catalogues.items())},
                   minimum_within_set_distance=minimum, readout_hsps=len(rows),
                   readout={"|".join(key): count for key, count in sorted(Counter((row["endpoint_state"], row["readout"]) for row in rows).items())})
    write_json(destination / "summary.json", summary)
    return summary
