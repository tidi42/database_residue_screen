import math
import re

from .contracts import ContractError
from .sequences import reverse_complement


LEGACY_FIELDS = ("qseqid", "sseqid", "staxids", "sscinames", "pident", "length", "qlen", "qcovs", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore", "stitle")
BLAST_FIELDS = LEGACY_FIELDS + ("qseq", "sseq", "btop", "sstrand", "qcovhsp")
INTEGER_FIELDS = ("length", "qlen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "slen")


def parse_legacy_hsp(values: list[str]) -> dict:
    if len(values) != len(LEGACY_FIELDS):
        raise ContractError(f"Expected {len(LEGACY_FIELDS)} legacy BLAST columns, observed {len(values)}")
    row = dict(zip(LEGACY_FIELDS, values))
    for name in INTEGER_FIELDS:
        row[name] = int(row[name])
    for name in ("pident", "qcovs", "evalue", "bitscore"):
        row[name] = float(row[name])
        if not math.isfinite(row[name]):
            raise ContractError(f"Non-finite HSP field: {name}")
    if not 0 <= row["pident"] <= 100 or not 0 <= row["qcovs"] <= 100 or row["mismatch"] < 0 or row["gapopen"] < 0:
        raise ContractError("Invalid HSP identity, coverage or error counts")
    return row


def parse_blast_hsp(values: list[str], queries: dict[str, str]) -> dict:
    if len(values) != len(BLAST_FIELDS):
        raise ContractError(f"Expected {len(BLAST_FIELDS)} BLAST columns, observed {len(values)}")
    try:
        row = parse_legacy_hsp(values[:len(LEGACY_FIELDS)])
        row.update(zip(BLAST_FIELDS[len(LEGACY_FIELDS):], values[len(LEGACY_FIELDS):]))
        row["qcovhsp"] = float(row["qcovhsp"])
    except ValueError as error:
        raise ContractError(f"Invalid numeric HSP field: {error}") from error
    if row["qseqid"] not in queries:
        raise ContractError(f"HSP query is not in the search manifest: {row['qseqid']}")
    query = queries[row["qseqid"]]
    if row["qlen"] != len(query):
        raise ContractError("HSP query length disagrees with the submitted FASTA")
    for prefix, total in (("q", len(query)), ("s", row["slen"])):
        start, end = row[prefix + "start"], row[prefix + "end"]
        if not 1 <= min(start, end) <= max(start, end) <= total:
            raise ContractError("HSP coordinates are outside the sequence bounds")
    aligned_query, aligned_subject = row["qseq"].upper(), row["sseq"].upper()
    if (len(aligned_query) != row["length"] or len(aligned_subject) != row["length"]
            or not row["length"] or set(aligned_query + aligned_subject) - set("ACGTRYSWKMBDHVN-")):
        raise ContractError("Invalid aligned sequence alphabet or alignment length")
    expected_query = query[min(row["qstart"], row["qend"]) - 1:max(row["qstart"], row["qend"])]
    if row["qstart"] > row["qend"]:
        expected_query = reverse_complement(expected_query)
    if aligned_query.replace("-", "") != expected_query:
        raise ContractError("Aligned query disagrees with the submitted FASTA")
    if len(aligned_subject.replace("-", "")) != abs(row["send"] - row["sstart"]) + 1:
        raise ContractError("Aligned subject disagrees with the subject span")
    if (row["sstrand"] not in {"plus", "minus"}
            or row["sstart"] < row["send"] and row["sstrand"] != "plus"
            or row["sstart"] > row["send"] and row["sstrand"] != "minus"):
        raise ContractError("Subject strand disagrees with HSP coordinates")
    identities = mismatches = 0
    for query_base, subject_base in zip(aligned_query, aligned_subject):
        if query_base == subject_base == "-":
            raise ContractError("Alignment contains a double-gap column")
        if "-" not in (query_base, subject_base):
            identities += query_base == subject_base
            mismatches += query_base != subject_base
    gapopen = sum(len(re.findall(r"-+", sequence)) for sequence in (aligned_query, aligned_subject))
    if (mismatches != row["mismatch"] or gapopen != row["gapopen"]
            or abs(100 * identities / row["length"] - row["pident"]) > 0.00051):
        raise ContractError("Identity, mismatch or gap counts disagree with the alignment")
    coverage = 100 * len(expected_query) / len(query)
    if (not math.isfinite(row["qcovhsp"]) or not 0 <= row["qcovhsp"] <= 100
            or abs(coverage - row["qcovhsp"]) > 0.5001):
        raise ContractError("Same-HSP query coverage disagrees with the alignment")
    tokens = re.findall(r"[0-9]+|[A-Z-]{2}", row["btop"])
    if not tokens or "".join(tokens) != row["btop"]:
        raise ContractError("Invalid BTOP encoding")
    offset = 0
    for token in tokens:
        if token.isdigit():
            stop = offset + int(token)
            if (stop > row["length"] or "-" in aligned_query[offset:stop]
                    or aligned_query[offset:stop] != aligned_subject[offset:stop]):
                raise ContractError("BTOP identities disagree with aligned sequences")
            offset = stop
        else:
            if (offset >= row["length"] or token != aligned_query[offset] + aligned_subject[offset]
                    or token[0] == token[1]):
                raise ContractError("BTOP differences disagree with aligned sequences")
            offset += 1
    if offset != row["length"]:
        raise ContractError("BTOP does not cover the complete alignment")
    row["qseq"], row["sseq"] = aligned_query, aligned_subject
    return row


def ungapped_geometry(row: dict) -> bool:
    return (row["qlen"] > 0 and row["gapopen"] == 0
            and 1 <= min(row["qstart"], row["qend"]) <= max(row["qstart"], row["qend"]) <= row["qlen"]
            and 1 <= min(row["sstart"], row["send"]) <= max(row["sstart"], row["send"]) <= row["slen"]
            and row["length"] == abs(row["qend"] - row["qstart"]) + 1 == abs(row["send"] - row["sstart"]) + 1)


def strict_hsp(row: dict) -> bool:
    return (ungapped_geometry(row) and row["pident"] == 100 and row["mismatch"] == 0
            and min(row["qstart"], row["qend"]) == 1 and max(row["qstart"], row["qend"]) == row["qlen"]
            and row["length"] == row["qlen"])


def tolerant_hsp(row: dict) -> bool:
    return (ungapped_geometry(row) and row["pident"] >= 95
            and 10 * (abs(row["qend"] - row["qstart"]) + 1) >= 9 * row["qlen"])
