import hashlib
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .contracts import ContractError


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def canonical(sequence: str) -> str:
    return min(sequence, reverse_complement(sequence))


def canonical_id(sequence: str) -> str:
    return "seq_" + hashlib.sha256(canonical(sequence).encode("ascii")).hexdigest()


def read_fasta(path: Path, allow_empty: bool = False, allow_ambiguity: bool = False) -> dict[str, str]:
    records = {}
    with path.open(encoding="utf-8") as handle:
        first = next((line for line in handle if line.strip()), "")
        if first and not first.startswith(">"):
            raise ContractError(f"{path.name}: content before the first FASTA header")
    with path.open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in records:
                raise ContractError(f"{path.name}: duplicate FASTA ID: {record.id}")
            if not re.fullmatch(r"[A-Za-z0-9_.:+-]+", record.id):
                raise ContractError(f"{path.name}: unsupported FASTA ID: {record.id}")
            sequence = str(record.seq).upper()
            alphabet = "ACGTRYSWKMBDHVN" if allow_ambiguity else "ACGT"
            if not sequence or set(sequence) - set(alphabet):
                raise ContractError(f"{path.name}: expected nonempty {alphabet} sequence: {record.id}")
            records[record.id] = sequence
    if not records and not allow_empty:
        raise ContractError(f"{path.name}: no sequences")
    return records


def accession_version(raw_id: str) -> str:
    parts = raw_id.split("|")
    known_prefixes = {"gb", "emb", "dbj", "ref", "tpg", "tpe", "tpd"}
    values = {parts[index + 1] for index, part in enumerate(parts[:-1])
              if part in known_prefixes and re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.[0-9]+", parts[index + 1])}
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.[0-9]+", raw_id):
        values.add(raw_id)
    if len(values) != 1:
        raise ContractError(f"Unresolved or ambiguous accession.version: {raw_id}")
    return values.pop()


def positional_eligibility(record: dict) -> tuple[bool, str]:
    if record["topology"] != "linear":
        return False, "topology_" + record["topology"]
    if record["record_type"] != "genomic":
        return False, "record_type_" + record["record_type"]
    return True, "eligible_linear_genomic"
