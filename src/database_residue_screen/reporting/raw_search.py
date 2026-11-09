"""Verified read access to the raw HSPs of a completed search directory (batch completion digests checked)."""
from pathlib import Path

from ..contracts import ContractError, read_json
from ..search import iter_hsps
from ..sequences import read_fasta
from .common import sha256


def completed_batches(search: Path):
    """Yield (batch directory, its queries) for every completed batch listed in complete.json."""
    search = Path(search)
    complete = read_json(search / "complete.json")
    if not complete.get("all_queries_completed", True):
        raise ContractError("Search does not cover every submitted query")
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        if directory.parent.resolve() != search.resolve() or sha256(directory / "complete.json") != batch["completion_sha256"]:
            raise ContractError("Changed raw-search completion: " + batch["directory"])
        evidence = read_json(directory / "complete.json")
        for name in ("queries.fasta", "blast.tsv"):
            if sha256(directory / name) != evidence["files"][name]:
                raise ContractError("Changed raw alignment input: " + str(directory / name))
        yield directory, read_fasta(directory / "queries.fasta")


def search_hsps(search: Path, wanted=None):
    """Yield (queries, HSP row) for every HSP of the wanted search sequences (all when ``wanted`` is None)."""
    for directory, queries in completed_batches(search):
        if wanted is not None and not set(queries) & wanted:
            continue
        for row in iter_hsps(directory / "blast.tsv", queries):
            if wanted is None or row["qseqid"] in wanted:
                yield queries, row
