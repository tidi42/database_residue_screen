"""Stage `exact_scan`: exhaustive exact occurrences of the search sequences on both strands in the database volumes."""
import shutil
from pathlib import Path

from ...contracts import ContractError, read_json
from ...design import digest
from ...exact_scan import run_exact_scan
from ...example import save_json
from .resume import copy_file, sources_with


def restore_volumes(context, destination: Path) -> dict | None:
    """Copy the completed volumes of earlier unfinished attempts with the same scan request into `destination`."""
    sources = sources_with(context, str(destination.relative_to(context.output)), "scan_request.json")
    if not sources:
        return None
    destination.mkdir(parents=True)
    copy_file(sources[0] / "scan_request.json", destination / "scan_request.json")
    copied, taken = [], set()
    for source in sources:
        for record_path in sorted(source.glob("*.complete.json")):
            record = read_json(record_path)
            if record["volume"] in taken:
                continue
            for name in [record_path.name, record["file"]] + ([record["accessions_file"]] if "accessions_file" in record else []):
                copy_file(source / name, destination / name)
            taken.add(record["volume"])
            copied.append(dict(volume=record["volume"], source=str(source), completion_sha256=digest(record_path)))
    return dict(source_attempts=[str(source) for source in sources], copied_volumes=copied)


def exact_scan(context):
    """Scan every volume (or the first `volume_limit` volumes of a dry run) for all search sequences of the minimum length."""
    section = context.config["exact_scan"]
    evidence = context.upstream["database"].output / "inventory" / "database_identity.json"
    identity = read_json(evidence)
    names = identity["runtime_metadata"]["volume_names"]
    volumes = names[:section["volume_limit"]] if section["volume_limit"] else None
    queries = context.upstream["design"].output / "search_panel" / "search_queries.fasta"
    destination = context.output / "scan"
    resumed = restore_volumes(context, destination) if section["resume"] else None

    def run():
        return run_exact_scan(queries, Path(identity["database_prefix"]), evidence, destination, context.tool("blastdbcmd"),
                              section["workers"], section["minimum_length"], resume=destination.exists(), volumes=volumes)

    try:
        report = run()
    except ContractError as error:
        if resumed is None or "identical request" not in str(error):
            raise
        shutil.rmtree(destination)
        resumed = dict(resumed, discarded="the earlier scan request differs from the request of this attempt")
        report = run()
    if resumed is not None:
        save_json(context.output / "scan.resume.json", resumed)
    return {key: report[key] for key in ("volume_scope", "scanned_volumes", "scanned_queries", "patterns", "records", "bases",
                                         "totals_match_scanned_volumes", "totals_match_database", "occurrences", "reused_volumes")}
