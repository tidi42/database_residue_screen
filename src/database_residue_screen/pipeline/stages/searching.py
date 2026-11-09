"""Resumable production searches shared by the `search` and `target_limit` stages."""
import shutil
from pathlib import Path

from ...contracts import ContractError, read_json
from ...design import digest
from ...example import save_json
from ...length_aware import SCORING_CONTRACT
from ...search import run_search, search_settings
from .resume import copy_directory, copy_file, sources_with


def check_scoring(section: dict) -> None:
    """The configured scoring must be the fixed scoring contract of the engine calibration."""
    differences = [f"{key} {section['scoring'].get(key)!r} != {value!r}" for key, value in SCORING_CONTRACT.items()
                   if section["scoring"].get(key) != value]
    if differences or set(section["scoring"]) != set(SCORING_CONTRACT):
        raise ContractError("Configured scoring differs from the scoring contract: " + "; ".join(differences or ["keys differ"]))


def database_inputs(context) -> tuple[Path, Path, Path | None]:
    """(database prefix, inventory evidence, acceptance artifact or None in working search mode) of the database stage."""
    database = context.upstream["database"].output
    evidence = database / "inventory" / "database_identity.json"
    acceptance = database / "acceptance" / "accepted.json"
    return Path(read_json(evidence)["database_prefix"]), evidence, acceptance


def restore_batches(context, destination: Path) -> dict | None:
    """Copy the completed batches of earlier unfinished attempts with the same search request into `destination`."""
    sources = sources_with(context, str(destination.relative_to(context.output)), "request.json")
    if not sources:
        return None
    destination.mkdir(parents=True)
    copy_file(sources[0] / "request.json", destination / "request.json")
    if (sources[0] / "database_evidence.json").is_file():
        copy_file(sources[0] / "database_evidence.json", destination / "database_evidence.json")
    copied, taken = [], set()
    for source in sources:
        for batch in sorted(source.glob("batch_*.attempt_*")):
            identifier = batch.name.split(".")[0]
            if batch.is_dir() and (batch / "complete.json").is_file() and identifier not in taken:
                copy_directory(batch, destination / batch.name)
                taken.add(identifier)
                copied.append(dict(batch=batch.name, source=str(batch), completion_sha256=digest(batch / "complete.json")))
    return dict(source_attempts=[str(source) for source in sources], copied_batches=copied)


def resumable_search(context, queries: Path, destination: Path, section: dict, target_limit: int) -> tuple[dict, dict | None]:
    """Run (or continue) one search with the configured settings at `target_limit`."""
    check_scoring(section)
    prefix, evidence, acceptance = database_inputs(context)
    acceptance = acceptance if section["mode"] == "production" else None
    resumed = restore_batches(context, destination) if section["resume"] else None

    def run():
        return run_search(queries, prefix, evidence, destination, search_settings(target_limit, section["evalue"]), section["threads"],
                          section["batch_size"], context.tool("blastn"), section.get("timeout_seconds"), section["mode"], acceptance,
                          None, section["threshold_mode"], None, None, section["informativeness_limit"], section["testable_chance_limit"])

    try:
        report = run()
    except ContractError as error:
        if resumed is None or "Existing search has different" not in str(error):
            raise
        shutil.rmtree(destination)
        resumed = dict(resumed, discarded="the earlier search request differs from the request of this attempt")
        context.log.warning("resume of %s discarded: %s", destination.name, resumed["discarded"])
        report = run()
    if resumed is not None:
        save_json(context.output / f"{destination.name}.resume.json", resumed)
        context.log.info("%s continued from %d completed batches of earlier attempts", destination.name, len(resumed["copied_batches"]))
    return report, resumed
