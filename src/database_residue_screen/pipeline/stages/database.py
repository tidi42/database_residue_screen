"""Stage `database`: inventory, acceptance (full integrity check, accession probes, taxonomy support) and volume audit."""
from pathlib import Path

from ...acceptance import accept_database
from ...contracts import ContractError
from ...database import inspect_database
from ...database_integrity import audit_integrity
from ...example import save_json

EXPECTED = (("sequences", "sequence_count"), ("bases", "base_count"), ("volumes", "volume_count"), ("build_date", "build_date"))


def database(context):
    """Inspect and hash the database, accept it for production searches and audit its volumes."""
    section = context.config["database"]
    prefix = Path(section["prefix"])
    blastdbcmd = context.tool("blastdbcmd")
    identity = inspect_database(prefix, context.output / "inventory", blastdbcmd, section["hash_all"])
    runtime = identity["runtime_metadata"]
    differences = [f"{name} {runtime[field]} != {section['expected'][name]}" for name, field in EXPECTED
                   if name in section["expected"] and runtime[field] != section["expected"][name]]
    if differences:
        raise ContractError("The database differs from the configured expectation: " + "; ".join(differences))
    probes = context.output / "probes.json"
    save_json(probes, dict(version="1.0.0", probes=section["probes"]))
    evidence = context.output / "inventory" / "database_identity.json"
    download = Path(section["download_evidence"]) if section.get("download_evidence") else None
    accepted = accept_database(prefix, evidence, probes, context.output / "acceptance", context.tool("blastn"), download,
                               section.get("check_timeout_seconds"), Path(section["taxonomy_support"]), blastdbcmd,
                               context.tool("blastdbcheck"))
    audit = audit_integrity(prefix, evidence, context.output / "integrity", tuple(row["accession_version"] for row in section["probes"]),
                            blastdbcmd)
    if not audit["consistent_volume_build"]:
        raise ContractError("Database volumes have inconsistent build dates")
    return dict(sequences=runtime["sequence_count"], bases=runtime["base_count"], volumes=runtime["volume_count"],
                build_date=runtime["build_date"], all_files_hashed=identity["all_database_files_hashed"], accepted=accepted["accepted"],
                accession_lookup_failures=audit["accession_lookup_failures"])
