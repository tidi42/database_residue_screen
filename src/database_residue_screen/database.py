import os
import re
import subprocess
from datetime import datetime
from pathlib import Path

from .contracts import ContractError, read_json
from .design import digest
from .example import save_json


def parse_database_info(text: str) -> dict:
    counts = re.search(r"([0-9,]+) sequences;\s*([0-9,]+) total bases", text)
    date = re.search(r"Date:\s*([A-Za-z]+\s+\d+,\s+\d{4})", text)
    version = re.search(r"BLASTDB Version:\s*(\d+)", text)
    if not counts or not date or not version or "Volumes:" not in text:
        raise ContractError("Cannot parse complete blastdbcmd database metadata")
    volumes = [line.strip() for line in text.split("Volumes:", 1)[1].splitlines() if line.strip()]
    if not volumes:
        raise ContractError("No database volumes reported")
    return dict(sequence_count=int(counts.group(1).replace(",", "")), base_count=int(counts.group(2).replace(",", "")),
                build_date=datetime.strptime(date.group(1), "%b %d, %Y").date().isoformat(),
                database_version=int(version.group(1)), volume_count=len(volumes), volume_names=[Path(value).name for value in volumes])


def metadata_conflicts(runtime: dict, metadata: dict) -> list[dict]:
    mappings = (("sequence_count", "number-of-sequences"), ("base_count", "number-of-letters"),
                ("database_version", "db-version"), ("volume_count", "number-of-volumes"))
    result = [dict(field=left, runtime=runtime.get(left), json_metadata=metadata.get(right))
              for left, right in mappings if runtime.get(left) != metadata.get(right)]
    metadata_date = str(metadata.get("last-updated", ""))[:10]
    if runtime["build_date"] != metadata_date:
        result.append(dict(field="build_date", runtime=runtime["build_date"], json_metadata=metadata_date))
    return result


def verify_inventory(report: dict, prefix: Path, hash_content: bool = True) -> dict:
    prefix = prefix.resolve()
    if Path(report["database_prefix"]).resolve() != prefix or report["missing_required_volume_files"] or not report["files"]:
        raise ContractError("Database evidence has the wrong prefix, missing volumes or no inventory")
    names = {record["name"] for record in report["files"]}
    observed_names = {path.name for path in prefix.parent.glob(prefix.name + ".*") if path.is_file() and path.suffix not in {".gz", ".md5"}}
    if len(names) != len(report["files"]) or names != observed_names:
        raise ContractError("Database file membership changed since the inventory")
    for record in report["files"]:
        path = prefix.parent / record["name"]
        if Path(record["name"]).name != record["name"] or path.is_symlink() or path.parent.resolve() != prefix.parent:
            raise ContractError("Unsafe database inventory path")
        observed = path.stat()
        if observed.st_size != record["bytes"] or observed.st_mtime_ns != record["mtime_ns"]:
            raise ContractError("Database file changed since the inventory: " + path.name)
        if hash_content and record["sha256"] and digest(path) != record["sha256"]:
            raise ContractError("Database content changed since the inventory: " + path.name)
    return report


def inspect_database(prefix: Path, destination: Path, executable: str = "blastdbcmd", hash_all: bool = False,
                     study_root: Path | None = None) -> dict:
    if study_root is not None and not destination.resolve().is_relative_to(study_root.resolve()):
        raise ContractError("Database evidence must remain inside the study directory")
    prefix = prefix.resolve()
    process = subprocess.run([executable, "-db", str(prefix), "-info"], capture_output=True, text=True, check=True,
                             env=dict(os.environ, LC_ALL="C"), timeout=120)
    version = subprocess.run([executable, "-version"], capture_output=True, text=True, check=True, timeout=30).stdout.strip()
    runtime = parse_database_info(process.stdout)
    json_path = prefix.with_suffix(".njs")
    metadata = read_json(json_path) if json_path.is_file() else {}
    conflicts = metadata_conflicts(runtime, metadata)
    candidates = sorted(path for path in prefix.parent.glob(prefix.name + ".*")
                        if path.suffix not in {".gz", ".md5"} and path.is_file())
    files = []
    for path in candidates:
        before = path.stat()
        checksum = digest(path) if hash_all or before.st_size <= 1048576 else None
        after = path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            raise ContractError("Database file changed while inspecting: " + path.name)
        files.append(dict(name=path.name, bytes=before.st_size, mtime_ns=before.st_mtime_ns, sha256=checksum))
    missing = [volume + extension for volume in runtime["volume_names"] for extension in (".nin", ".nsq", ".nhr")
               if not (prefix.parent / (volume + extension)).is_file()]
    report = dict(schema_version="1.0.0", database_prefix=str(prefix), command=[executable, "-db", str(prefix), "-info"],
                  executable_version=version, runtime_metadata=runtime, supplied_json_metadata=metadata,
                  metadata_conflicts=conflicts, metadata_consistent=not conflicts, missing_required_volume_files=missing,
                  files=files, total_bytes=sum(row["bytes"] for row in files),
                  all_database_files_hashed=bool(files) and not missing and all(row["sha256"] for row in files),
                  internal_integrity_check_performed=False, identity_verified=False,
                  meaning="Read-only inventory; partial hashes and metadata consistency are not a complete database-integrity test")
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / "database_identity.json", report)
    with (destination / "blastdbcmd_info.txt").open("x", encoding="utf-8") as handle:
        handle.write(process.stdout)
    save_json(destination / "database.json", dict(name=prefix.name, kind="blast", build_date=runtime["build_date"], blast_version=version,
                                                 identity_status="unverified", evidence_file=str((destination.resolve() / "database_identity.json").relative_to(study_root.resolve())) if study_root else "database_identity.json",
                                                 evidence_sha256=digest(destination / "database_identity.json"),
                                                 sequence_count=runtime["sequence_count"], base_count=runtime["base_count"]))
    return report
