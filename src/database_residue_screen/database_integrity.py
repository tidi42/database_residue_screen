import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .contracts import ContractError, read_json
from .database import parse_database_info
from .design import digest
from .example import save_json


def volume_evidence(prefix: Path, names: list[str], executable: str = "blastdbcmd", environment: dict | None = None) -> list[dict]:
    def inspect(name: str) -> dict:
        if Path(name).name != name:
            raise ContractError("Unsafe database volume name")
        command = [executable, "-db", str((prefix.parent / name).resolve()), "-info"]
        process = subprocess.run(command, capture_output=True, text=True, check=True, timeout=60, env=environment)
        return dict(volume=name, command=command, stdout=process.stdout, stderr=process.stderr, metadata=parse_database_info(process.stdout))

    with ThreadPoolExecutor(max_workers=4) as executor:
        return list(executor.map(inspect, names))


def assert_single_build(volumes: list[dict]) -> None:
    dates = {row["metadata"]["build_date"] for row in volumes}
    if not volumes or len(dates) != 1:
        raise ContractError("Database volumes have inconsistent build dates: " + ", ".join(sorted(dates)))


def audit_integrity(prefix: Path, evidence: Path, destination: Path) -> dict:
    recorded = read_json(evidence)
    if prefix.resolve() != Path(recorded["database_prefix"]).resolve():
        raise ContractError("Wrong database prefix for integrity evidence")
    volumes = volume_evidence(prefix, recorded["runtime_metadata"]["volume_names"])
    probes = []
    for accession in ("KC332287.1", "XM_069202981.1", "HM753700.1", "OY737857.1", "OV651819.1"):
        command = ["blastdbcmd", "-db", str(prefix.resolve()), "-entry", accession, "-outfmt", "%a\t%l\t%t"]
        process = subprocess.run(command, capture_output=True, text=True, check=False, timeout=60)
        returned = process.stdout.split("\t", 1)[0].strip()
        probes.append(dict(requested=accession, returned=returned, identity_matches=process.returncode == 0 and returned == accession,
                           command=command, exit_code=process.returncode, stdout=process.stdout, stderr=process.stderr))
    dates = Counter(row["metadata"]["build_date"] for row in volumes)
    report = dict(database_prefix=str(prefix.resolve()), source_evidence_sha256=digest(evidence),
                  volume_build_dates=dict(dates), volume_count=len(volumes), volumes=volumes, accession_probes=probes,
                  accession_lookup_failures=sum(not row["identity_matches"] for row in probes),
                  consistent_volume_build=len(dates) == 1,
                  research_search_usable=len(dates) == 1 and all(row["identity_matches"] for row in probes),
                  action="Preserve this database; obtain a coherent complete release in a new directory, verify it, then rerun real/control searches",
                  historical_scope="Retained HSP tables remain auditable observations but are not validated against this inconsistent live database",
                  new_core_nt_pilot_scope="Quarantined for scientific inference; successful execution and observed cap stability do not validate database identity")
    destination.mkdir(parents=True, exist_ok=False)
    save_json(destination / "integrity.json", report)
    return report