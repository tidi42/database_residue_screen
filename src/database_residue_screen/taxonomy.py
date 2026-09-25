import csv
import json
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

from .contracts import ContractError, read_json, read_table, table_schema
from .design import digest
from .eligibility import event_eligibility
from .example import save_json
from .outcomes import logical_memberships, project_real_events
from .panels import LOGICAL_SCHEMA
from .policy import load_policy
from .reanalysis import read_events, write_completed_analysis
from .search import STATUS_SCHEMA


TAXDUMP_FILES = ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")
TAXONOMY_SCHEMA = table_schema(("accession_version", "rank", "blast_taxids", "record_taxids", "resolved_taxids", "taxon_id", "taxon_name",
                               "status", "source_status", "context", "resolution_evidence", "snapshot_sha256"))
COMPOSITION_SCHEMA = table_schema(("group", "endpoint", "population", "rank", "category", "taxon_id", "percentage"),
                                  ("accessions", "denominator"))
UNRESOLVED_STATES = ("missing", "deleted", "unknown_id", "rank_unresolved", "multiple_ambiguous", "conflicting")


def dump_rows(path: Path, minimum_fields: int):
    with path.open(encoding="utf-8", newline="") as handle:
        for ordinal, values in enumerate(csv.reader(handle, delimiter="|", quoting=csv.QUOTE_NONE), 1):
            values = [value.strip() for value in values]
            if len(values) < minimum_fields or not values[0]:
                raise ContractError(f"Invalid taxonomy row: {path.name}:{ordinal}")
            yield values


def parse_taxids(value: str) -> tuple[int, ...]:
    if not value:
        return ()
    if not re.fullmatch(r"[0-9]+(?:;[0-9]+)*", value):
        raise ContractError("Invalid taxid list: " + value)
    return tuple(sorted({int(identifier) for identifier in value.split(";")}))


@dataclass
class Taxonomy:
    nodes: dict[int, tuple[int, str]]
    names: dict[int, str]
    merged: dict[int, int]
    deleted: set[int]

    @classmethod
    def from_directory(cls, directory: Path):
        nodes, names, merged, deleted = {}, {}, {}, set()
        for values in dump_rows(directory / "nodes.dmp", 3):
            identifier = int(values[0])
            if identifier in nodes:
                raise ContractError("Duplicate taxonomy node")
            nodes[identifier] = int(values[1]), values[2]
        for values in dump_rows(directory / "names.dmp", 4):
            if values[3] == "scientific name":
                identifier = int(values[0])
                if identifier in names and names[identifier] != values[1]:
                    raise ContractError("Conflicting scientific names")
                names[identifier] = values[1]
        for values in dump_rows(directory / "merged.dmp", 2):
            identifier = int(values[0])
            if identifier in merged:
                raise ContractError("Duplicate merged taxonomy ID")
            merged[identifier] = int(values[1])
        for values in dump_rows(directory / "delnodes.dmp", 1):
            deleted.add(int(values[0]))
        result = cls(nodes, names, merged, deleted)
        if not nodes or set(nodes) & (set(merged) | deleted) or set(merged) & deleted:
            raise ContractError("Taxonomy node, merged and deleted sets are inconsistent")
        checked = set()
        for identifier in nodes:
            chain, current = set(), identifier
            while current not in checked:
                if current not in nodes or current in chain:
                    raise ContractError("Taxonomy lineage has a cycle or missing parent")
                chain.add(current)
                parent = nodes[current][0]
                if current == 1 and parent == 1:
                    break
                current = parent
            checked.update(chain)
        for identifier in merged:
            if result.resolve(identifier)[1] in {"merge_cycle", "unknown_id"}:
                raise ContractError("Taxonomy merge chain has a cycle or invalid target")
        return result

    def resolve(self, identifier: int) -> tuple[int | None, str]:
        current, visited = identifier, set()
        while current in self.merged:
            if current in visited:
                return None, "merge_cycle"
            visited.add(current)
            current = self.merged[current]
        if current == 0:
            return None, "missing"
        if current in self.deleted:
            return None, "deleted"
        if current not in self.nodes:
            return None, "unknown_id"
        return current, "merged" if visited else "resolved"

    def lineage(self, identifier: int) -> list[int]:
        lineage, seen = [], set()
        while identifier not in seen:
            if identifier not in self.nodes:
                raise ContractError("Lineage contains an unknown node")
            lineage.append(identifier)
            seen.add(identifier)
            parent = self.nodes[identifier][0]
            if identifier == 1 and parent == 1:
                return lineage
            identifier = parent
        raise ContractError("Lineage cycle")

    def context(self, identifier: int) -> str:
        lineage = set(self.lineage(identifier))
        for ancestor, label in ((32630, "synthetic"), (81077, "vector"), (10239, "viral"), (131567, "cellular"),
                                (2, "cellular"), (2157, "cellular"), (2759, "cellular")):
            if ancestor in lineage:
                return label
        return "unknown"

    def assignment(self, raw: str, rank: str) -> dict:
        identifiers = parse_taxids(raw)
        resolutions = [dict(raw=identifier, resolved=self.resolve(identifier)[0], status=self.resolve(identifier)[1]) for identifier in identifiers]
        active = {row["resolved"] for row in resolutions if row["resolved"] is not None}
        contexts = {self.context(identifier) for identifier in active}
        result = dict(resolved_taxids=";".join(map(str, sorted(active))), taxon_id="", taxon_name="", status="missing",
                      context=next(iter(contexts)) if len(contexts) == 1 else "conflicting" if contexts else "unknown",
                      resolution_evidence=json.dumps(resolutions, sort_keys=True), all_resolved=bool(identifiers) and all(row["resolved"] is not None for row in resolutions))
        if not identifiers:
            return result
        if not result["all_resolved"]:
            result["status"] = "multiple_ambiguous" if len(identifiers) > 1 else resolutions[0]["status"]
            return result
        assignments = []
        for identifier in sorted(active):
            ranked = [ancestor for ancestor in self.lineage(identifier) if self.nodes[ancestor][1] == rank]
            assignments.append(ranked[0] if ranked and ranked[0] in self.names else None)
        if None in assignments:
            result["status"] = "rank_unresolved"
        elif len(set(assignments)) != 1:
            result["status"] = "multiple_ambiguous"
        else:
            assigned = assignments[0]
            result.update(taxon_id=str(assigned), taxon_name=self.names[assigned],
                          status="multiple_consensus" if len(identifiers) > 1 else "merged" if resolutions[0]["status"] == "merged" else "resolved")
        return result


def reconcile_taxids(taxonomy: Taxonomy, blast_taxids: str, record_taxids: str, rank: str) -> dict:
    blast, record = taxonomy.assignment(blast_taxids, rank), taxonomy.assignment(record_taxids, rank)
    blast_present, record_present = blast["status"] != "missing", record["status"] != "missing"
    if blast_present and record_present:
        agreed = blast["all_resolved"] and record["all_resolved"] and blast["resolved_taxids"] == record["resolved_taxids"]
        result = dict(record, source_status="agree" if agreed else "conflicting")
        if not agreed:
            result.update(status="conflicting", taxon_id="", taxon_name="", context="conflicting")
    elif record_present:
        result = dict(record, source_status="record_only")
    else:
        result = dict(blast, source_status="blast_only" if blast_present else "missing")
    result["resolution_evidence"] = json.dumps(dict(blast=blast, record=record), sort_keys=True)
    result.pop("all_resolved")
    return result


def pin_taxonomy(source: Path, destination: Path, source_url: str, retrieved_utc: str, release_date: str = "") -> dict:
    if urlsplit(source_url).scheme != "https":
        raise ContractError("Record an HTTPS taxonomy source URL")
    retrieved = datetime.fromisoformat(retrieved_utc.replace("Z", "+00:00"))
    if retrieved.tzinfo is None:
        raise ContractError("Taxonomy retrieval time must include its timezone")
    if release_date:
        date.fromisoformat(release_date)
    taxonomy = Taxonomy.from_directory(source)
    hashes = {name: digest(source / name) for name in TAXDUMP_FILES}
    destination.mkdir(parents=True, exist_ok=False)
    for name in TAXDUMP_FILES:
        shutil.copyfile(source / name, destination / name)
        if digest(destination / name) != hashes[name] or digest(source / name) != hashes[name]:
            raise ContractError("Taxonomy files changed during archival")
    report = dict(contract="pinned_ncbi_taxdump_v1", source_url=source_url, retrieved_utc=retrieved.isoformat(),
                  release_date=release_date or None, release_date_basis="provided_not_inferred" if release_date else "not_supplied",
                  archived_utc=datetime.now(timezone.utc).isoformat(), files=hashes, nodes=len(taxonomy.nodes),
                  identity="Exact archived file hashes; source and retrieval date are recorded declarations, not inferred from file mtime")
    save_json(destination / "taxonomy.json", report)
    return report


def load_taxonomy(directory: Path) -> tuple[Taxonomy, dict]:
    report = read_json(directory / "taxonomy.json")
    if report.get("contract") != "pinned_ncbi_taxdump_v1" or set(report["files"]) != set(TAXDUMP_FILES):
        raise ContractError("Incomplete pinned taxonomy contract")
    for name, checksum in report["files"].items():
        if digest(directory / name) != checksum:
            raise ContractError("Pinned taxonomy file changed: " + name)
    return Taxonomy.from_directory(directory), report


def taxonomy_rows(events: list[dict], taxonomy: Taxonomy, ranks: tuple[str, ...], snapshot_sha256: str) -> list[dict]:
    by_accession = defaultdict(lambda: {"blast": set(), "record": set()})
    for row in events:
        if row["accession_version"].startswith("unresolved:"):
            continue
        by_accession[row["accession_version"]]["blast"].add(row["blast_taxids"])
        by_accession[row["accession_version"]]["record"].add(row["record_taxids"])
    rows = []
    for accession, sources in sorted(by_accession.items()):
        values = {source: ";".join(map(str, sorted({identifier for raw in entries for identifier in parse_taxids(raw) if identifier != 0})))
                  for source, entries in sources.items()}
        for rank in ranks:
            conflict = False
            for entries in sources.values():
                identities = set()
                for raw in entries:
                    if parse_taxids(raw) not in {(), (0,)}:
                        identity = [taxonomy.resolve(identifier) for identifier in parse_taxids(raw)]
                        identities.add(tuple(sorted({resolved if resolved is not None else status for resolved, status in identity}, key=str)))
                conflict |= len(identities) > 1
            reconciled = reconcile_taxids(taxonomy, values["blast"], values["record"], rank)
            if conflict:
                reconciled.update(status="conflicting", source_status="within_source_conflict", taxon_id="", taxon_name="", context="conflicting")
            reconciled["resolution_evidence"] = json.dumps(dict(raw_sources={key: sorted(value) for key, value in sources.items()},
                                                                 reconciliation=json.loads(reconciled["resolution_evidence"])), sort_keys=True)
            rows.append(dict(accession_version=accession, rank=rank, blast_taxids=values["blast"], record_taxids=values["record"],
                             snapshot_sha256=snapshot_sha256, **reconciled))
    return rows


def taxonomy_composition(real_events: list[dict], assignments: list[dict], groups: list[str], policy: dict, ranks: tuple[str, ...], decisions: dict | None = None) -> list[dict]:
    lookup = {(row["accession_version"], row["rank"]): row for row in assignments}
    decisions = event_eligibility(real_events, policy["record_populations"]) if decisions is None else decisions
    populations = defaultdict(set)
    for row in real_events:
        for population in policy["record_populations"]:
            if decisions[row["accession_version"], population][0] != "included":
                continue
            for endpoint in ("strict", "tolerant"):
                if row[endpoint] == "yes":
                    populations[row["group"], endpoint, population].add(row["accession_version"])
                    populations["Global union", endpoint, population].add(row["accession_version"])
    rows = []
    for group in groups + ["Global union"]:
        for endpoint in ("strict", "tolerant"):
            for population in policy["record_populations"]:
                accessions = populations[group, endpoint, population]
                for rank in ranks:
                    counts = Counter({("unresolved:" + state, ""): 0 for state in UNRESOLVED_STATES})
                    for accession in accessions:
                        assignment = lookup.get((accession, rank))
                        if assignment is None:
                            key = "unresolved:missing", ""
                        elif assignment["taxon_id"]:
                            key = assignment["taxon_name"], assignment["taxon_id"]
                        else:
                            context = assignment["context"]
                            label = context + ":" if context in {"viral", "synthetic", "vector"} else "unresolved:"
                            key = label + assignment["status"], ""
                        counts[key] += 1
                    rows.extend(dict(group=group, endpoint=endpoint, population=population, rank=rank, category=category, taxon_id=taxon,
                                     accessions=count, denominator=len(accessions), percentage=str(100 * count / len(accessions)) if accessions else "not_estimated")
                                for (category, taxon), count in sorted(counts.items()))
    return rows


def annotate_taxonomy(analysis: Path, snapshot: Path, destination: Path, ranks: tuple[str, ...] = ("domain", "superkingdom", "genus")) -> dict:
    if not ranks or len(set(ranks)) != len(ranks) or any(not rank.strip() for rank in ranks):
        raise ContractError("Specify unique nonempty taxonomy ranks")
    source = read_json(analysis / "analysis.json")
    if source.get("result_contract") not in {"reference_screen_v1", "reference_screen_v2"}:
        raise ContractError("Taxonomy aggregation requires the new reference-screen result contract")
    for name, checksum in source["files"].items():
        if Path(name).name != name or digest(analysis / name) != checksum:
            raise ContractError("Source analysis artifact changed")
    if source.get("source_search") and (Path(source["source_search"]) / "quarantine.json").exists():
        raise ContractError("Source search is quarantined")
    taxonomy, evidence = load_taxonomy(snapshot)
    snapshot_hash = digest(snapshot / "taxonomy.json")
    events = read_events(analysis, source)
    assignments = taxonomy_rows(events, taxonomy, ranks, snapshot_hash)
    contexts = {row["accession_version"]: row["context"] for row in assignments}
    for row in events:
        context = contexts.get(row["accession_version"], "unknown")
        if context in {"synthetic", "vector"} and row["origin"] == "unknown":
            row["origin"] = context
        elif context in {"cellular", "viral", "conflicting"}:
            row["context"] = context if row["context"] in {"unknown", context} else "conflicting"
    logical = read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)
    statuses = read_table(analysis / "search_status.tsv", STATUS_SCHEMA)
    policy = load_policy(analysis / "analysis_policy.json")
    if policy != source["analysis_policy"]:
        raise ContractError("Analysis policy differs from its report")
    memberships = logical_memberships(logical, source["draws_per_parent"])
    real_events = project_real_events(events, memberships)
    composition = taxonomy_composition(real_events, assignments, sorted({row["group"] for row in logical if row["role"] == "real"}), policy, ranks,
                                       event_eligibility(events, policy["record_populations"]))
    provenance = dict(source, taxonomy_snapshot=dict(evidence, directory=str(snapshot.resolve()), manifest_sha256=snapshot_hash), taxonomy_ranks=list(ranks),
                      taxonomy_status="reconciled_with_explicit_unresolved_states", source_analysis_sha256=digest(analysis / "analysis.json"),
                      taxonomy_source_policy="Use sole available source; require active taxid-set agreement when both present; no title fallback")
    return write_completed_analysis(events, logical, statuses, policy, destination, provenance,
                                    {"accession_taxonomy.tsv": (assignments, TAXONOMY_SCHEMA), "taxonomic_composition.tsv": (composition, COMPOSITION_SCHEMA)})