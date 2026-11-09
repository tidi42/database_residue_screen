"""Eligibility report: accession-level eligibility flow, taxonomy context and record-field population.

Scopes (real strict events of every panel condition):
  testable_strict  events whose query length is strict-testable a priori (the positional and Figure 4A population)
  ungated_strict   every strict event, no testability gate (inventory sensitivity)
Reasons follow the code order record_type -> origin -> targeted -> RNA -> circular; each accession is counted once at its
first applicable reason. Taxonomy context separates cellular, viral, vector, synthetic and other artificial records.
Targeted/engineered fields are reported as populated or not; nothing is relabelled.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..analysis import ANNOTATED_EVENT_SCHEMA
from ..contracts import read_table
from ..eligibility import event_eligibility
from ..policy import applicability
from .common import read_tsv, write_json, write_tsv

FIELDS = ("targeted", "origin", "record_type", "molecule", "topology", "organelle")


def scoped_events(events: list[dict], policy: dict) -> dict[str, list[dict]]:
    strict = [row for row in events if row["strict"] == "yes"]
    return {"testable_strict": [row for row in strict if applicability(row["qlen"], "strict", policy) == "operational_retained_search"],
            "ungated_strict": strict}


def write_eligibility(destination: Path, *, analysis: Path, policy: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    events = read_table(analysis / "real_events.tsv", ANNOTATED_EVENT_SCHEMA)
    decisions = event_eligibility(events, ["all_records", "linear_genomic"])
    contexts = defaultdict(set)
    taxonomy = analysis / "accession_taxonomy.tsv"
    if taxonomy.is_file():
        for row in read_tsv(taxonomy):
            contexts[row["accession_version"]].add(row["context"])
    flow, context_rows, field_rows, summary = [], [], [], {}
    for scope, selected in scoped_events(events, policy).items():
        accessions = sorted({row["accession_version"] for row in selected})
        resolved = [accession for accession in accessions if not accession.startswith("unresolved:")]
        reasons = Counter(decisions[accession, "linear_genomic"] for accession in accessions)
        flow.append(dict(scope=scope, state="all", reason="distinct_accessions", accessions=len(accessions)))
        flow.append(dict(scope=scope, state="all", reason="resolved_identifiers", accessions=len(resolved)))
        flow.extend(dict(scope=scope, state=state, reason=reason, accessions=count) for (state, reason), count in sorted(reasons.items()))
        counts = Counter(";".join(sorted(contexts.get(accession, {"not_annotated"}))) for accession in resolved)
        context_rows.extend(dict(scope=scope, context=context, accessions=count) for context, count in sorted(counts.items()))
        first = {}
        for row in selected:
            first.setdefault(row["accession_version"], row)
        for field in FIELDS:
            values = Counter(first[accession].get(field, "") for accession in accessions)
            field_rows.extend(dict(scope=scope, field=field, value=value, accessions=count) for value, count in sorted(values.items()))
        summary[scope] = dict(events=len(selected), distinct_accessions=len(accessions), resolved_identifiers=len(resolved),
                              **{f"{state}|{reason}": count for (state, reason), count in sorted(reasons.items())})
    write_tsv(destination / "eligibility_flow.tsv", flow, ("scope", "state", "reason", "accessions"))
    write_tsv(destination / "taxonomy_context.tsv", context_rows, ("scope", "context", "accessions"))
    write_tsv(destination / "field_population.tsv", field_rows, ("scope", "field", "value", "accessions"))
    per_group = []
    for scope, selected in scoped_events(events, policy).items():
        by_group = defaultdict(set)
        for row in selected:
            if not row["accession_version"].startswith("unresolved:"):
                by_group[row["group"]].add(row["accession_version"])
        union = set().union(*by_group.values()) if by_group else set()
        for group, accessions in sorted(by_group.items()):
            per_group.append(dict(scope=scope, group=group, accessions=len(accessions)))
        per_group.append(dict(scope=scope, group="Global union", accessions=len(union)))
    write_tsv(destination / "accessions_by_group.tsv", per_group, ("scope", "group", "accessions"))
    write_json(destination / "summary.json", summary)
    return summary
