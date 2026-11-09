"""Positional report: terminal-event robustness and record type x taxonomy of the positional events.

The positional event population is the framework's: strict, a-priori strict-testable events of every panel condition in
eligible linear-genomic records. The sensitivity table is recomputed with the framework function (one event per
accession and per BioProject, 25-nt diagnostic window, INSDC targeted exclusion) from the archived ESummary and INSDC
snapshots and must equal the positional stage table when that is given. Events, accessions and projects are counted
separately; HSPs are not independent genomes.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..analysis import ANNOTATED_EVENT_SCHEMA
from ..contracts import ContractError, read_table
from ..eligibility import event_eligibility
from ..insdc import read_annotation
from ..policy import applicability
from ..positional import SENSITIVITY_SCHEMA, positional_sensitivity, primary_positional_events, projects_from_metadata
from .common import read_tsv, write_json, write_tsv


def write_positions(destination: Path, *, analysis: Path, metadata: Path, insdc: Path | None = None,
                    positional_table: Path | None = None) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    projects = projects_from_metadata(metadata)
    targeted = read_annotation(insdc) if insdc is not None else None
    rows = positional_sensitivity(analysis, projects, targeted)
    if positional_table is not None and read_table(positional_table, SENSITIVITY_SCHEMA) != rows:
        raise ContractError("The positional stage table differs from the recomputed positional sensitivity")
    write_tsv(destination / "positional_sensitivity.tsv", rows, list(SENSITIVITY_SCHEMA["properties"]))
    events, policy = primary_positional_events(analysis)
    all_events = read_table(analysis / "real_events.tsv", ANNOTATED_EVENT_SCHEMA)
    decisions = event_eligibility(all_events, ["linear_genomic"])
    contexts = defaultdict(set)
    if (analysis / "accession_taxonomy.tsv").is_file():
        for row in read_tsv(analysis / "accession_taxonomy.tsv"):
            contexts[row["accession_version"]].add(row["context"])
    strict = [row for row in all_events if row["strict"] == "yes"
              and applicability(row["qlen"], "strict", policy) == "operational_retained_search"]
    table, accessions = Counter(), defaultdict(set)
    for row in strict:
        key = (row["group"], row["condition"], row["record_type"], row["topology"],
               ";".join(sorted(contexts.get(row["accession_version"], {"not_annotated"}))), decisions[row["accession_version"], "linear_genomic"][1])
        table[key] += 1
        accessions[key].add(row["accession_version"])
    write_tsv(destination / "record_taxonomy.tsv",
              [dict(group=key[0], condition=key[1], record_type=key[2], topology=key[3], taxonomy_context=key[4], linear_genomic_reason=key[5],
                    events=count, accessions=len(accessions[key]),
                    bioprojects=len({projects.get(accession) or "record:" + accession for accession in accessions[key]}))
               for key, count in sorted(table.items())],
              ("group", "condition", "record_type", "topology", "taxonomy_context", "linear_genomic_reason", "events", "accessions", "bioprojects"))
    position = policy["position"]
    summary = dict(positional_events=len(events), accessions=len({row["accession_version"] for row in events}),
                   bioprojects=len({projects.get(row["accession_version"]) or "record:" + row["accession_version"] for row in events}),
                   accessions_without_project=len({row["accession_version"] for row in events if not projects.get(row["accession_version"])}),
                   insdc_targeted_variant=targeted is not None, stage_table_checked=positional_table is not None,
                   window=position["primary_window"], minimum_length=position["primary_minimum_length"],
                   headline=[{key: row[key] for key in ("variant", "group", "events", "terminal_events", "expected_n2", "fold_n2")}
                             for row in rows if row["window"] == position["primary_window"]
                             and row["minimum_length"] == position["primary_minimum_length"]])
    write_json(destination / "summary.json", summary)
    return summary
