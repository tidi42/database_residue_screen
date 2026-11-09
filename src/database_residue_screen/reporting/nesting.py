"""Nesting report: per-parent strict events versus distinct loci (descriptive diagnostic).

Nested queries are co-positive (the 24-nt ONT Native barcode lies inside the K40 and K46 contexts; Illumina bare
indexes inside the embedded adapters), so one physical site can yield several strict events. The event unit stays one
event per parent, accession, coordinates and strand (density and positional tables). This report counts, per record
population and group, how many of those events fall on distinct loci (overlapping strict intervals on one accession
merged, any parent or strand), with exactly the event gate of the analysis (strict, strict-testable length, eligible
record). Per-parent counts must equal the analysis density tables. Nothing is filtered or re-estimated.
"""
from collections import Counter, defaultdict
from pathlib import Path

from ..contracts import ContractError
from ..policy import applicability
from .common import ratio, read_tsv, write_json, write_tsv


def loci(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(item) for item in merged]


def write_nesting(destination: Path, *, analysis: Path, policy: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    included = {(row["sequence_id"], row["accession_version"], row["population"]) for row in read_tsv(analysis / "eligibility.tsv")
                if row["state"] == "included"}
    events = [row for row in read_tsv(analysis / "real_events.tsv") if row["strict"] == "yes"]
    group_rows, combination_rows, summary = [], [], {}
    for population in policy["record_populations"]:
        selected = [row for row in events if (row["query_id"], row["accession_version"], population) in included
                    and applicability(int(row["qlen"]), "strict", policy) == "operational_retained_search"]
        per_parent, intervals, conditions = defaultdict(set), defaultdict(list), defaultdict(set)
        for row in selected:
            start, end = sorted((int(row["sstart"]), int(row["send"])))
            key = row["group"], row["accession_version"]
            per_parent[key].add((row["parent_id"], start, end, row["canonical_strand"]))
            intervals[key].append((start, end))
            conditions[key].add((start, end, row["condition"]))
        density = {(row["group"], row["accession_version"]): int(row["strict_orientation_deduplicated_events"])
                   for row in read_tsv(analysis / f"{population}_density.tsv")}
        if density != {key: len(value) for key, value in per_parent.items()}:
            raise ContractError(f"Per-parent events differ from {population}_density.tsv; the event gate is not reproduced")
        by_group, combinations = defaultdict(Counter), defaultdict(Counter)
        for key, here in per_parent.items():
            merged = loci(intervals[key])
            counts = by_group[key[0]]
            counts["accessions"] += 1
            counts["per_parent_events"] += len(here)
            counts["distinct_loci"] += len(merged)
            counts["single_event_accessions"] += len(here) == 1
            counts["single_locus_accessions"] += len(merged) == 1
            for start, end in merged:
                inside = {condition for a, b, condition in conditions[key] if a <= end and start <= b}
                counts["loci_with_several_events"] += len({event for event in here if event[1] <= end and start <= event[2]}) > 1
                combinations[key[0]][";".join(sorted({part for condition in inside for part in condition.split(";")}))] += 1
        for group in sorted(by_group):
            counts = by_group[group]
            group_rows.append(dict(population=population, group=group, accessions=counts["accessions"],
                                   per_parent_events=counts["per_parent_events"], distinct_loci=counts["distinct_loci"],
                                   loci_with_several_events=counts["loci_with_several_events"],
                                   single_event_accessions=counts["single_event_accessions"],
                                   single_event_share=ratio(counts["single_event_accessions"], counts["accessions"]),
                                   single_locus_accessions=counts["single_locus_accessions"],
                                   single_locus_share=ratio(counts["single_locus_accessions"], counts["accessions"])))
            for combination, number in sorted(combinations[group].items(), key=lambda item: (-item[1], item[0])):
                combination_rows.append(dict(population=population, group=group, locus_conditions=combination, loci=number))
        summary[population] = {row["group"]: {key: row[key] for key in ("accessions", "per_parent_events", "distinct_loci", "single_event_share",
                                                                        "single_locus_share")}
                               for row in group_rows if row["population"] == population}
    write_tsv(destination / "events_and_loci.tsv", group_rows,
              ("population", "group", "accessions", "per_parent_events", "distinct_loci", "loci_with_several_events", "single_event_accessions",
               "single_event_share", "single_locus_accessions", "single_locus_share"))
    write_tsv(destination / "locus_condition_sets.tsv", combination_rows, ("population", "group", "locus_conditions", "loci"))
    write_json(destination / "summary.json", summary)
    return summary
