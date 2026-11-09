import shutil
from collections import Counter, defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, read_table, schema, write_table
from .controls import DRAW_SCHEMA
from .design import digest
from .example import save_fasta, save_json
from .policy import load_policy, validate_policy
from .sequences import canonical, canonical_id, read_fasta


LOGICAL_SCHEMA = dict(type="object", additionalProperties=False, properties={
    **{name: {"type": "string"} for name in ("logical_id", "role", "parent_id", "sequence_id", "orientation",
                                             "group", "condition", "construct_class", "family_id",
                                             "provenance_status", "architecture_status")},
    "draw_index": {"type": "integer"}, "length": {"type": "integer"}, "seed": {"type": "integer"},
})


def length_stratum(length: int) -> str:
    for bound in (7, 12, 20, 24, 32, 64, 128):
        if length < bound:
            return f"below_{bound:03d}"
    return "at_least_128"


def select_sample(metadata: list[dict], queries: dict[str, str], count: int) -> set[str]:
    strata = defaultdict(set)
    for row in metadata:
        sequence = queries[row["query_id"]]
        complexity = "low_complexity" if max(Counter(sequence).values()) / len(sequence) >= 0.8 or len({sequence[index:index + 2] for index in range(len(sequence) - 1)}) <= 3 else "other"
        key = row["group"], row["construct_class"], length_stratum(len(sequence)), complexity
        strata[key].add(canonical_id(sequence))
    by_group = defaultdict(list)
    for (group, construct, length, complexity), values in sorted(strata.items()):
        by_group[group].append((construct, length, complexity, min(values)))
    selected = set()
    while len(selected) < count and any(by_group.values()):
        for group in sorted(by_group):
            if by_group[group] and len(selected) < count:
                selected.add(by_group[group].pop(0)[3])
    return selected


def prepare_panel(query_fasta: Path, query_metadata: Path, controls_directory: Path, destination: Path,
                  draws: int, sample_parents: int = 0, policy_file: Path | None = None, *, policy: dict | None = None) -> dict:
    """Bind every physical query of one panel, its matched control draws and the analysis policy (from `policy_file`,
    or the mapping `policy`, else the defaults)."""
    if draws < 1 or sample_parents < 0:
        raise ContractError("Draw count must be positive; sample parent count must be nonnegative")
    if policy is not None and policy_file is not None:
        raise ContractError("Give the analysis policy as a file or as a mapping, not both")
    policy = validate_policy(dict(policy)) if policy is not None else load_policy(policy_file)
    queries = read_fasta(query_fasta)
    metadata = read_table(query_metadata, schema("queries"))
    if len(metadata) != len(queries) or {row["query_id"] for row in metadata} != set(queries):
        raise ContractError("Physical query FASTA and metadata membership disagree")
    if not metadata:
        raise ContractError("Panel contains no queries")
    control_qc = read_json(controls_directory / "control_qc.json")
    if (control_qc["source_sha256"] != digest(query_fasta) or not control_qc["all_invariants_valid"]
            or draws > control_qc["draws_per_parent"]):
        raise ContractError("Control QC does not cover the requested real inputs and draw count")
    for name, checksum in control_qc["files"].items():
        if name not in {"controls.fasta", "draws.tsv", "parents.tsv"} or digest(controls_directory / name) != checksum:
            raise ContractError("Control panel changed since its invariant audit")
    parents = (select_sample(metadata, queries, sample_parents) if sample_parents else
               {canonical_id(queries[row["query_id"]]) for row in metadata})
    control_sequences = read_fasta(controls_directory / "controls.fasta")
    draw_rows = read_table(controls_directory / "draws.tsv", DRAW_SCHEMA)
    search_sequences, logical = {}, []
    for row in metadata:
        sequence = queries[row["query_id"]]
        parent = canonical_id(sequence)
        if parent not in parents:
            continue
        search_sequences[parent] = canonical(sequence)
        logical.append(dict(logical_id=row["query_id"], role="real", parent_id=parent, sequence_id=parent,
                            orientation="+" if sequence == canonical(sequence) else "-", draw_index=0, length=len(sequence), seed=0,
                            **{key: row[key] for key in ("group", "condition", "construct_class", "family_id", "provenance_status", "architecture_status")}))
    seen_draws, seeds = set(), []
    for row in draw_rows:
        if row["parent_id"] not in parents or row["draw_index"] > draws:
            continue
        key = row["parent_id"], row["draw_index"]
        if key in seen_draws or row["draw_index"] < 1:
            raise ContractError("Control draw identity is duplicated or invalid")
        seen_draws.add(key)
        seeds.append(row["seed"])
        sequence = control_sequences[row["draw_id"]]
        sequence_id = canonical_id(sequence)
        if sequence_id != row["sequence_id"]:
            raise ContractError("Control sequence disagrees with its draw manifest")
        search_sequences[sequence_id] = canonical(sequence)
        logical.append(dict(logical_id=row["draw_id"], role="control", parent_id=row["parent_id"], sequence_id=sequence_id,
                            orientation="+" if sequence == canonical(sequence) else "-", draw_index=row["draw_index"], length=len(sequence), seed=row["seed"],
                            group="", condition="", construct_class="", family_id="",
                            provenance_status="generated", architecture_status="not_applicable"))
    if seen_draws != {(parent, index) for parent in parents for index in range(1, draws + 1)}:
        raise ContractError("One or more canonical parents lack the requested control draws")
    from .controls import null_space_status
    # Flag only proven singleton spaces; retain their declared sampling weights.
    policy = dict(policy, degenerate_parent_ids=sorted(parent for parent in parents
                       if null_space_status(search_sequences[parent]) == "proven_singleton_null"))
    destination.mkdir(parents=True, exist_ok=False)
    save_fasta(destination / "search_queries.fasta", dict(sorted(search_sequences.items())))
    write_table(destination / "logical_queries.tsv", logical, LOGICAL_SCHEMA)
    save_json(destination / "analysis_policy.json", policy)
    shutil.copyfile(query_fasta, destination / "physical_queries.fasta")
    shutil.copyfile(query_metadata, destination / "physical_queries.tsv")
    report = dict(stage="search_panel", scientific_approval_claimed=False, source_physical_queries=len(metadata),
                  source_queries_sha256=digest(query_fasta), source_metadata_sha256=digest(query_metadata),
                  physical_real_queries=sum(row["role"] == "real" for row in logical), canonical_parents=len(parents),
                  requested_sample_parent_limit=sample_parents, sample_selection="Round-robin groups; smallest canonical hash per construct/length/composition stratum; no hit-based selection",
                  sample_composition_rule="Low complexity if dominant base >=80% or <=3 observed dinucleotides; stratification only, no exclusion",
                  control_draws=len(seen_draws), draws_per_parent=draws, logical_queries=len(logical),
                  repeated_control_seed_values=len(seeds) - len(set(seeds)), control_algorithm=control_qc["algorithm"],
                  distinct_search_sequences=len(search_sequences), sequence_reuse="Identical and reverse-complement sequences searched once; logical parent and draw weights retained",
                  analysis_policy=policy, matching="strict_full_query_ungapped_and_identity95_span90_ungapped",
                  inference="descriptive_fixed_panel_no_hypothesis_tests", biology="Declared reference-panel matches, not physical residue authentication",
                  position_policy="Primary linear-genomic positions require explicit metadata; unknown records remain separate",
                  controls_qc_sha256=digest(controls_directory / "control_qc.json"),
                  files={name: digest(destination / name) for name in ("search_queries.fasta", "logical_queries.tsv", "physical_queries.fasta", "physical_queries.tsv", "analysis_policy.json")})
    save_json(destination / "panel.json", report)
    return report