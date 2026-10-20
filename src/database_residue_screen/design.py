import hashlib
from collections import Counter, defaultdict
from pathlib import Path

from .contracts import ContractError, read_json, read_table, schema
from .sequences import accession_version, canonical_id, read_fasta, reverse_complement


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def member(root: Path, name: str) -> Path:
    path = (root / name).resolve()
    if not name or Path(name).is_absolute() or not path.is_relative_to(root.resolve()):
        raise ContractError(f"Input must remain inside the study directory: {name}")
    if not path.is_file():
        raise ContractError(f"Missing input file: {name}")
    return path


def check_design(study_path: Path) -> dict:
    root = study_path.resolve().parent
    issues = []
    hashes = {study_path.name: digest(study_path)}

    def add(code, message, item="", severity="blocker"):
        issues.append(dict(code=code, severity=severity, item=item, message=message))

    def evidence(name, expected_hash, item):
        try:
            path = member(root, name)
            actual = digest(path)
            hashes[name] = actual
            if not expected_hash or actual != expected_hash:
                add("EVIDENCE_HASH", "Evidence checksum is missing or does not match", item)
        except (ContractError, OSError) as error:
            add("EVIDENCE_MISSING", str(error), item)

    try:
        study = read_json(study_path, schema("study"))
    except (ValueError, OSError) as error:
        return dict(schema_version="1.0.0", study_id="unreadable", scope="unknown", ready_for_freeze=False,
                    issues=[dict(code="STUDY_CONTRACT", severity="blocker", item=study_path.name, message=str(error))],
                    input_sha256=hashes, inventory=[], component_coverage=[], counts={})
    reference_screen = study.get("profile", "authenticated_residues") == "reference_screen"
    if study["decision_status"] != "approved":
        add("DESIGN_APPROVAL", "Scientific choices remain draft; review and approve them before freezing")
    if study["approval_file"]:
        try:
            path = member(root, study["approval_file"])
            if not path.read_text(encoding="utf-8").strip():
                raise ContractError("Empty approval record")
            hashes[study["approval_file"]] = digest(path)
        except (ContractError, OSError) as error:
            add("DESIGN_APPROVAL_EVIDENCE", str(error))
    else:
        add("DESIGN_APPROVAL_EVIDENCE", "A local decision/approval record is required")

    inputs = {}
    for name, relative in study["inputs"].items():
        try:
            path = member(root, relative)
            hashes[relative] = digest(path)
            if name in {"queries", "validation_subjects"}:
                inputs[name] = read_fasta(path, allow_empty=name == "validation_subjects", allow_ambiguity=name == "validation_subjects")
            elif name == "database":
                inputs[name] = read_json(path, schema("database"))
            else:
                inputs[name] = read_table(path, schema("queries" if name == "query_metadata" else name))
        except (ValueError, OSError) as error:
            add("INPUT_CONTRACT", str(error), name)

    search = study["search"]
    if search["real_target_limit"] != search["control_target_limit"]:
        add("UNEQUAL_RETRIEVAL", "Real and control target limits must be identical")
    if search["retention_ladder"] != sorted(search["retention_ladder"]) or search["real_target_limit"] not in search["retention_ladder"]:
        add("RETENTION_DESIGN", "An increasing retention ladder must include the starting target limit")
    if search["failure_policy"] != "abort_run_resume_retries_failed_batches":
        add("FAILURE_POLICY_UNIMPLEMENTED", "Runner aborts and retries on resume; selective matched-parent failure blocking is not implemented", severity="warning")
    controls = study["controls"]
    if not controls["algorithm_version"]:
        add("NULL_VERSION", "Pin the verified shuffle implementation version before freezing")
    if controls["draw_ladder"] != sorted(controls["draw_ladder"]):
        add("NULL_PRECISION_DESIGN", "Draw-count ladder must be increasing")
    if not study["metadata"]["taxonomy_snapshot"]:
        add("TAXONOMY_DESIGN", "Specify the dated taxonomy snapshot before interpreting taxonomic outputs")
    companion = study["companion"]
    if companion["reuse"] == "undecided":
        add("COMPANION_REUSE", "Decide shared-input reuse; matching build dates do not establish identity")
    elif companion["reuse"] == "shared":
        try:
            path = member(root, companion["reconciliation_file"])
            hashes[companion["reconciliation_file"]] = digest(path)
        except (ContractError, OSError) as error:
            add("COMPANION_RECONCILIATION", str(error))
    if companion["citation_status"] == "pending":
        add("COMPANION_CITATION", "Formal citation remains a publication requirement", severity="warning")

    database = inputs.get("database")
    if database:
        if database["identity_status"] == "unverified":
            add("DATABASE_IDENTITY", "Recorded build metadata is not verified database-file identity")
        if study["scope"] == "research" and database["identity_status"] == "synthetic":
            add("SYNTHETIC_DATABASE", "A synthetic example cannot establish a research database identity")
        if database["sequence_count"] <= 0 or database["base_count"] <= 0:
            add("DATABASE_COUNTS", "Database counts must be provided")
        evidence(database["evidence_file"], database["evidence_sha256"], "database")
        if database["kind"] == "blast" and database["evidence_file"]:
            try:
                identity = read_json(member(root, database["evidence_file"]))
                if identity.get("metadata_consistent") is not True:
                    add("DATABASE_METADATA_CONFLICT", "Database identity evidence has conflicting or missing metadata")
                if identity.get("all_database_files_hashed") is not True:
                    add("DATABASE_HASH_SCOPE", "Complete database content hashes are not recorded")
                observed = identity.get("runtime_metadata", {})
                if observed.get("sequence_count") != database["sequence_count"] or observed.get("base_count") != database["base_count"]:
                    add("DATABASE_EVIDENCE_COUNTS", "Declared counts disagree with captured runtime metadata")
            except (ValueError, OSError) as error:
                add("DATABASE_EVIDENCE_CONTRACT", str(error))

    sources = {}
    for source in inputs.get("sources", []):
        if source["source_id"] in sources:
            add("DUPLICATE_SOURCE", "Duplicate source ID", source["source_id"])
        sources[source["source_id"]] = source

    sequences = inputs.get("queries", {})
    metadata = inputs.get("query_metadata", [])
    query_map = {}
    for query in metadata:
        if query["query_id"] in query_map:
            add("DUPLICATE_QUERY_METADATA", "Duplicate physical query metadata", query["query_id"])
        query_map[query["query_id"]] = query
    if set(query_map) != set(sequences):
        add("QUERY_MEMBERSHIP", f"FASTA-only IDs: {len(set(sequences) - set(query_map))}; metadata-only IDs: {len(set(query_map) - set(sequences))}")
    inventory = []
    primary = []
    for query_id, query in query_map.items():
        if query_id not in sequences:
            continue
        record = dict(query, canonical_id=canonical_id(sequences[query_id]), length=len(sequences[query_id]))
        inventory.append(record)
        source = sources.get(query["source_id"])
        selected = query["role"] != "reference" if reference_screen else query["role"] == "primary"
        if source is None:
            add("QUERY_SOURCE_MISSING", "No source row for query", query_id, "blocker" if selected else "warning")
        if not selected:
            continue
        primary.append(record)
        if study["scope"] == "research" and (query["provenance_status"] == "synthetic" or query["architecture_status"] == "synthetic"):
            add("SYNTHETIC_PRIMARY", "Synthetic-example labels cannot authenticate a research query", query_id)
        if source:
            if source["verification"] == "unverified" or (study["scope"] == "research" and source["verification"] != "verified"):
                add("SOURCE_UNVERIFIED", "Reference basis is declared, not independently authenticated" if reference_screen else
                    "Primary query source is not verified for this scope", query_id, "warning" if reference_screen else "blocker")
            if not source["reference"] or not source["version"]:
                add("SOURCE_VERSION", "Available source/version detail must be reported without inventing proprietary information",
                    query_id, "warning" if reference_screen else "blocker")
            evidence(source["evidence_file"], source["evidence_sha256"], source["source_id"])
    if not primary:
        add("NO_PRIMARY_QUERIES", "No original reference query population has been declared" if reference_screen else
            "No authenticated primary query population has been declared")
    for group in study["analysis"]["primary_groups"]:
        if not any(query["group"] == group for query in primary):
            add("PRIMARY_GROUP_EMPTY", "No primary queries for declared group", group)
    undeclared = {query["group"] for query in primary} - set(study["analysis"]["primary_groups"])
    if undeclared:
        add("UNDECLARED_PRIMARY_GROUP", ", ".join(sorted(undeclared)))
    families = defaultdict(set)
    for query in primary:
        families[query["canonical_id"]].add(query["family_id"])
    for identifier, values in families.items():
        if len(values) > 1:
            add("CANONICAL_FAMILY_CONFLICT", "Identical canonical sequences have multiple family labels", identifier,
                "warning" if reference_screen else "blocker")

    components = defaultdict(list)
    seen_components = set()
    for component in ([] if reference_screen else inputs.get("components", [])):
        key = tuple(component.values())
        if key in seen_components:
            add("DUPLICATE_COMPONENT", "Duplicate component mapping", component["composite_id"])
        seen_components.add(key)
        composite = component["composite_id"]
        counterpart = component["component_id"]
        if composite not in sequences or counterpart not in sequences:
            add("COMPONENT_MISSING", "Composite or standalone counterpart is absent from FASTA", composite)
            continue
        if query_map.get(composite, {}).get("construct_class") != "composite":
            add("COMPONENT_PARENT_CLASS", "Mapped parent must be declared composite", composite)
        if component["end"] < component["start"] or component["end"] > len(sequences[composite]):
            add("COMPONENT_BOUNDS", "Invalid one-based inclusive component coordinates", composite)
            continue
        expected = sequences[counterpart]
        if component["strand"] == "-":
            expected = reverse_complement(expected)
        if sequences[composite][component["start"] - 1:component["end"]] != expected:
            add("COMPONENT_SEQUENCE", "Component sequence does not match its declared footprint", composite)
        if query_map.get(counterpart, {}).get("construct_class") != component["component_kind"]:
            add("COMPONENT_KIND", "Counterpart class disagrees with component kind", counterpart)
        components[composite].append(component)
    coverage = []
    for query in ([] if reference_screen else inventory):
        if query["construct_class"] != "composite":
            continue
        counts = Counter(position for component in components[query["query_id"]]
                         for position in range(component["start"], component["end"] + 1))
        kinds = {component["component_kind"] for component in components[query["query_id"]]}
        covered = len(counts)
        complete = covered == query["length"] and all(count == 1 for count in counts.values()) and {"barcode", "adapter"} <= kinds
        coverage.append(dict(query_id=query["query_id"], query_length=query["length"], covered_bases=covered,
                             component_count=len(components[query["query_id"]]), structurally_complete=complete,
                             physical_architecture_verified=query["architecture_status"] == "verified"))
        if not complete:
            add("COMPONENT_COVERAGE", "Composite must have a nonoverlapping full barcode/adapter decomposition with standalone counterparts",
                query["query_id"], "blocker" if query["role"] == "primary" else "warning")
        if query["role"] == "primary":
            for component in components[query["query_id"]]:
                counterpart = query_map[component["component_id"]]
                if counterpart["provenance_status"] == "unverified":
                    add("COMPONENT_PROVENANCE", "Primary component counterpart remains unverified", component["component_id"])

    subjects = inputs.get("subjects", [])
    subject_map = {}
    for subject in subjects:
        try:
            if accession_version(subject["raw_id"]) != subject["accession_version"]:
                add("SUBJECT_ID_MISMATCH", "Raw identifier does not map to declared accession.version", subject["raw_id"])
        except ContractError as error:
            add("SUBJECT_ID_AMBIGUOUS", str(error), subject["raw_id"])
        identifier = subject["accession_version"]
        if identifier in subject_map and subject_map[identifier] != subject:
            add("SUBJECT_METADATA_CONFLICT", "Conflicting metadata for one accession.version", identifier)
        subject_map[identifier] = subject

    cases = inputs.get("validation_cases", [])
    truth_sequences = inputs.get("validation_subjects", {})
    if not cases and not reference_screen:
        add("VALIDATION_DESIGN_EMPTY", "A defined independent evaluation truth set is required")
    case_ids = set()
    split_families = defaultdict(set)
    split_sources = defaultdict(set)
    split_sequences = defaultdict(set)
    split_subjects = defaultdict(set)
    for case in cases:
        if case["case_id"] in case_ids:
            add("DUPLICATE_TRUTH_CASE", "Duplicate validation case ID", case["case_id"])
        case_ids.add(case["case_id"])
        if case["query_id"] not in sequences or case["subject_id"] not in truth_sequences:
            add("TRUTH_INPUT_MISSING", "Truth query or subject is absent from the specified FASTA", case["case_id"])
            continue
        split_families[case["split"]].add(case["family_id"])
        split_sources[case["split"]].add(case["subject_source"])
        split_sequences[case["split"]].add(canonical_id(sequences[case["query_id"]]))
        split_subjects[case["split"]].add(canonical_id(truth_sequences[case["subject_id"]]))
        if case["family_id"] != query_map.get(case["query_id"], {}).get("family_id"):
            add("TRUTH_FAMILY_MISMATCH", "Truth family does not match query metadata", case["case_id"])
        query_sequence = sequences[case["query_id"]]
        exact = query_sequence in truth_sequences[case["subject_id"]] or reverse_complement(query_sequence) in truth_sequences[case["subject_id"]]
        if (case["truth"] == "strict_positive" and not exact) or (case["truth"] != "strict_positive" and exact):
            add("EXACT_TRUTH_CONTRADICTION", "Declared sequence truth contradicts independent exact matching", case["case_id"])
    evaluation = [case for case in cases if case["split"] == "evaluation"]
    if not reference_screen and not {"strict_positive", "negative"} <= {case["truth"] for case in evaluation}:
        add("EVALUATION_CLASSES", "Evaluation must contain strict positives and explicit negatives")
    missing_classes = set(study["validation"]["required_error_classes"]) - {case["error_class"] for case in evaluation}
    if missing_classes and not reference_screen:
        add("VALIDATION_STRATA", "Missing evaluation error classes: " + ", ".join(sorted(missing_classes)))
    if study["validation"]["tuning"] and (split_families["development"] & split_families["evaluation"] or split_sources["development"] & split_sources["evaluation"]):
        add("VALIDATION_LEAKAGE", "Development and evaluation share query families or subject sources")
    if study["validation"]["tuning"]:
        if not split_families["development"]:
            add("DEVELOPMENT_SPLIT_MISSING", "Tuning requires a separate development split")
        if split_sequences["development"] & split_sequences["evaluation"] or split_subjects["development"] & split_subjects["evaluation"]:
            add("VALIDATION_SEQUENCE_LEAKAGE", "Development and evaluation contain the same canonical query or subject sequence")

    counts = dict(physical_queries=len(sequences), metadata_rows=len(metadata), primary_queries=len(primary),
                  global_canonical_sequences=len({query["canonical_id"] for query in inventory}),
                  within_condition_canonical_entries=len({(query["group"], query["condition"], query["canonical_id"]) for query in inventory}),
                  composites=len(coverage), complete_component_decompositions=sum(row["structurally_complete"] for row in coverage),
                  validation_cases=len(cases), subject_metadata_records=len(subjects))
    return dict(schema_version="1.0.0", study_id=study["study_id"], scope=study["scope"],
                profile=study.get("profile", "authenticated_residues"),
                interpretation="Reference-sequence matches, not biological residue authentication" if reference_screen else "Authenticated-residue study design",
                ready_for_freeze=not any(issue["severity"] == "blocker" for issue in issues),
                scope_limit="G0 design completeness only; not search performance, completed searches or validated biological findings",
                input_sha256=hashes, issues=issues, counts=counts, inventory=inventory, component_coverage=coverage)