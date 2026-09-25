from datetime import datetime, timezone
from pathlib import Path

from .contracts import ContractError, read_table, schema
from .controls import generate_controls
from .design import digest
from .example import save_json
from .legacy import PANEL_FILES
from .panels import prepare_panel, select_original_queries
from .precision import write_pilot_spec
from .sequences import read_fasta


def prepare_reference_inputs(root: Path, imported_inputs: Path, destination: Path, draws: int = 100, seed: int = 42) -> dict:
    original, source_hashes = {}, {}
    for name in PANEL_FILES:
        path = root / "panels" / name
        source_hashes[str(path.resolve())] = digest(path)
        for identifier, sequence in read_fasta(path).items():
            if identifier in original and original[identifier] != sequence:
                raise ContractError("Original panel files disagree on query sequence: " + identifier)
            original[identifier] = sequence
    imported = read_fasta(imported_inputs / "queries.fasta")
    metadata = read_table(imported_inputs / "queries.tsv", schema("queries"))
    if len(metadata) != len(imported) or {row["query_id"] for row in metadata} != set(imported):
        raise ContractError("Imported query membership is inconsistent")
    selected = {row["query_id"]: imported[row["query_id"]] for row in metadata if row["role"] != "reference"}
    if selected != original:
        raise ContractError("Selected imported queries no longer match the original author panels; reconcile the inputs before generating controls")
    destination.mkdir(parents=True, exist_ok=False)
    originals = destination / "originals"
    selection = select_original_queries(imported_inputs / "queries.fasta", imported_inputs / "queries.tsv", originals)
    sources = read_table(imported_inputs / "sources.tsv", schema("sources"))
    save_json(destination / "reference_basis.json", dict(profile="reference_screen", source_panel_sha256=source_hashes,
                                                        imported_source_table=str((imported_inputs / "sources.tsv").resolve()),
                                                        imported_source_table_sha256=digest(imported_inputs / "sources.tsv"), sources=sources,
                                                        interpretation="Declared supplier references and author industry knowledge; no physical-junction authentication or new component comparison"))
    control_qc = generate_controls(originals / "queries.fasta", destination / "controls", draws, seed)
    full_panel = prepare_panel(originals / "queries.fasta", originals / "queries.tsv", destination / "controls", destination / "full-panel", draws)
    pilot = prepare_panel(originals / "queries.fasta", originals / "queries.tsv", destination / "controls", destination / "pilot-panel", draws, pilot_parents=16)
    write_pilot_spec(destination / "pilot-panel", destination / "pilot-spec")
    for name, checksum in source_hashes.items():
        if digest(Path(name)) != checksum:
            raise ContractError("Original panel changed during preparation")
    report = dict(stage="reference_inputs_prepared_not_searched", profile="reference_screen", final_freeze=False,
                  prepared_utc=datetime.now(timezone.utc).isoformat(), source_panel_sha256=source_hashes,
                  original_physical_queries=selection["original_physical_queries"], canonical_parents=selection["canonical_parents"],
                  excluded_reference_ids=selection["excluded_reference_ids"], total_control_draws=control_qc["total_draws"],
                  all_control_invariants_valid=control_qc["all_invariants_valid"], repeated_control_seed_values=control_qc["repeated_seed_values"],
                  distinct_full_search_sequences=full_panel["distinct_search_sequences"], pilot_parents=pilot["canonical_parents"],
                  pilot_search_sequences=pilot["distinct_search_sequences"], database_prefix=None, database_accepted=False,
                  searches_started=False, null_precision_established=False,
                  files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "preparation.json", report)
    return report