"""Bind executable settings to the declared study; other stages remain explicit."""
from pathlib import Path
from .contracts import ContractError, read_json, schema, validate_record
from .design import digest


def bind_study(path: Path, settings: dict, mode: str, maximum_evalue, informative_limit=.05, testable_limit=.01) -> dict:
    study = read_json(path)
    validate_record(study, schema("study"))
    declared = study["search"]
    expected = dict(task=settings["task"], word_size=settings["word_size"], evalue=settings["evalue"],
                    dust=settings["dust"] != "no", strand=settings["strand"],
                    real_target_limit=settings["max_target_seqs"], control_target_limit=settings["max_target_seqs"])
    if any(declared.get(key) != value for key, value in expected.items()):
        raise ContractError("Search request disagrees with study.json parameters")
    if declared["failure_policy"] != "abort_run_resume_retries_failed_batches":
        raise ContractError("Unsupported study failure_policy; runner aborts and retries on resume, never blocks a parent silently")
    if declared.get("threshold_mode", "fixed") != mode or declared.get("maximum_evalue") != maximum_evalue:
        raise ContractError("Length threshold policy disagrees with study.json")
    if (declared.get("search_space") != settings.get("searchsp") or
            declared.get("informativeness_limit", .05) != informative_limit or
            declared.get("testable_chance_limit", .01) != testable_limit):
        raise ContractError("Testability or fixed search space disagrees with study.json")
    return dict(sha256=digest(path), search=declared, analysis=study["analysis"],
                scope="Search settings enforced; analysis positions checked at analyse-search. Pilot/control budgets require their own specifications.")
