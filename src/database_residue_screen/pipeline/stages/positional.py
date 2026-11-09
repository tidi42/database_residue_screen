"""Stage `positional`: terminal-event robustness (one event per accession or BioProject, INSDC exclusion, windows)."""
import json

from ...contracts import write_table
from ...design import digest
from ...insdc import read_annotation
from ...positional import SENSITIVITY_SCHEMA, positional_sensitivity, projects_from_metadata


def positional(context):
    """Positional sensitivity rows of the final analysis with the archived BioProjects and INSDC annotation."""
    analysis = context.upstream["taxonomy"].output / "analysis"
    window = context.config["positional"]["diagnostic_window"]
    rows = positional_sensitivity(analysis, projects_from_metadata(context.upstream["metadata"].output / "metadata"),
                                  read_annotation(context.upstream["insdc"].output / "insdc"), window)
    write_table(context.output / "positional_sensitivity.tsv", rows, SENSITIVITY_SCHEMA)
    policy = json.loads((analysis / "analysis.json").read_text(encoding="utf-8"))["analysis_policy"]["position"]
    report = dict(kind="positional_sensitivity_v1", analysis_sha256=digest(analysis / "analysis.json"), position_policy=policy,
                  diagnostic_window=window, rows=len(rows), variants=sorted({row["variant"] for row in rows}),
                  files={"positional_sensitivity.tsv": digest(context.output / "positional_sensitivity.tsv")})
    (context.output / "positional.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return dict(rows=len(rows), variants=report["variants"], primary_window=policy["primary_window"],
                primary_minimum_length=policy["primary_minimum_length"], diagnostic_window=window)
