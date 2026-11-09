"""Stage `record_cases`: one evidence-complete, sequence-verified record case per declared category."""
from pathlib import Path

from ...contracts import read_table
from ...materials import SEGMENT_SCHEMA
from ...record_context import select_record_cases, selection_specification
from ...sequences import read_fasta


def record_cases(context):
    """Deterministic case selection with local flank retrieval and verification; failures are recorded."""
    panel = context.upstream["panel"].output
    report = select_record_cases(context.upstream["taxonomy"].output / "analysis", context.upstream["metadata"].output / "metadata",
                                 context.upstream["search"].output / "search", read_fasta(panel / "queries.fasta"),
                                 read_table(panel / "segments.tsv", SEGMENT_SCHEMA),
                                 context.upstream["database"].output / "acceptance" / "accepted.json", Path(context.tool("blastdbcmd")),
                                 context.output / "record_cases", selection_specification(context.config["record_cases"]))
    return dict(report["counts"], local_retrieval_attempts=report["local_retrieval_attempts"],
                categories={row["category"]: row["status"] for row in report["category_status"]})
