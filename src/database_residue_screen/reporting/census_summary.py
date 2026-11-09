"""Census summary: the descriptive flank census of retained strict hits (census stage), checked and condensed.

The census classifies the flanks of strict reference HSP windows of declared families (integrity-checked local
retrieval, no new search). This report re-checks the table digests recorded in census.json, copies the class summary and
counts windows per family and retrieval status. Unresolved windows are listed, never dropped.
"""
from collections import Counter
from pathlib import Path

from ..contracts import read_json
from .common import read_tsv, verified_files, write_json, write_tsv


def write_census(destination: Path, *, census: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    record = verified_files(census, read_json(census / "census.json"))
    summary_rows = read_tsv(census / "census_summary.tsv")
    write_tsv(destination / "census_summary.tsv", summary_rows,
              list(summary_rows[0]) if summary_rows else ("family", "rule", "census_class", "windows", "accessions"))
    windows = read_tsv(census / "census_windows.tsv") if (census / "census_windows.tsv").is_file() else []
    status = Counter((row.get("family", ""), row.get("status", "")) for row in windows)
    write_tsv(destination / "census_window_status.tsv", [dict(family=key[0], status=key[1], windows=count) for key, count in sorted(status.items())],
              ("family", "status", "windows"))
    summary = dict(windows=record.get("windows", len(windows)), status_counts=record.get("status_counts"),
                   integrity_complete=record.get("integrity_complete"), families=sorted({row.get("family", "") for row in summary_rows}))
    write_json(destination / "summary.json", summary)
    return summary
