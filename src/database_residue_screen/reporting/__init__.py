"""Publication analysis layer: central metrics, reports, figure sources and figures (see docs/reports.md).

Every function takes explicit input directories and writes one new output directory with a deterministic manifest.
Import the submodules directly (``from database_residue_screen.reporting.metrics import compute_metrics``); this
module imports nothing, so a stage's code closure holds only the reporting modules it uses.
"""
