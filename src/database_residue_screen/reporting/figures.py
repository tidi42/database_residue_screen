"""Figures: run the study's R figure build on one figure-source directory.

The build script (studies/paper3/figures/build_figures.sh) verifies the figure sources, renders every figure with the
configured Rscript, checks the artwork manifests, writes manuscript-numbered copies of the main figures and
figures_manifest.json (contract publication_figures_v1). PNG files, CSV tables and manifests are byte-deterministic;
PDFs carry a content digest that excludes cairo's embedded creation time.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

from ..contracts import ContractError, read_json
from .common import new_directory

CONTRACT = "publication_figures_v1"
SOURCE_CONTRACT = "publication_figure_sources_v1"
BUILD_SCRIPT = "build_figures.sh"


def build_figures(source: Path, destination: Path, *, scripts: Path, rscript: str, python: str | None = None, log=None) -> dict:
    source, destination, scripts = Path(source), Path(destination), Path(scripts)
    if read_json(source / "manifest.json").get("contract") != SOURCE_CONTRACT:
        raise ContractError(f"{source} is not a {SOURCE_CONTRACT} figure-source directory")
    script = scripts / BUILD_SCRIPT
    bash = shutil.which("bash")
    if not script.is_file() or bash is None:
        raise ContractError(f"Figure build needs bash and {script}")
    new_directory(destination)
    result = subprocess.run([bash, str(script), str(source), str(destination), str(rscript), python or sys.executable],
                            capture_output=True, text=True, env=dict(os.environ, LC_ALL="C.UTF-8"))
    if log is not None:
        for line in (result.stdout + result.stderr).splitlines():
            if line.strip():
                log(line)
    if result.returncode:
        raise ContractError(f"Figure build failed (exit {result.returncode}): " + "\n".join((result.stderr or result.stdout).splitlines()[-15:]))
    manifest = read_json(destination / "figures_manifest.json")
    if manifest.get("contract") != CONTRACT:
        raise ContractError("The figure build did not write a publication_figures_v1 manifest")
    return dict(contract=CONTRACT, figures={name: figure["manuscript_figure"] for name, figure in manifest["figures"].items()},
                files=len(manifest["files"]), R_version=manifest["R_version"])
