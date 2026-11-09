#!/usr/bin/env bash
# Render the publication figures from one figure-source directory (contract publication_figure_sources_v1, written by the
# figure_sources stage): plot_figures.R draws the R figures figure1-figure4 and figureS1-figureS3, plot_record_cases.R the
# record cases (R name figureS4). The script verifies the source manifest, renders, checks both artwork manifests, writes
# manuscript-numbered copies of the main figures and one figures_manifest.json.
#
# Manuscript numbering: the graphical abstract is Figure 1 (not rendered here), so figure1-figure4 are Figures 2-5;
# supplementary figures follow the order of first citation: figureS1 = Figure S1, figureS4 (record cases) = Figure S2,
# figureS3 = Figure S3, figureS2 = Figure S4. File names keep the R names; main figures are also copied to
# paper3_manuscript_FigureN.* (byte-identical).
#
# PNG files, CSV tables and manifests are byte-deterministic. cairo writes the current time into every PDF, so the manifest
# records a PDF content digest (decompressed content streams with the creation time removed) instead of the file digest.
# Usage: build_figures.sh SOURCE_DIR OUTPUT_DIR [RSCRIPT] [PYTHON]
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE=$(realpath "$1")
OUTPUT=$2
RSCRIPT=${3:-Rscript}
PYTHON=${4:-python3}
if [[ -e "$OUTPUT" && -n "$(ls -A "$OUTPUT")" ]]; then
  printf 'Refusing to write into non-empty %s; choose a new output directory.\n' "$OUTPUT" >&2
  exit 2
fi
"$PYTHON" - "$SOURCE" <<'EOF'
import hashlib, json, sys
from pathlib import Path
source = Path(sys.argv[1])
manifest = json.loads((source / "manifest.json").read_text())
assert manifest["contract"] == "publication_figure_sources_v1", manifest["contract"]
for name, checksum in manifest["files"].items():
    assert "/" not in name and hashlib.sha256((source / name).read_bytes()).hexdigest() == checksum, name
print(f"Verified figure source: {len(manifest['files'])} files")
EOF
mkdir -p "$OUTPUT"
export LC_ALL=C.UTF-8
"$RSCRIPT" "$HERE/plot_figures.R" "$SOURCE" "$HERE/theme_graphpad.R" "$OUTPUT" all
"$RSCRIPT" "$HERE/plot_record_cases.R" "$SOURCE" "$OUTPUT" compact
"$PYTHON" - "$SOURCE" "$OUTPUT" "$HERE" <<'EOF'
import hashlib, json, re, shutil, sys, zlib
from pathlib import Path

source, output, scripts = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sha256 = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pdf_content_sha256(path):
    """SHA-256 over every decompressed stream except cross-reference streams (byte offsets), creation time blanked."""
    data, digest = Path(path).read_bytes(), hashlib.sha256()
    for match in re.finditer(rb"(?<!end)stream\r?\n", data):
        dictionary = data[data.rfind(b" obj", 0, match.start()):match.start()]
        if b"/XRef" in dictionary:
            continue
        chunk = data[match.end():data.find(b"endstream", match.end())]
        try:
            chunk = zlib.decompress(chunk)
        except zlib.error:
            pass
        digest.update(re.sub(rb"/(CreationDate|ModDate) \([^)]*\)", rb"/\1 ()", chunk))
    return digest.hexdigest()


source_sha = sha256(source / "manifest.json")
artwork = json.loads((output / "paper3_figure_manifest.json").read_text())
cases = json.loads((output / "paper3_case_figure_manifest.json").read_text())
assert artwork["contract"] == "publication_figure_artwork_v1", artwork["contract"]
assert cases["contract"] == "publication_record_case_artwork_v1", cases["contract"]
for record in (artwork, cases):
    assert record["source_manifest_sha256"] == source_sha
    for name, checksum in record["files"].items():
        assert sha256(output / name) == checksum, name
for name, checksum in artwork["tables"].items():
    assert sha256(output / "tables" / name) == checksum, name
figures = dict(artwork["figures"])
figures["figureS4"] = dict(cases["figures"]["figureS4"], manuscript_figure="Figure S2")
assert sorted(figures) == ["figure1", "figure2", "figure3", "figure4", "figureS1", "figureS2", "figureS3", "figureS4"], sorted(figures)
assert all(figure["height_mm"] <= 240 for figure in figures.values())
copies = {}
for name, figure in sorted(figures.items()):
    number = figure["manuscript_figure"].replace(" ", "")
    for extension in ("png", "pdf"):
        assert (output / f"paper3_{name}.{extension}").is_file(), name
    if not number.startswith("FigureS"):
        for extension in ("png", "pdf"):
            original, copy = output / f"paper3_{name}.{extension}", output / f"paper3_manuscript_{number}.{extension}"
            assert not copy.exists(), copy
            shutil.copyfile(original, copy)
            copies[copy.name] = original.name
files = {}
for path in sorted(output.rglob("*")):
    if path.is_file():
        relative = str(path.relative_to(output))
        files[relative] = dict(content_sha256=pdf_content_sha256(path)) if path.suffix == ".pdf" else dict(sha256=sha256(path))
code = {path.name: sha256(path) for path in sorted(scripts.iterdir()) if path.is_file()}
manifest = dict(contract="publication_figures_v1", source_manifest_sha256=source_sha,
                artwork_manifest_sha256=sha256(output / "paper3_figure_manifest.json"),
                record_case_manifest_sha256=sha256(output / "paper3_case_figure_manifest.json"),
                R_version=artwork["R_version"], package_versions=artwork["package_versions"], code=code,
                figures={name: dict(manuscript_figure=figure["manuscript_figure"], height_mm=figure["height_mm"],
                                    width_mm=figure.get("width_mm", 180)) for name, figure in sorted(figures.items())},
                manuscript_copies=copies, graphical_abstract="Figure 1 (graphical abstract) is not rendered from the figure sources",
                pdf_digest="PDF entries carry content_sha256: decompressed content streams (cross-reference streams excluded) with the "
                           "embedded creation time removed",
                files=files)
(output / "figures_manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n", encoding="utf-8")
print(f"Verified artwork: {len(figures)} figures, {len(artwork['tables'])} CSV tables, {len(copies)} manuscript-numbered copies, "
      f"source {source_sha[:12]}")
EOF
