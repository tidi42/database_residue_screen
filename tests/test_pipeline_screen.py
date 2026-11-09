"""The `screen` pipeline through the study engine: pinned sources, panel, controls and design on synthetic sources."""
import hashlib
import json
import re
import shutil
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.constructs import DEFINITION_SCHEMA, build_panel, collect_sources
from database_residue_screen.contracts import ContractError, read_json, write_table
from database_residue_screen.materials import select_subset
from database_residue_screen.sequences import canonical_id, read_fasta
from database_residue_screen.study.runner import run_study
from test_constructs import definition_rows, write_sources

STUDY = Path(__file__).resolve().parents[1] / "studies" / "paper3" / "study.toml"
TOOLS_PRESENT = all(shutil.which(name) for name in ("blastn", "blastdbcmd", "blastdbcheck", "makeblastdb", "blast_formatter", "cc",
                                                     "pdftotext", "Rscript"))


def write_fixture(root: Path) -> tuple[Path, Path, dict[str, str]]:
    sources = root / "remote"
    write_sources(sources)
    rows, expected = definition_rows()
    definition = root / "panel_definition.tsv"
    write_table(definition, rows, DEFINITION_SCHEMA)
    manifest = [dict(name=path.name, url=path.resolve().as_uri(), retrieved="2026-11-09T00:00:00+00:00",
                     sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size, purpose="synthetic fixture",
                     redistribution="synthetic test data") for path in sorted(sources.iterdir())]
    (root / "sources.json").write_text(json.dumps(dict(kind="panel_sources", sources=manifest)), encoding="utf-8")
    return definition, root / "sources.json", expected


def study_config(root: Path, definition: Path, manifest: Path, query_ids: list[str]) -> Path:
    """The Paper 3 configuration with the synthetic panel in place of the Paper 3 panel and its materials."""
    text = STUDY.read_text(encoding="utf-8")
    panel = f'[panel]\ndefinition = "{definition}"\n\n'
    text = re.sub(r"\[panel\]\n.*?\n(?=\[subset\])", lambda match: panel, text, flags=re.S)
    text = text.replace('manifest = "sources.json"', f'manifest = "{manifest}"')
    text = re.sub(r"\[subset\]\nparents = 0\nquery_ids = \[\]", "[subset]\nparents = 0\nquery_ids = " + json.dumps(query_ids), text)
    text = text.replace("allow_source_drift = true", "allow_source_drift = false")
    path = root / "study.toml"
    path.write_text(text, encoding="utf-8")
    return path


class SubsetTests(unittest.TestCase):
    SEQUENCES = {"a": "ACGTACGTAA", "a_alias": "ACGTACGTAA", "b": "TTGACCAGGA", "c": "GGCATCAGTT"}

    def test_complete_panel_parents_and_query_ids(self):
        selected, record = select_subset(self.SEQUENCES)
        self.assertEqual((selected, record["mode"], record["selected_parents"]), (list(self.SEQUENCES), "complete_panel", 3))
        smallest = sorted({canonical_id(sequence) for sequence in self.SEQUENCES.values()})[0]
        selected, record = select_subset(self.SEQUENCES, parents=1)
        self.assertEqual(record["parent_ids"], [smallest])
        self.assertEqual({canonical_id(self.SEQUENCES[name]) for name in selected}, {smallest})
        selected, record = select_subset(self.SEQUENCES, query_ids=["c", "a_alias"])
        self.assertEqual((selected, record["mode"]), (["a_alias", "c"], "query_ids"))

    def test_refusals(self):
        for arguments in (dict(parents=1, query_ids=["a"]), dict(parents=-1), dict(parents=4), dict(query_ids=["absent"])):
            with self.subTest(arguments=arguments), self.assertRaises(ContractError):
                select_subset(self.SEQUENCES, **arguments)


class SourceDriftTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.definition, self.manifest, self.expected = write_fixture(self.root)

    def change_layout(self, directory: Path) -> str:
        """A byte change of an HTML source that leaves every derived sequence unchanged."""
        page = next(path for path in sorted(directory.iterdir()) if path.suffix == ".html")
        page.write_text(page.read_text(encoding="utf-8").replace("<p>", "<p class='x'>", 1), encoding="utf-8")
        return page.name

    def test_pinned_snapshot_and_refused_or_accepted_drift(self):
        local = self.root / "local"
        shutil.copytree(self.root / "remote", local)
        snapshot = collect_sources(self.manifest, local, self.root / "pinned", fetch=False)
        self.assertEqual(read_json(self.root / "pinned" / "collection.json")["drift"], [])
        self.assertEqual(snapshot["drift"], [])
        name = self.change_layout(local)
        with self.assertRaises(ContractError):
            collect_sources(self.manifest, local, self.root / "refused", fetch=False)
        accepted = collect_sources(self.manifest, local, self.root / "drift", fetch=False, allow_drift=True)
        self.assertEqual([row["name"] for row in accepted["drift"]], [name])
        snapshot = self.root / "drift"
        build = build_panel(self.definition, snapshot / "files", self.root / "panel", snapshot / "sources.json", shutil.which("pdftotext"))
        self.assertEqual(build["queries"], len(self.expected))
        self.assertEqual(read_fasta(self.root / "panel" / "queries.fasta"), self.expected)


@unittest.skipUnless(TOOLS_PRESENT, "BLAST+, cc, pdftotext or Rscript unavailable")
class ScreenPipelineTests(unittest.TestCase):
    def test_front_stages_through_the_engine(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            definition, manifest, expected = write_fixture(root)
            subset = ["BC1001", "BC1001_hairpin", "PRIMER"]
            config = study_config(root, definition, manifest, subset)
            environ = dict(CORE_NT_DB=str(root / "db"), CORE_NT_TAXDB=str(root), CORE_NT_DOWNLOAD_EVIDENCE=str(root),
                           BLAST_BIN=str(Path(shutil.which("blastn")).parent), BLAST_SOURCE_TARBALL=str(root / "src.tar.gz"),
                           TAXDUMP_DIR=str(root), PANEL_SOURCES=str(root / "local"), NCBI_EMAIL="someone@example.org",
                           PATH=str(Path(shutil.which("pdftotext")).parent) + ":" + str(Path(shutil.which("cc")).parent) + ":"
                           + str(Path(shutil.which("Rscript")).parent))
            run = root / "run"
            result = run_study(config, run, until="design", mode="working", environ=environ)
            self.assertEqual([row["stage"] for row in result["stages"] if row["status"] == "completed"],
                             ["sources", "panel", "controls", "design"])
            panel = next(run.glob("02_panel/attempt_*/output"))
            self.assertEqual(read_fasta(panel / "queries.fasta"), {name: expected[name] for name in subset})
            report = read_json(panel / "panel.json")
            self.assertEqual((report["build"]["queries"], report["selection"]["selected_queries"]), (len(expected), 3))
            self.assertTrue((root / "local").is_dir() and sorted(path.name for path in (root / "local").iterdir())
                            == sorted(path.name for path in (root / "remote").iterdir()))
            design = next(run.glob("04_design/attempt_*/output"))
            self.assertEqual(read_json(design / "search_panel" / "panel.json")["canonical_parents"],
                             len({canonical_id(expected[name]) for name in subset}))
            again = run_study(config, run, until="design", mode="working", environ=environ)
            self.assertTrue(all(row["action"] == "skipped" for row in again["stages"] if row["action"] != "not_selected"))


if __name__ == "__main__":
    unittest.main()
