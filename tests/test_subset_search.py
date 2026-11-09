import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_json, read_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.database import inspect_database
from database_residue_screen.design import digest
from database_residue_screen.example import create_example, save_fasta
from database_residue_screen.outcomes import OUTCOME_SCHEMA
from database_residue_screen.panels import prepare_panel
from database_residue_screen.analysis import analyse_search
from database_residue_screen.release import validate_analysis
from database_residue_screen.search import STATUS_SCHEMA, run_search, search_settings, validate_completed_batch
from database_residue_screen.sequences import read_fasta
from database_residue_screen.subset_search import assemble_subset, write_query_sets


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST tools unavailable")
class SubsetSearchTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        inputs = create_example(self.root / "study").parent / "inputs"
        generate_controls(inputs / "queries.fasta", self.root / "controls", 3)
        for draws in (1, 3):
            prepare_panel(inputs / "queries.fasta", inputs / "queries.tsv", self.root / "controls", self.root / f"panel{draws}", draws)
        self.prefix = self.root / "database"
        subprocess.run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(inputs / "validation_subjects.fasta"),
                        "-out", str(self.prefix)], capture_output=True, check=True)
        inspect_database(self.prefix, self.root / "evidence")
        self.evidence = self.root / "evidence/database_identity.json"
        self.source = self.root / "source"
        run_search(self.root / "panel3/search_queries.fasta", self.prefix, self.evidence, self.source, threads=1, batch_size=2, mode="working")

    def test_subset_matches_fresh_search_with_complete_provenance_and_zeros(self):
        before = {str(path): digest(path) for path in self.source.rglob("*") if path.is_file()}
        expected = read_fasta(self.root / "panel1/search_queries.fasta")
        result = assemble_subset(self.root / "panel1", [self.source], self.root / "assembled")
        self.assertEqual(result["completed_queries"], len(expected))
        run_search(self.root / "panel1/search_queries.fasta", self.prefix, self.evidence, self.root / "fresh", threads=2, batch_size=8, mode="working")
        for name in ("assembled", "fresh"):
            analyse_search(self.root / "panel1", self.root / name, self.root / (name + "_analysis"))
        self.assertEqual(read_table(self.root / "assembled_analysis/outcomes.tsv", OUTCOME_SCHEMA),
                         read_table(self.root / "fresh_analysis/outcomes.tsv", OUTCOME_SCHEMA))
        self.assertTrue(validate_analysis(self.root / "assembled_analysis")["valid"])
        self.assertEqual(before, {str(path): digest(path) for path in self.source.rglob("*") if path.is_file()})
        request = read_json(self.root / "assembled/request.json")
        batch = self.root / "assembled" / result["batches"][0]["directory"]
        (batch / "blast.tsv").write_text((batch / "blast.tsv").read_text() + "\n")
        with self.assertRaisesRegex(ContractError, "artifact changed"):
            validate_completed_batch(batch, request)

    def test_incomplete_source_never_becomes_a_zero_or_completed_panel(self):
        for path in self.source.glob("batch_*.attempt_*/complete.json"):
            path.rename(path.with_name("preserved_incomplete.json"))
        with self.assertRaisesRegex(ContractError, "without completed evidence"):
            assemble_subset(self.root / "panel1", [self.source], self.root / "invalid")

    def test_dedicated_reference_and_control_searches_cover_the_panel_once(self):
        sets = write_query_sets(self.root / "panel1", self.root / "sets")
        panel = read_fasta(self.root / "panel1/search_queries.fasta")
        references, controls = (read_fasta(self.root / "sets" / name) for name in ("references.fasta", "controls.fasta"))
        self.assertEqual((sets["counts"]["references"], sets["counts"]["controls"]), (len(references), len(controls)))
        self.assertEqual(dict(references, **controls), panel)
        self.assertFalse(references.keys() & controls.keys())
        sources = []
        for name in ("references", "controls"):
            run_search(self.root / "sets" / (name + ".fasta"), self.prefix, self.evidence, self.root / (name + "_search"), threads=2,
                       batch_size=8, mode="working")
            sources.append(self.root / (name + "_search"))
        result = assemble_subset(self.root / "panel1", sources, self.root / "assembled")
        rows = read_table(self.root / "assembled/status.tsv", STATUS_SCHEMA)
        self.assertEqual(len(rows), len({row["query_id"] for row in rows}))
        self.assertEqual({row["query_id"] for row in rows}, panel.keys())
        self.assertEqual(result["submitted_queries"], result["completed_queries"])
        run_search(self.root / "panel1/search_queries.fasta", self.prefix, self.evidence, self.root / "fresh", threads=1, mode="working")
        for name in ("assembled", "fresh"):
            analyse_search(self.root / "panel1", self.root / name, self.root / (name + "_analysis"))
        self.assertEqual(read_table(self.root / "assembled_analysis/outcomes.tsv", OUTCOME_SCHEMA),
                         read_table(self.root / "fresh_analysis/outcomes.tsv", OUTCOME_SCHEMA))

    def test_target_limit_arms_search_identical_batches(self):
        write_query_sets(self.root / "panel3", self.root / "sets")
        references = self.root / "sets/references.fasta"
        arms = {}
        for limit in (500, 2000):
            arms[limit] = run_search(references, self.prefix, self.evidence, self.root / f"references_{limit}", search_settings(limit),
                                     threads=1, batch_size=2, mode="working")
        batches = {limit: [read_json(self.root / f"references_{limit}" / batch["directory"] / "complete.json")["query_ids"]
                           for batch in report["batches"]] for limit, report in arms.items()}
        self.assertEqual(batches[500], batches[2000])
        self.assertGreater(len(batches[500]), 1)
        requests = {limit: read_json(self.root / f"references_{limit}/request.json") for limit in arms}
        self.assertEqual({key: value for key, value in requests[500]["settings"].items() if key != "max_target_seqs"},
                         {key: value for key, value in requests[2000]["settings"].items() if key != "max_target_seqs"})
        self.assertEqual(requests[500]["length_calibration"], requests[2000]["length_calibration"])
        self.assertEqual((requests[500]["queries_sha256"], requests[500]["batch_size"], requests[500]["threads"]),
                         (requests[2000]["queries_sha256"], requests[2000]["batch_size"], requests[2000]["threads"]))

    def test_rehashed_projection_edit_still_disagrees_with_original(self):
        result = assemble_subset(self.root / "panel1", [self.source], self.root / "assembled")
        request = read_json(self.root / "assembled/request.json")
        directory = self.root / "assembled" / result["batches"][0]["directory"]
        output = directory / "blast.tsv"
        output.write_text(output.read_text() + "# fabricated annotation\n")
        record = read_json(directory / "complete.json")
        record["files"]["blast.tsv"] = digest(output)
        record["output_bytes"] = output.stat().st_size
        (directory / "complete.json").write_text(json.dumps(record))
        with self.assertRaisesRegex(ContractError, "differ from completed source"):
            validate_completed_batch(directory, request)

    def test_cap_change_and_quarantined_source_are_rejected(self):
        run_search(self.root / "panel1/search_queries.fasta", self.prefix, self.evidence, self.root / "different",
                   search_settings(2000), threads=1, mode="working")
        with self.assertRaisesRegex(ContractError, "same database"):
            assemble_subset(self.root / "panel1", [self.source, self.root / "different"], self.root / "mixed")
        (self.source / "quarantine.json").write_text("{}")
        with self.assertRaisesRegex(ContractError, "quarantined"):
            assemble_subset(self.root / "panel1", [self.source], self.root / "quarantined")


if __name__ == "__main__":
    unittest.main()