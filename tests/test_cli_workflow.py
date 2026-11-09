import contextlib
import io
import json
import shutil
import subprocess
import tempfile
import unittest
import warnings
from pathlib import Path

from database_residue_screen.cli import main
from database_residue_screen.contracts import read_json, read_table
from database_residue_screen.outcomes import CONTRAST_SCHEMA, OUTCOME_SCHEMA
from database_residue_screen.panels import LOGICAL_SCHEMA


def invoke(arguments):
    output, errors = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        status = main(list(map(str, arguments)))
    text = output.getvalue()
    return status, json.loads(text[text.index("{"):]) if "{" in text else None, errors.getvalue()


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST+ unavailable")
class FreshWorkflowCliTests(unittest.TestCase):
    def test_example_controls_panel_search_and_analysis_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def run(*arguments):
                status, result, error = invoke(arguments)
                self.assertEqual(status, 0, error)
                return result

            run("example", root / "example")
            inputs = root / "example/inputs"
            controls = run("controls", inputs / "queries.fasta", "--output", root / "controls", "--draws", 3)
            self.assertTrue(controls["all_invariants_valid"])
            panel = run("prepare-panel", inputs / "queries.fasta", inputs / "queries.tsv", root / "controls",
                        "--output", root / "panel", "--draws", 2)
            self.assertEqual((panel["source_physical_queries"], panel["draws_per_parent"]), (3, 2))
            (root / "database").mkdir()
            subprocess.run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(inputs / "validation_subjects.fasta"),
                            "-out", str(root / "database/example")], capture_output=True, check=True)
            run("inspect-database", root / "database/example", "--output", root / "evidence")
            search = run("search", root / "panel/search_queries.fasta", root / "database/example",
                         "--evidence", root / "evidence/database_identity.json", "--output", root / "search",
                         "--mode", "working", "--threads", 1, "--batch-size", 2)
            self.assertTrue(search["all_queries_completed"])
            self.assertEqual(search["completed_queries"], panel["distinct_search_sequences"])
            analysis = run("analyse-search", root / "panel", root / "search", "--output", root / "analysis")
            self.assertEqual((analysis["result_contract"], analysis["event_contract"]), ("reference_screen_v2", "annotated_v1"))
            self.assertEqual(analysis["length_calibration"], read_json(root / "search/request.json")["length_calibration"])
            logical = read_table(root / "analysis/logical_queries.tsv", LOGICAL_SCHEMA)
            outcomes = read_table(root / "analysis/outcomes.tsv", OUTCOME_SCHEMA)
            self.assertEqual(len(outcomes), len(logical) * 2 * len(analysis["analysis_policy"]["record_populations"]))
            contrasts = read_table(root / "analysis/contrasts.tsv", CONTRAST_SCHEMA)
            union = {row["condition"] for row in contrasts if row["group"] == "Example"}
            self.assertEqual(union, {"all_conditions", "tag", "adapter", "composite"})
            self.assertTrue(any(row["group"] == "Global union" and row["planned_parents"] == 3 for row in contrasts))
            self.assertTrue(run("validate-analysis", root / "analysis")["valid"])
            run("precision-spec", root / "panel", "--output", root / "specification", "--draw-ladder", 1, 2, "--retention-ladder", 500, 2000)
            self.assertEqual(read_json(root / "specification/precision_specification.json")["draw_ladder"], [1, 2])
            precision = run("assess-precision", root / "analysis", root / "specification/precision_specification.json", "--output", root / "precision")
            self.assertIn("conditional_contrasts.tsv", precision["files"])


if __name__ == "__main__":
    unittest.main()
