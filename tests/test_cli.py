import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.cli import main


class CliTests(unittest.TestCase):
    def invoke(self, arguments):
        output = io.StringIO()
        errors = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            status = main(list(map(str, arguments)))
        return status, output.getvalue(), errors.getvalue()

    def test_documented_synthetic_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.assertEqual(self.invoke(["example", root / "example"])[0], 0)
            study = root / "example/study.json"
            self.assertEqual(self.invoke(["preflight", study, "--output", root / "report"])[0], 0)
            status, result, error = self.invoke(["validate-output", study, root / "report"])
            self.assertEqual(status, 0, error)
            self.assertTrue(json.loads(result)["valid"])
            self.assertTrue(json.loads(result)["design_ready"])

    def test_invalid_input_returns_nonzero_without_traceback(self):
        status, output, error = self.invoke(["preflight", "/does/not/exist/study.json"])
        self.assertEqual(status, 1)
        self.assertIn("Error:", error)
        self.assertNotIn("Traceback", error)

    def test_command_set_and_explicit_precision_ladders(self):
        import argparse
        from database_residue_screen.cli import parser
        commands = next(action for action in parser()._actions if isinstance(action, argparse._SubParsersAction)).choices
        self.assertEqual(set(commands), {
            "example", "preflight", "validate-output", "validate-table", "inspect-database", "accept-database",
            "audit-database-integrity", "controls", "prepare-panel", "precision-spec", "search", "analyse-search", "fetch-metadata",
            "annotate-events", "pin-taxonomy", "annotate-taxonomy", "fetch-insdc", "positional-sensitivity", "cap-review",
            "flank-census", "exact-scan", "calibrate-benchmark", "benchmark", "assess-precision", "validate-analysis",
            "figure-inputs", "fetch-sources", "build-panel", "run-study", "study-status"})
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser().parse_args(["precision-spec", "panel", "--output", "out"])
        args = parser().parse_args(["precision-spec", "panel", "--output", "out", "--draw-ladder", "5", "10", "15",
                                    "--retention-ladder", "500", "2000"])
        self.assertEqual((args.draw_ladder, args.retention_ladder, args.maximum_draws), ([5, 10, 15], [500, 2000], None))

    def test_blocked_research_design_returns_two(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.invoke(["example", root / "example"])
            study = root / "example/study.json"
            value = json.loads(study.read_text())
            value["scope"] = "research"
            study.write_text(json.dumps(value))
            self.assertEqual(self.invoke(["preflight", study, "--output", root / "report"])[0], 2)
            status, result, error = self.invoke(["validate-output", study, root / "report"])
            self.assertEqual(status, 0, error)
            self.assertFalse(json.loads(result)["design_ready"])


if __name__ == "__main__":
    unittest.main()