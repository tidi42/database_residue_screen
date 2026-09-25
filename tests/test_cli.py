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
            self.assertEqual(self.invoke(["freeze", study, "--output", root / "design.lock.json"])[0], 0)

    def test_invalid_input_returns_nonzero_without_traceback(self):
        status, output, error = self.invoke(["preflight", "/does/not/exist/study.json"])
        self.assertEqual(status, 2)
        self.assertIn("Error:", error)
        self.assertNotIn("Traceback", error)

    def test_blocked_research_does_not_write_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.invoke(["example", root / "example"])
            study = root / "example/study.json"
            value = json.loads(study.read_text())
            value["scope"] = "research"
            study.write_text(json.dumps(value))
            self.assertEqual(self.invoke(["preflight", study, "--output", root / "report"])[0], 2)
            self.assertEqual(self.invoke(["freeze", study, "--output", root / "bad.lock.json"])[0], 2)
            self.assertFalse((root / "bad.lock.json").exists())


if __name__ == "__main__":
    unittest.main()