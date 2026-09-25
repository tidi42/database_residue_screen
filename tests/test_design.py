import json
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import read_table, schema, write_table
from database_residue_screen.design import check_design
from database_residue_screen.example import create_example


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name) / "study"
        self.study_path = create_example(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def change_study(self, change):
        value = json.loads(self.study_path.read_text())
        change(value)
        self.study_path.write_text(json.dumps(value), encoding="utf-8")

    def codes(self):
        return {issue["code"] for issue in check_design(self.study_path)["issues"] if issue["severity"] == "blocker"}

    def change_table(self, name, change):
        path = self.root / "inputs" / (name + ".tsv")
        rows = read_table(path, schema("queries" if name == "queries" else name))
        change(rows)
        path.unlink()
        write_table(path, rows, schema(name))

    def test_complete_synthetic_example(self):
        result = check_design(self.study_path)
        self.assertTrue(result["ready_for_freeze"], result["issues"])
        self.assertEqual(result["scope"], "synthetic_example")
        self.assertEqual(result["counts"]["complete_component_decompositions"], 1)

    def test_reference_screen_accepts_declared_unverified_originals_without_components(self):
        self.change_study(lambda value: value.update(profile="reference_screen"))
        self.change_table("queries", lambda rows: [row.update(role="exploratory", provenance_status="unverified",
                                                             architecture_status="unverified") for row in rows])
        self.change_table("sources", lambda rows: rows[0].update(verification="unverified", version=""))
        self.change_table("components", lambda rows: rows.clear())
        self.change_table("validation_cases", lambda rows: rows.clear())
        result = check_design(self.study_path)
        self.assertTrue(result["ready_for_freeze"], result["issues"])
        self.assertEqual(result["counts"]["primary_queries"], 3)
        self.assertEqual(result["component_coverage"], [])
        self.change_study(lambda value: value["search"].update(control_target_limit=5))
        self.assertIn("UNEQUAL_RETRIEVAL", self.codes())
        self.change_table("sources", lambda rows: rows[0].update(evidence_sha256="0" * 64))
        self.assertIn("EVIDENCE_HASH", self.codes())

    def test_synthetic_labels_do_not_authenticate_research(self):
        self.change_study(lambda value: value.update(scope="research"))
        self.assertIn("SYNTHETIC_PRIMARY", self.codes())
        self.assertIn("SYNTHETIC_DATABASE", self.codes())

    def test_unequal_retrieval_blocks_freeze(self):
        self.change_study(lambda value: value["search"].update(control_target_limit=5))
        self.assertIn("UNEQUAL_RETRIEVAL", self.codes())

    def test_missing_input_not_zero(self):
        (self.root / "inputs/validation_subjects.fasta").unlink()
        self.assertIn("INPUT_CONTRACT", self.codes())
        self.assertIn("TRUTH_INPUT_MISSING", self.codes())

    def test_missing_components_blocks_primary_composite(self):
        self.change_table("components", lambda rows: rows.clear())
        self.assertIn("COMPONENT_COVERAGE", self.codes())

    def test_wrong_component_coordinates(self):
        self.change_table("components", lambda rows: rows[0].update(end=11))
        self.assertIn("COMPONENT_SEQUENCE", self.codes())

    def test_draft_cannot_freeze(self):
        self.change_study(lambda value: value.update(decision_status="draft"))
        self.assertIn("DESIGN_APPROVAL", self.codes())

    def test_evidence_digest_is_checked(self):
        self.change_table("sources", lambda rows: rows[0].update(evidence_sha256="0" * 64))
        self.assertIn("EVIDENCE_HASH", self.codes())

    def test_unknown_query_row(self):
        self.change_table("queries", lambda rows: rows.append(dict(rows[0], query_id="absent")))
        self.assertIn("QUERY_MEMBERSHIP", self.codes())

    def test_outside_study_path_rejected(self):
        self.change_study(lambda value: value["inputs"].update(queries="../outside.fasta"))
        self.assertIn("INPUT_CONTRACT", self.codes())

    def test_leakage_when_tuning(self):
        self.change_study(lambda value: value["validation"].update(tuning=True))
        self.change_table("validation_cases", lambda rows: rows.append(dict(rows[0], case_id="development", split="development")))
        self.assertIn("VALIDATION_LEAKAGE", self.codes())

    def test_negative_exact_match_is_not_a_valid_truth_label(self):
        self.change_table("validation_cases", lambda rows: rows[2].update(subject_id="EXAMPLE_POS.1"))
        self.assertIn("EXACT_TRUTH_CONTRADICTION", self.codes())


if __name__ == "__main__":
    unittest.main()