import json
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.artifacts import freeze_design, validate_outputs, write_preflight
from database_residue_screen.contracts import ContractError, schema, validate_record
from database_residue_screen.example import create_example


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.study = create_example(self.root / "study")
        self.output = self.root / "report"

    def tearDown(self):
        self.temporary.cleanup()

    def test_output_structure_and_contents(self):
        write_preflight(self.study, self.output)
        self.assertTrue(validate_outputs(self.study, self.output)["valid"])

    def test_output_tampering_fails(self):
        write_preflight(self.study, self.output)
        with (self.output / "query_inventory.tsv").open("a") as handle:
            handle.write("extra\n")
        with self.assertRaises(ContractError):
            validate_outputs(self.study, self.output)

    def test_modified_design_fails(self):
        write_preflight(self.study, self.output)
        with self.study.open("a") as handle:
            handle.write("\n")
        with self.assertRaises(ContractError):
            validate_outputs(self.study, self.output)

    def test_validator_code_checksum_is_part_of_contract(self):
        write_preflight(self.study, self.output)
        path = self.output / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest["software"]["code_sha256"] = "0" * 64
        path.write_text(json.dumps(manifest))
        with self.assertRaises(ContractError):
            validate_outputs(self.study, self.output)

    def test_extra_unmanifested_file_fails(self):
        write_preflight(self.study, self.output)
        (self.output / "old_results.tsv").write_text("stale")
        with self.assertRaises(ContractError):
            validate_outputs(self.study, self.output)

    def test_existing_outputs_are_not_overwritten(self):
        write_preflight(self.study, self.output)
        with self.assertRaises(FileExistsError):
            write_preflight(self.study, self.output)

    def test_synthetic_lock_cannot_claim_scientific_readiness(self):
        lock = freeze_design(self.study, self.root / "design.lock.json")
        self.assertEqual(lock["scope"], "synthetic_example")
        self.assertFalse(lock["production_search_ready"])
        self.assertFalse(lock["scientific_results_validated"])

    def test_failed_design_cannot_be_frozen(self):
        value = json.loads(self.study.read_text())
        value["decision_status"] = "draft"
        self.study.write_text(json.dumps(value))
        with self.assertRaises(ContractError):
            freeze_design(self.study, self.root / "design.lock.json")
        self.assertFalse((self.root / "design.lock.json").exists())

    def test_failed_search_cannot_be_reported_as_successful_zero(self):
        row = dict(query_id="tag", role="real", parent_id="", query_sha256="0" * 64, search_id="search", batch_id="batch",
                   status="success_zero", exit_code="1", retained_rows=0, parameters_sha256="0" * 64,
                   database_sha256="0" * 64, error="failed")
        with self.assertRaises(ContractError):
            validate_record(row, schema("search_status"))


if __name__ == "__main__":
    unittest.main()