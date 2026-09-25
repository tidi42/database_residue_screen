import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError
from database_residue_screen.example import create_example, save_fasta
from database_residue_screen.legacy import PANEL_FILES
from database_residue_screen.preparation import prepare_reference_inputs
from database_residue_screen.sequences import read_fasta


class ReferencePreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.inputs = create_example(self.root / "imported").parent / "inputs"
        (self.root / "panels").mkdir()
        sequences = read_fasta(self.inputs / "queries.fasta")
        for index, name in enumerate(PANEL_FILES):
            identifier = list(sequences)[index % len(sequences)]
            save_fasta(self.root / "panels" / name, {identifier: sequences[identifier]})

    def test_preparation_has_no_database_or_search_side_effect(self):
        report = prepare_reference_inputs(self.root, self.inputs, self.root / "prepared", draws=2)
        self.assertEqual(report["original_physical_queries"], 3)
        self.assertEqual(report["total_control_draws"], 6)
        self.assertTrue(report["all_control_invariants_valid"])
        self.assertIsNone(report["database_prefix"])
        self.assertFalse(report["searches_started"])
        self.assertFalse(report["final_freeze"])

    def test_changed_author_panel_is_not_silently_replaced_by_imported_data(self):
        (self.root / "panels" / PANEL_FILES[0]).write_text(">changed_reference\nACGTACGT\n")
        with self.assertRaisesRegex(ContractError, "no longer match"):
            prepare_reference_inputs(self.root, self.inputs, self.root / "blocked", draws=2)
        self.assertFalse((self.root / "blocked").exists())


if __name__ == "__main__":
    unittest.main()