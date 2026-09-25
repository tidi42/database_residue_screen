import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_table, schema, write_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.example import create_example, save_fasta
from database_residue_screen.panels import LOGICAL_SCHEMA, prepare_panel, select_original_queries
from database_residue_screen.sequences import read_fasta


class SearchPanelTests(unittest.TestCase):
    def test_original_selection_excludes_reference_alias_without_losing_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            study = create_example(root / "study")
            inputs = study.parent / "inputs"
            sequences = read_fasta(inputs / "queries.fasta")
            metadata = read_table(inputs / "queries.tsv", schema("queries"))
            sequences["added_reference"] = sequences["tag"]
            metadata.append(dict(metadata[0], query_id="added_reference", role="reference"))
            save_fasta(root / "queries.fasta", sequences)
            write_table(root / "queries.tsv", metadata, schema("queries"))
            generate_controls(root / "queries.fasta", root / "controls", 2)
            result = prepare_panel(root / "queries.fasta", root / "queries.tsv", root / "controls", root / "panel", 2)
            self.assertEqual(result["physical_real_queries"], 3)
            self.assertEqual(result["canonical_parents"], 3)
            self.assertEqual(result["excluded_reference_ids"], ["added_reference"])
            logical = read_table(root / "panel" / "logical_queries.tsv", LOGICAL_SCHEMA)
            self.assertNotIn("added_reference", {row["logical_id"] for row in logical})
            selected = select_original_queries(root / "queries.fasta", root / "queries.tsv", root / "originals")
            self.assertEqual(selected["original_physical_queries"], 3)
            original_sequences = read_fasta(root / "originals" / "queries.fasta")
            self.assertEqual(original_sequences, {name: sequence for name, sequence in sequences.items() if name != "added_reference"})
            inclusive = prepare_panel(root / "queries.fasta", root / "queries.tsv", root / "controls", root / "all", 2, selection="all")
            self.assertEqual(inclusive["physical_real_queries"], 4)

    def test_draw_identity_survives_sequence_reuse_and_pilot_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            study = create_example(root / "study")
            queries, metadata = study.parent / "inputs" / "queries.fasta", study.parent / "inputs" / "queries.tsv"
            generate_controls(queries, root / "controls", 3)
            full = prepare_panel(queries, metadata, root / "controls", root / "full", 2)
            self.assertEqual(full["canonical_parents"], 3)
            self.assertEqual(full["control_draws"], 6)
            logical = read_table(root / "full" / "logical_queries.tsv", LOGICAL_SCHEMA)
            searched = read_fasta(root / "full" / "search_queries.fasta")
            self.assertEqual({row["sequence_id"] for row in logical}, set(searched))
            pilot = prepare_panel(queries, metadata, root / "controls", root / "pilot", 2, 1)
            self.assertEqual(pilot["canonical_parents"], 1)
            self.assertEqual(pilot["control_draws"], 2)
            with self.assertRaises(ContractError):
                prepare_panel(queries, metadata, root / "controls", root / "too_many", 4)


if __name__ == "__main__":
    unittest.main()