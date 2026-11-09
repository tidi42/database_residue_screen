import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_table, schema, write_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.example import create_example, save_fasta
from database_residue_screen.panels import LOGICAL_SCHEMA, prepare_panel
from database_residue_screen.sequences import read_fasta


class SearchPanelTests(unittest.TestCase):
    def test_every_metadata_row_joins_the_panel_and_aliases_share_their_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            study = create_example(root / "study")
            inputs = study.parent / "inputs"
            sequences = read_fasta(inputs / "queries.fasta")
            metadata = read_table(inputs / "queries.tsv", schema("queries"))
            sequences["tag_alias"] = sequences["tag"]
            metadata.append(dict(metadata[0], query_id="tag_alias", condition="alias_condition"))
            save_fasta(root / "queries.fasta", sequences)
            write_table(root / "queries.tsv", metadata, schema("queries"))
            generate_controls(root / "queries.fasta", root / "controls", 2)
            result = prepare_panel(root / "queries.fasta", root / "queries.tsv", root / "controls", root / "panel", 2)
            self.assertEqual((result["source_physical_queries"], result["physical_real_queries"], result["canonical_parents"]), (4, 4, 3))
            logical = read_table(root / "panel" / "logical_queries.tsv", LOGICAL_SCHEMA)
            real = {row["logical_id"]: row for row in logical if row["role"] == "real"}
            self.assertEqual(set(real), {"tag", "adapter", "composite", "tag_alias"})
            self.assertEqual(real["tag_alias"]["parent_id"], real["tag"]["parent_id"])
            searched = read_fasta(root / "panel" / "search_queries.fasta")
            self.assertEqual({row["sequence_id"] for row in logical}, set(searched))
            self.assertEqual((result["control_draws"], result["distinct_search_sequences"]), (6, len(searched)))

    def test_draw_identity_survives_sequence_reuse_and_sample_selection(self):
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
            sample = prepare_panel(queries, metadata, root / "controls", root / "sample", 2, 1)
            self.assertEqual(sample["canonical_parents"], 1)
            self.assertEqual(sample["control_draws"], 2)
            with self.assertRaises(ContractError):
                prepare_panel(queries, metadata, root / "controls", root / "too_many", 4)


if __name__ == "__main__":
    unittest.main()