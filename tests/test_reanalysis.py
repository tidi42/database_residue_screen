import unittest
import tempfile
from pathlib import Path

from database_residue_screen.contracts import ContractError
from database_residue_screen.example import save_json
from database_residue_screen.reanalysis import analyse_search, event, orientation_key, position_summary
from database_residue_screen.sequences import reverse_complement


class EventAnalysisTests(unittest.TestCase):
    def test_quarantined_search_cannot_be_aggregated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_json(root / "quarantine.json", dict(reason="mixed database"))
            with self.assertRaisesRegex(ContractError, "quarantined"):
                analyse_search(root / "missing-panel", root, root / "output")
            self.assertFalse((root / "output").exists())

    def setUp(self):
        self.sequence = "ACGTTAGCTGATCGTACGTC"
        self.row = dict(qseqid="query", sseqid="ref|TEST0001.1|", staxids="0", stitle="record", pident=100,
                        length=20, qlen=20, gapopen=0, mismatch=0, qstart=1, qend=20, sstart=1, send=20, slen=100)

    def test_reverse_alias_deduplication_retains_opposite_biological_strands(self):
        first = event(self.row, {"group": "Example"}, self.sequence, "test", 1)
        reverse = event(dict(self.row, sstart=20, send=1), {"group": "Example"}, reverse_complement(self.sequence), "test", 2)
        self.assertEqual(orientation_key(first), orientation_key(reverse))
        opposite = event(dict(self.row, sstart=20, send=1), {"group": "Example"}, self.sequence, "test", 3)
        self.assertNotEqual(orientation_key(first), orientation_key(opposite))

    def test_inclusive_terminal_geometry_and_unknown_exclusion(self):
        first = event(self.row, {"group": "Example"}, self.sequence, "test", 1)
        observed = position_summary([first, first], "Example", 0, 1)
        self.assertEqual(observed["events"], 1)
        self.assertEqual(observed["terminal_events"], 1)
        self.assertAlmostEqual(observed["expected_n2"], 2 / 81)
        self.assertEqual(position_summary([first], "Example", 0, 1, primary=True)["events"], 0)
        first.update(topology="linear", record_type="genomic")
        self.assertEqual(position_summary([first], "Example", 0, 1, primary=True)["events"], 1)


if __name__ == "__main__":
    unittest.main()