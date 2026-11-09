import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError
from database_residue_screen.database import metadata_conflicts, parse_database_info
from database_residue_screen.sequences import read_fasta


class DatabaseTests(unittest.TestCase):
    def test_parse_runtime_metadata(self):
        text = "Database: example\n 12 sequences; 1,234 total bases\nDate: May 29, 2026 1:21 AM\nBLASTDB Version: 5\nVolumes:\n /data/example.00\n"
        metadata = parse_database_info(text)
        self.assertEqual(metadata["base_count"], 1234)
        self.assertEqual(metadata["build_date"], "2026-05-29")
        self.assertEqual(metadata["volume_count"], 1)

    def test_conflicting_source_metadata_not_accepted(self):
        runtime = dict(sequence_count=12, base_count=1234, build_date="2026-05-29", database_version=5, volume_count=1)
        metadata = {"number-of-sequences": 13, "number-of-letters": 1234, "last-updated": "2025-06-04T23:31:00", "db-version": 5, "number-of-volumes": 1}
        self.assertEqual({item["field"] for item in metadata_conflicts(runtime, metadata)}, {"sequence_count", "build_date"})

    def test_missing_runtime_metadata_fails(self):
        with self.assertRaises(ContractError):
            parse_database_info("Database: inaccessible")

    def test_ambiguous_bases_allowed_only_for_subjects(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "subject.fasta"
            path.write_text(">example\nNNNACGTR\n", encoding="utf-8")
            self.assertEqual(read_fasta(path, allow_ambiguity=True)["example"], "NNNACGTR")
            with self.assertRaises(ContractError):
                read_fasta(path)


if __name__ == "__main__":
    unittest.main()