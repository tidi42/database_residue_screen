import shutil
import subprocess
import tempfile
import unittest
import warnings
from pathlib import Path

from database_residue_screen.example import save_fasta
from database_residue_screen.cap_review import bits, review_saturated

QUERY = "ACGTTGCAAGGCTTACCGATAGCTTGACCATG"
OTHER = "TTGACCGATCAGGCATTCAGGACTTAGCCATGCAGTTAGC"


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST unavailable")
class CapReviewTests(unittest.TestCase):
    def test_saturated_exact_query_is_flagged_and_unsaturated_is_skipped(self):
        from database_residue_screen.database import inspect_database
        from database_residue_screen.search import run_search, search_settings
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records = {f"CAP{i:04d}.1": "N" * 20 + QUERY + "N" * 20 for i in range(6)}
            records["ONE0001.1"] = "N" * 20 + OTHER + "N" * 20
            save_fasta(root / "db.fa", records)
            subprocess.run(["makeblastdb", "-in", str(root / "db.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            inspect_database(root / "db", root / "evidence")
            save_fasta(root / "q.fa", {"capped": QUERY, "single": OTHER})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                run_search(root / "q.fa", root / "db", root / "evidence/database_identity.json", root / "search",
                           search_settings(target_limit=3), threads=1, mode="working")
            rows = review_saturated(root / "search", {"capped": {f"CAP{i:04d}.1" for i in range(6)}})
            self.assertEqual([row["query_id"] for row in rows], ["capped"])
            row = rows[0]
            self.assertEqual((row["retained_subjects"], row["strict_subjects"], row["exact_scan_accessions"]), (3, 3, 6))
            self.assertEqual(row["strict_status"], "cap_censored_possible")
            self.assertGreater(row["weakest_retained_bitscore"], row["strict_threshold_bits"] - 0.5)

    def test_bit_score_conversion(self):
        row = dict(lambda_value=1.37406312246009, kappa=0.710602795216363)
        self.assertAlmostEqual(bits(32, row), 63.93, places=1)


if __name__ == "__main__":
    unittest.main()
