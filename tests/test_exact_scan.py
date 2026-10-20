import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_json
from database_residue_screen.example import save_fasta
from database_residue_screen.exact_scan import read_occurrences, reverse_complement, run_exact_scan, scan_lines, scan_patterns

INNER = "AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT"          # 40 nt, nested in OUTER
OUTER = "ATTGCT" + INNER                                       # 46 nt
PALINDROME = "GAATTCGGATCCAAGCTTGCATGCAAGCTTGGATCCGAATTC"     # 42 nt, its own reverse complement
OTHER = "TTGACCGATCAGGCATTCAGGACTTAGCCATGCAG"               # 35 nt
SHORT = "ACGTACGTAC"                                           # below the minimum length


@unittest.skipUnless(shutil.which("blastdbcmd") and shutil.which("makeblastdb"), "BLAST unavailable")
class ExactScanTests(unittest.TestCase):
    def test_both_strands_nested_ends_palindrome_and_ambiguity(self):
        from database_residue_screen.database import inspect_database
        self.assertEqual(PALINDROME, reverse_complement(PALINDROME))
        records = {
            "REC0001.1": OUTER + "C" * 20 + reverse_complement(OTHER),       # nested pair at the record start, OTHER on '-' at the end
            "REC0002.1": "G" * 15 + PALINDROME + "T" * 10 + OTHER[:17] + "N" + OTHER[18:],  # palindrome once; N blocks OTHER
            "REC0003.1": "A" * 50,
        }
        queries = {"outer": OUTER, "inner": INNER, "pal": PALINDROME, "other": OTHER, "short": SHORT}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "db.fa", records)
            save_fasta(root / "q.fa", queries)
            subprocess.run(["makeblastdb", "-in", str(root / "db.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            inspect_database(root / "db", root / "evidence")
            report = run_exact_scan(root / "q.fa", root / "db", root / "evidence/database_identity.json", root / "scan",
                                    shutil.which("blastdbcmd"), workers=1, minimum_length=23)
            self.assertTrue(report["totals_match_database"])
            self.assertEqual(report["records"], 3)
            self.assertEqual(report["bases"], sum(map(len, records.values())))
            self.assertEqual(report["scanned_queries"], 4)
            found = {(row["query_id"], row["accession_version"], row["strand"], row["start"], row["end"]) for row in read_occurrences(root / "scan")}
            length = len(records["REC0001.1"])
            self.assertEqual(found, {("outer", "REC0001.1", "+", 1, 46), ("inner", "REC0001.1", "+", 7, 46),
                                     ("other", "REC0001.1", "-", length - len(OTHER) + 1, length),
                                     ("pal", "REC0002.1", "+", 16, 57)})
            # changing an output file invalidates the scan
            first = next(iter((root / "scan").glob("*.occurrences.tsv")))
            first.write_text(first.read_text() + "\n")
            with self.assertRaises(ContractError):
                read_occurrences(root / "scan")
            # totals that disagree with the evidence are not accepted
            evidence = read_json(root / "evidence/database_identity.json")
            evidence["runtime_metadata"]["base_count"] += 1
            (root / "bad.json").write_text(__import__("json").dumps(evidence))
            with self.assertRaises(ContractError):
                run_exact_scan(root / "q.fa", root / "db", root / "bad.json", root / "scan_bad", shutil.which("blastdbcmd"), workers=1)

    def test_scan_lines_reports_adjacent_and_overlapping_occurrences(self):
        from database_residue_screen.exact_scan import automaton
        machine = automaton(scan_patterns({"other": OTHER}, 23))
        sequence = OTHER + OTHER + reverse_complement(OTHER)
        rows = [row for row in scan_lines([f"X1.1 {sequence.lower()}\n"], machine, "v") if isinstance(row, dict)]
        self.assertEqual([(row["start"], row["strand"]) for row in rows], [(1, "+"), (36, "+"), (71, "-")])

    def test_ambiguous_patterns_are_refused(self):
        with self.assertRaises(ContractError):
            scan_patterns({"bad": "ACGTN" * 6}, 23)


if __name__ == "__main__":
    unittest.main()
