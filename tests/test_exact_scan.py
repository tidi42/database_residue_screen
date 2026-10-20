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

    def test_interrupted_scan_resumes_only_with_identical_request(self):
        from database_residue_screen.database import inspect_database
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "a.fa", {"VOLA0001.1": "C" * 10 + OUTER + "G" * 10})
            save_fasta(root / "b.fa", {"VOLB0001.1": "T" * 5 + reverse_complement(OTHER), "VOLB0002.1": "A" * 40})
            for name in ("a", "b"):
                subprocess.run(["makeblastdb", "-in", str(root / f"{name}.fa"), "-out", str(root / name), "-dbtype", "nucl", "-parse_seqids"],
                               check=True, capture_output=True)
            subprocess.run(["blastdb_aliastool", "-dblist", f"{root / 'a'} {root / 'b'}", "-dbtype", "nucl", "-out", str(root / "db"),
                            "-title", "two volumes"], check=True, capture_output=True)
            inspect_database(root / "db", root / "evidence")
            evidence = root / "evidence/database_identity.json"
            volumes = read_json(evidence)["runtime_metadata"]["volume_names"]
            self.assertEqual(len(volumes), 2)
            save_fasta(root / "q.fa", {"outer": OUTER, "inner": INNER, "other": OTHER})
            blastdbcmd = shutil.which("blastdbcmd")
            first = run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1)
            self.assertEqual(first["reused_volumes"], [])
            expected = sorted(map(str, read_occurrences(root / "scan")))
            with self.assertRaises(ContractError):
                run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1)
            with self.assertRaises(ContractError):
                run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1, resume=True)
            # Simulate an interruption after the first volume: no final report, second volume unfinished, stray partial file.
            (root / "scan/scan.json").unlink()
            second = Path(volumes[1]).name
            (root / f"scan/{second}.complete.json").unlink()
            (root / f"scan/{second}.occurrences.tsv.partial").write_text("truncated")
            with self.assertRaises(ContractError):
                run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1, minimum_length=30, resume=True)
            resumed = run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1, resume=True)
            self.assertEqual(resumed["reused_volumes"], [Path(volumes[0]).name])
            self.assertTrue(resumed["totals_match_database"])
            self.assertEqual(sorted(map(str, read_occurrences(root / "scan"))), expected)
            self.assertEqual(list((root / "scan").glob("*.partial")), [])
            # A changed completed volume is never reused.
            (root / "scan/scan.json").unlink()
            reused = root / f"scan/{Path(volumes[0]).name}.occurrences.tsv"
            reused.write_text(reused.read_text() + "\n")
            with self.assertRaises(ContractError):
                run_exact_scan(root / "q.fa", root / "db", evidence, root / "scan", blastdbcmd, workers=1, resume=True)


if __name__ == "__main__":
    unittest.main()
