import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.census import classify, orient_and_classify, retrieve
from database_residue_screen.example import save_fasta
from database_residue_screen.sequences import reverse_complement

NB = "CACAAAGACACCGACAACTTTCTT"
REAR50 = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
CATALOGUES = {"ont_barcodes": {NB}, "pacbio_barcodes": {"CACATATCAGAGTGCG"}}


class CensusClassificationTests(unittest.TestCase):
    def test_native_barcode_rules_and_precedence(self):
        self.assertEqual(classify("ont_native_barcode_24", "GG" + "ATTGCTAAGGTTAA", "CAGCACCT" + "A" * 10, CATALOGUES, False), "K46")
        self.assertEqual(classify("ont_native_barcode_24", "CCCCCC" + "AAGGTTAA", "CAGCACCT" + "A" * 10, CATALOGUES, False), "K40_only")
        self.assertEqual(classify("ont_native_barcode_24", "GCTTGGGTGTTTAACC", "GTTTTCGCATTTATCG", CATALOGUES, False), "rapid_context")
        self.assertEqual(classify("ont_native_barcode_24", "A" * 20, "CAGCACCT" + "C" * 10, CATALOGUES, False), "native_rear_only")
        self.assertEqual(classify("ont_native_barcode_24", "A" * 20, "C" * 20, CATALOGUES, False), "none")
        self.assertEqual(classify("ont_native_barcode_24", "A" * 5, "C" * 20, CATALOGUES, False), "insufficient_flank")
        one_off = "ATTGCTAAGGTTAT"
        self.assertEqual(classify("ont_native_barcode_24", one_off, "CAGCACCT", CATALOGUES, False, 0), "native_rear_only")
        self.assertEqual(classify("ont_native_barcode_24", one_off, "CAGCACCT", CATALOGUES, False, 1), "K46")

    def test_rapid_rear_smrtbell_and_tn5_rules(self):
        self.assertEqual(classify("ont_rapid_rear_50", "GCTTGGGTGTTTAACC" + NB, "", CATALOGUES, False), "rapid_front_and_barcode")
        self.assertEqual(classify("ont_rapid_rear_50", "A" * 16 + reverse_complement(NB), "", CATALOGUES, False), "catalogued_barcode_only")
        self.assertEqual(classify("ont_rapid_rear_50", "A" * 40, "", CATALOGUES, True), "phage_mu_record")
        self.assertEqual(classify("pacbio_smrtbell_45", "A" * 16, reverse_complement("CACATATCAGAGTGCG"), CATALOGUES, False), "barcode_adjacent")
        self.assertEqual(classify("illumina_tn5_19", "TCGTCGGCAGCGTC", "", CATALOGUES, False), "nextera_A")
        self.assertEqual(classify("illumina_tn5_19", "GTCTCGTGGGCTCGG", "", CATALOGUES, False), "nextera_B")


@unittest.skipUnless(shutil.which("blastdbcmd") and shutil.which("makeblastdb"), "BLAST unavailable")
class CensusRetrievalTests(unittest.TestCase):
    def test_batch_retrieval_orientation_integrity_and_unresolved(self):
        k46 = "TTTT" + "ATTGCTAAGGTTAA" + NB + "CAGCACCT" + "GGGG"
        records = {"CTX0001.1": "A" * 60 + k46 + "C" * 60, "CTX0002.1": "T" * 70 + reverse_complement(k46) + "G" * 3}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "db.fa", records)
            subprocess.run(["makeblastdb", "-in", str(root / "db.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            first = 60 + 4 + 14 + 1
            second_upper = 70 + len(k46) - 4 - 14
            requests = [dict(family="ont_native_barcode_24", physical_id="NB01_forward", query=NB, accession_version="CTX0001.1",
                             lower=first, upper=first + 23, record_length=len(records["CTX0001.1"]), strand="+", taxids="", testability="testable"),
                        dict(family="ont_native_barcode_24", physical_id="NB01_forward", query=NB, accession_version="CTX0002.1",
                             lower=second_upper - 23, upper=second_upper, record_length=len(records["CTX0002.1"]), strand="-", taxids="",
                             testability="testable"),
                        dict(family="ont_native_barcode_24", physical_id="NB01_forward", query=NB, accession_version="CTX0001.1",
                             lower=first + 1, upper=first + 24, record_length=len(records["CTX0001.1"]), strand="+", taxids="", testability="testable"),
                        dict(family="ont_native_barcode_24", physical_id="NB01_forward", query=NB, accession_version="unresolved:pdb|1MM8|B",
                             lower=1, upper=24, record_length=24, strand="+", taxids="", testability="testable")]
            retrieve(requests, root / "db", shutil.which("blastdbcmd"), root)
            rows = [orient_and_classify(item, CATALOGUES, None) for item in requests]
            self.assertEqual([row["status"] for row in rows], ["verified", "verified", "integrity_failed", "unresolved_identifier"])
            self.assertEqual([row["census_class"] for row in rows[:2]], ["K46", "K46"])
            self.assertEqual(len(rows[0]["upstream"]), 50)
            self.assertEqual(len(rows[1]["upstream"]), 14 + 4 + 3)
            self.assertEqual(len(rows[1]["downstream"]), 50)


if __name__ == "__main__":
    unittest.main()
