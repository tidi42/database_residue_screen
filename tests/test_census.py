import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.census import classify, load_elements, orient_and_classify, retrieve, validate_specification
from database_residue_screen.contracts import ContractError
from database_residue_screen.example import save_fasta
from database_residue_screen.sequences import reverse_complement

# Synthetic flank elements and barcodes; no vendor sequence.
FRONT = "GATCCTAGCATTGC"
REAR = "TGACCGTA"
BARCODE = "AGCTTACGGATCCAGTTGACGTCA"
OTHER = "CCGTAGTTCAGGACT"
ELEMENTS = {"front": FRONT, "rear": REAR, "other": OTHER}
CATALOGUES = {"barcodes": {BARCODE}}
FAMILY = {"name": "tag_family", "conditions": ["bare"], "orientation": "query",
          "required_flank": {"upstream": 14, "downstream": 8},
          "classes": [
              {"name": "front_and_rear", "require": [{"side": "upstream", "element": "front"}, {"side": "downstream", "element": "rear"}]},
              {"name": "front_offset", "require": [{"side": "upstream", "element": "front", "offset": 2}]},
              {"name": "rear_only", "require": [{"side": "downstream", "element": "rear"}]},
              {"name": "catalogued_barcode", "require": [{"side": "either", "catalogue": "barcodes", "length": 24}]},
              {"name": "taxon_record", "require": [{"taxon": 10}]},
          ]}


class CensusSpecificationTests(unittest.TestCase):
    def test_specification_rules(self):
        validate_specification({"flank": 50, "families": [FAMILY]})
        cases = [{"flank": 10, "families": [FAMILY]},
                 {"flank": 50, "families": [FAMILY, FAMILY]},
                 {"flank": 50, "families": [dict(FAMILY, conditions=[])]},
                 {"flank": 50, "families": [dict(FAMILY, classes=[dict(FAMILY["classes"][0], name="none")])]},
                 {"flank": 50, "families": [dict(FAMILY, classes=[{"name": "x", "require": [{"side": "upstream", "sequence": FRONT}]}])]}]
        for case in cases:
            with self.subTest(case=case), self.assertRaises(ContractError):
                validate_specification(case)

    def test_elements_must_be_unique_acgt(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "elements.json"
            path.write_text(json.dumps({"elements": [{"name": "front", "sequence": FRONT}]}), encoding="utf-8")
            self.assertEqual(load_elements(path), {"front": FRONT})
            for elements in ([{"name": "a", "sequence": "ACGN"}], [{"name": "a", "sequence": "ACG"}, {"name": "a", "sequence": "TTT"}]):
                path.write_text(json.dumps({"elements": elements}), encoding="utf-8")
                with self.subTest(elements=elements), self.assertRaises(ContractError):
                    load_elements(path)


class CensusClassificationTests(unittest.TestCase):
    def classify(self, up, down, lineage=lambda taxon: False, mismatches=0):
        return classify(FAMILY, up, down, ELEMENTS, CATALOGUES, lineage, mismatches)

    def test_ordered_rules_and_precedence(self):
        self.assertEqual(self.classify("AA" + FRONT, REAR + "A" * 10), "front_and_rear")
        self.assertEqual(self.classify(FRONT + "CC", "A" * 10), "front_offset")
        self.assertEqual(self.classify("A" * 20, REAR + "C" * 10), "rear_only")
        self.assertEqual(self.classify("A" * 20, BARCODE + "C" * 4), "catalogued_barcode")
        self.assertEqual(self.classify("C" * 4 + reverse_complement(BARCODE), "A" * 10), "catalogued_barcode")
        self.assertEqual(self.classify("A" * 20, "C" * 20, lambda taxon: taxon == 10), "taxon_record")
        self.assertEqual(self.classify("A" * 20, "C" * 20), "none")
        self.assertEqual(self.classify("A" * 5, "C" * 20), "insufficient_flank")

    def test_one_mismatch_secondary_classification(self):
        one_off = FRONT[:-1] + ("A" if FRONT[-1] != "A" else "C")
        self.assertEqual(self.classify(one_off, REAR, mismatches=0), "rear_only")
        self.assertEqual(self.classify(one_off, REAR, mismatches=1), "front_and_rear")


@unittest.skipUnless(shutil.which("blastdbcmd") and shutil.which("makeblastdb"), "BLAST unavailable")
class CensusRetrievalTests(unittest.TestCase):
    def test_batch_retrieval_orientation_integrity_and_unresolved(self):
        context = "TTTT" + FRONT + BARCODE + REAR + "GGGG"
        records = {"CTX0001.1": "A" * 60 + context + "C" * 60, "CTX0002.1": "T" * 70 + reverse_complement(context) + "G" * 3}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "db.fa", records)
            subprocess.run(["makeblastdb", "-in", str(root / "db.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            first = 60 + 4 + len(FRONT) + 1
            second_upper = 70 + len(context) - 4 - len(FRONT)
            base = dict(family="tag_family", physical_id="TAG01", query=BARCODE, taxids="", testability="testable")
            requests = [dict(base, accession_version="CTX0001.1", lower=first, upper=first + 23, record_length=len(records["CTX0001.1"]), strand="+"),
                        dict(base, accession_version="CTX0002.1", lower=second_upper - 23, upper=second_upper,
                             record_length=len(records["CTX0002.1"]), strand="-"),
                        dict(base, accession_version="CTX0001.1", lower=first + 1, upper=first + 24, record_length=len(records["CTX0001.1"]), strand="+"),
                        dict(base, accession_version="unresolved:pdb|1MM8|B", lower=1, upper=24, record_length=24, strand="+")]
            retrieve(requests, root / "db", shutil.which("blastdbcmd"), root)
            rows = [orient_and_classify(item, FAMILY, ELEMENTS, CATALOGUES, None) for item in requests]
            self.assertEqual([row["status"] for row in rows], ["verified", "verified", "integrity_failed", "unresolved_identifier"])
            self.assertEqual([row["census_class"] for row in rows[:2]], ["front_and_rear", "front_and_rear"])
            self.assertEqual(len(rows[0]["upstream"]), 50)
            self.assertEqual(len(rows[1]["upstream"]), len(FRONT) + 4 + 3)
            self.assertEqual(len(rows[1]["downstream"]), 50)
            reverse_family = dict(FAMILY, orientation="reverse_complement")
            reverse = orient_and_classify(dict(requests[0], status="retrieved", window=requests[0]["window"]), reverse_family, ELEMENTS,
                                          CATALOGUES, None)
            self.assertEqual(reverse["upstream"], reverse_complement(rows[0]["downstream"]))


if __name__ == "__main__":
    unittest.main()
