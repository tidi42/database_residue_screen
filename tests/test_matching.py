import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.matching import BLAST_FIELDS, parse_blast_hsp, strict_hsp, tolerant_hsp


class ProductionMatchingTests(unittest.TestCase):
    def setUp(self):
        self.queries = {"query": "ACGTTGCAACGTTGCAACGT"}
        self.row = dict(qseqid="query", sseqid="ref|TEST0001.1|", staxids="0", sscinames="N/A",
                        pident=100, length=20, qlen=20, qcovs=100, mismatch=0, gapopen=0,
                        qstart=1, qend=20, sstart=101, send=120, slen=500, evalue=0.01,
                        bitscore=40, stitle="test", qseq=self.queries["query"],
                        sseq=self.queries["query"], btop="20", sstrand="plus", qcovhsp=100)

    def parse(self, **changes):
        row = dict(self.row, **changes)
        return parse_blast_hsp([str(row[field]) for field in BLAST_FIELDS], self.queries)

    def test_exact_forward_and_reverse_subject(self):
        self.assertTrue(strict_hsp(self.parse()))
        self.assertTrue(strict_hsp(self.parse(sstart=120, send=101, sstrand="minus")))

    def test_subject_coverage_cannot_merge_partial_hsps(self):
        for start, end, sequence in ((1, 10, self.queries["query"][:10]),
                                     (11, 20, self.queries["query"][10:])):
            row = self.parse(length=10, qstart=start, qend=end, send=110,
                             qseq=sequence, sseq=sequence, btop="10", qcovhsp=50)
            self.assertFalse(strict_hsp(row))
            self.assertFalse(tolerant_hsp(row))

    def test_full_query_single_substitution_is_secondary_only(self):
        subject = "T" + self.queries["query"][1:]
        row = self.parse(sseq=subject, mismatch=1, pident=95, btop="AT19")
        self.assertFalse(strict_hsp(row))
        self.assertTrue(tolerant_hsp(row))

    def test_traceback_and_manifest_contradictions_are_rejected(self):
        for changes in ({"qseqid": "unknown"}, {"qlen": 21}, {"qseq": "T" * 20},
                        {"btop": "19"}, {"btop": "10000000"}, {"mismatch": 1},
                        {"sstrand": "minus"}, {"qcovhsp": "nan"}, {"qcovhsp": 99},
                        {"send": 121}, {"length": "bad"}, {"sseq": "?" * 20}):
            with self.subTest(changes=changes), self.assertRaises(ContractError):
                self.parse(**changes)

    def test_indel_alignment_is_parsed_but_not_accepted(self):
        query = self.queries["query"]
        row = self.parse(qseq=query[:10] + "-" + query[10:], sseq=query[:10] + "A" + query[10:],
                         btop="10-A10", pident=95.238, length=21, gapopen=1, send=121)
        self.assertFalse(strict_hsp(row))
        self.assertFalse(tolerant_hsp(row))


if __name__ == "__main__":
    unittest.main()