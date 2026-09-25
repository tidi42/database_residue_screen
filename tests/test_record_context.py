import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from database_residue_screen.contracts import ContractError
from database_residue_screen.record_context import category_flags, choose_cases, selection_specification, verified_flanks


class RecordContextTests(unittest.TestCase):
    def test_taxon_label_alone_cannot_select_synthetic_case(self):
        event = dict(accession_version="TEST001.1", identity_status="exact_version", slen=20000, distance_to_end=10)
        metadata = dict(accessionversion="TEST001.1", slen=20000, title="A deposited chromosome", taxid=32630)
        flags = category_flags(event, metadata, ("excluded", "synthetic"))
        self.assertFalse(flags["synthetic_or_vector"])
        self.assertFalse(flags["terminal"])
        self.assertTrue(category_flags(event, dict(metadata, title="Cloning vector pExample, complete sequence"), ("excluded", "synthetic"))["synthetic_or_vector"])

    def test_terminal_requires_linear_genomic_eligibility_and_exact_metadata(self):
        event = dict(accession_version="TEST001.1", identity_status="exact_version", slen=20000, distance_to_end=10)
        metadata = dict(accessionversion="TEST001.1", slen=20000, title="chromosome")
        self.assertTrue(category_flags(event, metadata, ("included", "linear_genomic"))["terminal"])
        self.assertFalse(category_flags(event, metadata, ("excluded", "circular"))["terminal"])
        self.assertFalse(category_flags(event, None, ("included", "linear_genomic"))["terminal"])
        self.assertFalse(category_flags(event, metadata, ("included", "linear_genomic"), completed=False)["terminal"])

    def test_minus_strand_and_clipped_flanks_are_verified(self):
        record = SeqRecord(Seq("NNAACCGGTTNN"), id="TEST001.1:1-12")
        hsp = dict(sstart=6, send=3, slen=12, sseq="GGTT")
        result = verified_flanks(record, "TEST001.1", hsp)
        self.assertEqual(result["hsp_subject_sequence"], "GGTT")
        self.assertEqual(result["display_left_flank"], "NNAACC")
        self.assertEqual(result["display_right_flank"], "NN")
        self.assertTrue(result["left_flank_clipped"])
        with self.assertRaises(ContractError):
            verified_flanks(SeqRecord(record.seq, id="OTHER001.1:1-12"), "TEST001.1", hsp)
        with self.assertRaises(ContractError):
            verified_flanks(record, "TEST001.1", dict(hsp, sseq="AAAA"))

    def test_selection_is_deterministic_and_failures_are_bounded(self):
        candidates = [dict(candidate_id=str(index), accession_version=f"TEST{index:03d}.1", parent_id="parent", lower=1, upper=20,
                           canonical_strand="+", terminal=True, synthetic_or_vector=False, internal=False, ambiguous=False) for index in range(6)]
        calls = []
        def failed(row, category, attempt):
            calls.append(row["candidate_id"])
            return dict(status="retrieval_or_verification_failed", error="missing")
        selected, attempts, categories = choose_cases(list(reversed(candidates)), failed, selection_specification())
        self.assertEqual(calls, ["0", "1", "2"])
        self.assertEqual(selected, [])
        self.assertEqual(sum(row["attempted"] for row in attempts), 3)
        self.assertEqual(categories[0]["status"], "unavailable_no_verified_candidate")


if __name__ == "__main__":
    unittest.main()