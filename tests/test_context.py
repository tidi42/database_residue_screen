import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from database_residue_screen.contracts import ContractError
from database_residue_screen.context import segment_coverage, traceback_segment_coverage, validate_flank_record


class ComponentCoverageTests(unittest.TestCase):
    def test_flank_identity_and_range_must_match_not_just_sequence(self):
        record = SeqRecord(Seq("ACGT"), id="TEST0001.1:10-13")
        self.assertEqual(validate_flank_record(record, "TEST0001.1", 10, 13), "ACGT")
        for identifier in ("OTHER0001.1:10-13", "TEST0001.1:11-14", "TEST0001.2:10-13", "TEST0001.1"):
            with self.subTest(identifier=identifier), self.assertRaises(ContractError):
                validate_flank_record(SeqRecord(Seq("ACGT"), id=identifier), "TEST0001.1", 10, 13)

    def test_tolerant_partial_barcode_is_not_declared_complete(self):
        event = dict(query_id="composite", qstart=1, qend=90, strict="no", source_row=1)
        component = dict(component_id="tag", component_kind="barcode", strand="-", start=81, end=100)
        observed = segment_coverage(event, component)
        self.assertEqual(observed["covered_bases"], 10)
        self.assertEqual(observed["complete_segment"], 0)
        self.assertEqual(observed["mismatch_allocation"], "unavailable_in_legacy_traceback")
        self.assertEqual(segment_coverage(dict(event, qstart=90, qend=1), component), observed)

    def test_strict_complete_segment_is_exact(self):
        event = dict(query_id="composite", qstart=100, qend=1, strict="yes", source_row=1)
        component = dict(component_id="tag", component_kind="barcode", strand="+", start=81, end=100)
        observed = segment_coverage(event, component)
        self.assertEqual(observed["complete_segment"], 1)
        self.assertEqual(observed["mismatch_allocation"], "all_bases_exact")

    def test_traceback_distinguishes_barcode_and_scaffold_mismatches(self):
        sequence = "ACGTACGTACGTACGTACGT"
        barcode = dict(start=1, end=10)
        hsp = dict(qstart=1, qend=20, qlen=20, qseq=sequence, sseq=sequence)
        barcode_error = dict(hsp, sseq=sequence[:2] + "A" + sequence[3:])
        scaffold_error = dict(hsp, sseq=sequence[:12] + "T" + sequence[13:])
        first = traceback_segment_coverage(barcode_error, barcode)
        second = traceback_segment_coverage(scaffold_error, barcode)
        self.assertEqual(first["covered_bases"], second["covered_bases"])
        self.assertEqual(first["mismatches"], 1)
        self.assertEqual(second["mismatches"], 0)
        self.assertEqual(first["matching_bases"], 9)
        self.assertEqual(second["matching_bases"], 10)

    def test_reverse_alias_and_clipped_barcode_use_physical_coordinates(self):
        hsp = dict(qstart=1, qend=18, qlen=20, qseq="ACGTACGTACGTACGTAC", sseq="ACGTACGTACGTACGTAC")
        result = traceback_segment_coverage(hsp, dict(start=1, end=10), alias_orientation="-")
        self.assertEqual((result["covered_start"], result["covered_end"]), (3, 10))
        self.assertEqual(result["covered_bases"], 8)
        self.assertEqual(result["support_state"], "partial_exact")
        empty = traceback_segment_coverage(hsp, dict(start=1, end=2), alias_orientation="-")
        self.assertEqual(empty["support_state"], "not_covered")
        reversed_hsp = dict(hsp, qstart=18, qend=1)
        self.assertEqual(traceback_segment_coverage(reversed_hsp, dict(start=1, end=10), "-"), result)

    def test_traceback_rejects_invalid_segments_or_gapped_hsps(self):
        hsp = dict(qstart=1, qend=4, qlen=4, qseq="ACGT", sseq="ACGT")
        with self.assertRaises(ContractError):
            traceback_segment_coverage(hsp, dict(start=0, end=4))
        with self.assertRaises(ContractError):
            traceback_segment_coverage(dict(hsp, sseq="AC-T"), dict(start=1, end=4))


if __name__ == "__main__":
    unittest.main()