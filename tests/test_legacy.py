import unittest

from database_residue_screen.legacy import query_classification, recover_components
from database_residue_screen.matching import strict_hsp, tolerant_hsp


class LegacyAdapterTests(unittest.TestCase):
    def test_c2_is_separate_from_smrtbell_condition(self):
        self.assertEqual(query_classification("C2_primer_forward", "pacbio_smrtbell_adapter"), ("primer", "pacbio_c2_primer"))
        self.assertEqual(query_classification("SMRTbell_blunt_adapter_forward", "pacbio_smrtbell_adapter"), ("adapter", "pacbio_smrtbell_adapter"))

    def test_missing_standalone_flank_is_derived_but_not_authenticated(self):
        sequences = {"NB01_forward": "ACGTAC", "NB01_forward_adapter_plus_barcode": "TTGGAACGTAC"}
        metadata = [dict(query_id="NB01_forward", group="ONT", condition="ont_barcode", construct_class="barcode",
                         family_id="bare-family", role="exploratory", source_id="legacy", provenance_status="unverified", architecture_status="not_applicable"),
                    dict(query_id="NB01_forward_adapter_plus_barcode", group="ONT", condition="ont_adapter_barcode", construct_class="composite",
                         family_id="composite-family", role="exploratory", source_id="legacy", provenance_status="unverified", architecture_status="unverified")]
        rows, added, unresolved = recover_components(sequences, metadata)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(added), 1)
        self.assertEqual(unresolved, [])
        self.assertEqual(sequences[added[0]], "TTGGA")
        self.assertEqual(metadata[-1]["role"], "reference")
        self.assertEqual(metadata[-1]["provenance_status"], "unverified")
        self.assertEqual(metadata[1]["family_id"], "bare-family")

    def test_subject_coverage_cannot_rescue_partial_hsp(self):
        row = dict(pident=100.0, mismatch=0, gapopen=0, qstart=1, qend=12, qlen=24,
                   length=12, sstart=1, send=12, slen=50, qcovs=100)
        self.assertFalse(strict_hsp(row))
        self.assertFalse(tolerant_hsp(row))
        row.update(qend=24, length=24, send=24)
        self.assertTrue(strict_hsp(row))
        row.update(pident=95.833, mismatch=1)
        self.assertFalse(strict_hsp(row))
        self.assertTrue(tolerant_hsp(row))
        row.update(gapopen=1)
        self.assertFalse(tolerant_hsp(row))


if __name__ == "__main__":
    unittest.main()