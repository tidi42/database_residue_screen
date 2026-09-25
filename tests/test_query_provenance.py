import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.query_provenance import component_map


class QuerySegmentTests(unittest.TestCase):
    def test_unique_barcode_and_recorded_scaffold_partition(self):
        query = dict(query_id="NB01_forward_adapter_plus_barcode", condition="ont_adapter_barcode", construct_class="composite")
        components = [dict(component_id="NB01_forward", component_kind="barcode", start=5, end=12, strand="+"),
                      dict(component_id="adapter", component_kind="adapter", start=1, end=4, strand="+")]
        segments, candidates = component_map(query, "GGGGACGTAGCA", components, {"NB01_forward": "ACGTAGCA", "adapter": "GGGG"})
        self.assertEqual([row["role"] for row in segments], ["scaffold", "barcode"])
        self.assertEqual(sum(len(row["component_sequence"]) for row in segments), 12)
        self.assertEqual(len(candidates), 1)
        self.assertTrue(all(not row["physical_architecture_authenticated"] for row in segments))

    def test_ambiguous_repeated_barcode_stays_unknown(self):
        query = dict(query_id="NB01_forward_adapter_plus_barcode", condition="ont_adapter_barcode", construct_class="composite")
        segments, candidates = component_map(query, "GGACGTAGCAACGTAGCATT", [], {"NB01_forward": "ACGTAGCA"})
        self.assertEqual(len(candidates), 2)
        self.assertEqual(segments[0]["role"], "unknown")
        self.assertEqual(segments[0]["assignment_status"], "barcode_boundary_ambiguous")

    def test_wrong_component_sequence_is_not_silently_repaired(self):
        query = dict(query_id="NB01_forward_adapter_plus_barcode", condition="ont_adapter_barcode", construct_class="composite")
        components = [dict(component_id="NB01_forward", component_kind="barcode", start=5, end=12, strand="+"),
                      dict(component_id="adapter", component_kind="adapter", start=1, end=4, strand="+")]
        with self.assertRaises(ContractError):
            component_map(query, "GGGGACGTAGCA", components, {"NB01_forward": "ACGTAGCA", "adapter": "CCCC"})

    def test_standalone_primer_is_not_relabelled_barcode(self):
        query = dict(query_id="primer", condition="primer", construct_class="primer")
        segments, candidates = component_map(query, "ACGTAGCA", [], {})
        self.assertEqual(segments[0]["role"], "scaffold")
        self.assertEqual(candidates, [])


if __name__ == "__main__":
    unittest.main()