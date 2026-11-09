import unittest

from database_residue_screen.eligibility import annotate_record, record_eligibility


class RecordEligibilityTests(unittest.TestCase):
    def setUp(self):
        self.event = dict(accession_version="TEST001.1", slen=20000, taxids="2", record_type="genomic", topology="linear")
        self.subject = dict(accession_version="TEST001.1", length=20000, taxids="3", record_type="genomic", topology="linear",
                            molecule="DNA", metadata_source="fixture", metadata_date="2026-09-15")

    def test_unknown_type_is_broad_match_not_a_genomic_negative(self):
        row = annotate_record(self.event, None)
        self.assertEqual(record_eligibility(row, "all_records")[0], "included")
        self.assertEqual(record_eligibility(row, "linear_genomic"), ("unresolved", "molecule_not_established"))

    def test_explicit_dna_is_required_and_taxid_evidence_is_not_overwritten(self):
        row = annotate_record(self.event, self.subject)
        self.assertEqual(record_eligibility(row, "linear_genomic")[0], "included")
        self.assertEqual(row["blast_taxids"], "2")
        self.assertEqual(row["record_taxids"], "3")
        row = annotate_record(self.event, {name: value for name, value in self.subject.items() if name != "molecule"})
        self.assertEqual(record_eligibility(row, "linear_genomic")[0], "unresolved")

    def test_known_exclusions_are_separate_from_unknowns(self):
        for change in (dict(topology="circular"), dict(molecule="RNA"), dict(record_type="transcript"),
                       dict(origin="synthetic"), dict(origin="vector"), dict(targeted="yes")):
            with self.subTest(change=change):
                row = annotate_record(self.event, dict(self.subject, **change))
                self.assertEqual(record_eligibility(row, "all_records")[0], "included")
                self.assertEqual(record_eligibility(row, "linear_genomic")[0], "excluded")
        for change in (dict(topology="unknown"), dict(molecule="conflicting"), dict(record_type="unknown")):
            row = annotate_record(self.event, dict(self.subject, **change))
            self.assertEqual(record_eligibility(row, "linear_genomic")[0], "unresolved")

    def test_identity_conflicts_never_become_resolved_zeros(self):
        for change in (dict(length=20001), dict(accession_version="TEST001.2")):
            row = annotate_record(self.event, dict(self.subject, **change))
            self.assertEqual(record_eligibility(row, "all_records")[0], "unresolved")
            self.assertEqual(record_eligibility(row, "linear_genomic")[0], "unresolved")
        row = annotate_record(dict(self.event, accession_version="unresolved:gi|123"), None)
        self.assertEqual(record_eligibility(row, "all_records")[0], "unresolved")


if __name__ == "__main__":
    unittest.main()