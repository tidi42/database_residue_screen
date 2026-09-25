import copy
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_json, read_table, schema, validate_record, write_table


class QueryContractTests(unittest.TestCase):
    def setUp(self):
        self.definition = schema("queries")
        self.record = dict(query_id="example", group="synthetic", condition="bare", construct_class="barcode",
                           family_id="family", role="primary", source_id="example-source",
                           provenance_status="synthetic", architecture_status="not_applicable")

    def test_primary_requires_verified_or_synthetic_provenance(self):
        record = dict(self.record, provenance_status="unverified")
        with self.assertRaises(ContractError):
            validate_record(record, self.definition)

    def test_exploratory_unknowns_remain_explicit(self):
        record = dict(self.record, role="exploratory", provenance_status="unverified", source_id="", family_id="")
        validate_record(record, self.definition)

    def test_primary_composite_requires_architecture_evidence(self):
        record = dict(self.record, construct_class="composite", architecture_status="unverified")
        with self.assertRaises(ContractError):
            validate_record(record, self.definition)

    def test_unknown_columns_are_rejected(self):
        with self.assertRaises(ContractError):
            validate_record(dict(self.record, inferred_genome="yes"), self.definition)

    def test_table_round_trip_and_refuse_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "queries.tsv"
            write_table(path, [self.record], self.definition)
            self.assertEqual(read_table(path, self.definition), [self.record])
            with self.assertRaises(FileExistsError):
                write_table(path, [self.record], self.definition)

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text('{"status": "failed", "status": "succeeded"}', encoding="utf-8")
            with self.assertRaises(ContractError):
                read_json(path)

    def test_nonfinite_table_number_rejected(self):
        definition = copy.deepcopy(self.definition)
        definition["properties"]["score"] = {"type": "number"}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.tsv"
            columns = list(definition["properties"])
            values = [str(dict(self.record, score="nan")[column]) for column in columns]
            path.write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n", encoding="utf-8")
            with self.assertRaises(ContractError):
                read_table(path, definition)


if __name__ == "__main__":
    unittest.main()