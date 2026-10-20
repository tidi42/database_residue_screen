import json
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError, write_table
from database_residue_screen.design import digest
from database_residue_screen.insdc import fetch_insdc_annotation, parse_genbank_headers, read_annotation, targeted_status
from database_residue_screen.policy import default_policy
from database_residue_screen.positional import one_event_per, positional_sensitivity, projects_from_metadata
from database_residue_screen.reanalysis import ANNOTATED_EVENT_SCHEMA

HEADERS = """LOCUS       AB000001                   1 bp    DNA     linear   BCT 01-FEB-2014
DEFINITION  Example 16S amplicon.
ACCESSION   AB000001 REGION: 1
VERSION     AB000001.2
KEYWORDS    TLS; Targeted Locus Study;
            AMPLICON.
FEATURES             Location/Qualifiers
     source          1
                     /organism="Example"
                     /PCR_primers="fwd_seq: acgt, rev_seq: tgca"
//
LOCUS       AB000002                   1 bp    DNA     linear   BCT 01-FEB-2014
DEFINITION  Example genome.
ACCESSION   AB000002 REGION: 1
VERSION     AB000002.1
KEYWORDS    .
FEATURES             Location/Qualifiers
     source          1
                     /organism="Example"
     gene            1
                     /PCR_primers="not a source qualifier"
//
"""


class InsdcTests(unittest.TestCase):
    def test_structured_fields_only(self):
        parsed = parse_genbank_headers(HEADERS)
        self.assertEqual(set(parsed), {"AB000001.2", "AB000002.1"})
        self.assertTrue(parsed["AB000001.2"]["pcr_primers"])
        self.assertFalse(parsed["AB000002.1"]["pcr_primers"])
        status, evidence = targeted_status(parsed["AB000001.2"], "")
        self.assertEqual(status, "yes")
        self.assertIn("KEYWORDS:AMPLICON", evidence)
        self.assertIn("KEYWORDS:TLS", evidence)
        self.assertEqual(targeted_status(parsed["AB000002.1"], "")[0], "no")
        self.assertEqual(targeted_status(parsed["AB000002.1"], "barcode"), ("yes", "esummary_tech:barcode"))
        self.assertEqual(targeted_status(None, "")[0], "unknown")

    def test_fetch_archives_versions_and_failures(self):
        calls = []

        def fake(chunk):
            calls.append(chunk)
            if "ZZ000009.1" in chunk:
                raise ContractError("down")
            return HEADERS, []

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = fetch_insdc_annotation(["AB000001.1", "AB000002.1", "unresolved:x", "ZZ000009.1"], {}, root / "insdc",
                                            batch_size=2, fetch=fake, delay=0)
            rows = read_annotation(root / "insdc")
            self.assertEqual(rows["AB000001.1"]["status"], "version_not_returned")
            self.assertEqual(rows["AB000002.1"]["targeted_insdc"], "no")
            self.assertEqual(rows["ZZ000009.1"]["status"], "retrieval_failed")
            self.assertEqual(report["requested"], 3)
            self.assertEqual(len(report["problems"]), 1)
            (root / "insdc/insdc_annotation.tsv").write_text("changed")
            with self.assertRaises(ContractError):
                read_annotation(root / "insdc")


def event(accession, start, slen, group="ONT", parent="p1", condition="keep", strand="+"):
    row = {name: "" for name in ANNOTATED_EVENT_SCHEMA["properties"]}
    row.update({name: 0 for name, field in ANNOTATED_EVENT_SCHEMA["properties"].items() if field["type"] in ("integer", "number")})
    row.update(query_id=parent, parent_id=parent, group=group, condition=condition, accession_version=accession, canonical_strand=strand,
               strict="yes", tolerant="yes", record_type="genomic", topology="linear", molecule="DNA", identity_status="exact_version",
               qlen=40, length=40, sstart=start, send=start + 39, slen=slen, distance_to_end=min(start - 1, slen - start - 39),
               origin="unknown", targeted="unknown", metadata_status="exact_version_record")
    return row


class PositionalTests(unittest.TestCase):
    def test_one_event_choice_is_position_independent(self):
        rows = [event("A.1", 1, 20000), event("A.1", 9000, 20000), event("B.1", 50, 20000)]
        first = one_event_per(rows, lambda row: row["accession_version"])
        second = one_event_per(list(reversed(rows)), lambda row: row["accession_version"])
        self.assertEqual(len(first), 2)
        self.assertEqual(sorted(map(str, first)), sorted(map(str, second)))

    def test_projects_and_variants(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "meta").mkdir()
            (root / "meta/esummary_0001.json").write_text(json.dumps({"result": {
                "uids": ["1", "2"], "1": {"accessionversion": "A.1", "projectid": "PRJ1"}, "2": {"accessionversion": "B.1", "projectid": "0"}}}))
            projects = projects_from_metadata(root / "meta")
            self.assertEqual(projects, {"A.1": "PRJ1", "B.1": ""})
            analysis = root / "analysis"
            analysis.mkdir()
            rows = [event("A.1", 1, 20000), event("A.1", 9000, 20000), event("B.1", 10, 20000),
                    event("C.1", 15, 20000, condition="old"), event("D.1", 10, 20000, group="PacBio", parent="p2")]
            write_table(analysis / "real_events.tsv", rows, ANNOTATED_EVENT_SCHEMA)
            policy = dict(default_policy(), condition_roles={"keep": dict(panel_version="v1", role="primary"),
                                                             "old": dict(panel_version="v1", role="superseded", superseded_by=["keep"])})
            (analysis / "analysis.json").write_text(json.dumps(dict(analysis_policy=policy,
                                                                    files={"real_events.tsv": digest(analysis / "real_events.tsv")})))
            result = positional_sensitivity(analysis, projects, {"B.1": dict(targeted_insdc="yes")})
            pick = {(r["variant"], r["group"], r["window"], r["minimum_length"]): r for r in result}
            self.assertEqual(pick["all_events", "ONT", 500, 10000]["events"], 3)
            self.assertEqual(pick["all_events", "ONT", 500, 10000]["terminal_events"], 2)
            self.assertEqual(pick["one_per_accession", "ONT", 500, 10000]["events"], 2)
            self.assertEqual(pick["excluding_targeted_insdc", "ONT", 500, 10000]["events"], 2)
            self.assertEqual(pick["all_events", "Global union", 25, 10000]["events"], 4)
            self.assertEqual(pick["all_events", "ONT", 25, 10000]["units_without_project"], 2)


if __name__ == "__main__":
    unittest.main()
