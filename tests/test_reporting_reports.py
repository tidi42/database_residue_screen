"""Publication reports on a synthetic study: every report, its checks and its refusal of inconsistent inputs."""
import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import reporting_fixture
from database_residue_screen.contracts import ContractError
from database_residue_screen.reporting.components import parse_segments, readout
from database_residue_screen.reporting.metrics import compute_metrics
from database_residue_screen.reporting.reports import build_reports
from database_residue_screen.reporting.target_limit import normalise_limit_row


def table(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def report_inputs(study: dict, metrics: Path) -> dict:
    return dict(metrics=metrics, analysis=study["taxonomy"], design=study["design"], search=study["search"], metadata=study["metadata"],
                panel_definition=study["panel"] / "panel_definition.tsv", insdc=study["insdc"],
                positional_table=study["positional"] / "positional_sensitivity.tsv", cap_review=study["cap_review"] / "cap_review.tsv",
                target_limit=study["target_limit"] / "target_limit.tsv", exact_scan=study["exact_scan"], scan_validation=study["scan_validation"],
                benchmarks={path.name.removeprefix("benchmark_"): path for path in sorted(study["benchmark"].glob("benchmark_*"))},
                inference=study["inference"], specification=study["inference"] / "specification.json", record_cases=study["record_cases"],
                census=study["census"], validation=study["validation"] / "validation.json")


class ReportsCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.study = reporting_fixture.build(cls.root / "study")
        cls.metrics = cls.root / "metrics"
        compute_metrics(cls.study["taxonomy"], cls.metrics, design=cls.study["design"], search=cls.study["search"])
        cls.inputs = report_inputs(cls.study, cls.metrics)
        cls.reports = cls.root / "reports"
        cls.result = build_reports(cls.reports, **cls.inputs)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def summary(self, name: str) -> dict:
        return json.loads((self.reports / name / "summary.json").read_text())


class ReportTests(ReportsCase):
    def test_every_report_is_written(self):
        self.assertEqual(self.result["skipped"], [])
        self.assertEqual(self.result["reports"], ["benchmark_misses", "cap_review", "census", "components", "eligibility", "exact_scan_validation",
                                                  "inference", "nesting", "positional", "record_cases", "target_limit"])
        self.assertEqual(self.result["validation"], dict(valid=True, result_contract="reference_screen_v2"))
        manifest = json.loads((self.reports / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "publication_reports_v1")
        self.assertIn("components/segment_readout.tsv", manifest["outputs"])

    def test_components(self):
        self.assertEqual(parse_segments("front:1-8;barcode:9-32=NB01;rear:33-40")[1], dict(role="barcode", start=9, end=32, element="NB01"))
        self.assertEqual(readout("ACGT", "ACGT", {"ACGT"}), "queried_exact")
        self.assertEqual(readout("ACGA", "ACGT", {"ACGT", "ACGA"}), "other_catalogued_exact")
        self.assertEqual(readout("ACGA", "ACGT", {"ACGT"}), "queried_1_mismatch")
        self.assertEqual(readout("TTTT", "ACGT", {"ACGT"}), "uncatalogued")
        summary = self.summary("components")
        self.assertEqual((summary["physical_queries"], summary["composite_queries"]), (17, 5))
        self.assertEqual(summary["minimum_within_set_distance"]["illumina_i7|10"], 1)
        rows = table(self.reports / "components" / "segment_readout.tsv")
        swaps = [row for row in rows if row["readout"] == "other_catalogued_exact"]
        self.assertEqual(len(swaps), 1)
        self.assertEqual((swaps[0]["physical_id"], swaps[0]["endpoint_state"]), ("ILLUMINA_I7ADP_A_forward", "tolerant_only"))
        sequences = reporting_fixture.read_fasta(self.study["design"] / "physical_queries.fasta")
        self.assertEqual(swaps[0]["observed_segment"], sequences["ILLUMINA_I7_B_forward"])
        inventory = table(self.reports / "components" / "component_inventory.tsv")
        self.assertEqual(len(inventory), 24)

    def test_cap_review_and_target_limit(self):
        cap = self.summary("cap_review")
        self.assertEqual((cap["capped_searches"], cap["reviewed_searches"], cap["unreviewed_capped_searches"]), (2, 2, 0))
        roles = {row["search_role"] for row in table(self.reports / "cap_review" / "cap_review_searches.tsv")}
        self.assertEqual(roles, {"reference", "control"})
        limit = self.summary("target_limit")
        self.assertEqual((limit["limits"], limit["parents"], limit["rows"], limit["production_parents_without_arms"]), ([500, 2000], 16, 32, 0))
        self.assertEqual((limit["detection_changed"], limit["production_differs"]), (0, 0))
        self.assertGreater(limit["breadth_changed"], 0)
        rows = table(self.reports / "target_limit" / "target_limit_parents.tsv")
        changed = [row for row in rows if row["breadth_agrees"] == "no"]
        self.assertTrue(all(int(row["retained_500"]) == int(row["retained_2000"]) - 1 for row in changed))
        summary = {(row["group"], row["endpoint"], row["population"], row["gate"]): row
                   for row in table(self.reports / "target_limit" / "target_limit_summary.tsv")}
        self.assertEqual(summary["Global union", "strict", "all_records", "all_parents"]["parents"], "16")
        self.assertEqual(normalise_limit_row(dict(max_target_seqs="500", parent_id="p", endpoint="strict", population="all_records",
                                                  completion="incomplete_not_zero", retained_accessions="0", unresolved_subjects="0"))["detection"], None)
        self.assertEqual(normalise_limit_row(dict(target_limit="limit_2000", sequence_id="p", endpoint="strict", population="all_records",
                                                  status="success_hits", retained_accessions="0", unresolved_subjects="2"))["breadth"], None)

    def test_eligibility_positional_and_nesting(self):
        eligibility = self.summary("eligibility")
        self.assertEqual(eligibility["testable_strict"]["resolved_identifiers"], 7)
        self.assertEqual(eligibility["testable_strict"]["unresolved|unresolved_identifier"], 1)
        self.assertGreater(eligibility["ungated_strict"]["distinct_accessions"], eligibility["testable_strict"]["distinct_accessions"])
        flow = table(self.reports / "eligibility" / "eligibility_flow.tsv")
        self.assertEqual({row["scope"] for row in flow}, {"testable_strict", "ungated_strict"})
        positional = self.summary("positional")
        self.assertTrue(positional["stage_table_checked"] and positional["insdc_targeted_variant"])
        self.assertEqual(positional["positional_events"], 8)
        variants = {row["variant"] for row in table(self.reports / "positional" / "positional_sensitivity.tsv")}
        self.assertEqual(variants, {"all_events", "one_per_accession", "one_per_bioproject", "excluding_targeted_insdc"})
        nesting = self.summary("nesting")["all_records"]["ONT"]
        # The 24-nt barcode lies inside the K40 context: two per-parent events on one locus.
        self.assertEqual((nesting["per_parent_events"], nesting["distinct_loci"]), (5, 4))

    def test_benchmark_misses_scan_inference_record_cases_and_census(self):
        misses = self.summary("benchmark_misses")
        self.assertEqual(sorted(misses), ["core_nt", "core_nt_div10", "core_nt_x10"])
        self.assertEqual(misses["core_nt"]["by_reason"], {"hsp_returned_not_qualifying": 1, "no_hsp_no_seed": 1})
        self.assertEqual(misses["core_nt"]["by_hsp_detail"]["hsp_span_below_90pct"], 1)
        scan = self.summary("exact_scan_validation")
        self.assertTrue(scan["accepted"] and scan["all_scan_only_explained"] and scan["stage_summary_checked"])
        self.assertEqual((scan["scan_only_target_cap"], scan["blast_only"]), (2, 0))
        inference = self.summary("inference")
        self.assertEqual((inference["draw_prefixes"], inference["holm_family_size"]), ([5, 10, 15], 8))
        self.assertAlmostEqual(inference["stage_alpha"], 0.05 / 3)
        slots = table(self.reports / "inference" / "conditional_slots.tsv")
        self.assertEqual([(row["group"], row["measure"]) for row in slots][:2], [("ONT", "detection"), ("ONT", "breadth")])
        self.assertTrue(all("post hoc" not in row["interpretation"].lower() for row in slots))
        cases = self.summary("record_cases")
        self.assertEqual((cases["cases"], cases["packets_checked"], cases["unavailable_categories"]), (2, 1, ["internal_or_ambiguous"]))
        rows = table(self.reports / "record_cases" / "record_case_table.tsv")
        self.assertTrue(all(row["displayed_window_equals_record"] == "True" and row["record_window_equals_hsp"] == "True" for row in rows))
        self.assertEqual([row["strict_testability"] for row in rows], ["testable", "censored_untestable"])
        self.assertEqual(self.summary("census")["windows"], 2)

    def test_outputs_are_deterministic(self):
        again = self.root / "reports-again"
        build_reports(again, **self.inputs)
        first = {path.relative_to(self.reports): path.read_bytes() for path in self.reports.rglob("*") if path.is_file()}
        second = {path.relative_to(again): path.read_bytes() for path in again.rglob("*") if path.is_file()}
        self.assertEqual(first, second)


class RefusalTests(ReportsCase):
    def copy(self, name: str) -> Path:
        target = self.root / "copies" / name
        shutil.rmtree(target, ignore_errors=True)
        shutil.copytree(self.study[name], target)
        return target

    def test_metrics_of_another_analysis_are_refused(self):
        other = self.root / "metrics-other"
        shutil.copytree(self.metrics, other)
        manifest = json.loads((other / "manifest.json").read_text())
        manifest["analysis_sha256"] = "0" * 64
        (other / "manifest.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ContractError, "different analysis"):
            build_reports(self.root / "refused-metrics", **dict(self.inputs, metrics=other))

    def test_changed_positional_stage_table_is_refused(self):
        positional = self.copy("positional")
        path = positional / "positional_sensitivity.tsv"
        lines = path.read_text().splitlines()
        values = lines[1].split("\t")
        values[5] = str(int(values[5]) + 1)
        path.write_text("\n".join([lines[0], "\t".join(values)] + lines[2:]) + "\n")
        with self.assertRaisesRegex(ContractError, "positional"):
            build_reports(self.root / "refused-positional", **dict(self.inputs, positional_table=path))

    def test_inconsistent_scan_validation_is_refused(self):
        validation = self.copy("scan_validation")
        summary = json.loads((validation / "scan_summary.json").read_text())
        summary["matched"] += 1
        (validation / "scan_summary.json").write_text(json.dumps(summary))
        with self.assertRaisesRegex(ContractError, "disagrees"):
            build_reports(self.root / "refused-scan", **dict(self.inputs, scan_validation=validation))

    def test_wrong_draw_prefixes_and_record_case_mismatch_are_refused(self):
        with self.assertRaisesRegex(ContractError, "draw prefixes"):
            build_reports(self.root / "refused-prefixes", **dict(self.inputs, draw_prefixes=(5, 10)))
        cases = self.copy("record_cases")
        path = cases / "record_cases.tsv"
        rows = table(path)
        rows[0]["display_left_flank"] = "A" + rows[0]["display_left_flank"][1:] if rows[0]["display_left_flank"][0] != "A" else \
            "C" + rows[0]["display_left_flank"][1:]
        reporting_fixture.write_rows(path, rows)
        with self.assertRaisesRegex(ContractError, "sequence level"):
            build_reports(self.root / "refused-cases", **dict(self.inputs, record_cases=cases))

    def test_missing_inputs_are_listed_as_skipped(self):
        minimal = {key: value for key, value in self.inputs.items() if key in {"metrics", "analysis", "design", "search", "metadata"}}
        result = build_reports(self.root / "minimal", **minimal)
        self.assertEqual(result["skipped"], ["components", "cap_review", "target_limit", "benchmark_misses", "exact_scan_validation",
                                             "inference", "record_cases", "census"])


if __name__ == "__main__":
    unittest.main()
