"""Figure sources and the R figure build on a synthetic study."""
import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import reporting_fixture
from test_reporting_reports import report_inputs
from database_residue_screen.contracts import ContractError
from database_residue_screen.reporting.figure_sources import CONTRACT, build_figure_sources
from database_residue_screen.reporting.figures import build_figures
from database_residue_screen.reporting.metrics import compute_metrics
from database_residue_screen.reporting.reports import build_reports

SCRIPTS = Path(__file__).resolve().parents[1] / "studies" / "paper3" / "figures"
RSCRIPT = shutil.which("Rscript")


def table(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def figure_inputs(study: dict, metrics: Path, reports: Path) -> dict:
    return dict(metrics=metrics, reports=reports, analysis=study["taxonomy"], design=study["design"], controls=study["controls"],
                search=study["search"], inference=study["inference"],
                benchmarks={path.name.removeprefix("benchmark_"): path for path in sorted(study["benchmark"].glob("benchmark_*"))},
                record_cases=study["record_cases"])


class FigureSourceCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.study = reporting_fixture.build(cls.root / "study")
        cls.metrics, cls.reports, cls.sources = cls.root / "metrics", cls.root / "reports", cls.root / "sources"
        compute_metrics(cls.study["taxonomy"], cls.metrics, design=cls.study["design"], search=cls.study["search"])
        build_reports(cls.reports, **report_inputs(cls.study, cls.metrics))
        cls.inputs = figure_inputs(cls.study, cls.metrics, cls.reports)
        cls.result = build_figure_sources(cls.sources, **cls.inputs)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()


class FigureSourceTests(FigureSourceCase):
    def test_manifest_and_tables(self):
        manifest = json.loads((self.sources / "manifest.json").read_text())
        self.assertEqual(manifest["contract"], CONTRACT)
        for name, checksum in manifest["files"].items():
            self.assertEqual(hashlib.sha256((self.sources / name).read_bytes()).hexdigest(), checksum, name)
        self.assertEqual(manifest["counts"], dict(physical_references=17, canonical_parents=16, conditions=11, controls=240, draws_per_parent=15,
                                                  distinct_searches=manifest["counts"]["distinct_searches"], testable_strict_accessions=7,
                                                  positional_events=8, cap_sequences=2))
        self.assertEqual(manifest["draw_prefixes"], [5, 10, 15])
        self.assertEqual(manifest["position_policy"]["primary_window"], 500)
        self.assertEqual(sorted(manifest["benchmarks"]), ["core_nt", "core_nt_div10", "core_nt_x10"])
        expected = {"detection.tsv", "paired_breadth.tsv", "eligibility_flow.tsv", "eligibility_flow_scopes.tsv", "position_events.tsv",
                    "position_sensitivity.tsv", "control_draw_qc.tsv", "control_parent_qc.tsv", "control_qc_summary.tsv", "query_properties.tsv",
                    "cap_sequences.tsv", "cap_memberships.tsv", "quality_marginal.tsv", "quality_joint.tsv", "quality_acceptance.tsv",
                    "precision.tsv", "precision_fixed_cohort.tsv", "matching_examples.tsv", "benchmark_metrics.tsv",
                    "short_query_detectability.tsv", "length_testability.tsv", "density.tsv", "logical_queries.tsv", "testability_classes.tsv",
                    "record_cases.tsv", "record_case_segments.tsv", "record_case_verification.tsv", "target_limit_summary.tsv",
                    "target_limit_parents.tsv"}
        self.assertEqual(set(manifest["files"]), expected)
        self.assertFalse({"pilot_common_summary.tsv", "u01_case_verification.tsv"} & set(manifest["files"]))

    def test_detection_rows_come_from_the_central_metrics(self):
        detection = table(self.sources / "detection.tsv")
        # 11 display classes, 4 platform unions and the global union; 2 populations x 2 rules x 2 cohorts each.
        self.assertEqual(len(detection), 16 * 8)
        self.assertEqual({row["block"] for row in detection}, {"Reference classes", "Platform unions"})
        metrics = {(row["label"], row["population"], row["endpoint"], row["cohort"]): row for row in table(self.metrics / "table1_source.tsv")}
        for row in detection:
            source = metrics[row["label"], row["population"], row["endpoint"], row["cohort"]]
            self.assertEqual((row["eligible_parents"], row["reference_positive"], row["control_positive"]),
                             (source["eligible_parents"], source["reference_positive"], source["control_positive"]))

    def test_figure_tables(self):
        breadth = table(self.sources / "paired_breadth.tsv")
        membership = [row for row in table(self.metrics / "cohort_membership.tsv") if (row["endpoint"], row["population"], row["measure"], row["state"])
                      == ("strict", "all_records", "breadth", "eligible")]
        self.assertEqual(len(breadth), sum(len(row["groups"].split(";")) for row in membership))
        self.assertTrue(all(row["reference_cap"] in {"True", "False"} for row in breadth))
        flow = table(self.sources / "eligibility_flow.tsv")
        self.assertNotIn("all", {row["state"] for row in flow})
        self.assertEqual(len(table(self.sources / "position_events.tsv")), 8)
        qc = table(self.sources / "control_qc_summary.tsv")[0]
        self.assertEqual((qc["parents"], qc["selected_draws"], qc["invariants_failed"]), ("16", "240", "0"))
        acceptance = {(row["group"], row["acceptance"]): int(row["hsp_rows"]) for row in table(self.sources / "quality_acceptance.tsv")}
        self.assertEqual(acceptance["PacBio", "Rejected"], 5)
        fixed = table(self.sources / "precision_fixed_cohort.tsv")
        self.assertEqual({row["draws_per_parent"] for row in fixed}, {"5", "10", "15"})
        examples = table(self.sources / "matching_examples.tsv")
        self.assertEqual([row["decision"] for row in examples], ["Strict", "Tolerant", "Tolerant", "Rejected", "Rejected", "Rejected"])
        self.assertEqual(len(table(self.sources / "record_case_verification.tsv")), 2)

    def test_deterministic_and_refusing_foreign_inputs(self):
        again = self.root / "sources-again"
        build_figure_sources(again, **self.inputs)
        self.assertEqual((self.sources / "manifest.json").read_bytes(), (again / "manifest.json").read_bytes())
        other = self.root / "reports-other"
        shutil.copytree(self.reports, other)
        manifest = json.loads((other / "manifest.json").read_text())
        manifest["analysis_sha256"] = "0" * 64
        (other / "manifest.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ContractError, "different analysis"):
            build_figure_sources(self.root / "sources-refused", **dict(self.inputs, reports=other))


@unittest.skipIf(RSCRIPT is None, "Rscript is not available")
class FigureBuildTests(FigureSourceCase):
    """Smoke test of the R figure build (skipped only when Rscript is unavailable)."""

    def test_build_renders_every_figure_deterministically(self):
        first, second = self.root / "figures", self.root / "figures-again"
        result = build_figures(self.sources, first, scripts=SCRIPTS, rscript=RSCRIPT)
        self.assertEqual(result["figures"], {"figure1": "Figure 2", "figure2": "Figure 3", "figure3": "Figure 4", "figure4": "Figure 5",
                                             "figureS1": "Figure S1", "figureS2": "Figure S4", "figureS3": "Figure S3", "figureS4": "Figure S2"})
        manifest = json.loads((first / "figures_manifest.json").read_text())
        self.assertEqual(manifest["contract"], "publication_figures_v1")
        for name in ("paper3_figure1.png", "paper3_figureS4.pdf", "paper3_manuscript_Figure5.png", "tables/table1_reference_classes.csv",
                     "tables/tableS7A_target_limit_summary.csv"):
            self.assertIn(name, manifest["files"])
        self.assertIn("content_sha256", manifest["files"]["paper3_figure1.pdf"])
        artwork = json.loads((first / "paper3_figure_manifest.json").read_text())
        self.assertTrue(all(figure["height_mm"] <= 240 for figure in artwork["figures"].values()))
        build_figures(self.sources, second, scripts=SCRIPTS, rscript=RSCRIPT)
        self.assertEqual((first / "figures_manifest.json").read_bytes(), (second / "figures_manifest.json").read_bytes())
        with self.assertRaises(ContractError):
            build_figures(self.sources, first, scripts=SCRIPTS, rscript=RSCRIPT)


if __name__ == "__main__":
    unittest.main()
