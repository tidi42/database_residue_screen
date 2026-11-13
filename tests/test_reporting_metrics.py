"""Central publication metrics on a synthetic study: cohorts, denominators, attrition, zero searches and sensitivities."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

import reporting_fixture
from database_residue_screen.contracts import ContractError, read_table
from database_residue_screen.outcomes import OUTCOME_SCHEMA
from database_residue_screen.reporting.metrics import (DEFAULT_CHANCE_LIMIT_VARIANTS, chance_limit_pairs, Cohorts, compute_metrics, limit_policy, parse_display_classes,
                                                       read_metrics, regate)
from database_residue_screen.sensitivity import ungated_outcomes


def table(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class MetricsCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.study = reporting_fixture.build(cls.root / "study")
        cls.metrics = cls.root / "metrics"
        cls.summary = compute_metrics(cls.study["taxonomy"], cls.metrics, design=cls.study["design"], search=cls.study["search"])
        cls.outcomes = read_table(cls.study["taxonomy"] / "outcomes.tsv", OUTCOME_SCHEMA)
        cls.policy = json.loads((cls.study["taxonomy"] / "analysis.json").read_text())["analysis_policy"]

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def rows(self, name: str, **selection) -> list[dict]:
        return [row for row in table(self.metrics / name) if all(row[key] == str(value) for key, value in selection.items())]

    def one(self, name: str, **selection) -> dict:
        rows = self.rows(name, **selection)
        self.assertEqual(len(rows), 1, selection)
        return rows[0]


class DenominatorTests(MetricsCase):
    def test_endpoint_native_and_common_cohorts(self):
        native = self.one("cohorts.tsv", scope="global", endpoint="strict", population="all_records", measure="detection", cohort="native")
        # 16 parents: 4 censored (10, 10, 16, 19 nt), 1 chance-limited (23 nt), 11 testable and resolved; 7 strict positive.
        self.assertEqual((native["planned_parents"], native["eligible_parents"], native["censored"], native["chance_limited"]), ("16", "11", "4", "1"))
        self.assertEqual((native["reference_positive"], native["control_draws"], native["control_positive"]), ("7", "165", "1"))
        linear = self.one("cohorts.tsv", scope="global", endpoint="strict", population="linear_genomic", measure="detection", cohort="native")
        # The unresolved record identifier of the TruSeq-like parent excludes it (neither 0 nor 1) in the linear-genomic population.
        self.assertEqual((linear["eligible_parents"], linear["unresolved"], linear["reference_positive"]), ("10", "1", "6"))
        tolerant = self.one("cohorts.tsv", scope="global", endpoint="tolerant", population="all_records", measure="detection", cohort="native")
        self.assertEqual((tolerant["eligible_parents"], tolerant["censored"]), ("9", "7"))
        common = self.one("cohorts.tsv", scope="global", endpoint="strict", population="all_records", measure="detection",
                          cohort="common_strict_tolerant")
        self.assertEqual((common["eligible_parents"], common["native_eligible_parents"], common["reference_positive"]), ("9", "11", "6"))
        breadth = self.one("cohorts.tsv", scope="global", endpoint="strict", population="all_records", measure="breadth", cohort="native")
        self.assertEqual((breadth["eligible_parents"], breadth["unresolved"]), ("10", "1"))

    def test_common_cohort_parents_are_resolved_under_both_rules_with_all_controls(self):
        membership = self.rows("cohort_membership.tsv", population="all_records", measure="detection")
        common = {row["parent_id"] for row in membership if row["common_strict_tolerant"] == "yes"}
        states = {(row["parent_id"], row["endpoint"]): row["state"] for row in membership}
        self.assertEqual(len(common), 9)
        self.assertTrue(all(states[parent, "strict"] == states[parent, "tolerant"] == "eligible" for parent in common))
        for parent in common:
            for endpoint in ("strict", "tolerant"):
                family = [row for row in self.outcomes if row["parent_id"] == parent and row["endpoint"] == endpoint
                          and row["population"] == "all_records"]
                self.assertEqual(sum(row["role"] == "control" for row in family), reporting_fixture.DRAWS)
                self.assertTrue(all(row["detection_resolution"] == "resolved" and row["completion"] == "successful" for row in family))

    def test_cohorts_agree_with_the_analysis_contrasts(self):
        contrasts = table(self.study["taxonomy"] / "contrasts.tsv")
        for row in contrasts:
            if row["condition"] != "all_conditions":
                continue
            scope = "global" if row["group"] == "Global union" else "group"
            mine = self.one("cohorts.tsv", scope=scope, group=row["group"], endpoint=row["endpoint"], population=row["population"],
                            measure=row["measure"], cohort="native")
            self.assertEqual((mine["planned_parents"], mine["eligible_parents"], mine["censored"], mine["chance_limited"], mine["unresolved"]),
                             (row["planned_parents"], row["eligible_parents"], row["censored_parents"], row["chance_limited_parents"],
                              row["unresolved_parents"]))
            if row["measure"] == "detection" and row["real_detection_rate"] != "not_estimated":
                self.assertAlmostEqual(float(mine["reference_rate"]), float(row["real_detection_rate"]))
            if row["measure"] == "breadth" and row["real_mean_breadth"] != "not_estimated":
                self.assertAlmostEqual(float(mine["reference_mean"]), float(row["real_mean_breadth"]))

    def test_table1_display_classes_and_conditions(self):
        classes = {(row["group"], row["label"]): row for row in table(self.metrics / "display_classes.tsv")}
        self.assertEqual(list(classes), [("ONT", "Native barcodes"), ("ONT", "K40 contexts"), ("ONT", "Adapters and primers"),
                                         ("Illumina", "Bare i7 indexes"), ("Illumina", "Embedded i7 indexes"), ("Illumina", "Adapters and primers"),
                                         ("BGI/DNBSEQ", "Adapters"), ("PacBio", "Barcodes"), ("PacBio", "Barcode-hairpin queries"),
                                         ("PacBio", "SMRTbell adapter"), ("PacBio", "C2 primer")])
        # Two queries (forward and reverse complement) form one canonical parent.
        bare_i7 = classes["Illumina", "Bare i7 indexes"]
        self.assertEqual((bare_i7["physical_ids"], bare_i7["canonical_parents"]), ("3", "2"))
        # The C2 primer shares its condition with the SMRTbell adapter and is separated by construct class.
        self.assertEqual(classes["PacBio", "C2 primer"]["members"], "pacbio_smrtbell_adapter:primer;pacbio_c2_primer")
        self.assertTrue(all(row["declared"] == "yes" for row in classes.values()))
        c2 = self.one("table1_source.tsv", scope="class", label="C2 primer", endpoint="strict", population="all_records", cohort="native")
        self.assertEqual((c2["planned_parents"], c2["eligible_parents"], c2["reference_positive"], c2["target_limited"]), ("1", "1", "0", "1"))
        condition = self.one("table1_source.tsv", scope="condition", label="pacbio_smrtbell_adapter", endpoint="strict",
                             population="all_records", cohort="native")
        self.assertEqual((condition["planned_parents"], condition["reference_positive"]), ("2", "1"))
        scopes = {(row["scope"], row["label"]) for row in table(self.metrics / "table1_source.tsv")}
        self.assertIn(("group", "BGI/DNBSEQ"), scopes)
        self.assertIn(("global", "Global union"), scopes)
        for row in table(self.metrics / "table1_source.tsv"):
            if row["cohort"] == "native":
                attrition = sum(int(row[name]) for name in ("eligible_parents", "censored", "chance_limited", "partially_censored", "unresolved",
                                                            "incomplete", "other_inapplicable"))
                self.assertEqual(attrition, int(row["planned_parents"]), row)

    def test_undeclared_conditions_get_their_own_class(self):
        classes = parse_display_classes([dict(group="ONT", label="Barcodes", members=["ont_barcode"])])
        destination = self.root / "metrics-undeclared"
        compute_metrics(self.study["taxonomy"], destination, display_classes=[dict(group="ONT", label="Barcodes", members=["ont_barcode"])])
        rows = table(destination / "display_classes.tsv")
        self.assertEqual(rows[0]["label"], "Barcodes")
        self.assertEqual(len(classes), 1)
        self.assertIn("pacbio_smrtbell_adapter (primer)", {row["label"] for row in rows if row["declared"] == "no"})
        with self.assertRaises(ContractError):
            parse_display_classes([("ONT", "Same", ("a",)), ("ONT", "Same", ("b",))])
        # The same label may name one class in each of two groups.
        self.assertEqual(len(parse_display_classes([("ONT", "Same", ("a",)), ("Illumina", "Same", ("b",))])), 2)

    def test_conversion_against_own_controls(self):
        row = self.one("conversion.tsv", cohort="testable", scope="global", population="all_records")
        self.assertEqual((row["common_parents"], row["strict_reference_positive"], row["tolerant_reference_positive"]), ("9", "6", "6"))
        self.assertEqual((row["strict_negative_references"], row["converted_references"]), ("3", "0"))
        self.assertEqual((row["converted_control_draws"], row["strict_negative_control_draws"]), ("1", "135"))
        framework = {(item["cohort"], item["group"], item["condition"], item["population"]): item
                     for item in table(self.study["taxonomy"] / "common_rule_contrasts.tsv")}
        names = {"testable": "primary_testable", "ungated": "ungated_diagnostic"}
        for mine in table(self.metrics / "conversion.tsv"):
            if mine["scope"] not in {"group", "global", "condition"}:
                continue
            condition = mine["label"] if mine["scope"] == "condition" else "all_conditions"
            other = framework[names[mine["cohort"]], mine["group"], condition, mine["population"]]
            for name in ("common_parents", "converted_references", "converted_control_draws", "own_converted_control_draws",
                         "strict_negative_control_draws", "references_without_strict_negative_controls"):
                self.assertEqual(mine[name], other[name], (mine["cohort"], mine["group"], condition, name))

    def test_breadth_lower_bound_adds_detection_resolved_parents(self):
        row = self.one("breadth_lower_bound.tsv", group="Illumina", endpoint="strict", population="all_records")
        self.assertEqual((row["breadth_resolved_parents"], row["unresolved_breadth_parents"]), ("2", "1"))
        self.assertIn("+1unresolved", row["unresolved_parents"])


class AttritionAndZeroSearchTests(MetricsCase):
    def test_attrition_chain(self):
        steps = {(row["step"], row["endpoint"], row["population"], row["measure"]): int(row["count"])
                 for row in self.rows("attrition.tsv", group="Global union")}
        self.assertEqual(steps["physical_ids", "", "", ""], 17)
        self.assertEqual(steps["canonical_parents", "", "", ""], 16)
        self.assertEqual(steps["logical_queries", "", "", ""], 17 + 16 * reporting_fixture.DRAWS)
        searches = steps["distinct_searches", "", "", ""]
        self.assertEqual(searches, steps["completed_with_hits", "", "", ""] + steps["completed_zero", "", "", ""] + steps["incomplete", "", "", ""])
        self.assertEqual((steps["incomplete", "", "", ""], steps["capped", "", "", ""]), (0, 2))
        self.assertEqual((steps["testability_testable", "strict", "", ""], steps["testability_censored", "strict", "", ""],
                          steps["testability_chance_limited", "strict", "", ""]), (11, 4, 1))
        self.assertEqual(steps["resolved", "strict", "linear_genomic", "detection"] + steps["unresolved", "strict", "linear_genomic", "detection"], 11)
        self.assertEqual(searches, self.summary["counts"]["distinct_searches"])

    def test_design_physical_ids_must_match(self):
        design = self.root / "design-copy"
        design.mkdir()
        rows = (self.study["design"] / "physical_queries.tsv").read_text().splitlines()
        (design / "physical_queries.tsv").write_text("\n".join(rows[:-1]) + "\n")
        with self.assertRaisesRegex(ContractError, "physical IDs"):
            compute_metrics(self.study["taxonomy"], self.root / "metrics-design", design=design)

    def test_zero_search_accounting(self):
        summary = self.summary["zero_searches"]
        self.assertEqual(summary["score_forced_lengths"], [10, 16])
        self.assertEqual(summary["capped"], 2)
        self.assertEqual(summary["incomplete"], 0)
        self.assertEqual(summary["distinct_searches"], summary["with_hits_uncapped"] + summary["capped"] + summary["score_forced_zero"]
                         + summary["ordinary_zero"] + summary["incomplete"])
        statuses = {row["query_id"]: row for row in table(self.study["taxonomy"] / "search_status.tsv")}
        lengths = {row["sequence_id"]: int(row["length"]) for row in table(self.study["taxonomy"] / "logical_queries.tsv")}
        forced = [sequence for sequence, length in lengths.items() if length in (10, 16)]
        self.assertEqual(summary["score_forced_zero"], len(forced))
        self.assertTrue(all(statuses[sequence]["status"] == "success_zero" for sequence in forced))
        by_length = self.rows("zero_searches.tsv", query_length=10)
        self.assertTrue(all(row["score_forced"] == "yes" and row["with_hits"] == "0" for row in by_length))
        units = {(row["unit"], row["length_class"], row["state"]): int(row["count"]) for row in table(self.metrics / "zero_search_units.tsv")}
        self.assertEqual(units["reference_physical_id", "score_forced", "score_forced_zero"], 4)
        self.assertEqual(self.summary["archived_search_attempts"], dict(searched=True, failed_attempts=0, failed_attempt_queries=0,
                                                                        incomplete_attempts=0, incomplete_attempt_queries=0))

    def test_archived_failed_attempts_are_counted(self):
        search = self.root / "search-with-failure"
        search.mkdir()
        attempt = search / "batch_00001.attempt_00001"
        attempt.mkdir()
        (attempt / "failed.json").write_text(json.dumps(dict(query_ids=["a", "b"], query_status="failed_not_zero")))
        from database_residue_screen.reporting.metrics import archived_attempts
        self.assertEqual(archived_attempts(search)["failed_attempt_queries"], 2)


class SensitivityTests(MetricsCase):
    def test_ungated_and_exact_score_panels(self):
        def global_row(panel, gating, population="all_records"):
            return self.one("sensitivity_panels.tsv", panel=panel, gating=gating, group="Global union", endpoint="strict", population=population)
        ungated = global_row("ungated_panel", "ungated")
        self.assertEqual((ungated["panel_parents"], ungated["eligible_parents"], ungated["reference_positive"]), ("16", "16", "9"))
        # Unresolved evidence stays excluded without the gate: it is never counted as zero.
        linear = global_row("ungated_panel", "ungated", "linear_genomic")
        self.assertEqual((linear["eligible_parents"], linear["unresolved"], linear["reference_positive"]), ("15", "1", "8"))
        exact = global_row("exact_score_panel", "testable")
        self.assertEqual((exact["panel_parents"], exact["eligible_parents"], exact["reference_positive"]), ("13", "11", "7"))
        exact_ungated = global_row("exact_score_panel", "ungated")
        self.assertEqual((exact_ungated["eligible_parents"], exact_ungated["reference_positive"]), ("13", "9"))
        self.assertEqual(self.summary["exact_score_panel"], dict(lengths=[19, 23, 24, 32, 33, 34, 35, 40, 45, 61], parents=13))

    def test_panels_agree_with_the_analysis_length_sensitivity(self):
        names = {"full_panel_inventory": ("ungated_panel", "ungated"), "baseline_exact_calibrated_length_panel": ("exact_score_panel", "ungated")}
        for row in table(self.study["taxonomy"] / "length_sensitivity.tsv"):
            if row["cohort"] not in names:
                continue
            panel, gating = names[row["cohort"]]
            mine = self.one("sensitivity_panels.tsv", panel=panel, gating=gating, group=row["group"], endpoint=row["endpoint"],
                            population=row["population"])
            self.assertEqual((mine["eligible_parents"], mine["reference_positive"]), (row["ungated_parents"], row["ungated_positive"]))

    def test_incomplete_evidence_is_never_zero(self):
        draws = reporting_fixture.DRAWS
        parent = next(row["parent_id"] for row in self.outcomes if row["logical_id"] == "ILLUMINA_I7ADP_B_forward")
        broken = [regate(dict(row, completion="incomplete_not_zero"), self.policy) if row["parent_id"] == parent else row for row in self.outcomes]
        base, gated, ungated = Cohorts(self.outcomes, draws), Cohorts(broken, draws), Cohorts(ungated_outcomes(broken, self.policy), draws)
        for cohorts in (gated, ungated):
            self.assertEqual(cohorts.states["strict", "all_records", "detection"][parent], "incomplete")
        parents = base.parents()
        before = base.summary("strict", "all_records", "detection", parents)
        after = gated.summary("strict", "all_records", "detection", parents)
        self.assertEqual((after["planned_parents"], after["eligible_parents"], after["incomplete"]),
                         (before["planned_parents"], before["eligible_parents"] - 1, 1))
        self.assertEqual(after["reference_positive"], before["reference_positive"] - 1)
        ungated_summary = ungated.summary("strict", "all_records", "detection", parents)
        self.assertEqual((ungated_summary["eligible_parents"], ungated_summary["incomplete"]), (15, 1))
        gated.check(broken)

    def test_chance_limit_pairs(self):
        self.assertEqual(DEFAULT_CHANCE_LIMIT_VARIANTS, (("stricter", 0.001, 0.01), ("looser", 0.05, 0.10)))
        policy = dict(detectability=[dict(testable_chance_limit=0.01, informativeness_limit=0.05)] * 2)
        self.assertEqual(chance_limit_pairs(policy), (("baseline", 0.01, 0.05), ("stricter", 0.001, 0.01), ("looser", 0.05, 0.10)))
        for broken in (dict(detectability=[dict(testable_chance_limit=0.01, informativeness_limit=0.05),
                                           dict(testable_chance_limit=0.02, informativeness_limit=0.05)]),):
            with self.assertRaises(ContractError):
                chance_limit_pairs(broken)
        with self.assertRaises(ContractError):
            chance_limit_pairs(policy, (("baseline", 0.001, 0.01),))
        with self.assertRaises(ContractError):
            chance_limit_pairs(policy, (("odd", 0.05, 0.01),))
        def pair(label, endpoint="strict"):
            return self.one("chance_limits.tsv", panel=label, group="Global union", endpoint=endpoint, population="all_records")
        baseline, stricter, looser = pair("baseline"), pair("stricter"), pair("looser")
        native = self.one("cohorts.tsv", scope="global", endpoint="strict", population="all_records", measure="detection", cohort="native")
        self.assertEqual((baseline["eligible_parents"], baseline["reference_positive"]), (native["eligible_parents"], native["reference_positive"]))
        # 0.001/0.01: 24 nt becomes chance-limited (bound 0.0077), 23 nt censored.
        self.assertEqual((stricter["eligible_parents"], stricter["reference_positive"], stricter["chance_limited"]), ("9", "6", "2"))
        # 0.05/0.10: the 23-nt query becomes testable (bound 0.031).
        self.assertEqual((looser["eligible_parents"], looser["reference_positive"], looser["chance_limited"]), ("12", "8", "0"))
        self.assertEqual(looser["testable_lengths"].split(";")[0], "23")
        self.assertEqual(pair("stricter", "tolerant")["eligible_parents"], "8")
        classes = {(row["query_length"], row["endpoint"]): row for row in table(self.metrics / "testability.tsv")}
        self.assertEqual((classes["24", "strict"]["class_baseline"], classes["24", "strict"]["class_stricter"],
                          classes["23", "strict"]["class_looser"]), ("testable", "chance_limited", "testable"))
        self.assertEqual((classes["10", "strict"]["physical_ids"], classes["10", "strict"]["canonical_parents"]), ("3", "2"))

    def test_baseline_regating_reproduces_the_outcomes(self):
        policy = limit_policy(self.policy, 0.01, 0.05)
        self.assertEqual([regate(row, policy) for row in self.outcomes], self.outcomes)
        looser = limit_policy(self.policy, 0.05, 0.10)
        changed = [row for row in self.outcomes if regate(row, looser) != row]
        self.assertTrue(changed)
        self.assertTrue(all(row["length"] == 23 and row["endpoint"] == "strict" for row in changed))


class OutputTests(MetricsCase):
    def test_outputs_are_deterministic_and_verified(self):
        again = self.root / "metrics-again"
        compute_metrics(self.study["taxonomy"], again, design=self.study["design"], search=self.study["search"])
        for path in sorted(self.metrics.iterdir()):
            self.assertEqual(path.read_bytes(), (again / path.name).read_bytes(), path.name)
        self.assertEqual(read_metrics(self.metrics)["contract"], "publication_metrics_v1")
        manifest = json.loads((self.metrics / "manifest.json").read_text())
        self.assertNotIn("created_utc", manifest)
        self.assertIn("table1_source.tsv", manifest["outputs"])

    def test_changed_inputs_and_existing_outputs_are_refused(self):
        with self.assertRaises(ContractError):
            compute_metrics(self.study["taxonomy"], self.metrics)
        tampered = self.root / "metrics-tampered"
        tampered.mkdir()
        (tampered / "x").write_text("x")
        with self.assertRaises(ContractError):
            compute_metrics(self.study["taxonomy"], tampered)
        path = self.metrics / "cohorts.tsv"
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"\n")
            with self.assertRaises(ContractError):
                read_metrics(self.metrics)
        finally:
            path.write_bytes(original)


if __name__ == "__main__":
    unittest.main()
