"""Testability reasons, runtime calibration, ungated sensitivity, common-rule and conditional contrasts."""
import json
import math
import shutil
import subprocess
import tempfile
import unittest
import warnings
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_json, read_table, schema, write_table
from database_residue_screen.design import digest
from database_residue_screen.example import save_fasta, save_json
from database_residue_screen.length_aware import score_calibration, verify_calibration_row
from database_residue_screen.outcomes import build_outcomes
from database_residue_screen.policy import default_policy, validate_policy
from database_residue_screen.precision import conditional_rows, precision_specification
from database_residue_screen.sensitivity import common_rule_rows, length_sensitivity, length_testability_rows
from test_length_aware import BASES, SEQUENCES, STATS


def calibrated(*lengths):
    return [dict(score_calibration(length, BASES, SEQUENCES, STATS), exact_witness_recovered=True) for length in lengths]


def family(parent, length, draws, group="G", condition="c"):
    rows = [dict(logical_id=parent, role="real", parent_id=parent, sequence_id=parent, group=group, condition=condition,
                 draw_index=0, length=length, seed=0)]
    rows += [dict(logical_id=f"{parent}d{index}", role="control", parent_id=parent, sequence_id=f"{parent}d{index}", group="",
                  condition="", draw_index=index, length=length, seed=sum(map(ord, parent)) * 1000 + index)
             for index in range(1, draws + 1)]
    return rows


def hit(sequence, accession="TEST001.1", strict=True):
    return dict(query_id=sequence, accession_version=accession, strict="yes" if strict else "no", tolerant="yes",
                molecule="DNA", record_type="genomic", topology="linear", identity_status="exact_version")


def statuses(logical, events):
    hits = {row["query_id"] for row in events}
    return [dict(query_id=row["sequence_id"], status="success_hits" if row["sequence_id"] in hits else "success_zero",
                 target_limit_reached="no") for row in logical]


class TestabilityReasonTests(unittest.TestCase):
    def test_reasons_separate_score_censoring_from_chance_saturation(self):
        row = {row["query_length"]: row for row in calibrated(16, 19, 23, 24, 26)}
        self.assertEqual(row[19]["strict_testability"], "censored_untestable")
        self.assertEqual(row[19]["strict_testability_reason"], "chance_bound_above_informativeness_limit")
        self.assertEqual(row[16]["strict_testability_reason"], "score_not_retainable+chance_bound_above_informativeness_limit")
        self.assertEqual(row[23]["strict_testability_reason"], "chance_bound_between_limits")
        self.assertEqual(row[24]["strict_testability_reason"], "chance_bound_within_testable_limit")
        self.assertEqual(row[24]["tolerant_testability_reason"], "score_not_retainable+chance_bound_above_informativeness_limit")
        # 26 nt: the weakest tolerant HSP is retained at E=10, but chance matches are expected.
        self.assertEqual(row[26]["tolerant_testability_reason"], "chance_bound_above_informativeness_limit")

    def test_tampered_calibration_rows_are_rejected(self):
        for field, value in (("strict_expected_chance_bound", 1e-9), ("strict_testable", False), ("exact_evalue", 1e-9),
                             ("strict_testability_reason", "chance_bound_between_limits"), ("minimum_tolerant_score", 99)):
            row = calibrated(24)[0]
            row[field] = value
            with self.assertRaises(ContractError, msg=field):
                verify_calibration_row(row)
        policy = default_policy()
        policy["detectability"] = calibrated(24)
        policy["detectability"][0]["strict_testable"] = False
        with self.assertRaises(ContractError):
            validate_policy(policy)

    def test_length_table_lists_both_rules_with_inputs(self):
        policy = default_policy()
        policy["detectability"] = calibrated(24, 8)
        rows = length_testability_rows(policy)
        self.assertEqual([(row["query_length"], row["endpoint"]) for row in rows], [(8, "strict"), (8, "tolerant"), (24, "strict"), (24, "tolerant")])
        tolerant24 = rows[3]
        self.assertEqual((tolerant24["minimum_qualifying_score"], tolerant24["testability"], tolerant24["score_retainable"],
                          tolerant24["partially_score_censored"]), (18, "censored_untestable", "no", "yes"))
        self.assertAlmostEqual(tolerant24["evalue_at_minimum_score"], 41.61, places=1)


class SensitivityAndCommonRuleTests(unittest.TestCase):
    def test_ungated_fractions_count_completed_untestable_searches_as_zero(self):
        logical = family("short8", 8, 2) + family("long32", 32, 2)
        events = [hit("long32")]
        policy = default_policy()
        policy["record_populations"] = ["all_records"]
        policy["detectability"] = calibrated(8, 32)
        outcomes, _, _, _ = build_outcomes(events, logical, statuses(logical, events), policy, 2)
        rows = {(row["cohort"], row["group"], row["endpoint"]): row for row in length_sensitivity(outcomes, policy, 2)}
        full = rows["full_panel_inventory", "Global union", "strict"]
        self.assertEqual((full["cohort_parents"], full["eligible_matched_parents"], full["censored_parents"], full["observed_positive"]), (2, 1, 1, 1))
        self.assertEqual(float(full["rate"]), 1.0)
        self.assertEqual((full["ungated_positive"], full["ungated_parents"]), (1, 2))
        self.assertEqual(float(full["ungated_rate"]), 0.5)
        calibrated_panel = rows["baseline_exact_calibrated_length_panel", "Global union", "strict"]
        self.assertEqual((calibrated_panel["cohort_parents"], calibrated_panel["ungated_parents"]), (1, 1))
        primary = rows["primary_testable_parents", "Global union", "strict"]
        self.assertEqual((primary["denominator"], primary["ungated_parents"]), (1, 1))

    def test_common_rule_conversion_uses_each_parents_own_controls(self):
        logical = family("p1", 40, 2) + family("p2", 40, 2) + family("p3", 40, 2)
        events = [hit("p1", strict=False), hit("p1d1", strict=False), hit("p2"), hit("p3d1", strict=False), hit("p3d2", strict=False)]
        policy = default_policy()
        policy["record_populations"] = ["all_records"]
        policy["detectability"] = calibrated(40)
        outcomes, _, _, _ = build_outcomes(events, logical, statuses(logical, events), policy, 2)
        row = next(row for row in common_rule_rows(outcomes, policy, 2)
                   if row["cohort"] == "primary_testable" and row["group"] == "Global union")
        self.assertEqual((row["common_parents"], row["strict_reference_positive"], row["tolerant_reference_positive"]), (3, 1, 2))
        self.assertEqual((row["strict_negative_references"], row["converted_references"], row["lost_strict_references"]), (2, 1, 0))
        self.assertEqual((row["strict_negative_control_draws"], row["converted_control_draws"]), (6, 3))
        self.assertEqual((row["own_strict_negative_control_draws"], row["own_converted_control_draws"]), (4, 3))
        self.assertAlmostEqual(float(row["expected_reference_conversions_under_own_controls"]), 1.5)
        self.assertEqual((float(row["reference_conversion_rate"]), float(row["own_control_conversion_rate"])), (0.5, 0.75))

    def test_untestable_tolerant_rule_leaves_only_the_ungated_diagnostic(self):
        # 24 nt: strict testable, tolerant chance-saturated; no primary common cohort exists.
        logical = family("p24", 24, 2)
        events = [hit("p24", strict=False)]
        policy = default_policy()
        policy["record_populations"] = ["all_records"]
        policy["detectability"] = calibrated(24)
        outcomes, _, _, _ = build_outcomes(events, logical, statuses(logical, events), policy, 2)
        rows = {row["cohort"]: row for row in common_rule_rows(outcomes, policy, 2) if row["group"] == "Global union"}
        self.assertEqual((rows["primary_testable"]["common_parents"], rows["primary_testable"]["strict_eligible_parents"]), (0, 1))
        self.assertEqual(rows["primary_testable"]["reference_conversion_rate"], "not_estimated")
        self.assertEqual(rows["ungated_diagnostic"]["converted_references"], 1)


class ConditionalContrastTests(unittest.TestCase):
    def outcomes(self, groups):
        """Synthetic all-record strict outcomes with fixed group counts and known conditional bounds."""
        rows, seed = [], 0
        for group, parents, positives, control_positive in groups:
            for index in range(parents):
                parent = f"{group}{index}"
                for draw in range(16):
                    seed += 1
                    positive = (index < positives) if draw == 0 else (index * 15 + draw - 1 < control_positive)
                    rows.append(dict(role="real" if draw == 0 else "control", parent_id=parent, group=group if draw == 0 else "",
                                     condition="c" if draw == 0 else "", endpoint="strict", population="all_records",
                                     applicability="operational_retained_search", completion="successful",
                                     detection_resolution="resolved", breadth_resolution="resolved", retrieval="retained_not_exhaustive",
                                     positive=int(positive), retained_accessions=int(positive), draw_index=draw, seed=seed, length=24))
        return rows

    def test_reproduces_known_detection_bounds(self):
        specification = precision_specification([5, 10, 15], [500, 2000], maximum_draws=100)
        groups = [("BGI", 2, 2, 0), ("Illumina", 796, 56, 17), ("ONT", 296, 70, 2), ("PacBio", 194, 1, 0)]
        rows = {(row["group"], row["measure"]): row for row in conditional_rows(self.outcomes(groups), specification, 2000)}
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["holm_family_size"] == 8 for row in rows.values()))
        alpha = 0.05 / 3
        bgi = rows["BGI", "detection"]
        self.assertAlmostEqual(float(bgi["interval_lower"]), (alpha / 2) ** (1 / 30), places=9)
        self.assertEqual((bgi["tail_bound"], bgi["holm_bound"]), ("0", "0"))
        ont = rows["ONT", "detection"]
        self.assertAlmostEqual(100 * float(ont["interval_lower"]), 23.4234, places=3)
        self.assertAlmostEqual(100 * float(ont["interval_upper"]), 23.6471, places=3)
        self.assertEqual(ont["holm_bound"][:5], "7.457")
        self.assertTrue(ont["holm_bound"].endswith("e-513"))
        self.assertTrue(rows["Illumina", "detection"]["holm_bound"].endswith("e-341"))
        pacbio = rows["PacBio", "detection"]
        expected_tail = 6 * math.exp(2910 * math.log(1 - 1 / 194))
        self.assertAlmostEqual(float(pacbio["tail_bound"]) / expected_tail, 1, places=9)
        self.assertAlmostEqual(float(pacbio["holm_bound"]) / (5 * expected_tail), 1, places=9)
        self.assertAlmostEqual(float(pacbio["interval_upper"]), 1 / 194, places=12)
        self.assertEqual(float(rows["PacBio", "breadth"]["holm_bound"]), 1.0)
        self.assertEqual(float(ont["observed_reference_total"]), 70)
        self.assertAlmostEqual(float(ont["expected_reference_total"]), 2 / 15)

    def test_subset_and_excluded_slots_stay_in_the_holm_family(self):
        specification = precision_specification([5, 10, 15], [500, 2000], maximum_draws=100)
        rows = self.outcomes([("A", 3, 3, 0), ("B", 3, 0, 0)])
        for row in rows:
            if row["parent_id"] == "B0":
                row.update(detection_resolution="unresolved", positive=-1)
        result = {(row["group"], row["measure"]): row for row in conditional_rows(rows, specification, 2000)}
        self.assertEqual(result["B", "detection"]["status"], "diagnostic_resolved_subset_no_coverage_claim")
        self.assertEqual(result["B", "detection"]["holm_bound"], "not_estimated")
        self.assertEqual(result["A", "detection"]["holm_family_size"], 4)
        self.assertLessEqual(float(result["A", "detection"]["holm_bound"]), 1)

    def test_extreme_bounds_are_not_printed_as_zero(self):
        from database_residue_screen.precision import format_log_probability
        self.assertTrue(format_log_probability(-2000).endswith("e-869"))
        self.assertEqual(float(format_log_probability(0.0)), 1.0)


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST unavailable")
class RuntimeCalibrationTests(unittest.TestCase):
    def strip_calibration(self, search):
        """Remove the runtime calibration and re-sign the completion records, as for a search from another tool."""
        from database_residue_screen.search import object_digest
        request = read_json(search / "request.json")
        for key in ("length_calibration", "calibration_code_sha256"):
            request.pop(key, None)
        (search / "request.json").unlink()
        save_json(search / "request.json", request)
        complete = read_json(search / "complete.json")
        complete["search_id"] = object_digest(request)
        for batch in complete["batches"]:
            path = search / batch["directory"] / "complete.json"
            record = dict(read_json(path), search_id=complete["search_id"])
            path.unlink()
            save_json(path, record)
            batch["completion_sha256"] = digest(path)
        (search / "complete.json").unlink()
        save_json(search / "complete.json", complete)

    def test_fresh_search_carries_its_calibration_and_uncalibrated_search_is_refused(self):
        from database_residue_screen.analysis import analyse_search
        from database_residue_screen.controls import generate_controls
        from database_residue_screen.database import inspect_database
        from database_residue_screen.panels import prepare_panel
        from database_residue_screen.release import validate_analysis
        from database_residue_screen.search import run_search
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "q.fa", {"short": "ACGTCAGT", "long": "ACGTTGCAAGGCTTACCGATAGCTTGACCATG"})
            rows = []
            for name in ("short", "long"):
                row = {key: "" for key in schema("queries")["properties"]}
                row.update(query_id=name, group="Custom", condition="bare", construct_class="barcode", family_id=name, source_id="synthetic",
                           provenance_status="synthetic", architecture_status="not_applicable")
                rows.append(row)
            write_table(root / "q.tsv", rows, schema("queries"))
            generate_controls(root / "q.fa", root / "controls", 2)
            prepare_panel(root / "q.fa", root / "q.tsv", root / "controls", root / "panel", 2)
            save_fasta(root / "s.fa", {"LONG0001.1": "N" * 30 + "ACGTTGCAAGGCTTACCGATAGCTTGACCATG" + "N" * 30})
            subprocess.run(["makeblastdb", "-in", str(root / "s.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            inspect_database(root / "db", root / "evidence")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                for name in ("search", "uncalibrated"):
                    run_search(root / "panel/search_queries.fasta", root / "db", root / "evidence/database_identity.json", root / name,
                               threads=1, mode="working")
            analyse_search(root / "panel", root / "search", root / "analysis")
            self.assertTrue(validate_analysis(root / "analysis")["valid"])
            report = read_json(root / "analysis/analysis.json")
            self.assertEqual(report["length_calibration_source"], dict(kind="runtime_request"))
            self.assertEqual(report["length_calibration"], read_json(root / "search/request.json")["length_calibration"])
            self.assertTrue((root / "analysis/length_testability.tsv").is_file())
            # Precision assessment also writes the conditional contrasts for the largest declared draw prefix.
            from database_residue_screen.precision import CONDITIONAL_SCHEMA, assess_precision
            specification = dict(precision_specification([1, 2], [500, 2000]), panel_sha256=digest(root / "panel/panel.json"))
            save_json(root / "precision_specification.json", specification)
            precision = assess_precision(root / "analysis", root / "precision_specification.json", root / "precision")
            self.assertIn("conditional_contrasts.tsv", precision["files"])
            conditional = read_table(root / "precision/conditional_contrasts.tsv", CONDITIONAL_SCHEMA)
            self.assertEqual({(row["group"], row["measure"]) for row in conditional}, {("Custom", "detection"), ("Custom", "breadth")})
            # A search without its runtime calibration is never analysed.
            self.strip_calibration(root / "uncalibrated")
            with self.assertRaisesRegex(ContractError, "runtime length calibration"):
                analyse_search(root / "panel", root / "uncalibrated", root / "refused_analysis")
            self.assertFalse((root / "refused_analysis").exists())
            # Changing the calibration recorded by the search afterwards invalidates the result release.
            request = read_json(root / "search/request.json")
            request["length_calibration"]["lengths"][0]["evalue"] = 1000.0
            (root / "search/request.json").write_text(json.dumps(request))
            with self.assertRaises(ContractError):
                validate_analysis(root / "analysis")


class EventTableGateTests(unittest.TestCase):
    def test_event_tables_exclude_untestable_and_intermediate_lengths(self):
        from database_residue_screen.analysis import derive_event_tables, POSITION_SCHEMA, DENSITY_SCHEMA
        policy = default_policy()
        policy['detectability'] = calibrated(19, 23, 24)
        events = [dict(group='G', parent_id=str(length), query_id=str(length), qlen=length, length=length,
                       accession_version=f'TEST{length}.1', strict='yes', tolerant='yes', sstart=1, send=length,
                       slen=20000, distance_to_end=0, canonical_strand='+', metadata_status='exact_version',
                       topology='linear', record_type='genomic') for length in (19, 23, 24)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = derive_event_tables(events, ['G'], root, '', policy)
            density = read_table(root/'density.tsv', DENSITY_SCHEMA)
            positions = read_table(root/'positions.tsv', POSITION_SCHEMA)
            self.assertEqual([row['accession_version'] for row in density], ['TEST24.1'])
            self.assertTrue(all(row['events'] == int(row['minimum_length'] <= 20000) for row in positions))
            self.assertEqual(report['strict_distinct_accessions'], 1)
        # Filtering the derived views must preserve the complete evidence inventory.
        self.assertEqual(len(events), 3)

    def test_taxonomy_unions_apply_native_endpoint_testability(self):
        from database_residue_screen.taxonomy import taxonomy_composition
        policy = default_policy()
        policy['record_populations'] = ['all_records']
        policy['detectability'] = calibrated(19, 24, 32)
        events = [dict(hit(str(length), accession=f'TEST{length}.1'), group='G', qlen=length) for length in (19, 24, 32)]
        assignments = [dict(accession_version=row['accession_version'], rank='genus', taxon_id='1', taxon_name='example',
                            context='cellular', status='resolved') for row in events]
        rows = taxonomy_composition(events, assignments, ['G'], policy, ('genus',))
        denominators = {(row['group'], row['endpoint']): row['denominator'] for row in rows}
        self.assertEqual(denominators['G', 'strict'], 2)
        self.assertEqual(denominators['G', 'tolerant'], 1)
        self.assertEqual(denominators['Global union', 'strict'], 2)

    def test_missing_own_control_conversion_is_not_a_partial_expected_total(self):
        logical = family('p1', 40, 2) + family('p2', 40, 2)
        events = [hit('p1', strict=False), hit('p1d1'), hit('p1d2'), hit('p2d1', strict=False)]
        policy = default_policy()
        policy['record_populations'] = ['all_records']
        policy['detectability'] = calibrated(40)
        outcomes, _, _, _ = build_outcomes(events, logical, statuses(logical, events), policy, 2)
        row = next(row for row in common_rule_rows(outcomes, policy, 2)
                   if row['cohort'] == 'primary_testable' and row['group'] == 'Global union')
        self.assertEqual(row['references_without_strict_negative_controls'], 1)
        self.assertEqual(row['expected_reference_conversions_under_own_controls'], 'not_estimated')


if __name__ == "__main__":
    unittest.main()
