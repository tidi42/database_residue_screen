import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.outcomes import build_outcomes
from database_residue_screen.policy import default_policy
import math

from database_residue_screen.precision import (grouped_half_width, mean_interval, precision_specification, precision_rows, relative_entropy,
                                             relative_entropy_limits, required_units, validate_precision_specification)
import test_outcomes as outcome_fixtures


class ControlPrecisionTests(unittest.TestCase):
    def test_seed_reuse_widens_bound_and_zero_outcomes_do_not_mean_zero_uncertainty(self):
        independent = grouped_half_width([1, 2, 3, 4], 0.05)
        shared = grouped_half_width([1, 1, 2, 2], 0.05)
        self.assertGreater(shared, independent)
        self.assertGreater(independent, 0)
        with self.assertRaises(ContractError):
            grouped_half_width([0], 0.05)

    def test_relative_entropy_limits_match_closed_forms(self):
        alpha = 0.05 / 3
        for units in (30, 1470, 13860):
            with self.subTest(units=units):
                lower, upper = relative_entropy_limits(0.0, units, alpha)
                self.assertEqual(lower, 0.0)
                self.assertAlmostEqual(upper, 1 - (alpha / 2) ** (1 / units), places=12)
                lower, upper = relative_entropy_limits(1.0, units, alpha)
                self.assertAlmostEqual(lower, (alpha / 2) ** (1 / units), places=12)
                self.assertEqual(upper, 1.0)
        lower, upper = relative_entropy_limits(0.3, 500, alpha)
        level = math.log(2 / alpha) / 500
        self.assertAlmostEqual(relative_entropy(0.3, lower), level, places=9)
        self.assertAlmostEqual(relative_entropy(0.3, upper), level, places=9)
        self.assertEqual(relative_entropy(0.0, 1.0), math.inf)
        self.assertEqual(relative_entropy(1.0, 1.0), 0.0)

    def test_relative_entropy_interval_is_never_wider_than_additive_hoeffding(self):
        alpha = 0.05 / 3
        for values in ([0] * 100, [0] * 99 + [1], [0, 1] * 50, [0] * 98 + [2000, 0]):
            bound = max(1, max(values))
            seeds = list(range(1, len(values) + 1))
            lower, upper, half_width = mean_interval(values, seeds, alpha, bound)
            estimate = sum(values) / len(values)
            self.assertLessEqual(half_width, grouped_half_width(seeds, alpha, bound) + 1e-12)
            self.assertLessEqual(lower, estimate)
            self.assertGreaterEqual(upper, estimate)
            self.assertTrue(0 <= lower and upper <= bound)

    def test_unequal_seed_streams_use_the_grouped_additive_bound(self):
        alpha = 0.05
        values, seeds = [0, 0, 0, 1], [1, 1, 2, 3]
        lower, upper, _ = mean_interval(values, seeds, alpha)
        radius = grouped_half_width(seeds, alpha)
        self.assertEqual((lower, upper), (max(0.0, 0.25 - radius), min(1.0, 0.25 + radius)))
        lower, upper, half_width = mean_interval([0, 0, 1, 1], [1, 1, 2, 2], alpha)
        self.assertAlmostEqual(upper, relative_entropy_limits(0.5, 2, alpha)[1])

    def test_required_units_reach_the_target_and_one_fewer_does_not(self):
        alpha = 0.05 / 3
        for estimate, target, bound in ((0.0, 0.01, 1), (0.0003, 0.01, 1), (0.0, 1.0, 2000), (0.5, 0.05, 1)):
            with self.subTest(estimate=estimate, target=target):
                units = required_units(estimate, target, alpha, bound)
                for count, met in ((units, True), (units - 1, False)):
                    lower, upper = relative_entropy_limits(estimate / bound, count, alpha)
                    width = bound * max(upper - estimate / bound, estimate / bound - lower)
                    self.assertEqual(width <= target, met)
        self.assertEqual(required_units(0.0, 0.01, alpha), math.ceil(math.log(2 / alpha) / -math.log(1 - 0.01)))

    def test_draw_prefixes_keep_parent_weights_and_unsearched_prefixes_explicit(self):
        fixture = outcome_fixtures.MatchedOutcomeTests()
        fixture.setUp()
        outcomes, _, _, _ = build_outcomes(fixture.events, fixture.logical, fixture.statuses, default_policy(), 2)
        spec = precision_specification([1, 2, 3], [500, 2000])
        self.assertEqual(spec["maximum_draws"], 3)
        rows, parents = precision_rows(outcomes, spec, 500, 2)
        selected = next(row for row in rows if row["group"] == "Global union" and row["endpoint"] == "strict" and
                        row["population"] == "all_records" and row["draws_per_parent"] == 2 and row["measure"] == "detection")
        self.assertEqual(float(selected["estimate"]), 0.75)
        self.assertEqual(selected["control_draws"], 4)
        self.assertGreater(float(selected["half_width"]), 0)
        self.assertTrue(all(row["status"] == "not_searched" for row in rows if row["draws_per_parent"] == 3))
        zero = next(row for row in parents if row["parent_id"] == "parent2" and row["endpoint"] == "strict" and row["population"] == "linear_genomic" and row["draws"] == 2)
        self.assertGreater(float(zero["zero_success_upper"]), 0)
        self.assertEqual(zero["distinct_sequences"], 1)
        self.assertEqual(zero["draws"], 2)

    def test_outcome_dependent_attrition_does_not_certify_coverage(self):
        fixture = outcome_fixtures.MatchedOutcomeTests()
        fixture.setUp()
        fixture.events[-1].update(molecule="unknown", record_type="genomic")
        outcomes, _, _, _ = build_outcomes(fixture.events, fixture.logical, fixture.statuses, default_policy(), 2)
        spec = precision_specification([2], [500, 2000])
        rows, _ = precision_rows(outcomes, spec, 500, 2)
        selected = next(row for row in rows if row["group"] == "Global union" and row["endpoint"] == "strict" and
                        row["population"] == "linear_genomic" and row["measure"] == "detection")
        self.assertEqual(selected["precision_met"], "not_estimated")
        self.assertEqual(selected["status"], "diagnostic_resolved_subset_no_coverage_claim")

    def test_ladders_are_explicit_increasing_and_within_budget(self):
        spec = precision_specification([5, 10, 15], [500, 2000])
        self.assertEqual((spec["draw_ladder"], spec["retention_ladder"], spec["maximum_draws"]), ([5, 10, 15], [500, 2000], 15))
        for draws, retention, maximum in (([10, 5], [500, 2000], None), ([5, 10], [2000, 500], None), ([5, 10], [500], None),
                                          ([], [500, 2000], None), ([5, 10, 15], [500, 2000], 10)):
            with self.subTest(draws=draws, retention=retention, maximum=maximum), self.assertRaises(ContractError):
                precision_specification(draws, retention, maximum)
        with self.assertRaises(ContractError):
            validate_precision_specification(dict(spec, stopping="stop_when_significant"))

    def test_positive_with_unresolved_breadth_has_measure_specific_coverage(self):
        fixture = outcome_fixtures.MatchedOutcomeTests()
        fixture.setUp()
        fixture.events.append(fixture.event("parent1", "unresolved:extra", "DNA", "genomic"))
        outcomes, _, _, _ = build_outcomes(fixture.events, fixture.logical, fixture.statuses, default_policy(), 2)
        spec = precision_specification([2], [500, 2000])
        rows, parents = precision_rows(outcomes, spec, 2000, 2)
        selected = {row["measure"]: row for row in rows if row["group"] == "Global union" and row["endpoint"] == "strict" and row["population"] == "all_records"}
        self.assertEqual(selected["detection"]["matched_parents"], 2)
        self.assertEqual(selected["detection"]["status"], "conditional_bound")
        self.assertEqual(selected["breadth"]["matched_parents"], 1)
        self.assertEqual(selected["breadth"]["status"], "diagnostic_resolved_subset_no_coverage_claim")
        selected_parents = {row["measure"]: row for row in parents if row["parent_id"] == "parent1" and row["endpoint"] == "strict" and row["population"] == "all_records"}
        self.assertEqual(selected_parents["detection"]["unresolved_draws"], 0)
        self.assertEqual(selected_parents["breadth"]["unresolved_draws"], 1)
        self.assertEqual(selected_parents["breadth"]["mean_breadth"], "not_estimated")


if __name__ == "__main__":
    unittest.main()