import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.outcomes import build_outcomes, historical_v1_results, matched_population, project_real_events
from database_residue_screen.policy import default_policy


class MatchedOutcomeTests(unittest.TestCase):
    def setUp(self):
        self.logical = [self.query("alias_a", "real", "parent1", "parent1", 0, "A"),
                        self.query("alias_b", "real", "parent1", "parent1", 0, "B"),
                        self.query("original2", "real", "parent2", "parent2", 0, "A"),
                        self.query("draw11", "control", "parent1", "parent2", 1),
                        self.query("draw12", "control", "parent1", "parent1", 2),
                        self.query("draw21", "control", "parent2", "shuffled", 1),
                        self.query("draw22", "control", "parent2", "shuffled", 2)]
        self.statuses = [dict(query_id=name, status="success_zero" if name == "parent2" else "success_hits", target_limit_reached="no")
                         for name in ("parent1", "parent2", "shuffled")]
        self.events = [self.event("parent1", "TEST001.1", "DNA", "genomic"), self.event("parent1", "TEST002.1", "RNA", "transcript"),
                       self.event("shuffled", "TEST003.1", "RNA", "transcript")]

    def query(self, identifier, role, parent, sequence, draw, group=""):
        return dict(logical_id=identifier, role=role, parent_id=parent, sequence_id=sequence, group=group,
                    condition="condition" if group else "", draw_index=draw, length=12,
                    seed=0 if role == "real" else draw + (100 if parent == "parent1" else 200))

    def event(self, sequence, accession, molecule, kind):
        return dict(query_id=sequence, accession_version=accession, strict="yes", tolerant="yes", molecule=molecule,
                    record_type=kind, topology="linear", identity_status="exact_version")

    def test_aliases_and_repeated_draws_have_correct_matched_weights(self):
        outcomes, contrasts, memberships, eligibility = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        global_rows = {row["population"]: row for row in contrasts if row["group"] == "Global union" and row["endpoint"] == "strict" and row["measure"] == "detection"}
        self.assertEqual(global_rows["all_records"]["eligible_parents"], 2)
        self.assertEqual(float(global_rows["all_records"]["parent_weighted_control_rate"]), 0.75)
        self.assertEqual(float(global_rows["linear_genomic"]["parent_weighted_control_rate"]), 0.25)
        self.assertEqual(float(global_rows["linear_genomic"]["real_detection_rate"]), 0.5)
        self.assertEqual(len([row for row in outcomes if row["role"] == "control" and row["endpoint"] == "strict" and row["population"] == "all_records"]), 4)
        real_events = project_real_events(self.events, memberships)
        self.assertEqual(len(real_events), 4)
        self.assertNotIn("shuffled", {row["query_id"] for row in real_events})
        control_memberships = [row for row in memberships if row["sequence_id"] == "parent1" and row["role"] == "control"]
        self.assertEqual({row["group"] for row in control_memberships}, {"A", "B"})

    def test_missing_metadata_excludes_matched_parent_not_control_draws_selectively(self):
        self.events[-1].update(molecule="unknown", record_type="genomic")
        outcomes, contrasts, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        row = next(row for row in contrasts if row["group"] == "Global union" and row["endpoint"] == "strict" and row["population"] == "linear_genomic")
        self.assertEqual(row["planned_parents"], 2)
        self.assertEqual(row["eligible_parents"], 1)
        self.assertEqual(row["unresolved_parents"], 1)
        self.assertEqual(row["control_draws"], 2)
        unresolved = next(row for row in outcomes if row["logical_id"] == "draw21" and row["endpoint"] == "strict" and row["population"] == "linear_genomic")
        self.assertEqual(unresolved["positive"], -1)
        self.assertEqual(unresolved["resolution"], "unresolved")

    def test_positive_with_unresolved_breadth_is_retained_but_not_fully_resolved(self):
        self.events.append(self.event("parent1", "unresolved:gi|123", "DNA", "genomic"))
        outcomes, contrasts, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        original = next(row for row in outcomes if row["logical_id"] == "alias_a" and row["endpoint"] == "strict" and row["population"] == "all_records")
        self.assertEqual(original["positive"], 1)
        self.assertEqual(original["resolution"], "positive_with_unresolved_breadth")
        self.assertEqual(original["detection_resolution"], "resolved")
        self.assertEqual(original["breadth_resolution"], "unresolved")
        detection = matched_population(outcomes, "Global union", "all_conditions", "strict", "all_records", 2, measure="detection")
        breadth = matched_population(outcomes, "Global union", "all_conditions", "strict", "all_records", 2, measure="breadth")
        self.assertEqual(set(detection[2]), {"parent1", "parent2"})
        self.assertEqual(set(breadth[2]), {"parent2"})
        self.assertEqual(breadth[3]["unresolved"], 1)
        comparison = {row["measure"]: row for row in contrasts if row["group"] == "Global union" and row["endpoint"] == "strict" and row["population"] == "all_records"}
        self.assertEqual(comparison["detection"]["eligible_parents"], 2)
        self.assertEqual(comparison["breadth"]["eligible_parents"], 1)
        self.assertEqual(comparison["detection"]["control_draws"], 4)
        self.assertEqual(comparison["breadth"]["control_draws"], 2)
        legacy_outcomes, legacy_contrasts = historical_v1_results(outcomes, 2)
        legacy = next(row for row in legacy_contrasts if row["group"] == "Global union" and row["endpoint"] == "strict" and row["population"] == "all_records")
        self.assertEqual(legacy["eligible_parents"], 1)
        self.assertEqual(float(legacy["real_detection_rate"]), 0)
        self.assertNotIn("detection_resolution", legacy_outcomes[0])
        self.assertNotIn("measure", legacy)

    def test_missing_status_or_duplicate_draw_is_not_zero(self):
        with self.assertRaises(ContractError):
            build_outcomes(self.events, self.logical, self.statuses[:-1], default_policy(), 2)
        self.logical[-1]["draw_index"] = 1
        with self.assertRaises(ContractError):
            build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)

    def test_conflicting_accession_metadata_across_queries_remains_unresolved(self):
        self.events.append(self.event("shuffled", "TEST001.1", "RNA", "transcript"))
        outcomes, _, _, eligibility = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        conflicts = [row for row in eligibility if row["accession_version"] == "TEST001.1" and row["population"] == "linear_genomic"]
        self.assertEqual(len(conflicts), 2)
        self.assertTrue(all(row["state"] == "unresolved" for row in conflicts))
        original = next(row for row in outcomes if row["logical_id"] == "alias_a" and row["endpoint"] == "strict" and row["population"] == "linear_genomic")
        self.assertEqual(original["positive"], -1)

    def test_seed_reuse_does_not_claim_independent_monte_carlo_precision(self):
        _, contrasts, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2, seed_reuse=1)
        self.assertTrue(all(row["mc_detection_se"] == "not_estimated" for row in contrasts))

    def test_partial_positive_from_incomplete_search_is_excluded_for_both_measures(self):
        self.statuses[0]["status"] = "failed"
        outcomes, _, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        for measure in ("detection", "breadth"):
            _, _, eligible, attrition = matched_population(outcomes, "Global union", "all_conditions", "strict", "all_records", 2, measure=measure)
            self.assertNotIn("parent1", eligible)
            self.assertEqual(attrition["incomplete"], 1)
        original = next(row for row in outcomes if row["logical_id"] == "alias_a" and row["endpoint"] == "strict" and row["population"] == "all_records")
        self.assertEqual(original["positive"], 1)
        self.assertEqual(original["detection_resolution"], "not_completed")

    def test_cap_flag_does_not_convert_retained_zero_to_absence_or_unknown(self):
        self.statuses[1]["target_limit_reached"] = "yes"
        outcomes, _, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        original = next(row for row in outcomes if row["logical_id"] == "original2" and row["endpoint"] == "strict" and row["population"] == "all_records")
        self.assertEqual(original["positive"], 0)
        self.assertEqual(original["detection_resolution"], "resolved")
        self.assertEqual(original["retrieval"], "target_limit_reached")

    def test_broad_positive_does_not_resolve_unknown_genomic_eligibility(self):
        self.events[0]["molecule"] = "unknown"
        outcomes, _, _, _ = build_outcomes(self.events, self.logical, self.statuses, default_policy(), 2)
        original = {row["population"]: row for row in outcomes if row["logical_id"] == "alias_a" and row["endpoint"] == "strict"}
        self.assertEqual(original["all_records"]["detection_resolution"], "resolved")
        self.assertEqual(original["linear_genomic"]["detection_resolution"], "unresolved")
        self.assertEqual(original["linear_genomic"]["positive"], -1)


if __name__ == "__main__":
    unittest.main()