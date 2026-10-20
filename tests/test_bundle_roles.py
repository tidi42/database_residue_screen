import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.outcomes import CONTRAST_V1_SCHEMA, OUTCOME_V1_SCHEMA, build_outcomes, matched_population
from database_residue_screen.policy import default_policy, primary_event, validate_policy
from database_residue_screen.sensitivity import common_rule_rows, detection_summaries

ROLES = {"keep": dict(panel_version="v1", role="primary"),
         "old": dict(panel_version="v1", role="superseded", superseded_by=["new"]),
         "new": dict(panel_version="v2", role="primary"),
         "alt": dict(panel_version="v2", role="alternative")}


class BundleRoleTests(unittest.TestCase):
    def setUp(self):
        self.policy = dict(default_policy(), condition_roles=ROLES)
        self.logical, self.statuses, self.events = [], [], []
        for parent, condition in (("p_keep", "keep"), ("p_old", "old"), ("p_new", "new"), ("p_alt", "alt")):
            self.logical.append(dict(logical_id=parent, role="real", parent_id=parent, sequence_id=parent, group="ONT",
                                     condition=condition, draw_index=0, length=40, seed=0))
            self.statuses.append(dict(query_id=parent, status="success_hits", target_limit_reached="no"))
            self.events.append(dict(query_id=parent, accession_version=f"{parent}.1", strict="yes", tolerant="yes", molecule="DNA",
                                    record_type="genomic", topology="linear", identity_status="exact_version"))
            for draw in (1, 2):
                control = f"{parent}_d{draw}"
                self.logical.append(dict(logical_id=control, role="control", parent_id=parent, sequence_id=control, group="",
                                         condition="", draw_index=draw, length=40, seed=draw + 10 * len(self.logical)))
                self.statuses.append(dict(query_id=control, status="success_zero", target_limit_reached="no"))

    def build(self):
        return build_outcomes(self.events, self.logical, self.statuses, self.policy, 2)

    def test_unions_use_primary_roles_and_conditions_keep_flagged_estimates(self):
        outcomes, contrasts, _, _ = self.build()
        def row(group, condition):
            return next(r for r in contrasts if r["group"] == group and r["condition"] == condition and r["endpoint"] == "strict"
                        and r["population"] == "all_records" and r["measure"] == "detection")
        union = row("Global union", "all_conditions")
        self.assertEqual((union["planned_parents"], union["real_positive"]), (2, 2))
        self.assertEqual((union["condition_role"], union["panel_version"]), ("primary_union", "v1;v2"))
        self.assertEqual(row("ONT", "all_conditions")["planned_parents"], 2)
        for condition, role, version in (("old", "superseded", "v1"), ("alt", "alternative", "v2"), ("new", "primary", "v2")):
            flagged = row("ONT", condition)
            self.assertEqual((flagged["condition_role"], flagged["panel_version"], flagged["eligible_parents"]), (role, version, 1))
        controls = {r["condition_role"] for r in outcomes if r["role"] == "control" and r["parent_id"] == "p_old"}
        self.assertEqual(controls, {"superseded"})

    def test_historical_inventory_selects_panel_version_not_role(self):
        outcomes, _, _, _ = self.build()
        v1 = detection_summaries(outcomes, 2, roles=None, panel_versions=("v1",))
        primary = detection_summaries(outcomes, 2)
        self.assertEqual(v1["Global union", "strict", "all_records"]["planned_parents"], 2)
        self.assertEqual(primary["Global union", "strict", "all_records"]["planned_parents"], 2)
        _, _, eligible, _ = matched_population(outcomes, "Global union", "all_conditions", "strict", "all_records", 2,
                                               measure="detection", roles=None)
        self.assertEqual(set(eligible), {"p_keep", "p_old", "p_new", "p_alt"})

    def test_common_rule_rows_carry_roles(self):
        outcomes, _, _, _ = self.build()
        rows = common_rule_rows(outcomes, self.policy, 2)
        union = next(r for r in rows if r["cohort"] == "primary_testable" and r["group"] == "Global union" and r["population"] == "all_records")
        self.assertEqual((union["condition_role"], union["common_parents"]), ("primary_union", 2))
        self.assertIn("alternative", {r["condition_role"] for r in rows})

    def test_v1_exports_exclude_bundle_columns(self):
        self.assertNotIn("condition_role", OUTCOME_V1_SCHEMA["properties"])
        self.assertNotIn("panel_version", CONTRAST_V1_SCHEMA["properties"])

    def test_policy_without_roles_keeps_historical_behaviour(self):
        self.policy = default_policy()
        outcomes, contrasts, _, _ = self.build()
        self.assertEqual({(r["panel_version"], r["condition_role"]) for r in outcomes}, {("v1", "primary")})
        union = next(r for r in contrasts if r["group"] == "Global union" and r["endpoint"] == "strict" and r["population"] == "all_records")
        self.assertEqual(union["planned_parents"], 4)

    def test_undeclared_condition_and_conflicting_parent_roles_fail(self):
        self.policy = dict(default_policy(), condition_roles={k: v for k, v in ROLES.items() if k != "alt"})
        with self.assertRaises(ContractError):
            self.build()
        self.policy = dict(default_policy(), condition_roles=ROLES)
        self.logical.append(dict(self.logical[0], logical_id="alias_old", condition="old"))
        with self.assertRaises(ContractError):
            self.build()

    def test_policy_validation_of_supersession(self):
        for broken in ({"old": dict(panel_version="v1", role="superseded")},
                       {"old": dict(panel_version="v1", role="superseded", superseded_by=["alt"]), "alt": ROLES["alt"]},
                       {"keep": dict(panel_version="v1", role="primary", superseded_by=["keep"])}):
            with self.assertRaises(ContractError):
                validate_policy(dict(default_policy(), condition_roles=broken))
        validate_policy(dict(default_policy(), condition_roles=ROLES))

    def test_primary_event_gate(self):
        self.assertTrue(primary_event(dict(condition="keep"), self.policy))
        self.assertFalse(primary_event(dict(condition="old"), self.policy))
        self.assertFalse(primary_event(dict(condition="alt"), self.policy))
        self.assertTrue(primary_event(dict(condition="anything"), default_policy()))


if __name__ == "__main__":
    unittest.main()
