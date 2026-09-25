import unittest

from database_residue_screen.contracts import ContractError
from database_residue_screen.policy import applicability, default_policy, validate_policy


class AnalysisPolicyTests(unittest.TestCase):
    def test_short_original_queries_are_visible_without_fixture_cutoffs(self):
        policy = validate_policy(default_policy())
        for endpoint in ("strict", "tolerant"):
            self.assertEqual(applicability(6, endpoint, policy), "below_blast_word_size")
            for length in (7, 12, 20, 24, 32):
                self.assertEqual(applicability(length, endpoint, policy), "operational_retained_search")

    def test_explicit_sensitivity_is_not_an_implicit_benchmark_dependency(self):
        policy = default_policy()
        policy["minimum_query_length"]["strict"] = 24
        validate_policy(policy)
        self.assertEqual(applicability(20, "strict", policy), "excluded_by_declared_length_policy")
        self.assertEqual(applicability(20, "tolerant", policy), "operational_retained_search")

    def test_unknown_fields_and_invalid_primary_settings_are_rejected(self):
        policy = default_policy()
        policy["assume_unknown_linear"] = True
        with self.assertRaises(ContractError):
            validate_policy(policy)
        policy = default_policy()
        policy["position"]["primary_window"] = 999
        with self.assertRaises(ContractError):
            validate_policy(policy)


if __name__ == "__main__":
    unittest.main()