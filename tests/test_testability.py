"""Protocol-v2 classification, independent of observed reference/control output."""
import unittest
from database_residue_screen.length_aware import chance_bound, score_calibration
# Aliased so that pytest does not collect the imported classifier as a test function.
from database_residue_screen.length_aware import testability as classify_testability
from database_residue_screen.contracts import ContractError
from database_residue_screen.policy import applicability, default_policy, validate_policy
from database_residue_screen.outcomes import build_outcomes
from test_length_aware import BASES, SEQUENCES, STATS

class TestabilityTests(unittest.TestCase):
    def test_declared_band_boundaries_and_invalid_limits(self):
        self.assertEqual(classify_testability(True,.01), 'testable')
        self.assertEqual(classify_testability(True,.05), 'chance_limited')
        self.assertEqual(classify_testability(True,.050001), 'censored_untestable')
        self.assertEqual(classify_testability(False,0), 'censored_untestable')
        for limits in ((.01,.05),(float('inf'),.01),(0,0)):
            with self.assertRaises(ContractError): classify_testability(True,0,*limits)

    def test_chance_model_and_fixed_space_do_not_hide_database_growth(self):
        self.assertAlmostEqual(chance_bound(20,1000000,'strict'),2e6/4**20)
        # A 19-nt tolerant query permits exact spans 18 and 19, at two and one positions.
        self.assertAlmostEqual(chance_bound(19,1000000,'tolerant'),2e6*(2/4**18+1/4**19))
        small=score_calibration(24,1000000,10,STATS,searchsp=1000)
        large=score_calibration(24,BASES,SEQUENCES,STATS,searchsp=1000)
        self.assertEqual(small['exact_evalue'],large['exact_evalue'])
        self.assertEqual(small['strict_testability'],'testable')
        self.assertEqual(large['strict_testability'],'testable')
        self.assertEqual(large['tolerant_testability'],'censored_untestable')
        self.assertGreater(large['tolerant_expected_chance_bound'],small['tolerant_expected_chance_bound'])

    def test_diagnostic_retention_never_promotes_chance_saturation(self):
        for length in (8,10,16,19,20,21,22):
            production=score_calibration(length,BASES,SEQUENCES,STATS)
            diagnostic=score_calibration(length,BASES,SEQUENCES,STATS,mode='adaptive')
            self.assertEqual(production['evalue'],10)
            self.assertTrue(diagnostic['strict_score_retainable'])
            self.assertEqual(diagnostic['strict_testability'],'censored_untestable')
            self.assertTrue(diagnostic['diagnostic_chance_saturated_output'])
        for length in range(19,26):
            self.assertTrue(score_calibration(length,BASES,SEQUENCES,STATS)['partially_score_censored'])

    def test_chance_limited_positive_and_zero_are_both_non_estimates(self):
        # Strict 23 nt lies between the declared .01/.05 thresholds.
        policy=default_policy();policy['record_populations']=['all_records']
        policy['detectability']=[dict(score_calibration(23,BASES,SEQUENCES,STATS),exact_witness_recovered=True)]
        self.assertEqual(applicability(23,'strict',policy),'chance_limited')
        logical=[dict(logical_id=str(i),parent_id='p',sequence_id='p' if i==0 else 'c',role='real' if i==0 else 'control',
                      group='g' if i==0 else '',condition='x' if i==0 else '',draw_index=i,length=23,seed=i) for i in (0,1)]
        statuses=[dict(query_id=q,status='success_zero',target_limit_reached='no') for q in ('p','c')]
        outcomes,contrasts,_,_=build_outcomes([],logical,statuses,policy,1)
        self.assertTrue(all(o['positive']==-1 for o in outcomes))
        strict=next(c for c in contrasts if c['endpoint']=='strict')
        self.assertEqual(strict['chance_limited_parents'],1)
        self.assertEqual(strict['eligible_parents'],0)
        # Tampered class cannot bypass the denominator gate.
        policy['detectability'][0]['strict_testability']='testable'
        with self.assertRaises(ContractError): validate_policy(policy)

    def test_observed_positive_retains_evidence_without_entering_primary(self):
        from test_outcomes import MatchedOutcomeTests
        case=MatchedOutcomeTests();case.setUp()
        for row in case.logical: row['length']=23
        policy=default_policy()
        policy['detectability']=[dict(score_calibration(23,BASES,SEQUENCES,STATS),exact_witness_recovered=True)]
        outcomes,contrasts,_,_=build_outcomes(case.events,case.logical,case.statuses,policy,2)
        self.assertTrue(any(o['retained_accessions']>0 for o in outcomes))
        self.assertTrue(all(o['positive']==-1 for o in outcomes))
        self.assertTrue(all(c['eligible_parents']==0 for c in contrasts))
