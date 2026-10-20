import io
import json
import math
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from database_residue_screen.contracts import ContractError, read_json
from database_residue_screen.controls import null_space_status
from database_residue_screen.database import inspect_database
from database_residue_screen.example import create_example, save_fasta, save_json
from database_residue_screen.length_aware import engine_statistics, score_calibration, minimum_tolerant_score
from database_residue_screen.metadata import summary_record, summary_with_retry
from database_residue_screen.outcomes import build_outcomes
from database_residue_screen.policy import default_policy, applicability
from database_residue_screen.search import run_search, search_settings
from database_residue_screen.sensitivity import length_sensitivity
from database_residue_screen.study_binding import bind_study

STATS = dict(kappa=.710602795216363, lambda_value=1.37406312246009, entropy=1.30724660390929)
BASES, SEQUENCES = 1079005865156, 130787235


class LengthStatisticsTests(unittest.TestCase):
    def test_ncbi_reference_values_and_database_dependence(self):
        for length, adjustment, expected in ((8,6,25785793.2039874),(10,8,1651121.60360146),(16,14,433.485524603448)):
            row = score_calibration(length,BASES,SEQUENCES,STATS,mode='adaptive')
            self.assertEqual(row['length_adjustment'],adjustment)
            self.assertAlmostEqual(row['exact_evalue']/expected,1,places=10)
            self.assertGreater(row['evalue'],row['tolerant_evalue'])
            self.assertTrue(row['strict_testable'])
        small = score_calibration(8,1000,2,STATS,mode='fixed')
        big = score_calibration(8,BASES,SEQUENCES,STATS,mode='fixed')
        self.assertTrue(small['strict_testable'])
        self.assertFalse(big['strict_testable'])

    def test_minimal_tolerant_scores_are_integer_endpoint_boundaries(self):
        for length in range(4,160):
            independently = min(span-4*errors for span in range(1,length+1) for errors in range(span+1)
                                if 10*span >= 9*length and 20*(span-errors) >= 19*span)
            self.assertEqual(minimum_tolerant_score(length),independently)
        self.assertEqual([minimum_tolerant_score(l) for l in (19,20,21,22)],[18,16,16,16])

    def test_full_and_partial_censoring_and_random_budget(self):
        policy=default_policy()
        rows=[score_calibration(l,BASES,SEQUENCES,STATS,mode='fixed') for l in (8,16,19,20,21,22)]
        policy['detectability']=[dict(row,exact_witness_recovered=True) for row in rows]
        self.assertEqual(applicability(16,'strict',policy),'censored_untestable')
        self.assertEqual(applicability(19,'strict',policy),'censored_untestable')
        for length in (19,20,21,22):
            self.assertEqual(applicability(length,'tolerant',policy),'censored_untestable')
            capped=score_calibration(length,BASES,SEQUENCES,STATS,maximum_evalue=10)
            self.assertFalse(capped['tolerant_fully_testable'])
        self.assertFalse(score_calibration(3,BASES,SEQUENCES,STATS)['strict_testable'])

    def test_degeneracy_is_proven_not_inferred_from_draw_repetition(self):
        self.assertEqual(null_space_status('AAAAAAAA'),'proven_singleton_null')
        self.assertEqual(null_space_status('AACAGACAA'),'multiple_sequences_proven')
        self.assertIn('budget',null_space_status('AACAGACAA',budget=1))

    def test_censored_outcome_is_not_zero_or_in_denominator(self):
        logical=[]; statuses=[]
        for length in (8,32):
            parent=f'p{length}'
            for index in (0,1,2):
                name=parent if index==0 else f'{parent}d{index}'
                logical.append(dict(logical_id=name,parent_id=parent,sequence_id=name,role='real' if index==0 else 'control',
                                    group='Other' if index==0 else '',condition='x' if index==0 else '',draw_index=index,
                                    length=length,seed=0 if index==0 else length+index))
                statuses.append(dict(query_id=name,status='success_zero',target_limit_reached='no'))
        policy=default_policy()
        policy['record_populations']=['all_records']
        policy['detectability']=[dict(score_calibration(l,BASES,SEQUENCES,STATS,mode='fixed'),exact_witness_recovered=True) for l in (8,32)]
        outcomes,contrasts,_,_=build_outcomes([],logical,statuses,policy,2)
        short=[r for r in outcomes if r['length']==8]
        self.assertTrue(all(r['positive']==-1 and r['resolution']=='censored_untestable' for r in short))
        strict=next(r for r in contrasts if r['group']=='Global union' and r['endpoint']=='strict' and r['measure']=='detection')
        self.assertEqual((strict['planned_parents'],strict['eligible_parents'],strict['censored_parents']),(2,1,1))
        self.assertAlmostEqual(float(strict['zero_control_probability_upper95']),1-.05**.5)
        sensitivity=length_sensitivity(outcomes,policy,2)
        self.assertEqual({r['cohort_parents'] for r in sensitivity if r['cohort']=='baseline_exact_calibrated_length_panel'},{1})

    def test_zero_bound_not_given_after_control_success(self):
        from test_outcomes import MatchedOutcomeTests
        case=MatchedOutcomeTests();case.setUp()
        _,contrasts,_,_=build_outcomes(case.events,case.logical,case.statuses,default_policy(),2)
        self.assertTrue(all(r['zero_control_probability_upper95']=='not_estimated' for r in contrasts if r['measure']=='detection'))

    def test_study_binding_rejects_settings_drift_and_false_failure_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            study=create_example(Path(directory)/'study')
            value=read_json(study)
            value['search']['real_target_limit']=value['search']['control_target_limit']=500
            save_json(Path(directory)/'bound.json',value)
            bound=Path(directory)/'bound.json'
            self.assertEqual(bind_study(bound,search_settings(), 'fixed',None)['search']['evalue'],10)
            with self.assertRaisesRegex(ContractError,'disagrees'):
                bind_study(bound,search_settings(evalue=100),'fixed',None)
            value['search']['failure_policy']='retry_then_block_matched_parent'
            bound.write_text(json.dumps(value))
            with self.assertRaisesRegex(ContractError,'failure_policy'):
                bind_study(bound,search_settings(),'fixed',None)

    def test_summary_retries_and_keeps_failures(self):
        payload={'result':{'uids':[]}}
        with patch('database_residue_screen.metadata.urlopen',side_effect=[OSError('temporary'),io.BytesIO(json.dumps(payload).encode())]), patch('database_residue_screen.metadata.time.sleep') as sleep:
            result,failures=summary_with_retry(__import__('urllib.request',fromlist=['Request']).Request('https://example.test/'))
        self.assertEqual(result,payload)
        self.assertEqual(failures[0]['attempt'],1)
        sleep.assert_called_once_with(.4)
        record=dict(accessionversion='TEST0001.1',slen=100,biomol='genomic',moltype='dna',taxid=81077)
        self.assertEqual(summary_record(record,{'TEST0001.1'},'2026-10-20')['origin'],'artificial_other')


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),'BLAST unavailable')
class LengthEngineTests(unittest.TestCase):
    def test_runtime_statistics_and_minimum_word_sizes(self):
        observed,witnesses=engine_statistics(shutil.which('blastn'),[4,5,6,7,8,10,16,19,20,21,22,111])
        for key,value in STATS.items():
            self.assertAlmostEqual(observed[key],value)
        self.assertTrue(all(r['exact_witness_recovered'] for r in witnesses))

    def test_censoring_survives_analysis_release_and_sensitivity_export(self):
        from database_residue_screen.contracts import schema, write_table, read_table
        from database_residue_screen.controls import generate_controls
        from database_residue_screen.panels import prepare_panel
        from database_residue_screen.reanalysis import analyse_search
        from database_residue_screen.release import validate_analysis
        from database_residue_screen.outcomes import OUTCOME_SCHEMA
        from database_residue_screen.sensitivity import SENSITIVITY_SCHEMA
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            save_fasta(root/'q.fa',{'short':'ACGTCAGT'})
            row={k:'' for k in schema('queries')['properties']}
            row.update(query_id='short',group='Custom',condition='bare',construct_class='barcode',family_id='f',role='exploratory',
                       provenance_status='unverified',architecture_status='not_applicable')
            write_table(root/'q.tsv',[row],schema('queries'))
            generate_controls(root/'q.fa',root/'controls',2)
            prepare_panel(root/'q.fa',root/'q.tsv',root/'controls',root/'panel',2)
            save_fasta(root/'s.fa',{'SHORT001.1':'N'*30+'ACGTCAGT'+'N'*30})
            subprocess.run(['makeblastdb','-in',str(root/'s.fa'),'-out',str(root/'db'),'-dbtype','nucl','-parse_seqids'],check=True,capture_output=True)
            inspect_database(root/'db',root/'evidence')
            settings=search_settings(search_space=score_calibration(8,BASES,SEQUENCES,STATS)['effective_search_space'])
            with self.assertWarns(RuntimeWarning):
                run_search(root/'panel/search_queries.fasta',root/'db',root/'evidence/database_identity.json',root/'search',settings,
                           threads=1,mode='working',threshold_mode='fixed')
            analyse_search(root/'panel',root/'search',root/'analysis')
            self.assertTrue(validate_analysis(root/'analysis')['valid'])
            rows=read_table(root/'analysis/outcomes.tsv',OUTCOME_SCHEMA)
            self.assertTrue(all(row['resolution']=='censored_untestable' and row['positive']==-1 for row in rows))
            rows=read_table(root/'analysis/length_sensitivity.tsv',SENSITIVITY_SCHEMA)
            self.assertTrue(all(row['denominator']==0 and row['rate']=='not_estimated' for row in rows))

    def test_real_search_recovers_short_exact_and_records_fixed_censoring(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            save_fasta(root/'subjects.fa',{'SHORT001.1':'N'*30+'ACGTCAGT'+'N'*30})
            save_fasta(root/'q.fa',{'short':'ACGTCAGT'})
            subprocess.run(['makeblastdb','-in',str(root/'subjects.fa'),'-out',str(root/'db'),'-dbtype','nucl','-parse_seqids'],capture_output=True,check=True)
            inspect_database(root/'db',root/'evidence')
            settings=search_settings(search_space=score_calibration(8,BASES,SEQUENCES,STATS)['effective_search_space'])
            with self.assertWarns(RuntimeWarning):
                old=run_search(root/'q.fa',root/'db',root/'evidence/database_identity.json',root/'old',settings,threads=1,mode='working',threshold_mode='fixed')
            new=run_search(root/'q.fa',root/'db',root/'evidence/database_identity.json',root/'new',settings,threads=1,mode='working',threshold_mode='adaptive')
            self.assertEqual((old['success_zero'],new['success_zero']),(1,0))
            self.assertFalse(read_json(root/'old/request.json')['length_calibration']['lengths'][0]['strict_testable'])
            self.assertTrue(read_json(root/'new/request.json')['length_calibration']['lengths'][0]['strict_testable'])
