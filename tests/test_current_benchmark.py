import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from database_residue_screen.benchmark_calibration import archive_search, compile_length_adapter, xml_statistics
from database_residue_screen.benchmark_engine import direct_search
from database_residue_screen.benchmark_truth import precompute_truth
from database_residue_screen.current_benchmark import analyse_arm, prepare_benchmark
from database_residue_screen.contracts import ContractError, write_table
from database_residue_screen.example import save_fasta, save_json
from database_residue_screen.panels import LOGICAL_SCHEMA
from database_residue_screen.search import search_settings
from database_residue_screen.sequences import canonical_id


class CurrentBenchmarkTests(unittest.TestCase):
    def test_source_adapter_rejects_unpinned_archive(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.tar.gz';p.write_bytes(b'not NCBI')
            with self.assertRaisesRegex(ContractError,'differs'):
                compile_length_adapter(p,Path(d)/'out')

    def test_prepare_accepts_foreign_panel_and_database_dimensions(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);panel=root/'panel';cal=root/'calibration';panel.mkdir();cal.mkdir()
            sequence='AGCCTGATCAGT';parent=canonical_id(sequence)
            save_fasta(panel/'search_queries.fasta',{parent:sequence})
            row={key:('' if prop['type']=='string' else 0) for key,prop in LOGICAL_SCHEMA['properties'].items()}
            row.update(role='real',logical_id='foreign',parent_id=parent,sequence_id=parent,length=12,group='Foreign',condition='custom')
            write_table(panel/'logical_queries.tsv',[row],LOGICAL_SCHEMA)
            save_json(panel/'panel.json',dict(files={}))
            save_json(cal/'calibration.json',dict(files={},status='engine_validated_version_matched',route='direct_per_query_XML_and_format7',
                       binaries={},search_spaces=[dict(query_length=12,database_bases=12345,database_sequences=7,effective_search_space=100,exact_evalue=.1)]))
            spec=prepare_benchmark(panel,cal,root/'benchmark')
            self.assertEqual(spec['original_parents'],1)
            self.assertEqual(spec['actual_lengths'],[12])
            self.assertEqual(spec['boundary_templates'],8)
            self.assertEqual(spec['production_database'],dict(bases=12345,sequences=7))


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'),'BLAST unavailable')
class BenchmarkEngineTests(unittest.TestCase):
    def test_direct_and_archive_xml_calibration_and_independent_aggregation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            queries={'q':'ACGTGCTATCGAGTGCATCG'}
            subjects={'FIX0001.1':'N'*30+queries['q']+'N'*30,'FIX0002.1':'N'*60}
            save_fasta(root/'q.fa',queries);save_fasta(root/'s.fa',subjects)
            subprocess.run(['makeblastdb','-in',str(root/'s.fa'),'-out',str(root/'db'),'-dbtype','nucl','-parse_seqids'],check=True,capture_output=True)
            binaries={name:{'path':shutil.which(name)} for name in ('blastn','blast_formatter')}
            direct=direct_search(root/'q.fa',queries,root/'db',root/'direct',binaries,search_settings())
            self.assertTrue(direct['xml_tabular_hsp_multisets_equal'])
            archive_search(root/'q.fa',root/'db',root/'archive',binaries,search_settings())
            self.assertEqual(xml_statistics(root/'direct/blast.xml'),xml_statistics(root/'archive/blast.xml'))
            truth,_=precompute_truth(queries,subjects)
            pairs=[dict(query_id='q',subject_id=s,strict=int(('q',s) in truth['strict']),tolerant=int(('q',s) in truth['tolerant'])) for s in subjects]
            rows=analyse_arm(root/'direct',queries,subjects,pairs,{},dict(q={'complexity':'other'}),dict(shard_id='tiny',panel='A'),'local_fixture')
            self.assertEqual([r['strict_classification'] for r in rows],['recovered','correct_nonpositive'])
