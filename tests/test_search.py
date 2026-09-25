import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from database_residue_screen.contracts import ContractError, read_json, read_table
from database_residue_screen.database import inspect_database
from database_residue_screen.design import digest
from database_residue_screen.example import save_fasta
from database_residue_screen.search import STATUS_SCHEMA, iter_hsps, run_search, search_settings, validate_completed_batch


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST+ integration tools unavailable")
class SearchIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        sequence = "ACGTATCGATTCGAGCTGCTAGCTAGGCTTCGATCGGTATCGAATCG"
        save_fasta(self.root / "subjects.fasta", {"TEST0001.1": "GGG" + sequence + "TTT"})
        self.prefix = self.root / "database"
        subprocess.run(["makeblastdb", "-in", str(self.root / "subjects.fasta"), "-dbtype", "nucl", "-parse_seqids", "-out", str(self.prefix)],
                       capture_output=True, check=True)
        inspect_database(self.prefix, self.root / "evidence")
        self.evidence = self.root / "evidence" / "database_identity.json"
        self.queries = self.root / "queries.fasta"
        save_fasta(self.queries, {"positive": sequence, "zero": "C" * 40})

    def test_real_blast_completion_zero_resume_and_batching(self):
        first = run_search(self.queries, self.prefix, self.evidence, self.root / "together", threads=1, mode="working")
        resumed = run_search(self.queries, self.prefix, self.evidence, self.root / "together", threads=1, mode="working")
        self.assertEqual(first, resumed)
        self.assertEqual(first["success_zero"], 1)
        run_search(self.queries, self.prefix, self.evidence, self.root / "separate", threads=1, batch_size=1, mode="working")
        together = read_table(self.root / "together" / "status.tsv", STATUS_SCHEMA)
        separate = read_table(self.root / "separate" / "status.tsv", STATUS_SCHEMA)
        for row in together + separate:
            row.pop("batch_id")
        self.assertEqual(together, separate)
        self.assertEqual(together[0]["strict_targets"], 1)
        with self.assertRaises(ContractError):
            run_search(self.queries, self.prefix, self.evidence, self.root / "together", search_settings(2000), threads=1, mode="working")

    def test_truncated_or_appended_output_is_rejected(self):
        run_search(self.queries, self.prefix, self.evidence, self.root / "run", threads=1, mode="working")
        path = self.root / "run" / "batch_00001.attempt_00001" / "blast.tsv"
        original = path.read_text()
        from database_residue_screen.sequences import read_fasta

        for text in (original + original, "\n".join(original.splitlines()[:-1])):
            path.write_text(text)
            with self.assertRaises(ContractError):
                list(iter_hsps(path, read_fasta(self.queries)))
        with self.assertRaises(ContractError):
            run_search(self.queries, self.prefix, self.evidence, self.root / "run", threads=1, mode="working")

    def test_failed_query_is_not_zero_and_can_be_retried(self):
        real_run = subprocess.run

        def fail_batch(command, **kwargs):
            if "-query" in command:
                return subprocess.CompletedProcess(command, 2)
            return real_run(command, **kwargs)

        with patch("database_residue_screen.search.subprocess.run", side_effect=fail_batch), self.assertRaises(ContractError):
            run_search(self.queries, self.prefix, self.evidence, self.root / "retry", threads=1, mode="working")
        self.assertFalse((self.root / "retry" / "complete.json").exists())
        result = run_search(self.queries, self.prefix, self.evidence, self.root / "retry", threads=1, mode="working")
        self.assertTrue(result["all_queries_completed"])
        self.assertTrue((self.root / "retry" / "batch_00001.attempt_00001" / "failed.json").is_file())
        self.assertTrue((self.root / "retry" / "batch_00001.attempt_00002" / "complete.json").is_file())

    def test_final_report_interruption_is_resumable(self):
        expected = run_search(self.queries, self.prefix, self.evidence, self.root / "interrupted", threads=1, mode="working")
        (self.root / "interrupted" / "complete.json").rename(self.root / "interrupted" / "complete.pending.json")
        self.assertEqual(run_search(self.queries, self.prefix, self.evidence, self.root / "interrupted", threads=1, mode="working"), expected)

    def test_timeout_is_incomplete_not_zero(self):
        real_run = subprocess.run

        def timeout_batch(command, **kwargs):
            if "-query" in command:
                raise subprocess.TimeoutExpired(command, 0.01)
            return real_run(command, **kwargs)

        with patch("database_residue_screen.search.subprocess.run", side_effect=timeout_batch), self.assertRaises(ContractError):
            run_search(self.queries, self.prefix, self.evidence, self.root / "timeout", threads=1, timeout_seconds=0.01, mode="working")
        self.assertTrue((self.root / "timeout" / "batch_00001.attempt_00001" / "incomplete.json").is_file())
        self.assertFalse((self.root / "timeout" / "complete.json").exists())

    def test_production_requires_acceptance_before_any_query_execution(self):
        with self.assertRaisesRegex(ContractError, "requires --acceptance"):
            run_search(self.queries, self.prefix, self.evidence, self.root / "blocked")
        self.assertFalse((self.root / "blocked").exists())

    def test_thread_change_reuses_only_completed_batches_with_original_provenance(self):
        actual_run = subprocess.run
        submitted = []

        def interrupt_second_batch(command, **kwargs):
            if "-query" in command:
                submitted.append(command)
                if len(submitted) == 2:
                    return subprocess.CompletedProcess(command, 2)
            return actual_run(command, **kwargs)

        source = self.root / "original"
        with patch("database_residue_screen.search.subprocess.run", side_effect=interrupt_second_batch), self.assertRaises(ContractError):
            run_search(self.queries, self.prefix, self.evidence, source, threads=1, batch_size=1, mode="working")
        original = source / "batch_00001.attempt_00001/complete.json"
        original_checksum = digest(original)
        submitted.clear()

        def record_execution(command, **kwargs):
            if "-query" in command:
                submitted.append(command)
            return actual_run(command, **kwargs)

        destination = self.root / "continued"
        with patch("database_residue_screen.search.subprocess.run", side_effect=record_execution):
            result = run_search(self.queries, self.prefix, self.evidence, destination, threads=2, batch_size=1,
                                mode="working", reuse_search=source)
        self.assertTrue(result["all_queries_completed"])
        self.assertEqual(len(submitted), 1)
        self.assertEqual(submitted[0][submitted[0].index("-num_threads") + 1], "2")
        self.assertEqual(digest(original), original_checksum)
        request = read_json(destination / "request.json")
        reused = validate_completed_batch(destination / "batch_00001.attempt_00001", request)
        self.assertEqual(reused["reused_from"]["threads"], 1)
        self.assertEqual(reused["reused_from"]["completion_sha256"], original_checksum)
        self.assertEqual(run_search(self.queries, self.prefix, self.evidence, destination, threads=2, batch_size=1,
                                    mode="working", reuse_search=source), result)
        with self.assertRaisesRegex(ContractError, "identical queries"):
            run_search(self.queries, self.prefix, self.evidence, self.root / "wrong_cap", search_settings(2000),
                       threads=2, batch_size=1, mode="working", reuse_search=source)
        (destination / "batch_00001.attempt_00001/source_complete.json").write_text('{}')
        with self.assertRaises((ContractError, KeyError)):
            validate_completed_batch(destination / "batch_00001.attempt_00001", request)


if __name__ == "__main__":
    unittest.main()