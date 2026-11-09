import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from database_residue_screen.acceptance import accept_database, probe_result, validate_acceptance, verify_download_evidence
from database_residue_screen.contracts import ContractError, read_json
from database_residue_screen.database import inspect_database, metadata_conflicts, parse_database_info
from database_residue_screen.design import digest
from database_residue_screen.example import save_fasta, save_json


class DownloadEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.evidence = Path(self.temporary.name) / "evidence"
        self.evidence.mkdir()
        self.runtime = dict(sequence_count=1, base_count=100, database_version=5, volume_count=1, build_date="2026-09-12")
        self.metadata = {"version": "1.1", "dbname": "core_nt", "number-of-sequences": 1, "number-of-letters": 100,
                         "number-of-volumes": 1, "last-updated": "2026-09-12T00:00:00",
                         "files": ["ftp://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.00.tar.gz"]}
        save_json(self.evidence / "DOWNLOAD_VERIFIED.json", dict(status="all_archives_checksum_verified", archives=1,
                                                                remote_release_metadata_unchanged=True))
        save_json(self.evidence / "last_exit.json", dict(exit_code=0))
        (self.evidence / "archives.txt").write_text("core_nt.00.tar.gz\n")
        (self.evidence / "archives.md5").write_text("0" * 32 + "  core_nt.00.tar.gz\n")
        (self.evidence / "download_core_nt_snapshot.sh").write_text("fixture launcher\n")
        (self.evidence / "archive_checksum_verification.log").write_text("core_nt.00.tar.gz: OK\n")
        self.pin_metadata()

    def pin_metadata(self):
        (self.evidence / "core_nt-nucl-metadata.json").write_text(json.dumps(self.metadata))
        (self.evidence / "metadata-after-download.json").write_text(json.dumps(self.metadata))
        names = ("core_nt-nucl-metadata.json", "archives.txt", "archives.md5", "download_core_nt_snapshot.sh")
        (self.evidence / "request.sha256").write_text("".join(digest(self.evidence / name) + "  evidence/" + name + "\n" for name in names))

    def test_download_manifest_may_omit_database_version_without_inventing_it(self):
        result = verify_download_evidence(self.evidence, self.runtime)
        self.assertEqual(result["release_metadata"], self.metadata)
        self.assertNotIn("db-version", result["release_metadata"])
        self.assertEqual(metadata_conflicts(self.runtime, self.metadata),
                         [dict(field="database_version", runtime=5, json_metadata=None)])

    def test_explicit_database_version_must_match(self):
        self.metadata["db-version"] = 5
        self.pin_metadata()
        self.assertEqual(verify_download_evidence(self.evidence, self.runtime)["release_metadata"]["db-version"], 5)
        for version in (4, None):
            with self.subTest(version=version):
                self.metadata["db-version"] = version
                self.pin_metadata()
                with self.assertRaisesRegex(ContractError, "release metadata"):
                    verify_download_evidence(self.evidence, self.runtime)

    def test_release_identity_fields_remain_required_and_must_match(self):
        original = dict(self.metadata)
        for field, changed in (("number-of-sequences", 2), ("number-of-letters", 101),
                               ("number-of-volumes", 2), ("last-updated", "2026-09-13T00:00:00")):
            for missing in (False, True):
                with self.subTest(field=field, missing=missing):
                    self.metadata = dict(original)
                    if missing:
                        del self.metadata[field]
                    else:
                        self.metadata[field] = changed
                    self.pin_metadata()
                    with self.assertRaisesRegex(ContractError, "release metadata"):
                        verify_download_evidence(self.evidence, self.runtime)


class AccessionProbeTests(unittest.TestCase):
    def test_wrong_accession_or_length_is_never_accepted(self):
        requested = dict(accession_version="TEST0001.1", expected_length=47, allow_absent=True, basis="fixture")
        self.assertFalse(probe_result(requested, subprocess.CompletedProcess([], 0, "OTHER0001.1\t47\n", ""))["passed"])
        self.assertFalse(probe_result(requested, subprocess.CompletedProcess([], 0, "TEST0001.1\t48\n", ""))["passed"])
        self.assertTrue(probe_result(requested, subprocess.CompletedProcess([], 0, "TEST0001.1\t47\n", ""))["passed"])
        self.assertEqual(probe_result(requested, subprocess.CompletedProcess([], 1, "", "Entry not found: TEST0001.1"))["status"], "absent_recorded")
        self.assertFalse(probe_result(requested, subprocess.CompletedProcess([], 1, "", "Database not found"))["passed"])


@unittest.skipUnless(all(shutil.which(name) for name in ("blastn", "makeblastdb", "blastdbcmd", "blastdbcheck")), "BLAST tools unavailable")
class DatabaseAcceptanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        sequence = "ACGTATCGATTCGAGCTGCTAGCTAGGCTTCGATCGGTATCGAATCG"
        save_fasta(self.root / "subjects.fasta", {"TEST0001.1": sequence})
        self.prefix = self.root / "temporary_db"
        subprocess.run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(self.root / "subjects.fasta"), "-out", str(self.prefix)], capture_output=True, check=True)
        runtime = parse_database_info(subprocess.run(["blastdbcmd", "-db", str(self.prefix), "-info"], capture_output=True, text=True, check=True).stdout)
        if not self.prefix.with_suffix(".njs").exists():
            save_json(self.prefix.with_suffix(".njs"), {"number-of-sequences": runtime["sequence_count"], "number-of-letters": runtime["base_count"],
                      "db-version": runtime["database_version"], "number-of-volumes": runtime["volume_count"], "last-updated": runtime["build_date"]})
        inspect_database(self.prefix, self.root / "inventory", hash_all=True)
        self.evidence = self.root / "inventory/database_identity.json"
        self.probes = self.root / "probes.json"
        save_json(self.probes, dict(version="1.0.0", probes=[dict(accession_version="TEST0001.1", expected_length=len(sequence), allow_absent=False, basis="Temporary fixture")]))
        self.support = self.root / "simulated_taxonomy_support"
        self.support.mkdir()
        for name in ("taxdb.bti", "taxdb.btd"):
            (self.support / name).write_bytes(b"Simulated support for the mocked integrity checker; not a real taxdb")

    def simulate_checker(self, returncode):
        actual_run = subprocess.run

        def run(command, **kwargs):
            if "-full" in command and Path(command[0]).name == "blastdbcheck":
                kwargs["stdout"].write("Simulated full-check result for gate tests only\n")
                return subprocess.CompletedProcess(command, returncode)
            return actual_run(command, **kwargs)

        return patch("database_residue_screen.acceptance.subprocess.run", side_effect=run)

    def test_acceptance_identity_guards_with_simulated_checker_and_real_probes(self):
        with self.simulate_checker(0):
            result = accept_database(self.prefix, self.evidence, self.probes, self.root / "accepted", taxonomy_support=self.support)
        self.assertTrue(result["accepted"])
        accepted = self.root / "accepted/accepted.json"
        checks = read_json(self.root / "accepted/checks.json")
        self.assertTrue(all("-target_only" in row["command"] for row in checks["accession_probes"]))
        self.assertTrue(validate_acceptance(accepted, self.prefix, self.evidence, "blastn")["full_integrity_passed"])
        self.prefix.with_suffix(".extra").write_bytes(b"unexpected")
        with self.assertRaisesRegex(ContractError, "membership"):
            validate_acceptance(accepted, self.prefix, self.evidence, "blastn")

    def test_failed_checker_never_writes_acceptance(self):
        with self.simulate_checker(1), self.assertRaisesRegex(ContractError, "integrity check"):
            accept_database(self.prefix, self.evidence, self.probes, self.root / "failed", taxonomy_support=self.support)
        self.assertFalse((self.root / "failed/accepted.json").exists())
        self.assertTrue((self.root / "failed/failed.json").exists())

    def test_unhashed_inventory_and_missing_support_are_rejected(self):
        with self.assertRaisesRegex(ContractError, "taxonomy support"):
            accept_database(self.prefix, self.evidence, self.probes, self.root / "no_support")
        inventory = read_json(self.evidence)
        inventory["all_database_files_hashed"] = False
        self.evidence.write_text(json.dumps(inventory))
        with self.assertRaisesRegex(ContractError, "complete content hashes"):
            accept_database(self.prefix, self.evidence, self.probes, self.root / "unhashed", taxonomy_support=self.support)

    def test_altered_acceptance_logs_and_executable_are_rejected(self):
        with self.simulate_checker(0):
            accept_database(self.prefix, self.evidence, self.probes, self.root / "accepted", taxonomy_support=self.support)
        accepted = self.root / "accepted/accepted.json"
        with self.assertRaisesRegex(ContractError, "executable"):
            validate_acceptance(accepted, self.prefix, self.evidence, "blastdbcmd")
        (accepted.parent / "integrity_stdout.txt").write_text("altered")
        with self.assertRaisesRegex(ContractError, "evidence changed"):
            validate_acceptance(accepted, self.prefix, self.evidence, "blastn")


if __name__ == "__main__":
    unittest.main()