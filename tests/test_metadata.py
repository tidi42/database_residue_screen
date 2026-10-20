import io
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from database_residue_screen.contracts import ContractError, read_json, read_table, write_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.database import inspect_database
from database_residue_screen.design import digest
from database_residue_screen.example import create_example
from database_residue_screen.metadata import annotate_events, fetch_metadata, summary_record
from database_residue_screen.outcomes import CONTRAST_SCHEMA, OUTCOME_SCHEMA
from database_residue_screen.panels import prepare_panel
from database_residue_screen.reanalysis import analyse_search, read_events
from database_residue_screen.release import export_figure_inputs, validate_analysis
from database_residue_screen.search import run_search
from database_residue_screen.sequences import read_fasta
from database_residue_screen.taxonomy import COMPOSITION_SCHEMA, annotate_taxonomy, pin_taxonomy
from test_taxonomy import make_taxdump


class MetadataEvidenceTests(unittest.TestCase):
    def test_exact_version_and_unknown_topology_are_not_inferred(self):
        record = dict(accessionversion="TEST0001.1", slen=100, biomol="genomic", taxid=123)
        result = summary_record(record, {"TEST0001.1"}, "2026-09-15")
        self.assertEqual(result["topology"], "unknown")
        self.assertEqual(result["record_type"], "genomic")
        self.assertEqual(result["molecule"], "unknown")
        self.assertEqual(summary_record(dict(record, moltype="dna"), {"TEST0001.1"}, "2026-09-15")["molecule"], "DNA")
        self.assertIsNone(summary_record(dict(record, accessionversion="TEST0001.2"), {"TEST0001.1"}, "2026-09-15"))

    def test_transcripts_and_synthetic_records_are_separate(self):
        record = dict(accessionversion="TEST0001.1", slen=100, biomol="mRNA", topology="linear", taxid=123)
        self.assertEqual(summary_record(record, {"TEST0001.1"}, "2026-09-15")["record_type"], "transcript")
        self.assertEqual(summary_record(dict(record, taxid=32630), {"TEST0001.1"}, "2026-09-15")["record_type"], "synthetic")


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST+ unavailable")
class AnnotatedSearchIntegrationTests(unittest.TestCase):
    def test_annotation_recalculates_matched_outcomes_and_preserves_broad_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = create_example(root / "study").parent / "inputs"
            generate_controls(inputs / "queries.fasta", root / "controls", 2)
            prepare_panel(inputs / "queries.fasta", inputs / "queries.tsv", root / "controls", root / "panel", 2)
            prefix = root / "temporary_database"
            subprocess.run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(inputs / "validation_subjects.fasta"), "-out", str(prefix)],
                           capture_output=True, check=True)
            inspect_database(prefix, root / "database_evidence")
            run_search(root / "panel/search_queries.fasta", prefix, root / "database_evidence/database_identity.json", root / "search", threads=1, mode="working")
            source = analyse_search(root / "panel", root / "search", root / "analysis")
            run_search(root / "panel/search_queries.fasta", prefix, root / "database_evidence/database_identity.json", root / "reused_search",
                       threads=2, mode="working", reuse_search=root / "search")
            analyse_search(root / "panel", root / "reused_search", root / "reused_analysis")
            self.assertEqual(read_table(root / "analysis/outcomes.tsv", OUTCOME_SCHEMA),
                             read_table(root / "reused_analysis/outcomes.tsv", OUTCOME_SCHEMA))
            self.assertEqual(read_table(root / "analysis/contrasts.tsv", CONTRAST_SCHEMA),
                             read_table(root / "reused_analysis/contrasts.tsv", CONTRAST_SCHEMA))
            self.assertTrue(validate_analysis(root / "reused_analysis")["valid"])
            events = read_events(root / "analysis", source)
            requested = sorted({row["accession_version"] for row in events})
            sequences = read_fasta(inputs / "validation_subjects.fasta")
            payload = dict(result=dict(uids=requested, **{name: dict(accessionversion=name, slen=len(sequences[name]), moltype="dna",
                                                                    biomol="genomic", topology="linear", taxid=2) for name in requested}))
            with patch("database_residue_screen.metadata.urlopen", return_value=io.BytesIO(json.dumps(payload).encode())):
                fetch_metadata(root / "analysis", root / "metadata")
            with patch("database_residue_screen.metadata.urlopen", side_effect=AssertionError("Cached exact versions must not be fetched again")):
                cached = fetch_metadata(root / "analysis", root / "cached_metadata", reuse_metadata=root / "metadata")
            self.assertEqual(cached["matched_exact_version_records"], len(requested))
            self.assertEqual(cached["reused_metadata"]["matching_records"], len(requested))
            # A hash-verified chained cache may be reused; a changed upstream response may not.
            with patch("database_residue_screen.metadata.urlopen", side_effect=AssertionError("Chained cache must not fetch")):
                chained = fetch_metadata(root / "analysis", root / "chained_metadata", reuse_metadata=root / "cached_metadata")
            self.assertEqual(chained["matched_exact_version_records"], len(requested))
            upstream = next((root / "metadata").glob("esummary_*.json"))
            original = upstream.read_bytes()
            upstream.write_bytes(original + b" ")
            with self.assertRaises(ContractError):
                fetch_metadata(root / "analysis", root / "tampered_chain", reuse_metadata=root / "cached_metadata")
            upstream.write_bytes(original)
            annotate_events(root / "analysis", root / "cached_metadata", root / "cached_annotated")
            self.assertTrue(validate_analysis(root / "cached_annotated")["valid"])
            result = annotate_events(root / "analysis", root / "metadata", root / "annotated")
            before = read_table(root / "analysis/outcomes.tsv", OUTCOME_SCHEMA)
            after = read_table(root / "annotated/outcomes.tsv", OUTCOME_SCHEMA)
            self.assertEqual([row for row in before if row["population"] == "all_records"],
                             [row for row in after if row["population"] == "all_records"])
            self.assertTrue(any(row["resolution"] == "unresolved" for row in before if row["population"] == "linear_genomic"))
            self.assertTrue(all(row["resolution"] in {"resolved_zero", "resolved_positive"} for row in after))
            contrasts = read_table(root / "annotated/contrasts.tsv", CONTRAST_SCHEMA)
            self.assertTrue(all(row["eligible_parents"] == row["planned_parents"] for row in contrasts))
            self.assertEqual(result["result_contract"], "reference_screen_v2")
            self.assertTrue(all(row["record_taxids"] == "2" for row in read_events(root / "annotated", result)))
            make_taxdump(root / "taxdump")
            pin_taxonomy(root / "taxdump", root / "taxonomy", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/", "2026-09-15T12:00:00Z")
            reconciled = annotate_taxonomy(root / "annotated", root / "taxonomy", root / "reconciled", ("superkingdom", "genus"))
            self.assertEqual(reconciled["taxonomy_status"], "reconciled_with_explicit_unresolved_states")
            composition = read_table(root / "reconciled/taxonomic_composition.tsv", COMPOSITION_SCHEMA)
            pooled = [row for row in composition if row["group"] == "Global union" and row["endpoint"] == "strict"
                      and row["population"] == "all_records" and row["rank"] == "superkingdom"]
            self.assertEqual(sum(row["accessions"] for row in pooled), len(requested))
            self.assertEqual(next(row["accessions"] for row in pooled if row["taxon_id"] == "2"), len(requested))
            self.assertTrue(validate_analysis(root / "reconciled")["valid"])
            with self.assertRaisesRegex(ContractError, "accepted database search"):
                export_figure_inputs(root / "reconciled", root / "blocked_figures")
            self.assertFalse((root / "blocked_figures").exists())
            altered_scope = dict(reconciled, search_mode="production", database_acceptance=dict(full_integrity_passed=True))
            (root / "reconciled/analysis.json").write_text(json.dumps(altered_scope))
            with self.assertRaisesRegex(ContractError, "provenance disagrees"):
                validate_analysis(root / "reconciled", require_production=True)
            (root / "reconciled/analysis.json").write_text(json.dumps(reconciled))
            changed_payload = dict(result=dict(uids=requested, **{name: dict(accessionversion=name.rsplit(".", 1)[0] + ".2", slen=len(sequences[name]),
                                                                            moltype="dna", biomol="genomic", topology="linear", taxid=2) for name in requested}))
            with patch("database_residue_screen.metadata.urlopen", return_value=io.BytesIO(json.dumps(changed_payload).encode())):
                mismatch = fetch_metadata(root / "analysis", root / "version_mismatch")
            self.assertEqual({row["reason"] for row in mismatch["problems"]}, {"version_mismatch"})
            annotate_events(root / "analysis", root / "version_mismatch", root / "unresolved")
            unresolved_outcomes = read_table(root / "unresolved/outcomes.tsv", OUTCOME_SCHEMA)
            self.assertTrue(any(row["resolution"] == "unresolved" for row in unresolved_outcomes if row["population"] == "all_records"))
            changed_contrasts = read_table(root / "annotated/contrasts.tsv", CONTRAST_SCHEMA)
            changed_contrasts[0]["real_detection_rate"] = "0.123456"
            (root / "annotated/contrasts.tsv").unlink()
            write_table(root / "annotated/contrasts.tsv", changed_contrasts, CONTRAST_SCHEMA)
            changed_manifest = read_json(root / "annotated/analysis.json")
            changed_manifest["files"]["contrasts.tsv"] = digest(root / "annotated/contrasts.tsv")
            (root / "annotated/analysis.json").write_text(json.dumps(changed_manifest))
            with self.assertRaisesRegex(ContractError, "contrasts.tsv disagrees"):
                validate_analysis(root / "annotated")


if __name__ == "__main__":
    unittest.main()