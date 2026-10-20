import shutil
import subprocess
import tempfile
import unittest
import warnings
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_json, read_table, schema, write_table
from database_residue_screen.design import digest
from database_residue_screen.example import save_fasta, save_json
from database_residue_screen.outcomes import CONTRAST_SCHEMA, OUTCOME_SCHEMA
from database_residue_screen.policy import default_policy
from database_residue_screen.sequences import canonical, read_fasta

LONG = "ACGTTGCAAGGCTTACCGATAGCTTGACCATG"
NEW = "TTGACCGATCAGGCATTCAGGACTTAGCCATGCAGTTAGC"
ROLES = {"bare": dict(panel_version="v1", role="primary"),
         "old_context": dict(panel_version="v1", role="superseded", superseded_by=["new_context"]),
         "new_context": dict(panel_version="v2", role="primary")}


def strip_calibration(search: Path) -> None:
    """Make a runtime-calibrated working search look like an 8db62fa-era fixed-threshold search."""
    from database_residue_screen.search import object_digest
    request = read_json(search / "request.json")
    for key in ("length_calibration", "calibration_code_sha256", "failure_policy", "study_binding", "database_runtime_metadata"):
        request.pop(key, None)
    (search / "request.json").unlink()
    save_json(search / "request.json", request)
    complete = read_json(search / "complete.json")
    complete["search_id"] = object_digest(request)
    for batch in complete["batches"]:
        path = search / batch["directory"] / "complete.json"
        record = dict(read_json(path), search_id=complete["search_id"])
        path.unlink()
        save_json(path, record)
        batch["completion_sha256"] = digest(path)
    (search / "complete.json").unlink()
    save_json(search / "complete.json", complete)


@unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST unavailable")
class BundleAssemblyTests(unittest.TestCase):
    def metadata(self, path, entries):
        rows = []
        for name, condition in entries:
            row = {key: "" for key in schema("queries")["properties"]}
            row.update(query_id=name, group="Custom", condition=condition, construct_class="barcode", family_id=name, role="exploratory",
                       provenance_status="unverified", architecture_status="not_applicable")
            rows.append(row)
        write_table(path, rows, schema("queries"))

    def search(self, root, queries, output, **settings):
        from database_residue_screen.search import run_search, search_settings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            run_search(queries, root / "db", root / "evidence/database_identity.json", output, search_settings(**settings) if settings else None,
                       threads=1, mode="working")

    def test_bundle_combines_migrated_and_new_searches_with_roles(self):
        from database_residue_screen.bundle import assemble_bundle
        from database_residue_screen.controls import generate_controls
        from database_residue_screen.database import inspect_database
        from database_residue_screen.panels import prepare_panel
        from database_residue_screen.reanalysis import analyse_search
        from database_residue_screen.release import validate_analysis
        from database_residue_screen.testability_migration import calibrate_existing_search
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            save_fasta(root / "s.fa", {"LONG0001.1": "N" * 30 + LONG + "N" * 30 + NEW + "N" * 30})
            subprocess.run(["makeblastdb", "-in", str(root / "s.fa"), "-out", str(root / "db"), "-dbtype", "nucl", "-parse_seqids"],
                           check=True, capture_output=True)
            inspect_database(root / "db", root / "evidence")
            # historical panel and fixed-threshold search without runtime calibration
            save_fasta(root / "h.fa", {"short": "ACGTCAGT", "long": LONG, "old": LONG[:20] + "GATTACAGATTACA"})
            self.metadata(root / "h.tsv", [("short", "bare"), ("long", "bare"), ("old", "old_context")])
            generate_controls(root / "h.fa", root / "h_controls", 2)
            prepare_panel(root / "h.fa", root / "h.tsv", root / "h_controls", root / "h_panel", 2)
            self.search(root, root / "h_panel/search_queries.fasta", root / "historical")
            strip_calibration(root / "historical")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                calibrate_existing_search(root / "historical", root / "h_panel/search_queries.fasta", root / "migration")
            # bundle panel: historical queries plus the corrected construct, with declared roles
            union = dict(read_fasta(root / "h.fa"), new=NEW)
            save_fasta(root / "u.fa", union)
            self.metadata(root / "u.tsv", [("short", "bare"), ("long", "bare"), ("old", "old_context"), ("new", "new_context")])
            save_json(root / "policy.json", dict(default_policy(), condition_roles=ROLES))
            generate_controls(root / "u.fa", root / "u_controls", 2)
            prepare_panel(root / "u.fa", root / "u.tsv", root / "u_controls", root / "panel", 2, policy_file=root / "policy.json")
            historical_ids = set(read_fasta(root / "h_panel/search_queries.fasta"))
            bundle_queries = read_fasta(root / "panel/search_queries.fasta")
            self.assertTrue(historical_ids < set(bundle_queries))
            bridge = [name for name, sequence in bundle_queries.items() if sequence == canonical(LONG)]
            self.assertEqual(len(bridge), 1)
            new = {name: sequence for name, sequence in bundle_queries.items() if name not in historical_ids or name in bridge}
            save_fasta(root / "new.fa", new)
            self.search(root, root / "new.fa", root / "new_search")
            report = assemble_bundle(root / "panel", root / "historical", root / "migration/calibration_migration.json",
                                     [root / "new_search"], root / "bundle")
            self.assertEqual(report["completed_queries"], len(bundle_queries))
            request = read_json(root / "bundle/request.json")
            self.assertEqual(request["assembly_contract"], "verified_bundle_projection_v1")
            per_length = request["calibration_provenance"]["per_length"]
            self.assertEqual(per_length["8"], ["historical_migration"])
            self.assertEqual(per_length[str(len(LONG))][0], "historical_migration")
            self.assertEqual(len(per_length[str(len(LONG))]), 2)
            self.assertEqual(len(per_length[str(len(NEW))]), 1)
            analyse_search(root / "panel", root / "bundle", root / "analysis")
            self.assertTrue(validate_analysis(root / "analysis")["valid"])
            self.assertEqual(read_json(root / "analysis/analysis.json")["length_calibration_source"]["kind"], "bundle_assembly")
            outcomes = read_table(root / "analysis/outcomes.tsv", OUTCOME_SCHEMA)
            self.assertEqual({row["condition_role"] for row in outcomes}, {"primary", "superseded"})
            contrasts = read_table(root / "analysis/contrasts.tsv", CONTRAST_SCHEMA)
            union_row = next(row for row in contrasts if row["group"] == "Global union" and row["endpoint"] == "strict"
                             and row["population"] == "all_records" and row["measure"] == "detection")
            self.assertEqual(union_row["planned_parents"], 3)
            self.assertEqual(union_row["panel_version"], "v1;v2")
            # a source with different scientific settings is refused
            self.search(root, root / "new.fa", root / "other_settings", target_limit=10)
            with self.assertRaises(ContractError):
                assemble_bundle(root / "panel", root / "historical", root / "migration/calibration_migration.json",
                                [root / "other_settings"], root / "refused")


if __name__ == "__main__":
    unittest.main()
