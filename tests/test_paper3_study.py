"""The Paper 3 study configurations resolve against the `screen` pipeline and reference only declared panel queries."""
import json
import re
import unittest
from pathlib import Path

from database_residue_screen.census import specification_references, validate_specification
from database_residue_screen.constructs import load_manifest, read_definition
from database_residue_screen.pipeline.screen import PIPELINE
from database_residue_screen.record_context import selection_specification
from database_residue_screen.study.config import load_config, read_config_file

STUDY = Path(__file__).resolve().parents[1] / "studies" / "paper3"
ENVIRONMENT = dict(CORE_NT_DB="/db/core_nt", CORE_NT_TAXDB="/db", CORE_NT_DOWNLOAD_EVIDENCE="/db/evidence", BLAST_BIN="/blast/bin",
                   BLAST_SOURCE_TARBALL="/blast/src.tar.gz", TAXDUMP_DIR="/taxdump", PANEL_SOURCES="/sources",
                   NCBI_EMAIL="someone@example.org")
DRY_RUN_KEYS = {("study", "id"), ("study", "title"), ("subset", "query_ids"), ("exact_scan", "volume_limit"),
                ("benchmark", "run_sizes"), ("benchmark", "shard_limit"), ("validation", "require_production")}


def flatten(value, prefix=()):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten(item, prefix + (key,))
    else:
        yield prefix, value


class Paper3StudyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.configs = {name: load_config(STUDY / name, PIPELINE.sections(), environ=ENVIRONMENT)
                       for name in ("study.toml", "study_dryrun.toml")}
        cls.entries = {entry["query_id"]: entry for entry in read_definition(STUDY / "panel_definition.tsv")}

    def test_both_configurations_resolve(self):
        for name, config in self.configs.items():
            with self.subTest(name=name):
                self.assertEqual(config.values["study"]["pipeline"], "screen")
                self.assertEqual(config.values["tools"]["blastn"]["command"], "/blast/bin/blastn")
                self.assertEqual(Path(config.values["panel"]["definition"]), STUDY / "panel_definition.tsv")

    def test_dry_run_differs_only_in_its_declared_keys(self):
        full, dry = (dict(flatten(read_config_file(STUDY / name).data)) for name in ("study.toml", "study_dryrun.toml"))
        self.assertEqual(set(full), set(dry))
        self.assertEqual({key for key in full if full[key] != dry[key]}, DRY_RUN_KEYS)
        self.assertEqual(self.configs["study.toml"].values["subset"], dict(parents=0, query_ids=[]))
        self.assertLessEqual(set(self.configs["study_dryrun.toml"].values["subset"]["query_ids"]), set(self.entries))
        families = {self.entries[name]["family"] for name in self.configs["study_dryrun.toml"].values["subset"]["query_ids"]}
        self.assertEqual(families, {entry["family"] for entry in self.entries.values()})

    def test_configurations_hold_no_sequences(self):
        for name in ("study.toml", "study_dryrun.toml", "sources.json", "README.md"):
            with self.subTest(name=name):
                self.assertIsNone(re.search(r"[ACGT]{12,}", (STUDY / name).read_text(encoding="utf-8")))

    def test_material_rules_name_panel_queries_and_sources(self):
        panel = self.configs["study.toml"].values["panel"]
        sources = set(load_manifest(STUDY / "sources.json"))
        names = [rule["query"] for rule in panel["census_elements"] + panel["benchmark_templates"]]
        names += [rule["strip_suffix"]["query"] for rule in panel["census_elements"] if "strip_suffix" in rule]
        self.assertLessEqual(set(names), set(self.entries))
        self.assertLessEqual({rule["source"] for rule in panel["census_catalogue"]}, sources)
        self.assertLessEqual({rule["source"] for rule in panel["census_catalogue_checks"] if "source" in rule}, sources)
        families = {entry["family"] for entry in self.entries.values()}
        self.assertLessEqual({rule["family"] for rule in panel["census_catalogue_checks"] if "family" in rule}, families)

    def test_census_families_use_declared_members_elements_and_catalogues(self):
        data = self.configs["study.toml"].values
        census = validate_specification(json.loads(json.dumps(data["census"])))
        used = specification_references(census)
        self.assertEqual(used["elements"], {rule["name"] for rule in data["panel"]["census_elements"]})
        self.assertEqual(used["catalogues"], {rule["set"] for rule in data["panel"]["census_catalogue"]})
        conditions = {entry["condition"] for entry in self.entries.values()}
        for family in census["families"]:
            self.assertLessEqual(set(family.get("conditions", ())), conditions)
            self.assertLessEqual(set(family.get("query_ids", ())), set(self.entries))

    def test_settings_shared_by_stages_agree(self):
        data = self.configs["study.toml"].values
        self.assertIn(data["search"]["target_limit"], data["target_limit"]["arms"])
        self.assertEqual(sorted(data["inference"]["retention_ladder"]), sorted(data["target_limit"]["arms"]))
        self.assertEqual(data["reports"]["draw_prefixes"], data["inference"]["draw_ladder"])
        self.assertEqual(max(data["inference"]["draw_ladder"]), data["controls"]["draws"])
        self.assertEqual(data["inference"]["maximum_draws"], data["controls"]["draws"])
        self.assertIn(data["analysis_policy"]["primary_window"], data["analysis_policy"]["windows"])
        self.assertEqual(data["record_cases"]["terminal_window"], data["analysis_policy"]["primary_window"])
        self.assertEqual(data["record_cases"]["terminal_minimum_length"], data["analysis_policy"]["primary_minimum_length"])
        selection_specification(data["record_cases"])


if __name__ == "__main__":
    unittest.main()
