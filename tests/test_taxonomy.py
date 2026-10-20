import tempfile
import unittest
from pathlib import Path

from database_residue_screen.contracts import ContractError
from database_residue_screen.policy import default_policy
from database_residue_screen.taxonomy import Taxonomy, load_taxonomy, pin_taxonomy, reconcile_taxids, taxonomy_composition, taxonomy_rows


def make_taxdump(directory: Path) -> None:
    directory.mkdir()
    nodes = {1: (1, "no rank"), 2: (1, "superkingdom"), 10: (2, "genus"), 11: (10, "species"), 12: (10, "species"),
             20: (2, "genus"), 21: (20, "species"), 30: (2, "family"), 31: (30, "species"),
             10239: (1, "acellular root"), 50: (10239, "species"), 32630: (1, "no rank"), 60: (32630, "species")}
    names = {identifier: "taxon_" + str(identifier) for identifier in nodes}
    (directory / "nodes.dmp").write_text("".join(f"{identifier}\t|\t{parent}\t|\t{rank}\t|\n" for identifier, (parent, rank) in nodes.items()))
    (directory / "names.dmp").write_text("".join(f"{identifier}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for identifier, name in names.items()))
    (directory / "merged.dmp").write_text("100\t|\t101\t|\n101\t|\t11\t|\n")
    (directory / "delnodes.dmp").write_text("200\t|\n")


class TaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        make_taxdump(self.root / "source")
        self.taxonomy = Taxonomy.from_directory(self.root / "source")

    def test_active_merged_deleted_missing_and_exact_ranks(self):
        self.assertEqual(self.taxonomy.resolve(100), (11, "merged"))
        self.assertEqual(self.taxonomy.assignment("100", "genus")["taxon_id"], "10")
        self.assertEqual(self.taxonomy.assignment("200", "genus")["status"], "deleted")
        self.assertEqual(self.taxonomy.assignment("999", "genus")["status"], "unknown_id")
        self.assertEqual(self.taxonomy.assignment("0", "genus")["status"], "missing")
        self.assertEqual(self.taxonomy.assignment("31", "genus")["status"], "rank_unresolved")
        self.assertEqual(self.taxonomy.assignment("11", "domain")["status"], "rank_unresolved")
        self.assertEqual(self.taxonomy.assignment("11", "superkingdom")["taxon_id"], "2")

    def test_multiple_ids_require_unanimous_rank_assignment(self):
        self.assertEqual(self.taxonomy.assignment("11;12", "genus")["status"], "multiple_consensus")
        self.assertEqual(self.taxonomy.assignment("11;21", "genus")["status"], "multiple_ambiguous")
        self.assertEqual(self.taxonomy.assignment("11;200", "genus")["status"], "multiple_ambiguous")
        self.assertEqual(self.taxonomy.assignment("11;31", "genus")["status"], "rank_unresolved")

    def test_sources_are_reconciled_not_overwritten(self):
        self.assertEqual(reconcile_taxids(self.taxonomy, "100", "11", "genus")["source_status"], "agree")
        conflict = reconcile_taxids(self.taxonomy, "11", "12", "genus")
        self.assertEqual(conflict["status"], "conflicting")
        self.assertEqual(conflict["taxon_id"], "")
        self.assertEqual(reconcile_taxids(self.taxonomy, "", "11", "genus")["source_status"], "record_only")

    def test_artificial_sequences_are_not_all_vectors(self):
        source = self.root / "source"
        with (source / "nodes.dmp").open("a") as handle:
            for identifier, parent in ((81077,1),(29278,81077),(90,29278),(91,81077)):
                handle.write(f"{identifier}\t|\t{parent}\t|\tno rank\t|\n")
        with (source / "names.dmp").open("a") as handle:
            for identifier in (81077,29278,90,91):
                handle.write(f"{identifier}\t|\ttaxon_{identifier}\t|\t\t|\tscientific name\t|\n")
        taxonomy = Taxonomy.from_directory(source)
        self.assertEqual(taxonomy.context(90), "vector")
        self.assertEqual(taxonomy.context(91), "artificial_other")
        self.assertEqual(taxonomy.context(81077), "artificial_other")

    def test_viruses_and_synthetic_nodes_are_not_forced_into_cellular_domains(self):
        self.assertEqual(self.taxonomy.assignment("50", "superkingdom")["context"], "viral")
        self.assertEqual(self.taxonomy.assignment("50", "superkingdom")["taxon_id"], "")
        self.assertEqual(self.taxonomy.assignment("60", "genus")["context"], "synthetic")

    def test_pinned_hashes_and_separate_release_retrieval_dates(self):
        report = pin_taxonomy(self.root / "source", self.root / "snapshot", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/", "2026-09-15T12:00:00Z")
        self.assertIsNone(report["release_date"])
        self.assertEqual(load_taxonomy(self.root / "snapshot")[0].resolve(100), (11, "merged"))
        (self.root / "snapshot/merged.dmp").write_text("100\t|\t21\t|\n")
        with self.assertRaises(ContractError):
            load_taxonomy(self.root / "snapshot")

    def test_merge_cycles_and_missing_parents_are_rejected(self):
        (self.root / "source/merged.dmp").write_text("100\t|\t101\t|\n101\t|\t100\t|\n")
        with self.assertRaises(ContractError):
            Taxonomy.from_directory(self.root / "source")
        (self.root / "source/merged.dmp").write_text("")
        (self.root / "source/nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t999\t|\tgenus\t|\n")
        with self.assertRaises(ContractError):
            Taxonomy.from_directory(self.root / "source")

    def test_pooled_denominator_is_an_accession_union_including_unresolved_ranks(self):
        common = dict(strict="yes", tolerant="yes", molecule="DNA", record_type="genomic", topology="linear", identity_status="exact_version")
        events = [dict(common, accession_version="TEST001.1", group="A", blast_taxids="11", record_taxids="11"),
                  dict(common, accession_version="TEST002.1", group="A", blast_taxids="12", record_taxids="12"),
                  dict(common, accession_version="TEST002.1", group="B", blast_taxids="12", record_taxids="12"),
                  dict(common, accession_version="TEST003.1", group="B", blast_taxids="50", record_taxids="50")]
        assignments = taxonomy_rows(events, self.taxonomy, ("genus",), "0" * 64)
        summary = taxonomy_composition(events, assignments, ["A", "B"], default_policy(), ("genus",))
        repeated = taxonomy_composition(events + events, assignments, ["A", "B"], default_policy(), ("genus",))
        self.assertEqual(summary, repeated)
        pooled = [row for row in summary if row["group"] == "Global union" and row["endpoint"] == "strict" and row["population"] == "all_records"]
        self.assertEqual({row["denominator"] for row in pooled}, {3})
        self.assertEqual(sum(row["accessions"] for row in pooled), 3)
        self.assertEqual(next(row["accessions"] for row in pooled if row["taxon_id"] == "10"), 2)
        self.assertEqual(next(row["accessions"] for row in pooled if row["category"] == "viral:rank_unresolved"), 1)

    def test_conflicting_hsp_taxids_for_one_accession_are_not_merged_into_consensus(self):
        rows = taxonomy_rows([dict(accession_version="TEST001.1", blast_taxids="11", record_taxids=""),
                              dict(accession_version="TEST001.1", blast_taxids="12", record_taxids="")], self.taxonomy, ("genus",), "0" * 64)
        self.assertEqual(rows[0]["source_status"], "within_source_conflict")
        self.assertEqual(rows[0]["taxon_id"], "")


if __name__ == "__main__":
    unittest.main()