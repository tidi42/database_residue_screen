import itertools
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from database_residue_screen.contracts import ContractError, read_table
from database_residue_screen.controls import DRAW_SCHEMA, check_draw, dinucleotides, generate_controls, shuffle_draw
from database_residue_screen.example import save_fasta
from database_residue_screen.sequences import read_fasta, reverse_complement


class ControlGenerationTests(unittest.TestCase):
    def test_reseeding_and_all_invariants(self):
        for sequence in ("A", "AC", "A" * 40, "ACGTACACGGTTACGATCGATCCG", "AACACC"):
            initial = shuffle_draw(sequence, 12345)
            shuffle_draw(sequence, 987)
            self.assertEqual(shuffle_draw(sequence, 12345), initial)
            check_draw(sequence, initial)

    def test_invalid_control_is_rejected(self):
        with self.assertRaises(ContractError):
            check_draw("AACCGGTT", "ACGTACGT")
        for sequence in ("", "AN", "acgt"):
            with self.assertRaises(ContractError):
                shuffle_draw(sequence, 42)

    def test_enumerable_support_and_frequency_sanity(self):
        parent = "AACACC"
        support = {"".join(parts) for parts in set(itertools.permutations(parent))
                   if parts[0] == parent[0] and parts[-1] == parent[-1]
                   and dinucleotides("".join(parts)) == dinucleotides(parent)}
        self.assertGreater(len(support), 1)
        observed = Counter(shuffle_draw(parent, seed) for seed in range(1, 6001))
        self.assertEqual(set(observed), support)
        expected = 6000 / len(support)
        self.assertTrue(all(abs(count - expected) < 0.15 * expected for count in observed.values()), observed)

    def test_order_aliases_and_draw_extension_do_not_change_draws(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "queries.fasta"
            sequence = "AACACCTAGGTCACGT"
            save_fasta(source, {"first": sequence, "reverse": reverse_complement(sequence), "second": "T" * 12})
            short = generate_controls(source, root / "short", 3)
            extended = generate_controls(source, root / "extended", 7)
            reversed_source = root / "reversed.fasta"
            save_fasta(reversed_source, dict(reversed(list(read_fasta(source).items()))))
            generate_controls(reversed_source, root / "reordered", 3)
            first = read_fasta(root / "short" / "controls.fasta")
            longer = read_fasta(root / "extended" / "controls.fasta")
            self.assertEqual(first, {name: longer[name] for name in first})
            self.assertEqual(first, read_fasta(root / "reordered" / "controls.fasta"))
            self.assertEqual(short["canonical_parents"], 2)
            self.assertEqual(extended["total_draws"], 14)
            rows = read_table(root / "short" / "draws.tsv", DRAW_SCHEMA)
            self.assertEqual(len({row["draw_id"] for row in rows}), 6)
            self.assertGreaterEqual(short["parent_identical_draws"], 3)


if __name__ == "__main__":
    unittest.main()