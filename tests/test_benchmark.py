import random
import unittest

from database_residue_screen.benchmark_truth import independent_hsp_truth, precompute_truth, sequence_witnesses
from database_residue_screen.benchmark_fixtures import make_panels, template_sequence, variants
from database_residue_screen.matching import BLAST_FIELDS, parse_blast_hsp, strict_hsp, tolerant_hsp
from database_residue_screen.sequences import reverse_complement


def exact_reference(query, subject):
    return query in subject or reverse_complement(query) in subject


def gapfree_reference(query, subject):
    return sequence_witnesses(query, subject)["tolerant"] is not None


class IndependentReferenceTests(unittest.TestCase):
    def test_exact_both_strands(self):
        query = "ACGGTCAAGCTATCCGACTT"
        for subject, expected in (("N" + query + "N", True), ("N" + reverse_complement(query) + "N", True), ("N" + query[:-1] + "N", False)):
            self.assertEqual(sequence_witnesses(query, subject)["strict"] is not None, expected)

    def test_gapfree_boundaries_and_substitution(self):
        query = "ACGGTCAAGCTATCCGACTT"
        self.assertTrue(gapfree_reference(query, query[:18]))
        self.assertFalse(gapfree_reference(query, query[:17]))
        self.assertTrue(gapfree_reference(query, query[:10] + "A" + query[11:]))
        self.assertFalse(gapfree_reference(query, "T" * 20))

    def test_truth_contains_undeclared_cross_pairs_and_is_immutable(self):
        queries = {"first": "ACGGTCAAGCTATCCGACTT", "alias": "ACGGTCAAGCTATCCGACTT"}
        subjects = {"declared": "TT" + queries["first"] + "CC", "extra": queries["first"][:18]}
        truth, witnesses = precompute_truth(queries, subjects)
        self.assertIn(("alias", "extra"), truth["tolerant"])
        self.assertNotIn(("alias", "extra"), truth["strict"])
        self.assertIsInstance(truth["tolerant"], frozenset)
        self.assertEqual((truth, witnesses), precompute_truth(dict(reversed(list(queries.items()))), subjects))

    def test_diagonal_oracle_matches_naive_enumeration(self):
        generator = random.Random(20260922)
        for length in (8, 10, 16, 19, 20, 21, 24, 40, 61):
            for replicate in range(10):
                query = "".join(generator.choices("ACGT", k=length))
                subject = "".join(generator.choices("ACGT", k=length + 4)) if replicate % 2 else "AA" + query + "CC"
                if replicate % 3 == 0:
                    middle = len(subject) // 2
                    subject = subject[:middle] + next(base for base in "ACGT" if base != subject[middle]) + subject[middle + 1:]
                expected = any(20 * sum(first != second for first, second in zip(oriented[query_offset:query_offset + span], subject[subject_offset:subject_offset + span])) <= span
                               for oriented in (query, reverse_complement(query))
                               for span in range((9 * length + 9) // 10, length + 1)
                               for query_offset in range(length - span + 1)
                               for subject_offset in range(len(subject) - span + 1))
                actual = sequence_witnesses(query, subject)
                self.assertEqual(actual["tolerant"] is not None, expected, (length, replicate))
                self.assertEqual(actual["strict"] is not None, exact_reference(query, subject))

    def test_panel_fixtures_are_bounded_and_deterministic(self):
        originals = {"parent_a": "ACGTACGT", "parent_b": "AGGTAACGTT"}
        first = make_panels(originals)
        self.assertEqual(first, make_panels(originals))
        self.assertNotEqual(first, make_panels(originals, seed=7))
        self.assertEqual(sum(len(row["queries"]) for row in first[0] if row["panel"] == "A"), 2)
        self.assertEqual(sum(len(row["queries"]) for row in first[0] if row["panel"] == "B"), 16)
        self.assertEqual(sum(len(row["queries"]) for row in first[0] if row["panel"] == "C"), 2)
        self.assertFalse(any(row["panel"] == "D" for row in first[0]))
        # Extra lengths only add panel D templates; the other fixtures are unchanged.
        extended = make_panels(originals, (22,))
        self.assertEqual([row for row in extended[0] if row["panel"] != "D"], first[0])
        self.assertEqual([row["shard_id"] for row in extended[0] if row["panel"] == "D"], ["D_022"])
        self.assertTrue(all(len(row["subjects"]) < 2000 for row in first[0]))
        self.assertEqual(len(variants("ACGGTCAAGCTATCCGACTT")), 20)
        for shard in first[0]:
            for case in shard["cases"]:
                if case["error_class"].startswith("negative_"):
                    self.assertIsNone(sequence_witnesses(shard["queries"][case["query_id"]], shard["subjects"][case["subject_id"]])["tolerant"])

    def test_declared_templates_appear_at_every_extra_length(self):
        originals = {"parent_a": "ACGTACGT"}
        template = dict(name="tmpl", sequence="ACGTTGCAAGGCTTACCG", repeat_from=5)
        shards, _ = make_panels(originals, (19, 22), templates=[template])
        for length in (19, 22):
            shard = next(row for row in shards if row["shard_id"] == f"D_{length:03d}")
            self.assertEqual(shard["queries"][f"extra_{length:03d}_tmpl"], template_sequence(template, length))
            self.assertEqual(len(shard["queries"][f"extra_{length:03d}_tmpl"]), length)

    def test_direct_predicates_agree_even_when_short_hsps_are_not_retrieved(self):
        query = "ACGGTCAAGCTATCCGACTT"
        mismatch = query[:10] + "A" + query[11:]
        gapped_subject = query[:10] + "A" + query[10:]
        examples = [
            (query[:8], query[:8], query[:8], query[:8], 1, 8, 1, 8, True, True),
            (query, query[:18], query[:18], query[:18], 1, 18, 1, 18, False, True),
            (query, query[:17], query[:17], query[:17], 1, 17, 1, 17, False, False),
            (query, mismatch, query, mismatch, 1, 20, 1, 20, False, True),
            (query, reverse_complement(query), reverse_complement(query), reverse_complement(query), 20, 1, 1, 20, True, True),
            (query, reverse_complement(query), query, query, 1, 20, 20, 1, True, True),
            (query, gapped_subject, query[:10] + "-" + query[10:], gapped_subject, 1, 20, 1, 21, False, False),
        ]
        for sequence, subject, aligned_query, aligned_subject, query_start, query_end, subject_start, subject_end, strict, tolerant in examples:
            tokens, exact_run, exact, mismatches = [], 0, 0, 0
            for query_base, subject_base in zip(aligned_query, aligned_subject):
                if query_base == subject_base:
                    exact_run += 1
                    exact += 1
                else:
                    if exact_run:
                        tokens.append(str(exact_run))
                    exact_run = 0
                    tokens.append(query_base + subject_base)
                    mismatches += "-" not in (query_base, subject_base)
            if exact_run:
                tokens.append(str(exact_run))
            span = abs(query_end - query_start) + 1
            values = dict(qseqid="query", sseqid="FIX000001.1", staxids="0", sscinames="", pident=100 * exact / len(aligned_query),
                          length=len(aligned_query), qlen=len(sequence), qcovs=100 * span / len(sequence), mismatch=mismatches,
                          gapopen=int("-" in aligned_query + aligned_subject), qstart=query_start, qend=query_end,
                          sstart=subject_start, send=subject_end, slen=len(subject), evalue=1, bitscore=1, stitle="fixture",
                          qseq=aligned_query, sseq=aligned_subject, btop="".join(tokens),
                          sstrand="plus" if subject_end >= subject_start else "minus", qcovhsp=100 * span / len(sequence))
            row = parse_blast_hsp([str(values[field]) for field in BLAST_FIELDS], {"query": sequence})
            expected = dict(strict=strict, tolerant=tolerant)
            self.assertEqual(independent_hsp_truth(row, sequence, subject), expected)
            self.assertEqual(dict(strict=strict_hsp(row), tolerant=tolerant_hsp(row)), expected)

    def test_single_base_subject_span_is_accepted_in_either_orientation(self):
        # BLAST+ 2.17.0 reports a one-base minus-strand HSP with sstart == send and sstrand "plus" (XML hit frame -1).
        query, subject = "ACGTTGCA" * 3, "GGGGAGGGG"
        values = dict(qseqid="query", sseqid="FIX000001.1", staxids="0", sscinames="", pident=100, length=1,
                      qlen=len(query), qcovs=4, mismatch=0, gapopen=0, qstart=4, qend=4, sstart=5, send=5,
                      slen=len(subject), evalue=1e11, bitscore=2.5, stitle="fixture", qseq="T", sseq="T", btop="1",
                      sstrand="plus", qcovhsp=4)
        row = parse_blast_hsp([str(values[field]) for field in BLAST_FIELDS], {"query": query})
        self.assertEqual(independent_hsp_truth(row, query, subject), dict(strict=False, tolerant=False))
        self.assertEqual(dict(strict=strict_hsp(row), tolerant=tolerant_hsp(row)), dict(strict=False, tolerant=False))
        plus = dict(row, sseq="A", qseq="A", qstart=1, qend=1)
        self.assertEqual(independent_hsp_truth(plus, query, subject), dict(strict=False, tolerant=False))
        with self.assertRaises(ValueError):
            independent_hsp_truth(dict(row, sseq="C"), query, subject)


if __name__ == "__main__":
    unittest.main()