"""Planted-sequence benchmark fixtures derived from the real parents of one panel.

Panel A plants every parent exactly and reverse-complemented in internal and terminal contexts, plus a negative
without a shared seed. Panel B tests boundary templates (GC-balanced, periodic, homopolymer) at the panel's
lengths with all 20 boundary variants. Panel C plants the tolerant-boundary and combined clip/substitution
variants of every parent. Panel D repeats the boundary templates, with the combined variants and the declared
panel-derived templates (for example a primer-like template, supplied at run time as data), at the declared extra
lengths. Panels A/B and C/D use separate seeded random streams, so the extra lengths never change the A/B fixtures.
"""
import random
from collections import Counter, defaultdict

from .benchmark_truth import reverse_complement, sequence_witnesses
from .contracts import ContractError

FIXTURE_SEED = 20260922
VARIANT_SEED = 20261020
REAL_VARIANTS = ("exact", "reverse", "substitution_at_center", "substitution_over_center", "clip_both_at", "clip_both_below")


def _random_sequence(generator, length):
    return "".join(generator.choices("ACGT", k=length))


def _context(sequence, context, index, generator):
    if context == "internal":
        return _random_sequence(generator, 100) + sequence + _random_sequence(generator, 100), 101
    flank = _random_sequence(generator, 200)
    return (sequence + flank, 1) if index % 2 == 0 else (flank + sequence, 201)


def _complexity(sequence):
    return "low_complexity" if max(Counter(sequence).values()) / len(sequence) >= 0.8 or len({sequence[index:index + 2] for index in range(len(sequence) - 1)}) <= 3 else "other"


def _negative(query, generator, kind, context, index):
    seeds = {oriented[offset:offset + 7] for oriented in (query, reverse_complement(query)) for offset in range(len(query) - 6)}
    for attempt in range(1, 101):
        sequence = _random_sequence(generator, len(query))
        if kind == "negative_seed":
            sequence = query[:7] + sequence[7:]
        subject, start = _context(sequence, context, index, generator)
        has_seed = any(seed in subject for seed in seeds)
        if kind == "negative_no_seed" and has_seed or kind == "negative_seed" and not has_seed:
            continue
        if sequence_witnesses(query, subject)["tolerant"] is None:
            return subject, start, attempt
    return None, None, 100


def _mutate(sequence, count, placement):
    positions = list(range(len(sequence))) if placement == "end" else sorted(range(len(sequence)), key=lambda index: (abs(index - len(sequence) // 2), index))
    result = list(sequence)
    for position in positions[:count]:
        result[position] = next(base for base in "ACGT" if base != result[position])
    return "".join(result)


def variants(sequence):
    length = len(sequence)
    span = (9 * length + 9) // 10
    middle = length // 2
    change = next(base for base in "ACGT" if base != sequence[middle])
    result = {"exact": sequence, "reverse": reverse_complement(sequence)}
    for label, count in (("at", length // 20), ("over", length // 20 + 1)):
        for placement in ("center", "end"):
            result[f"substitution_{label}_{placement}"] = _mutate(sequence, count, placement)
    for label, retained in (("at", span), ("below", span - 1)):
        removed = length - retained
        result[f"clip_left_{label}"] = sequence[removed:]
        result[f"clip_right_{label}"] = sequence[:retained]
        result[f"clip_both_{label}"] = sequence[removed // 2:removed // 2 + retained]
    result.update(insert_internal=sequence[:middle] + change + sequence[middle:],
                  delete_internal=sequence[:middle] + sequence[middle + 1:],
                  insert_terminal=change + sequence, delete_terminal=sequence[1:],
                  split_partial=sequence[:middle] + "N" * (length + 1) + sequence[middle:],
                  negative_no_seed=None, negative_seed=None,
                  reverse_substitution=reverse_complement(_mutate(sequence, length // 20 + 1, "center")))
    if len(result) != 20:
        raise ContractError("Boundary fixture variant count changed")
    return result


def template_sequence(template: dict, length: int) -> str:
    """A declared template at `length`: its sequence, extended by repeating the part from `repeat_from` (one-based)."""
    sequence, tail = template["sequence"], template["sequence"][template["repeat_from"] - 1:]
    if not tail:
        raise ContractError(f"Template {template['name']} has nothing to repeat")
    while len(sequence) < length:
        sequence += tail
    return sequence[:length]


def make_panels(originals: dict[str, str], extra_lengths: tuple[int, ...] = (), seed: int = FIXTURE_SEED,
                variant_seed: int = VARIANT_SEED, templates=()) -> tuple[list, list]:
    """Shards of panels A-D for the given real parents; extra lengths add panel D boundary templates and the declared
    `templates` (dicts with name, sequence and repeat_from)."""
    generator = random.Random(seed)
    grouped = defaultdict(list)
    for query_id, sequence in sorted(originals.items()):
        grouped[len(sequence)].append((query_id, sequence))
    shards, exceptions = [], []
    subject_index = 0

    def add_subject(shard, query_id, sequence, variant, context, index, complexity, template=None):
        nonlocal subject_index
        attempt = 1
        if template is None:
            subject, start, attempt = _negative(sequence, generator, variant, context, index)
        else:
            subject, start = _context(template, context, index, generator)
        if subject is None:
            exceptions.append(dict(shard_id=shard["shard_id"], query_id=query_id, error_class=variant, context=context,
                                   attempts=attempt, status="negative_unavailable_after_predeclared_attempts"))
            return
        subject_index += 1
        subject_id = f"FIX{subject_index:07d}.1"
        shard["subjects"][subject_id] = subject
        shard["cases"].append(dict(query_id=query_id, subject_id=subject_id, error_class=variant, context=context,
                                   length=len(sequence), complexity=complexity, insert_start=start,
                                   insert_end=start + (len(template) if template is not None else len(sequence)) - 1,
                                   attempts=attempt, designated_orientation="-" if variant.startswith("reverse") else "+"))

    for length, parents in sorted(grouped.items()):
        for offset in range(0, len(parents), 16):
            shard = dict(shard_id=f"A_{length:03d}_{offset // 16:03d}", panel="A", length=length, queries=dict(parents[offset:offset + 16]), subjects={}, cases=[])
            for index, (query_id, sequence) in enumerate(shard["queries"].items(), offset):
                for context in ("internal", "terminal"):
                    for variant, template in (("exact", sequence), ("reverse", reverse_complement(sequence))):
                        add_subject(shard, query_id, sequence, variant, context, index, _complexity(sequence), template)
                add_subject(shard, query_id, sequence, "negative_no_seed", "internal", index, _complexity(sequence))
            shards.append(shard)
    for length in sorted(grouped):
        queries, compositions = {}, {}
        for fraction in (0.25, 0.50, 0.75):
            for replicate in range(2):
                query_id = f"boundary_{length:03d}_gc{int(100 * fraction)}_{replicate}"
                gc_count = round(length * fraction)
                bases = generator.choices("GC", k=gc_count) + generator.choices("AT", k=length - gc_count)
                generator.shuffle(bases)
                queries[query_id] = "".join(bases)
                compositions[query_id] = f"gc_{int(100 * fraction)}"
        for label, sequence in (("periodic", ("ACG" * length)[:length]), ("homopolymer", "A" * length)):
            query_id = f"boundary_{length:03d}_{label}"
            queries[query_id], compositions[query_id] = sequence, label
        shard = dict(shard_id=f"B_{length:03d}", panel="B", length=length, queries=queries, subjects={}, cases=[])
        for index, (query_id, sequence) in enumerate(queries.items()):
            for variant, template in variants(sequence).items():
                for context in ("internal", "terminal"):
                    add_subject(shard, query_id, sequence, variant, context, index, compositions[query_id], template)
        shards.append(shard)
    variant_shards, variant_exceptions = make_variant_panels(originals, extra_lengths, variant_seed, templates)
    shards, exceptions = shards + variant_shards, exceptions + variant_exceptions
    if any(len(shard["subjects"]) >= 2000 or shard["panel"] == "A" and (len(shard["queries"]) > 16 or len(shard["subjects"]) > 80) for shard in shards):
        raise ContractError("Fixture shard exceeds the declared query or target budget")
    return shards, exceptions


def combined_variants(sequence):
    """Clipping to the minimum tolerant span combined with substitutions at and just over the identity limit."""
    length = len(sequence)
    span = (9 * length + 9) // 10
    removed = length - span
    clipped = sequence[removed // 2:removed // 2 + span]
    return {"clip_both_at_substitution_at": _mutate(clipped, span // 20, "center"),
            "clip_both_at_substitution_over": _mutate(clipped, span // 20 + 1, "center"),
            "reverse_clip_both_at_substitution_at": reverse_complement(_mutate(clipped, span // 20, "center"))}


def make_variant_panels(originals: dict[str, str], extra_lengths: tuple[int, ...], seed: int = VARIANT_SEED,
                        templates=()) -> tuple[list, list]:
    """Panel C: tolerant and combined variants of every real parent. Panel D: synthetic and declared templates at extra
    boundary lengths."""
    generator = random.Random(seed)
    shards, exceptions = [], []
    counter = 0

    def add(shard, query_id, sequence, variant, context, index, complexity, template):
        nonlocal counter
        attempt = 1
        if template is None:
            subject, start, attempt = _negative(sequence, generator, variant, context, index)
        else:
            subject, start = _context(template, context, index, generator)
        if subject is None:
            exceptions.append(dict(shard_id=shard["shard_id"], query_id=query_id, error_class=variant, context=context,
                                   attempts=attempt, status="negative_unavailable_after_predeclared_attempts"))
            return
        counter += 1
        subject_id = f"VAR{counter:07d}.1"
        shard["subjects"][subject_id] = subject
        shard["cases"].append(dict(query_id=query_id, subject_id=subject_id, error_class=variant, context=context,
                                   length=len(sequence), complexity=complexity, insert_start=start,
                                   insert_end=start + (len(template) if template is not None else len(sequence)) - 1,
                                   attempts=attempt, designated_orientation="-" if variant.startswith("reverse") else "+"))

    grouped = defaultdict(list)
    for query_id, sequence in sorted(originals.items()):
        grouped[len(sequence)].append((query_id, sequence))
    for length, parents in sorted(grouped.items()):
        for offset in range(0, len(parents), 8):
            shard = dict(shard_id=f"C_{length:03d}_{offset // 8:03d}", panel="C", length=length,
                         queries=dict(parents[offset:offset + 8]), subjects={}, cases=[])
            for index, (query_id, sequence) in enumerate(shard["queries"].items(), offset):
                variant_templates = {name: value for name, value in variants(sequence).items() if name in REAL_VARIANTS}
                variant_templates.update(combined_variants(sequence))
                for name, template in variant_templates.items():
                    for context in ("internal", "terminal"):
                        add(shard, query_id, sequence, name, context, index, _complexity(sequence), template)
                add(shard, query_id, sequence, "negative_no_seed", "internal", index, _complexity(sequence), None)
            shards.append(shard)
    for length in sorted(set(extra_lengths)):
        queries, compositions = {}, {}
        for fraction in (0.25, 0.50, 0.75):
            for replicate in range(2):
                query_id = f"extra_{length:03d}_gc{int(100 * fraction)}_{replicate}"
                gc_count = round(length * fraction)
                bases = generator.choices("GC", k=gc_count) + generator.choices("AT", k=length - gc_count)
                generator.shuffle(bases)
                queries[query_id] = "".join(bases)
                compositions[query_id] = f"gc_{int(100 * fraction)}"
        named = [("periodic", ("ACG" * length)[:length], "low_complexity"), ("homopolymer", "A" * length, "low_complexity")]
        named += [(template["name"], template_sequence(template, length), template["name"]) for template in templates]
        for label, sequence, composition in named:
            query_id = f"extra_{length:03d}_{label}"
            if query_id in queries:
                raise ContractError("Duplicate benchmark template name: " + label)
            queries[query_id], compositions[query_id] = sequence, composition
        shard = dict(shard_id=f"D_{length:03d}", panel="D", length=length, queries=queries, subjects={}, cases=[])
        for index, (query_id, sequence) in enumerate(queries.items()):
            variant_templates = dict(variants(sequence))
            variant_templates.update(combined_variants(sequence))
            for name, template in variant_templates.items():
                for context in ("internal", "terminal"):
                    add(shard, query_id, sequence, name, context, index, compositions[query_id], template)
        shards.append(shard)
    return shards, exceptions
