import hashlib
import json
import platform
import sys
import threading
from collections import Counter, defaultdict
from importlib.metadata import distribution, version
from pathlib import Path

from .contracts import ContractError, write_table
from .design import digest
from .example import save_fasta, save_json
from .sequences import canonical, canonical_id, read_fasta


UPSTREAM_COMMIT = "b54df34f60246ff236ec58bbccb40af3cea608ed"
SHUFFLE_LOCK = threading.Lock()
DRAW_SCHEMA = dict(type="object", additionalProperties=False, properties={
    "draw_id": {"type": "string"}, "parent_id": {"type": "string"},
    "draw_index": {"type": "integer"}, "seed": {"type": "integer"},
    "sequence_id": {"type": "string"}, "parent_identical": {"type": "string", "enum": ["yes", "no"]},
})
PARENT_SCHEMA = dict(type="object", additionalProperties=False, properties={
    "parent_id": {"type": "string"}, "physical_ids_json": {"type": "string"},
    "length": {"type": "integer"}, "draws": {"type": "integer"},
    "distinct_sequences": {"type": "integer"}, "distinct_canonical_sequences": {"type": "integer"},
    "parent_identical_draws": {"type": "integer"},
})


def dinucleotides(sequence: str) -> Counter:
    return Counter(sequence[index:index + 2] for index in range(len(sequence) - 1))


def check_draw(parent: str, draw: str) -> None:
    if (not parent or set(parent + draw) - set("ACGT") or len(parent) != len(draw)
            or Counter(parent) != Counter(draw) or dinucleotides(parent) != dinucleotides(draw)
            or parent[0] != draw[0] or parent[-1] != draw[-1]):
        raise ContractError("Control failed length, alphabet, base, dinucleotide or endpoint invariants")


def draw_seed(seed: int, parent_id: str, draw_index: int) -> int:
    value = f"ushuffle-k2-v1:{seed}:{parent_id}:{draw_index}".encode("ascii")
    return 1 + int.from_bytes(hashlib.sha256(value).digest()[:8], "big") % 2147483647


def shuffle_draw(sequence: str, seed: int) -> str:
    import ushuffle

    if not sequence or set(sequence) - set("ACGT") or not 0 < seed <= 2147483647:
        raise ContractError("uShuffle requires a nonempty ACGT sequence and a positive signed 32-bit seed")
    with SHUFFLE_LOCK:
        ushuffle.set_seed(seed)
        result = ushuffle.shuffle(sequence.encode("ascii"), 2).decode("ascii")
    check_draw(sequence, result)
    return result


def algorithm_record() -> dict:
    import ushuffle

    installed = version("ushuffle")
    if installed != "1.1.3":
        raise ContractError(f"Expected the pinned uShuffle 1.1.3 build, observed {installed}")
    direct_url = distribution("ushuffle").read_text("direct_url.json")
    origin = json.loads(direct_url) if direct_url else {}
    if UPSTREAM_COMMIT not in origin.get("url", "") and origin.get("vcs_info", {}).get("commit_id") != UPSTREAM_COMMIT:
        raise ContractError("uShuffle install does not identify the pinned upstream commit")
    return dict(name="ushuffle", version=installed, upstream_commit=UPSTREAM_COMMIT, direct_url=origin,
                extension_sha256=digest(Path(ushuffle.__file__)), python=sys.version,
                platform=platform.platform(), libc=list(platform.libc_ver()), k=2,
                sampling="with_replacement", endpoints="first_and_last_base_fixed",
                seed_derivation="sha256_ushuffle-k2-v1_to_positive_int32_per_parent_draw",
                sampling_caveat="Published algorithm; finite PRNG and modulo arithmetic. Enumerable tests are diagnostics, not a proof of uniformity.")


def generate_controls(query_fasta: Path, destination: Path, draws: int, seed: int = 42) -> dict:
    if draws < 1 or seed < 0:
        raise ContractError("Draw count must be positive and master seed must be nonnegative")
    queries = read_fasta(query_fasta)
    algorithm = algorithm_record()
    parents = {}
    aliases = defaultdict(list)
    for query_id, sequence in queries.items():
        parent_id = canonical_id(sequence)
        parents[parent_id] = canonical(sequence)
        aliases[parent_id].append(query_id)
    sequences, draw_rows, parent_rows = {}, [], []
    seed_counts = Counter()
    for parent_id, sequence in sorted(parents.items()):
        parent_draws = []
        for draw_index in range(1, draws + 1):
            identifier = "ctrl_" + parent_id.removeprefix("seq_") + f"_d{draw_index:06d}"
            effective_seed = draw_seed(seed, parent_id, draw_index)
            generated = shuffle_draw(sequence, effective_seed)
            sequences[identifier] = generated
            parent_draws.append(generated)
            seed_counts[effective_seed] += 1
            draw_rows.append(dict(draw_id=identifier, parent_id=parent_id, draw_index=draw_index, seed=effective_seed,
                                  sequence_id=canonical_id(generated), parent_identical="yes" if generated == sequence else "no"))
        parent_rows.append(dict(parent_id=parent_id, physical_ids_json=json.dumps(sorted(aliases[parent_id])), length=len(sequence),
                                draws=draws, distinct_sequences=len(set(parent_draws)),
                                distinct_canonical_sequences=len({canonical(draw) for draw in parent_draws}),
                                parent_identical_draws=parent_draws.count(sequence)))
    destination.mkdir(parents=True, exist_ok=False)
    save_fasta(destination / "controls.fasta", sequences)
    write_table(destination / "draws.tsv", draw_rows, DRAW_SCHEMA)
    write_table(destination / "parents.tsv", parent_rows, PARENT_SCHEMA)
    report = dict(stage="control_generation_qc", research_freeze=False, source_sha256=digest(query_fasta),
                  algorithm=algorithm, master_seed=seed, canonical_parents=len(parents), physical_queries=len(queries),
                  draws_per_parent=draws, total_draws=len(draw_rows), all_invariants_valid=True,
                  repeated_draws_within_parent=sum(row["draws"] - row["distinct_sequences"] for row in parent_rows),
                  parent_identical_draws=sum(row["parent_identical_draws"] for row in parent_rows),
                  repeated_seed_values=sum(count - 1 for count in seed_counts.values()),
                  zero_diversity_parents=sum(row["distinct_sequences"] == 1 for row in parent_rows),
                  files={name: digest(destination / name) for name in ("controls.fasta", "draws.tsv", "parents.tsv")})
    save_json(destination / "control_qc.json", report)
    return report