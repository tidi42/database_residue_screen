import math
import random
import xml.etree.ElementTree as ElementTree
from collections import Counter
from pathlib import Path

from .benchmark_calibration import (ACTUAL_LENGTHS, compile_length_adapter, effective_search_space,
                                    executable_identities, execute, xml_statistics)
from .contracts import ContractError
from .design import digest
from .example import save_fasta, save_json
from .matching import BLAST_FIELDS
from .search import iter_hsps, search_settings, summarise_batch
from .sequences import accession_version, reverse_complement


def direct_search(query_file: Path, queries: dict, database: Path, destination: Path, binaries: dict, settings: dict, timeout: int = 3600) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    commands = {}
    common = [binaries["blastn"]["path"], "-query", str(query_file.resolve()), "-db", str(database.resolve()), "-num_threads", "1"]
    for name, value in settings.items():
        if value is not None:
            common.extend(["-" + name, str(value)])
    for label, specification, filename in (("tabular", "7 " + " ".join(BLAST_FIELDS), "blast.tsv"), ("xml", "5", "blast.xml")):
        commands[label] = execute(common + ["-outfmt", specification, "-out", str((destination / filename).resolve())], timeout)
    rows = list(iter_hsps(destination / "blast.tsv", queries))
    tabular = Counter((row["qseqid"], accession_version(row["sseqid"]), row["qstart"], row["qend"], row["sstart"], row["send"], row["qseq"], row["sseq"]) for row in rows)
    xml_hsps = Counter()
    for iteration in ElementTree.parse(destination / "blast.xml").getroot().findall(".//Iteration"):
        query_id = iteration.findtext("Iteration_query-def").split()[0]
        for hit in iteration.findall(".//Hit"):
            subject_id = accession_version(hit.findtext("Hit_id")) or accession_version(hit.findtext("Hit_def").split()[0])
            for hsp in hit.findall(".//Hsp"):
                xml_hsps[query_id, subject_id, int(hsp.findtext("Hsp_query-from")), int(hsp.findtext("Hsp_query-to")),
                         int(hsp.findtext("Hsp_hit-from")), int(hsp.findtext("Hsp_hit-to")), hsp.findtext("Hsp_qseq"), hsp.findtext("Hsp_hseq")] += 1
    if tabular != xml_hsps:
        raise ContractError("Direct XML/format-7 searches disagree on exact HSP content or multiplicity")
    statuses = summarise_batch(destination / "blast.tsv", queries, destination.name, settings["max_target_seqs"])
    report = dict(mode="working", production_database=False, commands=commands, settings=settings, query_sha256=digest(query_file),
                  all_queries_completed=True, xml_tabular_hsp_multisets_equal=True, statuses=statuses,
                  files={path.name: digest(path) for path in destination.iterdir() if path.is_file()})
    save_json(destination / "complete.json", report)
    return report


def calibrate_direct(archive: Path, binary_directory: Path, destination: Path, bases: int, sequences: int) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    binaries = executable_identities(binary_directory)
    adapter = compile_length_adapter(archive, destination / "adapter")
    generator = random.Random(20260922)
    queries = {f"length_{length}": "".join(generator.choices("ACGT", k=length)) for length in ACTUAL_LENGTHS}
    queries.update({f"poly_{length}": "A" * length for length in (8, 10, 16, 24, 111)})
    save_fasta(destination / "queries.fasta", queries)
    validations, parameter_rows = [], []
    library = destination / "adapter/length_adapter.so"
    for label, flank, repeats in (("tiny", 30, 1), ("larger", 250, 3)):
        fixture = destination / label
        fixture.mkdir()
        subjects = {f"CAL{index:05d}.1": "N" * flank + sequence + "N" * flank + reverse_complement(sequence) + "N" * flank
                    for index, sequence in enumerate(list(queries.values()) * repeats, 1)}
        save_fasta(fixture / "subjects.fasta", subjects)
        save_json(fixture / "database_build.json", execute([binaries["makeblastdb"]["path"], "-in", str(fixture / "subjects.fasta"),
                                                           "-dbtype", "nucl", "-parse_seqids", "-out", str(fixture / "database")]))
        for query_id, query in queries.items():
            query_file, search = fixture / (query_id + ".fasta"), fixture / query_id
            save_fasta(query_file, {query_id: query})
            direct_search(query_file, {query_id: query}, fixture / "database", search, binaries, search_settings(2000))
            stats = xml_statistics(search / "blast.xml")
            if len(stats) != 1 or stats[0]["query_id"] != query_id:
                raise ContractError("Direct calibration XML query membership changed")
            row = stats[0]
            actual_bases, actual_sequences = sum(map(len, subjects.values())), len(subjects)
            if (row["database_bases"], row["database_sequences"]) != (actual_bases, actual_sequences):
                raise ContractError("Calibration database statistics disagree with frozen subjects")
            calculated = effective_search_space(library, len(query), actual_bases, actual_sequences, row)
            if calculated["effective_search_space"] != row["effective_search_space"]:
                raise ContractError("NCBI adapter disagrees with direct per-query engine search space")
            for hsp in row["hsps"]:
                expected = row["kappa"] * row["effective_search_space"] * math.exp(-row["lambda_value"] * hsp["score"])
                if not math.isclose(expected, hsp["evalue"], rel_tol=1e-5, abs_tol=1e-300):
                    raise ContractError("Direct engine E-value disagrees with its reported statistical parameters")
            validations.append(dict(fixture=label, query_id=query_id, **calculated, hsps=len(row["hsps"]), status="engine_agreement"))
            parameter_rows.append(row)
    reference = parameter_rows[0]
    if any(any(not math.isclose(row[field], reference[field], rel_tol=1e-10) for field in ("kappa", "lambda_value", "entropy")) for row in parameter_rows):
        raise ContractError("Per-query statistical parameters differ; length-only calibration is insufficient")
    spaces = [effective_search_space(library, length, bases, sequences, reference) for length in ACTUAL_LENGTHS]
    if max(row["exact_evalue"] for row in spaces) >= 1e12:
        raise ContractError("Predeclared diagnostic E-value fails to cover exact witnesses")
    report = dict(status="engine_validated_version_matched", route="direct_per_query_XML_and_format7", binaries=binaries, adapter=adapter,
                  validation=validations, search_spaces=spaces, diagnostic_evalue=1e12, calibration_queries=len(queries), calibration_databases=2,
                  production_searches=0, limitation="Statistical scale calibration, not production retrieval completeness",
                  files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "calibration.json", report)
    return report