import ctypes
import math
import os
import random
import subprocess
import tarfile
import xml.etree.ElementTree as ElementTree
from pathlib import Path

from .contracts import ContractError
from .design import digest
from .example import save_fasta, save_json
from .matching import BLAST_FIELDS
from .search import iter_hsps, search_settings
from .sequences import reverse_complement


NCBI_SOURCE_SHA256 = "502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147"
NCBI_SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-src.tar.gz"


def execute(command: list[str], timeout: int = 3600) -> dict:
    process = subprocess.run(command, capture_output=True, text=True, timeout=timeout,
                             env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"))
    if process.returncode:
        raise ContractError(f"Fixture command failed ({process.returncode}): {command[0]}: {process.stderr}")
    return dict(command=command, exit_code=process.returncode, stdout=process.stdout, stderr=process.stderr)


def executable_identities(binary_directory: Path) -> dict:
    identities = {}
    for name in ("blastn", "makeblastdb", "blast_formatter", "blastdbcmd"):
        binary = (binary_directory / name).resolve()
        version = execute([str(binary), "-version"], 30)["stdout"].strip()
        if "2.17.0" not in version:
            raise ContractError("Current-version benchmark requires BLAST+ 2.17.0: " + name)
        identities[name] = dict(path=str(binary), version=version, sha256=digest(binary))
    return identities


def compile_length_adapter(archive: Path, destination: Path) -> dict:
    if digest(archive) != NCBI_SOURCE_SHA256:
        raise ContractError("NCBI statistical source archive differs from pinned 2.17.0 release")
    destination.mkdir(parents=True, exist_ok=False)
    source_name = "ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/blast_stat.c"
    setup_name = "ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/blast_setup.c"
    with tarfile.open(archive, "r:gz") as package:
        source = package.extractfile(source_name).read().decode("utf-8")
        setup = package.extractfile(setup_name).read()
    start = source.index("Int4\nBLAST_ComputeLengthAdjustment(")
    end = source.index("\n}\n", start) + 3
    routine = source[start:end]
    preamble = "#include <stdint.h>\n#include <math.h>\ntypedef int32_t Int4;\ntypedef int64_t Int8;\ntypedef int Boolean;\n#define FALSE 0\n#define TRUE 1\n#define MAX(first, second) ((first) > (second) ? (first) : (second))\n"
    (destination / "blast_stat.c").write_text(source, encoding="utf-8")
    (destination / "blast_setup.c").write_bytes(setup)
    (destination / "length_adapter.c").write_text(preamble + routine, encoding="utf-8")
    library = destination / "length_adapter.so"
    build = execute(["cc", "-shared", "-fPIC", "-O2", "-o", str(library), str(destination / "length_adapter.c"), "-lm"], 60)
    report = dict(source_url=NCBI_SOURCE_URL, archive_sha256=digest(archive), source_member=source_name,
                  extraction="Unmodified BLAST_ComputeLengthAdjustment; integer typedefs and MAX macro only",
                  scoring="reward=1 penalty=-3 gapopen=5 gapextend=2: ungapped K/Lambda/H; alpha/Lambda=1/H; beta=0",
                  strand_convention="Per nucleotide context, identical plus/minus blocks; no factor of two in query length",
                  compiler=execute(["cc", "--version"], 30)["stdout"], build=build,
                  files={path.name: digest(path) for path in destination.iterdir()})
    save_json(destination / "adapter.json", report)
    return report


def effective_search_space(library: Path, length: int, bases: int, sequences: int, statistics: dict) -> dict:
    if not 0 < length < 2 ** 31 or not 0 < sequences < 2 ** 31 or not 0 < bases < 2 ** 63:
        raise ContractError("Statistical inputs exceed the pinned engine's integer domain")
    routine = ctypes.CDLL(str(library.resolve())).BLAST_ComputeLengthAdjustment
    routine.argtypes = [ctypes.c_double] * 4 + [ctypes.c_int32, ctypes.c_int64, ctypes.c_int32, ctypes.POINTER(ctypes.c_int32)]
    routine.restype = ctypes.c_int32
    adjustment = ctypes.c_int32()
    result = routine(statistics["kappa"], math.log(statistics["kappa"]), 1 / statistics["entropy"], 0,
                     length, bases, sequences, ctypes.byref(adjustment))
    effective = max(1, bases - sequences * adjustment.value) * (length - adjustment.value)
    return dict(query_length=length, database_bases=bases, database_sequences=sequences,
                length_adjustment=adjustment.value, routine_status=result, effective_search_space=effective,
                kappa=statistics["kappa"], lambda_value=statistics["lambda_value"], entropy=statistics["entropy"],
                alpha_over_lambda=1 / statistics["entropy"], beta=0,
                exact_score=length, exact_evalue=statistics["kappa"] * effective * math.exp(-statistics["lambda_value"] * length))


def archive_search(query_file: Path, database: Path, destination: Path, binaries: dict, settings: dict, timeout: int = 3600) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    command = [binaries["blastn"]["path"], "-query", str(query_file.resolve()), "-db", str(database.resolve()),
               "-num_threads", "1", "-outfmt", "11", "-out", str((destination / "blast.asn").resolve())]
    for name, value in settings.items():
        if value is not None:
            command.extend(["-" + name, str(value)])
    search = execute(command, timeout)
    formats = {}
    for label, specification, filename in (("tabular", "7 " + " ".join(BLAST_FIELDS), "blast.tsv"), ("xml", "5", "blast.xml")):
        formats[label] = execute([binaries["blast_formatter"]["path"], "-archive", str((destination / "blast.asn").resolve()),
                                  "-outfmt", specification, "-out", str((destination / filename).resolve())], timeout)
    report = dict(mode="working", production_database=False, search=search, formats=formats, settings=settings,
                  query_sha256=digest(query_file), files={path.name: digest(path) for path in destination.iterdir()})
    save_json(destination / "request.json", report)
    return report


def xml_statistics(path: Path) -> list[dict]:
    root = ElementTree.parse(path).getroot()
    rows = []
    for iteration in root.findall(".//Iteration"):
        stats = iteration.find("Iteration_stat/Statistics")
        if stats is None:
            raise ContractError("BLAST XML lacks effective-search-space statistics")
        rows.append(dict(query_id=iteration.findtext("Iteration_query-def").split()[0],
                         query_length=int(iteration.findtext("Iteration_query-len")),
                         database_bases=int(stats.findtext("Statistics_db-len")),
                         database_sequences=int(stats.findtext("Statistics_db-num")),
                         effective_search_space=int(stats.findtext("Statistics_eff-space")),
                         kappa=float(stats.findtext("Statistics_kappa")), lambda_value=float(stats.findtext("Statistics_lambda")),
                         entropy=float(stats.findtext("Statistics_entropy")),
                         hsps=[dict(score=float(hsp.findtext("Hsp_score")), evalue=float(hsp.findtext("Hsp_evalue")),
                                    bitscore=float(hsp.findtext("Hsp_bit-score"))) for hsp in iteration.findall(".//Hsp")]))
    return rows


def calibrate(archive: Path, binary_directory: Path, destination: Path, bases: int, sequences: int, lengths: tuple[int, ...]) -> dict:
    if not lengths or any(length < 4 for length in lengths):
        raise ContractError("Calibration requires explicit query lengths >=4")
    destination.mkdir(parents=True, exist_ok=False)
    binaries = executable_identities(binary_directory)
    adapter = compile_length_adapter(archive, destination / "adapter")
    generator = random.Random(20260922)
    queries = {f"length_{length}": "".join(generator.choices("ACGT", k=length)) for length in sorted(set(lengths))}
    queries.update({f"poly_{length}": "A" * length for length in sorted(set(lengths))})
    save_fasta(destination / "queries.fasta", queries)
    validations, parameter_rows = [], []
    library = destination / "adapter/length_adapter.so"
    for label, flank, repeats in (("tiny", 30, 1), ("larger", 250, 3)):
        fixture = destination / label
        fixture.mkdir()
        subjects = {f"CAL{index:05d}.1": "N" * flank + sequence + "N" * flank + reverse_complement(sequence) + "N" * flank
                    for index, sequence in enumerate(list(queries.values()) * repeats, 1)}
        save_fasta(fixture / "subjects.fasta", subjects)
        build = execute([binaries["makeblastdb"]["path"], "-in", str(fixture / "subjects.fasta"), "-dbtype", "nucl", "-parse_seqids", "-out", str(fixture / "database")])
        save_json(fixture / "database_build.json", build)
        tabular, stats = [], []
        for length in sorted(set(lengths)):
            same_length = {name: sequence for name, sequence in queries.items() if len(sequence) == length}
            query_file = fixture / f"queries_{length}.fasta"
            search = fixture / f"search_{length}"
            save_fasta(query_file, same_length)
            archive_search(query_file, fixture / "database", search, binaries, search_settings(2000))
            tabular.extend(iter_hsps(search / "blast.tsv", same_length))
            stats.extend(xml_statistics(search / "blast.xml"))
        if {row["query_id"] for row in stats} != set(queries):
            raise ContractError("Calibration XML query membership changed")
        for row in stats:
            actual_bases, actual_sequences = sum(map(len, subjects.values())), len(subjects)
            if (row["database_bases"], row["database_sequences"]) != (actual_bases, actual_sequences):
                raise ContractError("Calibration database statistics disagree with known fixture")
            calculated = effective_search_space(library, row["query_length"], actual_bases, actual_sequences, row)
            if calculated["effective_search_space"] != row["effective_search_space"]:
                raise ContractError("NCBI length adapter disagrees with emitted effective search space")
            raw = [hsp for hsp in tabular if hsp["qseqid"] == row["query_id"]]
            if len(raw) != len(row["hsps"]):
                raise ContractError("Calibration archive XML and format-7 HSP counts disagree")
            tabular_values = sorted((hsp["evalue"], hsp["bitscore"]) for hsp in raw)
            xml_values = sorted((hsp["evalue"], hsp["bitscore"]) for hsp in row["hsps"])
            for formatted, precise in zip(tabular_values, xml_values):
                bit_score_agrees = (formatted[1] == math.floor(precise[1]) if formatted[1] >= 100
                                    else math.isclose(formatted[1], precise[1], abs_tol=0.051))
                if not math.isclose(formatted[0], precise[0], rel_tol=0.051, abs_tol=1e-300) or not bit_score_agrees:
                    raise ContractError("Calibration archive XML and tabular numerical output disagree")
            for hsp in row["hsps"]:
                expected = row["kappa"] * row["effective_search_space"] * math.exp(-row["lambda_value"] * hsp["score"])
                if not math.isclose(expected, hsp["evalue"], rel_tol=1e-5, abs_tol=1e-300):
                    raise ContractError("NCBI statistical parameters do not reproduce emitted HSP E-values")
            validations.append(dict(fixture=label, query_id=row["query_id"], **calculated, hsps=len(raw), status="engine_agreement"))
            parameter_rows.append(row)
    reference = parameter_rows[0]
    if any(any(not math.isclose(row[field], reference[field], rel_tol=1e-10) for field in ("kappa", "lambda_value", "entropy")) for row in parameter_rows):
        raise ContractError("Query-specific statistical parameters differ; per-length calibration is insufficient")
    spaces = [effective_search_space(library, length, bases, sequences, reference) for length in sorted(set(lengths))]
    if max(row["exact_evalue"] for row in spaces) >= 1e12:
        raise ContractError("Predeclared diagnostic E-value does not cover exact query witnesses")
    report = dict(status="engine_validated_version_matched", binaries=binaries, adapter=adapter,
                  validation=validations, search_spaces=spaces, diagnostic_evalue=1e12,
                  calibration_queries=len(queries), calibration_databases=2, production_searches=0,
                  limitation="Statistical scale calibration only, not a production database search or a completeness guarantee",
                  files={str(path.relative_to(destination)): digest(path) for path in destination.rglob("*") if path.is_file()})
    save_json(destination / "calibration.json", report)
    return report