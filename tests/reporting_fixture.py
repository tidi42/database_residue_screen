"""Synthetic study outputs for the reporting tests.

A small panel of random sequences (no vendor data) is pushed through the framework's own control generator, design
step and analysis writer; search output, metadata, INSDC, scan, benchmark, target-limit and record-case directories are
written in the stage output layouts. ``build(root)`` returns the directory of every upstream stage.
"""
import json
import random
from pathlib import Path

from database_residue_screen.analysis import write_completed_analysis, event
from database_residue_screen.cap_review import REVIEW_SCHEMA, review_saturated
from database_residue_screen.contracts import read_table, table_schema, write_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.design import digest
from database_residue_screen.eligibility import annotate_record, annotation_defaults
from database_residue_screen.example import save_fasta, save_json
from database_residue_screen.insdc import ANNOTATION_SCHEMA
from database_residue_screen.length_aware import score_calibration
from database_residue_screen.matching import BLAST_FIELDS, tolerant_hsp
from database_residue_screen.outcomes import OUTCOME_SCHEMA
from database_residue_screen.panels import LOGICAL_SCHEMA, prepare_panel
from database_residue_screen.positional import SENSITIVITY_SCHEMA, positional_sensitivity, projects_from_metadata
from database_residue_screen.precision import assess_precision, precision_specification
from database_residue_screen.search import STATUS_SCHEMA, iter_hsps, object_digest, summarise_batch
from database_residue_screen.sequences import canonical, canonical_id, read_fasta, reverse_complement
from database_residue_screen.taxonomy import TAXONOMY_SCHEMA

BASES, SEQUENCES = 1079005865156, 130787235
STATS = dict(kappa=.710602795216363, lambda_value=1.37406312246009, entropy=1.30724660390929)
DRAWS = 15
TARGETS = 5
RANDOM = random.Random(20261109)


def random_sequence(length: int) -> str:
    return "".join(RANDOM.choice("ACGT") for _ in range(length))


def calibrated(lengths, bases=BASES) -> list[dict]:
    return [dict(score_calibration(length, bases, SEQUENCES, STATS), exact_witness_recovered=True) for length in sorted(set(lengths))]


def panel_definition() -> tuple[list[dict], dict[str, str]]:
    nb1, nb2, front, rear = random_sequence(24), random_sequence(24), random_sequence(8), random_sequence(8)
    i7a = random_sequence(10)
    i7b = i7a[:4] + {"A": "C", "C": "G", "G": "T", "T": "A"}[i7a[4]] + i7a[5:]
    i7front, bc01, hairpin = random_sequence(24), random_sequence(16), random_sequence(45)
    sequences = {
        "ONT_NB01_forward": nb1, "ONT_NB02_forward": nb2,
        "ONT_K40_NB01_forward": front + nb1 + rear, "ONT_K40_NB02_forward": front + nb2 + rear,
        "ONT_PCR_adapter": random_sequence(23),
        "ILLUMINA_I7_A_forward": i7a, "ILLUMINA_I7_A_reverse": reverse_complement(i7a), "ILLUMINA_I7_B_forward": i7b,
        "ILLUMINA_I7ADP_A_forward": i7front + i7a, "ILLUMINA_I7ADP_B_forward": i7front + i7b,
        "ILLUMINA_NEXTERA_motif": random_sequence(19), "ILLUMINA_TRUSEQ_adapter": random_sequence(33),
        "BGI_filter_adapter": random_sequence(32),
        "PACBIO_BC01_forward": bc01, "PACBIO_BC01_hairpin_forward": bc01 + hairpin, "PACBIO_SMRTbell_forward": hairpin,
        "PACBIO_C2_primer_forward": "A" * 18 + random_sequence(17),
    }
    rows = [
        ("ONT_NB01_forward", "ONT", "ont_barcode", "ont_native_barcode_24", "barcode", "barcode:1-24=NB01", ""),
        ("ONT_NB02_forward", "ONT", "ont_barcode", "ont_native_barcode_24", "barcode", "barcode:1-24=NB02", ""),
        ("ONT_K40_NB01_forward", "ONT", "ont_native_barcode_context_k40", "ont_native_barcode_context_k40", "composite",
         "front:1-8;barcode:9-32=NB01;rear:33-40", ""),
        ("ONT_K40_NB02_forward", "ONT", "ont_native_barcode_context_k40", "ont_native_barcode_context_k40", "composite",
         "front:1-8;barcode:9-32=NB02;rear:33-40", ""),
        ("ONT_PCR_adapter", "ONT", "ont_pcr_adapter", "ont_pcr_adapter", "adapter", "oligo:1-23", ""),
        ("ILLUMINA_I7_A_forward", "Illumina", "illumina_i7_index", "illumina_i7_bare", "barcode", "index:1-10=A", "i7"),
        ("ILLUMINA_I7_A_reverse", "Illumina", "illumina_i7_index", "illumina_i7_bare", "barcode", "index:1-10=A", "i7"),
        ("ILLUMINA_I7_B_forward", "Illumina", "illumina_i7_index", "illumina_i7_bare", "barcode", "index:1-10=B", "i7"),
        ("ILLUMINA_I7ADP_A_forward", "Illumina", "illumina_i7_index_adapter", "illumina_i7_embedded", "composite", "front:1-24;index:25-34=A", "i7"),
        ("ILLUMINA_I7ADP_B_forward", "Illumina", "illumina_i7_index_adapter", "illumina_i7_embedded", "composite", "front:1-24;index:25-34=B", "i7"),
        ("ILLUMINA_NEXTERA_motif", "Illumina", "illumina_nextera", "illumina_nextera", "adapter", "oligo:1-19", ""),
        ("ILLUMINA_TRUSEQ_adapter", "Illumina", "illumina_truseq", "illumina_truseq", "adapter", "oligo:1-33", ""),
        ("BGI_filter_adapter", "BGI/DNBSEQ", "illumina_bgi", "mgi_dnbseq_filter_adapter", "adapter", "oligo:1-32", ""),
        ("PACBIO_BC01_forward", "PacBio", "pacbio_barcode", "pacbio_barcode_16", "barcode", "barcode:1-16=BC01", ""),
        ("PACBIO_BC01_hairpin_forward", "PacBio", "pacbio_barcode_adapter", "pacbio_barcode_plus_hairpin_61", "composite",
         "barcode:1-16=BC01;hairpin:17-61", ""),
        ("PACBIO_SMRTbell_forward", "PacBio", "pacbio_smrtbell_adapter", "pacbio_smrtbell_adapter", "adapter", "oligo:1-45", ""),
        ("PACBIO_C2_primer_forward", "PacBio", "pacbio_smrtbell_adapter", "pacbio_smrtbell_adapter", "primer", "oligo:1-35", ""),
    ]
    definition = [dict(query_id=query_id, group=group, condition=condition, family=family, construct_class=construct, length=len(sequences[query_id]),
                       segments=segments, index_side=side) for query_id, group, condition, family, construct, segments, side in rows]
    return definition, sequences


RECORDS = {
    "CP000001.1": dict(length=50000, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA1", organism="Escherichia coli", title="Escherichia coli strain X chromosome", taxids="562"),
    "CP000002.1": dict(length=30000, topology="circular", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA1", organism="Escherichia coli", title="Escherichia coli strain X plasmid p1", taxids="562"),
    "MN000003.1": dict(length=1500, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="", organism="Bacillus subtilis", title="Bacillus subtilis 16S rRNA gene", taxids="1423"),
    "AB000004.1": dict(length=40000, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA2", organism="Homo sapiens", title="Homo sapiens genomic clone", taxids="9606"),
    "OQ000005.1": dict(length=6000, topology="linear", record_type="genomic", molecule="DNA", origin="vector", context="unknown",
                       project="", organism="Cloning vector pTEST", title="Cloning vector pTEST, complete sequence", taxids="29278"),
    "LC000006.1": dict(length=2000, topology="linear", record_type="transcript", molecule="RNA", origin="unknown", context="cellular",
                       project="", organism="Mus musculus", title="Mus musculus mRNA", taxids="10090"),
    "OX000007.1": dict(length=120000, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA3", organism="Danio rerio", title="Danio rerio genome assembly, chromosome 7", taxids="7955"),
    "KX000008.1": dict(length=20000, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA4", organism="Arabidopsis thaliana", title="Arabidopsis thaliana amplicon", taxids="3702"),
    "CT000009.1": dict(length=10000, topology="linear", record_type="genomic", molecule="DNA", origin="unknown", context="cellular",
                       project="PRJNA5", organism="Saccharomyces cerevisiae", title="Saccharomyces cerevisiae chromosome", taxids="4932"),
}
TARGETED = {"KX000008.1"}


def btop(query: str, subject: str) -> str:
    parts, run = [], 0
    for a, b in zip(query, subject):
        if a == b:
            run += 1
        else:
            if run:
                parts.append(str(run))
            parts.append(a + b)
            run = 0
    if run:
        parts.append(str(run))
    return "".join(parts)


def hsp(query_id: str, query: str, subject_id: str, slen: int, start: int, plus: bool = True, qstart: int = 1, qend: int | None = None,
        mismatches=(), subject: str | None = None, taxid: str = "0", name: str = "N/A") -> str:
    """One ungapped tabular HSP line; ``start`` is the lower subject coordinate, mismatches are 0-based alignment offsets
    (or ``subject`` gives the aligned subject bases in query orientation)."""
    qend = qend or len(query)
    aligned_query = query[qstart - 1:qend]
    aligned_subject = subject or "".join(({"A": "C", "C": "G", "G": "T", "T": "A"}[base] if index in mismatches else base)
                                         for index, base in enumerate(aligned_query))
    differences = sum(a != b for a, b in zip(aligned_query, aligned_subject))
    length = len(aligned_query)
    sstart, send = (start, start + length - 1) if plus else (start + length - 1, start)
    identity = 100 * (length - differences) / length
    values = dict(qseqid=query_id, sseqid=subject_id, staxids=taxid, sscinames=name, pident=f"{identity:.3f}", length=length, qlen=len(query),
                  qcovs=f"{100 * length / len(query):.0f}", mismatch=differences, gapopen=0, qstart=qstart, qend=qend, sstart=sstart,
                  send=send, slen=slen, evalue="1e-05", bitscore=f"{2 * (length - 4 * differences):.1f}", stitle="synthetic record",
                  qseq=aligned_query, sseq=aligned_subject, btop=btop(aligned_query, aligned_subject),
                  sstrand="plus" if plus else "minus", qcovhsp=f"{100 * length / len(query):.0f}")
    return "\t".join(str(values[name]) for name in BLAST_FIELDS)


def blast_output(queries: dict[str, str], lines: dict[str, list[str]]) -> str:
    text = ["# BLASTN 2.17.0+\n"]
    for query_id in queries:
        text.append(f"# Query: {query_id}\n# Database: synthetic\n# Fields: {', '.join(BLAST_FIELDS)}\n")
        text.append(f"# {len(lines.get(query_id, []))} hits found\n")
        text.extend(line + "\n" for line in lines.get(query_id, []))
    text.append(f"# BLAST processed {len(queries)} queries\n")
    return "".join(text)


def design_hits(sequences: dict[str, str], search: dict[str, str], logical: list[dict]) -> dict[str, list[str]]:
    """Synthetic HSPs per searched sequence (canonical orientation)."""
    parent = {query_id: canonical_id(sequence) for query_id, sequence in sequences.items()}
    controls = {(row["parent_id"], row["draw_index"]): row["sequence_id"] for row in logical if row["role"] == "control"}
    lines = {}

    def add(sequence_id, accession, start, plus=True, **options):
        record = RECORDS.get(accession, dict(length=300 if accession.startswith("pdb") else 100000, taxids="0", organism="N/A"))
        lines.setdefault(sequence_id, []).append(hsp(sequence_id, search[sequence_id], accession, record["length"], start, plus,
                                                     taxid=record.get("taxids", "0") or "0", name=record.get("organism", "N/A"), **options))
    nb1, k40, bgi = parent["ONT_NB01_forward"], parent["ONT_K40_NB01_forward"], parent["BGI_filter_adapter"]
    add(nb1, "CP000001.1", 120)
    add(nb1, "CP000002.1", 5000, plus=False)
    add(nb1, "OX000007.1", 60000)
    add(nb1, "MN000003.1", 10)
    add(k40, "CP000001.1", 112)
    add(parent["ONT_NB02_forward"], "OX000007.1", 70000, mismatches=(12,))
    add(parent["ONT_PCR_adapter"], "OX000007.1", 80000)
    adapter_a, adapter_b = parent["ILLUMINA_I7ADP_A_forward"], sequences["ILLUMINA_I7ADP_B_forward"]
    add(adapter_a, "OX000007.1", 90000)
    # The B index read in an A adapter: one substitution, a scaffold-preserving swap to another catalogued index.
    swapped = adapter_b if search[adapter_a] == sequences["ILLUMINA_I7ADP_A_forward"] else reverse_complement(adapter_b)
    add(adapter_a, "KX000008.1", 100, subject=swapped)
    add(parent["ILLUMINA_I7ADP_B_forward"], "KX000008.1", 5000)
    add(parent["ILLUMINA_NEXTERA_motif"], "CP000001.1", 20000)
    add(parent["ILLUMINA_NEXTERA_motif"], "OQ000005.1", 1, plus=False)
    add(parent["ILLUMINA_TRUSEQ_adapter"], "LC000006.1", 900)
    add(parent["ILLUMINA_TRUSEQ_adapter"], "pdb|4FCY|C", 100)
    add(bgi, "AB000004.1", 39950, plus=False)
    add(bgi, "CP000001.1", 30000, mismatches=(15,))
    add(parent["PACBIO_SMRTbell_forward"], "OX000007.1", 110000)
    c2 = parent["PACBIO_C2_primer_forward"]
    for number in range(1, TARGETS + 1):
        add(c2, f"CAP0000{number}.1", 1000 * number, qstart=1, qend=20)
        add(controls[c2, 2], f"CAP0000{number}.1", 2000 * number, qstart=1, qend=20)
    add(controls[nb1, 3], "CT000009.1", 500)
    add(controls[bgi, 1], "CP000001.1", 40000, mismatches=(10,))
    add(controls[parent["ILLUMINA_TRUSEQ_adapter"], 4], "OX000007.1", 3000, qstart=5, qend=20)
    return lines


def write_search(root: Path, design: Path, policy: dict, hits: dict[str, list[str]]) -> tuple[Path, list[dict], dict]:
    search = root / "search"
    search.mkdir()
    queries = read_fasta(design / "search_queries.fasta")
    lengths = sorted({len(sequence) for sequence in queries.values()})
    settings = dict(task="blastn-short", word_size=7, evalue=10.0, dust="no", soft_masking="false", strand="both", reward=1, penalty=-3,
                    gapopen=5, gapextend=2, max_target_seqs=TARGETS, searchsp=None)
    request = dict(settings=settings, mode="working", queries_sha256=digest(design / "search_queries.fasta"),
                   length_calibration=dict(lengths=policy["detectability"]))
    save_json(search / "request.json", request)
    search_id = object_digest(request)
    statuses, batches = [], []
    for number, length in enumerate(lengths, 1):
        batch_queries = {name: sequence for name, sequence in sorted(queries.items()) if len(sequence) == length}
        directory = search / f"batch_{number:05d}.attempt_00001"
        directory.mkdir()
        save_fasta(directory / "queries.fasta", batch_queries)
        (directory / "blast.tsv").write_text(blast_output(batch_queries, hits), encoding="utf-8")
        for name in ("stderr.txt", "stdout.txt"):
            (directory / name).write_text("", encoding="utf-8")
        save_json(directory / "command.json", dict(command=["blastn"], search_id=search_id))
        rows = summarise_batch(directory / "blast.tsv", batch_queries, f"batch_{number:05d}", TARGETS)
        write_table(directory / "status.tsv", rows, STATUS_SCHEMA)
        save_json(directory / "complete.json", dict(search_id=search_id, batch_id=f"batch_{number:05d}", query_ids=list(batch_queries),
                                                    files={name: digest(directory / name) for name in ("queries.fasta", "blast.tsv", "stderr.txt",
                                                                                                      "stdout.txt", "status.tsv", "command.json")}))
        statuses.extend(rows)
        batches.append(dict(directory=directory.name, completion_sha256=digest(directory / "complete.json")))
    statuses.sort(key=lambda row: row["query_id"])
    write_table(search / "status.tsv", statuses, STATUS_SCHEMA)
    save_json(search / "complete.json", dict(search_id=search_id, all_queries_completed=True, batches=batches))
    return search, statuses, request


def subject_metadata(accession: str) -> dict:
    record = RECORDS[accession]
    return dict(raw_id=accession, accession_version=accession, length=record["length"], topology=record["topology"],
                record_type=record["record_type"], taxids=record["taxids"], assembly_id="", project_id=record["project"],
                metadata_source="synthetic ESummary", metadata_date="2026-11-01", molecule=record["molecule"], origin=record["origin"],
                targeted="unknown", organelle="no", context=record["context"])


def write_analysis(root: Path, design: Path, search: Path, statuses: list[dict], policy: dict, request: dict) -> Path:
    queries = read_fasta(design / "search_queries.fasta")
    logical = read_table(design / "logical_queries.tsv", LOGICAL_SCHEMA)
    complete = json.loads((search / "complete.json").read_text())
    events = []
    for batch in complete["batches"]:
        directory = search / batch["directory"]
        batch_queries = read_fasta(directory / "queries.fasta")
        for ordinal, row in enumerate(iter_hsps(directory / "blast.tsv", batch_queries), 1):
            if not tolerant_hsp(row):
                continue
            accepted = event(row, {}, queries[row["qseqid"]], directory.name + "/blast.tsv", ordinal)
            accepted.update(annotation_defaults(accepted))
            accepted["parent_id"] = ""
            subject = subject_metadata(accepted["accession_version"]) if accepted["accession_version"] in RECORDS else None
            events.append(annotate_record(accepted, subject))
    taxonomy = []
    for accession in sorted({row["accession_version"] for row in events if row["accession_version"] in RECORDS}):
        record = RECORDS[accession]
        taxonomy.append(dict(accession_version=accession, rank="domain", blast_taxids=record["taxids"], record_taxids=record["taxids"],
                             resolved_taxids=record["taxids"], taxon_id="2", taxon_name="Synthetic", status="resolved", source_status="agree",
                             context="vector" if record["origin"] == "vector" else record["context"], resolution_evidence="{}",
                             snapshot_sha256="0" * 64))
    definition = json.loads((design / "panel.json").read_text())
    provenance = dict(panel_sha256=digest(design / "panel.json"), search_completion_sha256=digest(search / "complete.json"),
                      search_mode="working", search_id=complete["search_id"], source_search=str(search.resolve()), draws_per_parent=DRAWS,
                      repeated_control_seed_values=definition["repeated_control_seed_values"], control_algorithm=definition["control_algorithm"],
                      search_settings=request["settings"], length_calibration=request["length_calibration"],
                      length_calibration_source=dict(kind="runtime_request"), study_binding=None,
                      search_request_sha256=digest(search / "request.json"), database_acceptance=None)
    destination = root / "taxonomy"
    write_completed_analysis(events, logical, statuses, policy, destination, provenance,
                             {"accession_taxonomy.tsv": (taxonomy, TAXONOMY_SCHEMA)})
    return destination


def write_metadata(root: Path) -> Path:
    metadata = root / "metadata"
    metadata.mkdir()
    result = {"uids": [str(index) for index in range(1, len(RECORDS) + 1)]}
    for index, (accession, record) in enumerate(sorted(RECORDS.items()), 1):
        result[str(index)] = dict(accessionversion=accession, projectid=record["project"], tech="barcode" if accession in TARGETED else "",
                                  title=record["title"], slen=record["length"])
    save_json(metadata / "esummary_0001.json", dict(result=result))
    save_json(metadata / "metadata.json", dict(kind="synthetic"))
    return metadata


def write_insdc(root: Path) -> Path:
    insdc = root / "insdc"
    insdc.mkdir()
    rows = [dict(accession_version=accession, status="retrieved", keywords="BARCODE" if accession in TARGETED else "", pcr_primers="",
                 esummary_tech="", targeted_insdc="yes" if accession in TARGETED else "no", evidence="synthetic", archive_file="efetch_0000.txt")
            for accession in sorted(RECORDS)]
    write_table(insdc / "insdc_annotation.tsv", rows, ANNOTATION_SCHEMA)
    (insdc / "insdc.json").write_text(json.dumps(dict(annotation_sha256=digest(insdc / "insdc_annotation.tsv"))) + "\n", encoding="utf-8")
    return insdc


def write_positional(root: Path, analysis: Path, metadata: Path, insdc: Path) -> Path:
    from database_residue_screen.insdc import read_annotation
    positional = root / "positional"
    positional.mkdir()
    rows = positional_sensitivity(analysis, projects_from_metadata(metadata), read_annotation(insdc))
    write_table(positional / "positional_sensitivity.tsv", rows, SENSITIVITY_SCHEMA)
    return positional


def write_scan(root: Path, analysis: Path, statuses: list[dict]) -> tuple[Path, Path]:
    exact, validation = root / "exact_scan", root / "scan_validation"
    exact.mkdir()
    validation.mkdir()
    lengths = {row["sequence_id"]: row["length"] for row in read_table(analysis / "logical_queries.tsv", LOGICAL_SCHEMA)}
    lines = (analysis / "events.tsv").read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    strict = {}
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        if row["strict"] == "yes":
            strict[row["query_id"]] = strict.get(row["query_id"], 0) + 1
    rows, discordant = [], []
    for status in statuses:
        length = lengths[status["query_id"]]
        if length < 23:
            continue
        count = strict.get(status["query_id"], 0)
        capped = status["target_limit_reached"] == "yes"
        rows.append(dict(query_id=status["query_id"], target_limit_reached="yes" if capped else "no", blast_strict_hsps=count,
                         scan_occurrences=count + capped, matched=count, blast_only=0, scan_only=int(capped), scan_only_target_cap=int(capped),
                         scan_only_overlap=0, scan_only_unexplained=0, blast_strict_accessions=count, scan_accessions=count + capped,
                         accessions_equal="no" if capped else "yes", length=length))
        if capped:
            discordant.append(dict(query_id=status["query_id"], accession_version="CAP00099.1", start=1, end=length, strand="+",
                                   side="scan_only", classification="target_cap"))
    write_rows(validation / "scan_concordance_by_query.tsv", rows)
    write_rows(validation / "scan_discordant_occurrences.tsv", discordant,
               ["query_id", "accession_version", "start", "end", "strand", "side", "classification"])
    totals = {name: sum(row[name] for row in rows) for name in ("blast_strict_hsps", "scan_occurrences", "matched", "blast_only", "scan_only",
                                                                "scan_only_target_cap", "scan_only_overlap", "scan_only_unexplained")}
    save_json(validation / "scan_summary.json", dict(scanned_search_sequences=len(rows), minimum_length=23,
                                                     capped_searches=sum(row["target_limit_reached"] == "yes" for row in rows),
                                                     uncapped_searches_with_unequal_accessions=0, accepted=True, **totals))
    save_json(exact / "scan.json", dict(kind="exact_occurrence_scan_v1", records=SEQUENCES, bases=BASES, expected_records=SEQUENCES,
                                        expected_bases=BASES, totals_match_database=True, minimum_length=23, scanned_queries=len(rows),
                                        deflines=SEQUENCES, occurrences=totals["scan_occurrences"], volumes=[]))
    return exact, validation


def write_rows(path: Path, rows: list[dict], fields=None) -> None:
    fields = list(fields or rows[0])
    path.write_text("\t".join(fields) + "\n" + "".join("\t".join(str(row[name]) for name in fields) + "\n" for row in rows), encoding="utf-8")


def write_target_limit(root: Path, analysis: Path) -> Path:
    directory = root / "target_limit"
    directory.mkdir()
    outcomes = [row for row in read_table(analysis / "outcomes.tsv", OUTCOME_SCHEMA) if row["role"] == "real"]
    seen, rows = set(), []
    for row in sorted(outcomes, key=lambda row: (row["parent_id"], row["endpoint"], row["population"])):
        key = row["parent_id"], row["endpoint"], row["population"]
        if row["population"] != "all_records":
            continue
        if key in seen:
            continue
        seen.add(key)
        for limit in (500, 2000):
            retained = row["retained_accessions"]
            if limit == 500 and retained >= 3:
                retained -= 1
            completion = ("success_hits" if retained or row["unresolved_subjects"] else "success_zero") if row["completion"] == "successful" \
                else "failed"
            rows.append(dict(search=f"limit_{limit}", target_limit=limit, parent_id=row["parent_id"], endpoint=row["endpoint"],
                             population=row["population"], completion=completion, retained_accessions=retained,
                             unresolved_subjects=row["unresolved_subjects"],
                             target_limit_reached="yes" if row["retrieval"] == "target_limit_reached" else "no"))
    write_rows(directory / "target_limit.tsv", rows)
    return directory


def write_cap_review(root: Path, search: Path) -> Path:
    directory = root / "cap_review"
    directory.mkdir()
    write_table(directory / "cap_review.tsv", review_saturated(search), REVIEW_SCHEMA)
    return directory


def write_inference(root: Path, analysis: Path, design: Path) -> Path:
    specification = precision_specification([5, 10, 15], [500, 2000], threads=60)
    specification["panel_sha256"] = digest(design / "panel.json")
    spec_dir = root / "inference_spec"
    spec_dir.mkdir()
    save_json(spec_dir / "specification.json", specification)
    destination = root / "inference"
    assess_precision(analysis, spec_dir / "specification.json", destination)
    (destination / "specification.json").write_text((spec_dir / "specification.json").read_text())
    return destination


def write_benchmarks(root: Path, sequences: dict[str, str]) -> Path:
    directory = root / "benchmark"
    directory.mkdir()
    parents = sorted({canonical(sequence) for sequence in sequences.values()})
    for database, factor in (("core_nt_div10", 0.1), ("core_nt", 1.0), ("core_nt_x10", 10.0)):
        folder = directory / f"benchmark_{database}"
        (folder / "fixtures" / "A_001").mkdir(parents=True)
        (folder / "searches" / "A_001" / "permissive_diagnostic").mkdir(parents=True)
        bases = int(BASES * factor)
        calibration = {row["query_length"]: row for row in calibrated([len(parent) for parent in parents], bases)}
        spikes = []
        for parent in parents:
            row = calibration[len(parent)]
            for orientation in ("+", "-"):
                recovered = int(row["exact_evalue"] <= 10)
                spikes.append(dict(parent_id=canonical_id(parent), aliases="synthetic", groups="synthetic", conditions="synthetic",
                                   length=len(parent), context="internal", orientation=orientation, subject_id=f"FIX{len(spikes) + 1:07d}.1",
                                   local_fixture_recovered=1, production_calibrated_recovered=recovered, permissive_diagnostic_recovered=1,
                                   local_fixture_raw_hsps=1, production_calibrated_raw_hsps=recovered, permissive_diagnostic_raw_hsps=1,
                                   exact_score=len(parent), calibrated_exact_evalue=row["exact_evalue"],
                                   effective_search_space=row["effective_search_space"],
                                   classification="recovered" if recovered else "evalue_censoring_confirmed",
                                   zero_interpretation="synthetic"))
        metrics = []
        lengths = sorted({len(parent) for parent in parents})
        for arm in ("local_fixture", "production_calibrated", "permissive_diagnostic"):
            for endpoint in ("strict", "tolerant"):
                per_length = []
                for length in lengths:
                    positive = sum(len(parent) == length for parent in parents) * 2
                    recovered = positive if arm != "production_calibrated" or calibration[length]["exact_evalue"] <= 10 else 0
                    per_length.append((length, positive, recovered))
                positive_total = sum(item[1] for item in per_length)
                recovered_total = sum(item[2] for item in per_length)
                metrics.append(dict(search_space=arm, endpoint=endpoint, population="all_lengths", error_class="all", length=0,
                                    positive_pairs=positive_total, negative_pairs=40, true_positive=recovered_total,
                                    false_negative=positive_total - recovered_total, false_positive=0, true_negative=40,
                                    recall=recovered_total / positive_total, precision=1.0 if recovered_total else "not_estimated"))
                for length, positive, recovered in per_length:
                    metrics.append(dict(search_space=arm, endpoint=endpoint, population="length", error_class="all", length=length,
                                        positive_pairs=positive, negative_pairs=2, true_positive=recovered, false_negative=positive - recovered,
                                        false_positive=0, true_negative=2, recall=recovered / positive,
                                        precision=1.0 if recovered else "not_estimated"))
        query, subject = random_sequence(40), random_sequence(200)
        subject = subject[:50] + query[:30] + subject[80:]
        missed = random_sequence(40)
        subjects = {"FIX0000001.1": subject, "FIX0000002.1": random_sequence(60)}
        save_fasta(folder / "fixtures" / "A_001" / "queries.fasta", {"q1": query, "q2": missed})
        save_fasta(folder / "fixtures" / "A_001" / "subjects.fasta", subjects)
        pair_fields = ("shard_id", "panel", "arm", "query_id", "subject_id", "length", "complexity", "error_class", "context", "raw_hsps",
                       "strict_truth", "strict_independent", "strict_accepted", "strict_classification", "tolerant_truth", "tolerant_independent",
                       "tolerant_accepted", "tolerant_classification")
        pairs = [dict(shard_id="A_001", panel="A", arm="permissive_diagnostic", query_id="q1", subject_id="FIX0000001.1", length=40,
                      complexity="other", error_class="clip_both_below", context="internal", raw_hsps=1, strict_truth=0, strict_independent=0,
                      strict_accepted=0, strict_classification="correct_nonpositive", tolerant_truth=1, tolerant_independent=0, tolerant_accepted=0,
                      tolerant_classification="retrieval_or_representation_failure"),
                 dict(shard_id="A_001", panel="A", arm="permissive_diagnostic", query_id="q2", subject_id="FIX0000002.1", length=40,
                      complexity="other", error_class="cross_pair", context="cross_pair", raw_hsps=0, strict_truth=0, strict_independent=0,
                      strict_accepted=0, strict_classification="correct_nonpositive", tolerant_truth=1, tolerant_independent=0, tolerant_accepted=0,
                      tolerant_classification="retrieval_or_representation_failure"),
                 dict(shard_id="A_001", panel="A", arm="permissive_diagnostic", query_id="q1", subject_id="FIX0000002.1", length=40,
                      complexity="other", error_class="exact", context="internal", raw_hsps=1, strict_truth=1, strict_independent=1,
                      strict_accepted=1, strict_classification="recovered", tolerant_truth=1, tolerant_independent=1, tolerant_accepted=1,
                      tolerant_classification="recovered")]
        witnesses = [dict(shard_id="A_001", query_id="q1", subject_id="FIX0000001.1", endpoint="tolerant", qstart=1, qend=37, sstart=51,
                          send=87, strand="+", span=37, matches=35, mismatches=2),
                     dict(shard_id="A_001", query_id="q2", subject_id="FIX0000002.1", endpoint="tolerant", qstart=1, qend=37, sstart=1, send=37,
                          strand="+", span=37, matches=33, mismatches=4)]
        for path, rows in ((folder / "benchmark_pair_results.tsv", pairs), (folder / "truth_witnesses.tsv", witnesses),
                           (folder / "metrics.tsv", metrics), (folder / "short_query_detectability.tsv", spikes)):
            write_rows(path, rows)
        line = hsp("q1", query, "FIX0000001.1", len(subject), 51, qstart=1, qend=30)
        (folder / "searches" / "A_001" / "permissive_diagnostic" / "blast.tsv").write_text(
            blast_output({"q1": query, "q2": missed}, {"q1": [line]}), encoding="utf-8")
        save_json(folder / "benchmark.json", dict(contract="panel_benchmark_synthetic", status="complete", predicate_disagreements=0,
                                                  independent_pair_universe=sum(row["positive_pairs"] + row["negative_pairs"] for row in metrics
                                                                                if row["population"] == "all_lengths" and row["search_space"] == "local_fixture"
                                                                                and row["endpoint"] == "strict")))
    return directory


def write_record_cases(root: Path, analysis: Path, design: Path, sequences: dict[str, str]) -> Path:
    directory = root / "record_cases"
    (directory / "cases" / "terminal_attempt_1").mkdir(parents=True)
    search = read_fasta(design / "search_queries.fasta")
    cases, segments = [], []
    for number, (category, query_id, accession, start, plus) in enumerate((("terminal", "BGI_filter_adapter", "AB000004.1", 39950, False),
                                                                           ("synthetic_or_vector", "ILLUMINA_NEXTERA_motif", "OQ000005.1", 1, False)), 1):
        parent = canonical_id(sequences[query_id])
        searched = search[parent]
        record = RECORDS[accession]
        lower, upper = start, start + len(searched) - 1
        flank_start, flank_end = max(1, lower - 50), min(record["length"], upper + 50)
        left, right = random_sequence(lower - flank_start), random_sequence(flank_end - upper)
        match = searched if plus else reverse_complement(searched)
        forward = left + match + right
        display_left, display_right = (left, right) if plus else (reverse_complement(right), reverse_complement(left))
        if number == 1:
            (directory / "cases" / "terminal_attempt_1" / "retrieved.fasta").write_text(f">{accession}:{flank_start}-{flank_end}\n{forward}\n")
        cases.append(dict(case_id=f"case_{number:02d}", category=category, selection_subcategory=category, candidate_id=f"candidate_{number}",
                          accession_version=accession, query_id=query_id, parent_id=parent, physical_aliases=query_id,
                          group="BGI/DNBSEQ" if query_id.startswith("BGI") else "Illumina",
                          condition="illumina_bgi" if query_id.startswith("BGI") else "illumina_nextera", query_sequence=sequences[query_id],
                          alias_orientation="+" if sequences[query_id] == searched else "-", qstart=1, qend=len(searched),
                          sstart=lower if plus else upper, send=upper if plus else lower, record_length=record["length"],
                          distance_to_end=min(lower - 1, record["length"] - upper), qseq=searched, sseq=searched,
                          display_segment_roles=";".join(["scaffold"] * len(searched)), lower=lower, upper=upper, flank_start=flank_start,
                          flank_end=flank_end, forward_sequence=forward, display_left_flank=display_left, display_right_flank=display_right,
                          display_orientation="+" if plus else "-", left_flank_clipped=flank_start == 1 and lower - 1 < 50,
                          right_flank_clipped=flank_end == record["length"] and record["length"] - upper < 50, molecule=record["molecule"],
                          topology=record["topology"], record_type=record["record_type"], origin=record["origin"], context=record["context"],
                          taxids=record["taxids"], organism=record["organism"], title=record["title"], linear_eligibility="included",
                          evidence_file=f"cases/{category}_attempt_1/case_evidence.json"))
        segments.append(dict(case_id=f"case_{number:02d}", query_id=query_id, segment_id=f"{query_id}:segment_01", start=1, end=len(searched),
                             role="scaffold", assignment_status="declared", covered_bases=len(searched), matching_bases=len(searched),
                             mismatches=0, component_length=len(searched), coverage_fraction=1.0, complete_segment=1, support_state="complete_exact",
                             covered_start=1, covered_end=len(searched), mismatch_positions=""))
    write_rows(directory / "record_cases.tsv", cases)
    write_rows(directory / "record_case_segments.tsv", segments)
    status = [dict(category="terminal", selection_subcategory="terminal", eligible_candidates=3, attempted_accessions=1, selected_candidate="candidate_1",
                   status="selected"),
              dict(category="synthetic_or_vector", selection_subcategory="synthetic_or_vector", eligible_candidates=1, attempted_accessions=1,
                   selected_candidate="candidate_2", status="selected"),
              dict(category="internal_or_ambiguous", selection_subcategory="internal", eligible_candidates=0, attempted_accessions=0,
                   selected_candidate="", status="unavailable_no_verified_candidate")]
    write_rows(directory / "case_category_status.tsv", status)
    return directory


def write_census(root: Path) -> Path:
    census = root / "census"
    census.mkdir()
    write_rows(census / "census_windows.tsv", [
        dict(family="ont_native_barcode_24", physical_id="ONT_NB01_forward", accession_version="CP000001.1", strand="+", status="retrieved",
             census_class="none"),
        dict(family="ont_native_barcode_24", physical_id="ONT_NB01_forward", accession_version="OX000007.1", strand="+",
             status="insufficient_flank", census_class="insufficient_flank")])
    write_rows(census / "census_summary.tsv", [
        dict(family="ont_native_barcode_24", rule="exact", census_class="none", windows=1, accessions=1, family_windows=2, family_accessions=2),
        dict(family="ont_native_barcode_24", rule="exact", census_class="insufficient_flank", windows=1, accessions=1, family_windows=2,
             family_accessions=2)])
    save_json(census / "census.json", dict(windows=2, status_counts=dict(retrieved=1, insufficient_flank=1), integrity_complete=True,
                                           files={name: digest(census / name) for name in ("census_windows.tsv", "census_summary.tsv")}))
    return census


def build(root: Path) -> dict[str, Path]:
    """Write every upstream stage output of a small synthetic study below ``root`` (same content on every call)."""
    RANDOM.seed(20261109)
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    definition, sequences = panel_definition()
    panel = root / "panel"
    panel.mkdir()
    write_rows(panel / "panel_definition.tsv", definition)
    save_fasta(panel / "queries.fasta", sequences)
    metadata_fields = ("query_id", "group", "condition", "construct_class", "family_id", "source_id", "provenance_status",
                       "architecture_status")
    (panel / "queries.tsv").write_text("\t".join(metadata_fields) + "\n" + "".join(
        "\t".join((row["query_id"], row["group"], row["condition"], row["construct_class"], canonical_id(sequences[row["query_id"]]),
                   "synthetic", "synthetic", "synthetic" if row["construct_class"] == "composite" else "not_applicable")) + "\n"
        for row in definition), encoding="utf-8")
    controls = root / "controls"
    generate_controls(panel / "queries.fasta", controls, DRAWS)
    design = root / "design"
    prepare_panel(panel / "queries.fasta", panel / "queries.tsv", controls, design, DRAWS)
    queries = read_fasta(design / "search_queries.fasta")
    policy = json.loads((design / "analysis_policy.json").read_text())
    policy["detectability"] = calibrated(len(sequence) for sequence in queries.values())
    logical = read_table(design / "logical_queries.tsv", LOGICAL_SCHEMA)
    search, statuses, request = write_search(root, design, policy, design_hits(sequences, queries, logical))
    analysis = write_analysis(root, design, search, statuses, policy, request)
    metadata = write_metadata(root)
    insdc = write_insdc(root)
    exact, validation = write_scan(root, analysis, statuses)
    outputs = dict(panel=panel, controls=controls, design=design, search=search, taxonomy=analysis, metadata=metadata, insdc=insdc,
                   positional=write_positional(root, analysis, metadata, insdc), exact_scan=exact, scan_validation=validation,
                   cap_review=write_cap_review(root, search), target_limit=write_target_limit(root, analysis),
                   inference=write_inference(root, analysis, design), benchmark=write_benchmarks(root, sequences),
                   record_cases=write_record_cases(root, analysis, design, sequences), census=write_census(root))
    validation_dir = root / "validation"
    validation_dir.mkdir()
    save_json(validation_dir / "validation.json", dict(valid=True, result_contract="reference_screen_v2"))
    outputs["validation"] = validation_dir
    return outputs
