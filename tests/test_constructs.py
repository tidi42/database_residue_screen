import contextlib
import csv
import hashlib
import io
import json
import random
import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from database_residue_screen.cli import main
from database_residue_screen.constructs import (DEFINITION_SCHEMA, RULES, Library, build_panel, dorado_kit, fetch_sources, locate_oligo,
                                                metadata_schema, parse_bracket_parameters, parse_dorado, parse_illumina_formulas,
                                                parse_illumina_indexes, rule_documented_oligo)
from database_residue_screen.contracts import ContractError, read_table, schema, write_table
from database_residue_screen.controls import generate_controls
from database_residue_screen.panels import prepare_panel
from database_residue_screen.sequences import canonical_id, read_fasta, reverse_complement

RANDOM = random.Random(20261109)


def made_up(length: int) -> str:
    return "".join(RANDOM.choice("ACGT") for _ in range(length))


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


# Synthetic stand-ins for every source format; all sequences are random and unrelated to any vendor product.
S = {name: made_up(length) for name, length in dict(
    nb1=24, nb2=24, front14=14, rear8=8, bc1=24, rbk2=24, doc_rb2=24, rapid_front=10, rapid_rear=20, p1=22, p2=26, p3=12,
    i7=10, i5=10, i7b=10, i5b=10, ud_i7_front=12, ud_i7_rear=9, ud_i5_front=13, ud_i5_rear=8, nx_i7_front=11, nx_i7_rear=7,
    nx_i5_front=14, nx_i5_rear=9, hn701=8, ehns501=8, trim=19, tn5a=20, tn5b=21, pcr1=20, pcr2=21, pe1=24, pe2=18, ra5=26,
    mgi1_head=14, mgi1_tail=9, mgi2=30, b1=16, b2=16, b3=16, hairpin=30, primer=20).items()}

S["ra5"] = S["ra5"][:-1] + "T"  # the RNA fixture needs at least one U

SPECIAL_GLYPHS = {"【": ("uni3010", 0x3010), "】": ("uni3011", 0x3011), "′": ("uni2032", 0x2032)}


def write_pdf(path: Path, pages: list[list[str]]) -> None:
    """Minimal text PDF (Courier, one line per string) that pdftotext extracts, including 【 】 and ′."""
    codes = {char: 0x80 + index for index, char in enumerate(SPECIAL_GLYPHS)}

    def encode(line):
        return b"".join(bytes([codes[char]]) if char in codes else
                        char.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)").encode("ascii") for char in line)

    objects = []

    def add(body: bytes) -> int:
        objects.append(body)
        return len(objects)

    cmap = (f"/CIDInit /ProcSet findresource begin 12 dict begin begincmap /CMapName /Fixture-UCS def /CMapType 2 def\n"
            f"1 begincodespacerange <00> <FF> endcodespacerange\n{len(codes)} beginbfchar\n"
            + "".join(f"<{codes[char]:02X}> <{unicode:04X}>\n" for char, (_, unicode) in SPECIAL_GLYPHS.items())
            + "endbfchar endcmap CMapName currentdict /CMap defineresource pop end end").encode("ascii")
    to_unicode = add(b"<< /Length %d >>\nstream\n" % len(cmap) + cmap + b"\nendstream")
    differences = " ".join(f"{codes[char]} /{name}" for char, (name, _) in SPECIAL_GLYPHS.items())
    font = add(f"<< /Type /Font /Subtype /Type1 /BaseFont /Courier /Encoding << /Type /Encoding /BaseEncoding /WinAnsiEncoding "
               f"/Differences [{differences}] >> /ToUnicode {to_unicode} 0 R >>".encode("ascii"))
    contents = []
    for lines in pages:
        stream = b"BT /F1 9 Tf 11 TL 36 806 Td\n" + b"".join(b"(" + encode(line) + b") Tj T*\n" for line in lines) + b"ET"
        contents.append(add(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream"))
    parent = len(objects) + len(pages) + 1
    kids = [add(f"<< /Type /Page /Parent {parent} 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 {font} 0 R >> >> "
                f"/Contents {content} 0 R >>".encode("ascii")) for content in contents]
    add(f"<< /Type /Pages /Kids [{' '.join(f'{kid} 0 R' for kid in kids)}] /Count {len(kids)} >>".encode("ascii"))
    catalog = add(f"<< /Type /Catalog /Pages {parent} 0 R >>".encode("ascii"))
    data, offsets = bytearray(b"%PDF-1.4\n"), []
    for number, body in enumerate(objects, 1):
        offsets.append(len(data))
        data += b"%d 0 obj\n" % number + body + b"\nendobj\n"
    xref = len(data)
    data += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objects) + 1) + b"".join(b"%010d 00000 n \n" % offset for offset in offsets)
    data += b"trailer\n<< /Size %d /Root %d 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objects) + 1, catalog, xref)
    path.write_bytes(bytes(data))


def write_xlsx(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Workbook whose header uses shared strings and whose data rows use inline strings."""
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    shared = "".join(f"<si><t>{escape(value)}</t></si>" for value in header)
    cells = ['<row r="1">' + "".join(f'<c r="{column}1" t="s"><v>{index}</v></c>' for index, column in enumerate("ABCDEF")) + "</row>"]
    cells += [f'<row r="{number}">' + "".join(f'<c r="{column}{number}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'
                                             for column, value in zip("ABCDEF", row)) + "</row>" for number, row in enumerate(rows, 2)]
    with zipfile.ZipFile(path, "w") as book:
        book.writestr("xl/workbook.xml", f'<workbook xmlns="{main_ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/'
                                         'relationships"><sheets><sheet name="adapters" sheetId="4" r:id="rId9"/></sheets></workbook>')
        book.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                                                    '<Relationship Id="rId9" Target="worksheets/sheet7.xml"/></Relationships>')
        book.writestr("xl/sharedStrings.xml", f'<sst xmlns="{main_ns}">{shared}</sst>')
        book.writestr("xl/worksheets/sheet7.xml", f'<worksheet xmlns="{main_ns}"><sheetData>{"".join(cells)}</sheetData></worksheet>')


def html(body: str) -> str:
    return f"<html><head><title>Fixture</title></head><body>{body}</body></html>\n"


def write_sources(directory: Path) -> None:
    rc = reverse_complement
    directory.mkdir(parents=True)
    (directory / "ont.html").write_text(html(
        "<h2>Native barcode sequences</h2><p>Barcode flanking sequences:</p>"
        f"<p>Forward sequence: 5' - {S['front14'][-8:]} - barcode - {S['rear8']} - 3'</p>"
        f"<p>Reverse sequence: 5' - {rc(S['rear8'][:7])} - barcode - {rc(S['front14'])} - 3'</p>"
        "<table><tr><td>Component</td><td>Forward sequence</td><td>Reverse sequence</td></tr>"
        + "".join(f"<tr><td>NB0{n}</td><td>{S[f'nb{n}']}</td><td>{rc(S[f'nb{n}'])}</td></tr>" for n in (1, 2)) + "</table>"
        "<h2>Rapid barcoding kits</h2><p>Rapid Barcoding Kit 96 V14 (SQK-RBK114.96): RB01-96</p>"
        f"<p>Barcode flanking sequence: 5' - {S['rapid_front']} - barcode - {S['rapid_rear']} - 3'</p>"
        f"<table><tr><td>BC01 / RB01</td><td>{S['bc1']}</td></tr><tr><td>BC02 / RB02</td><td>{S['doc_rb2']}</td></tr></table>"
        f"<p>Ligation Adapter bottom strand: 5'-{S['p1']}{S['p3']}-3'</p>"), encoding="utf-8")
    (directory / "dorado.cpp").write_text(
        f'const std::string NB_1st_FRONT = "{S["front14"]}";\nconst std::string NB_1st_REAR = "{S["rear8"]}";\n'
        f'const std::string NB_2nd_FRONT = "{S["front14"]}";\nconst std::string NB_2nd_REAR = "{S["rear8"][:7]}";\n'
        f'const std::string RBK4_kit14_FRONT = "{S["rapid_front"]}";\nconst std::string RBK4_kit14_REAR = "{S["rapid_rear"]}";\n'
        'const std::vector<std::string> NB_1_96 = {"NB01", "NB02"};\n'
        '// RBK_1_96 is the generic list except for position 2.\nconst std::vector<std::string> RBK_1_96 = {"BC01", "RBK02"};\n'
        'const KitInfo kit_nb96 = {\n        "NB96", true, true, false, false, false, NB_1st_FRONT, "",\n'
        '        "", NB_1st_REAR, NB_2nd_FRONT, "", "", NB_2nd_REAR, NB_1_96, NB_1_96,\n        {}, {}, DEFAULT_PARAMS,\n};\n'
        'const std::unordered_map<std::string, KitInfo> kit_info_map = {\n        {"SQK-NBD114-96", kit_nb96},\n'
        '        {"SQK-RBK114-96",\n         {\n                 "RBK096_kit14", false, false, false, false, false,\n'
        '                 RBK4_kit14_FRONT, "", "", RBK4_kit14_REAR, "", "", "", "",\n'
        '                 RBK_1_96, {}, {}, {}, RBK114_PARAMS,\n         }},\n};\n'
        'const std::unordered_map<std::string, std::string> barcodes = {\n'
        f'        {{"NB01", "{S["nb1"]}"}}, {{"NB02", "{S["nb2"]}"}},\n        {{"BC01", "{S["bc1"]}"}}, {{"RBK02", "{S["rbk2"]}"}},\n}};\n',
        encoding="utf-8")
    (directory / "porechop.py").write_text(
        "class Adapter(object):\n    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):\n"
        "        self.name = name\n\n\n"
        f"ADAPTERS = [Adapter('Kit A',\n                    start_sequence=('Kit_A_top', '{S['p1']}'),\n"
        f"                    end_sequence=('Kit_A_bottom', '{S['p2']}')),\n"
        f"            Adapter('Barcode 1', start_sequence=('BC_computed', 'AC' + 'GT'))]\n", encoding="utf-8")
    (directory / "ud.html").write_text(html(
        f"<p>Index 1 (i7) Adapters</p><p>{S['ud_i7_front']}[i7]{S['ud_i7_rear']}</p>"
        f"<p>Index 2 (i5) Adapters</p><p>{S['ud_i5_front']}[i5]{S['ud_i5_rear']}</p><table>"
        f"<tr><td>UDP0001</td><td>{S['i7']}</td><td>{rc(S['i7'])}</td><td>{S['i5']}</td><td>{S['i5']}</td><td>{rc(S['i5'])}</td></tr>"
        f"<tr><td>UDP0002V3</td><td>{S['i7b']}</td><td>{rc(S['i7b'])}</td><td>{S['i5b']}</td><td>{S['i5b']}</td><td>{rc(S['i5b'])}</td>"
        "</tr><tr><td>A</td><td>UDP0001</td><td>UDP0002V3</td></tr></table>"), encoding="utf-8")
    (directory / "nextera.html").write_text(html(
        "<h1>Sequences for Illumina Prep and Nextera Kits</h1><h2>Adapter Trimming</h2>"
        f"<p>The following sequence is used for Read 1 and Read 2 adapter trimming.</p><p>{S['trim']}</p>"
        f"<h2>Nextera Transposase Adapters</h2><p>Read 1</p><p>5′ {S['tn5a']}</p><p>Read 2</p><p>5′ {S['tn5b']}</p>"
        f"<h2>PCR Primers</h2><p>Index 1 Read</p><p>5′ {S['nx_i7_front']}[i7]{S['nx_i7_rear']}</p>"
        f"<p>Index 2 Read</p><p>5′ {S['nx_i5_front']}[i5]{S['nx_i5_rear']}</p>"), encoding="utf-8")
    write_pdf(directory / "vendor.pdf", [
        ["Contents", "   Legacy Kit Z (Obsolete)        7", "Nextera DNA Indexes",
         "    i7 Index Name      Bases in Adapter       i7 Bases for Sample Sheet",
         f"    [H/N]701           {S['hn701']}               {rc(S['hn701'])}",
         f" [E/H/N/S]501       {S['ehns501']}       {S['ehns501']}       {rc(S['ehns501'])}"],
        ["Legacy Kit Z", "(Obsolete)", "PCR Primers", f"5′ {S['pcr1']}", f"5′ {S['pcr2']}", "", "Paired End (Obsolete)", "PE Adapters",
         f"5′ P-{S['pe1']}", f"5′ {S['pe2']}", "RNA 5′ Adapter (RA5)", f"5′ {S['ra5'].replace('T', 'U')}"]])
    write_pdf(directory / "manual.pdf", [["Front matter"], [
        "        Name                 Default", f"                     【{S['mgi1_head']}", "        【adapter1】            adapter 1 sequence",
        f"                     {S['mgi1_tail']}】", f"        【adapter2】   【{S['mgi2']}】", "        【Phred】      【2】"]])
    with zipfile.ZipFile(directory / "official.zip", "w") as archive:
        archive.writestr("barcodes.fasta", f">bc1001\n{S['b1']}\n")
    (directory / "mirror.fasta").write_text(f">bc1001\n{S['b1']}\n>bc1002\n{S['b2'][:8]}\n{S['b2'][8:]}\n", encoding="utf-8")
    (directory / "catalogue.fa").write_text(f">gnl|uv|NGB90001.1:1-30 Synthetic hairpin\n{S['hairpin']}\n"
                                            f">gnl|uv|NGB90002.1:1-20 Synthetic primer\n{S['primer'].lower()}\n", encoding="utf-8")
    write_xlsx(directory / "layout.xlsx", ["barcode ID", "Left Barcode", "Right Barcode", "Adapter Sequence", "Oligo Modification",
                                           "Barcoded Adapter Sequence "],
               [[name, rc(barcode), barcode, S["hairpin"], "/5Phos/", "/5Phos/" + barcode + S["hairpin"] + rc(barcode)]
                for name, barcode in (("_0001", S["b1"]), ("_0003", S["b3"]))])
    (directory / "reference.pdf").write_bytes(b"reference-only evidence file\n")


def definition_rows() -> tuple[list[dict], dict[str, str]]:
    rows, expected = [], {}
    ont, dorado, pacbio = "ont=ont.html", "dorado=dorado.cpp", "official=official.zip;mirror=mirror.fasta"

    def add(query_id, group, condition, rule, element, sources, sequence, segments, verification="verified", construct_class="adapter",
            transformation="none", index_side=""):
        expected[query_id] = sequence
        rows.append(dict(query_id=query_id, group=group, condition=condition, family=condition, construct_class=construct_class,
                         length=len(sequence), sequence_sha256=sha256(sequence), kit_or_set="synthetic fixture set",
                         orientation="as documented", segments=segments, transformation=transformation, rule=rule, element=element,
                         sources=sources, verification=verification, index_name=element if index_side else "", index_side=index_side))

    add("NB01_bare", "ONT", "ont_barcode", "ont_native_barcode", "NB01", f"{ont};{dorado}", S["nb1"], "barcode:1-24=NB01",
        construct_class="barcode")
    add("K40_NB01", "ONT", "k40", "ont_native_context_k40", "NB01", f"{ont};{dorado}", S["front14"][-8:] + S["nb1"] + S["rear8"],
        "front:1-8;barcode:9-32=NB01;rear:33-40", construct_class="composite")
    add("K46_NB02", "ONT", "k46", "ont_native_context_k46", "NB02", f"{dorado};{ont}", S["front14"] + S["nb2"] + S["rear8"],
        "front:1-14;barcode:15-38=NB02;rear:39-46", construct_class="composite")
    for position, name in ((1, "BC01"), (2, "RBK02")):
        add(f"RAPID_RB0{position}", "ONT", "rapid", "ont_rapid_context_90", f"RB0{position}", f"{dorado};{ont}",
            S["rapid_front"] + S["bc1" if position == 1 else "rbk2"] + S["rapid_rear"], f"front:1-10;barcode:11-34={name};rear:35-54",
            construct_class="composite")
    add("LIG_top", "ONT", "ligation", "porechop_adapter", "Kit_A_top", f"porechop=porechop.py;{ont}", S["p1"], "oligo:1-22")
    add("LIG_bottom", "ONT", "ligation", "porechop_adapter", "Kit_A_bottom", f"porechop=porechop.py;{ont}", S["p2"], "oligo:1-26",
        verification="verified_secondary_source")
    add("UDP0001_i7", "Illumina", "i7_index", "illumina_index", "UDP0001", "table=ud.html", S["i7"], "index:1-10=UDP0001",
        construct_class="barcode", index_side="i7")
    add("UDP0002V3_i5", "Illumina", "i5_index", "illumina_index", "UDP0002V3", "table=ud.html", S["i5b"], "index:1-10=UDP0002V3",
        construct_class="barcode", index_side="i5")
    add("UDP0001_i7_adapter", "Illumina", "i7_index_adapter", "illumina_index_adapter", "UDP0001", "table=ud.html;formula=ud.html",
        S["ud_i7_front"] + S["i7"] + S["ud_i7_rear"], "front:1-12;index:13-22=UDP0001;rear:23-31", construct_class="composite",
        index_side="i7")
    add("HN701_i7", "Illumina", "i7_index", "illumina_index", "HN701", "table=vendor.pdf", S["hn701"], "index:1-8=HN701",
        construct_class="barcode", index_side="i7")
    add("EHNS501_i5_adapter", "Illumina", "i5_index_adapter", "illumina_index_adapter", "EHNS501", "table=vendor.pdf;formula=nextera.html",
        S["nx_i5_front"] + S["ehns501"] + S["nx_i5_rear"], "front:1-14;index:15-22=EHNS501;rear:23-31", construct_class="composite",
        index_side="i5")
    add("I7_UNIVERSAL", "Illumina", "truseq", "illumina_formula_segment", "i7:front", "formula=ud.html", S["ud_i7_front"], "oligo:1-12")
    add("TN5_A", "Illumina", "nextera", "documented_oligo", "Nextera Transposase Adapters > Read 1", "document=nextera.html", S["tn5a"],
        "oligo:1-20")
    add("TRIM", "Illumina", "nextera", "documented_oligo",
        "Sequences for Illumina Prep and Nextera Kits > The following sequence is used for Read 1 and Read 2 adapter trimming.",
        "document=nextera.html", S["trim"], "oligo:1-19")
    add("PCR2", "Illumina", "truseq", "documented_oligo", "Legacy Kit Z (Obsolete) > PCR Primers #2", "document=vendor.pdf", S["pcr2"],
        "oligo:1-21")
    add("PE_ATAIL", "Illumina", "truseq", "documented_oligo", "Paired End (Obsolete) > PE Adapters #1", "document=vendor.pdf",
        "A" + S["pe1"], "a_tail:1-1;oligo:2-25", transformation="A_tail_base")
    add("RA5", "Illumina", "smallrna", "documented_oligo", "RNA 5′ Adapter (RA5)", "document=vendor.pdf", S["ra5"], "oligo:1-26",
        transformation="U_to_T")
    add("MGI_1", "BGI/DNBSEQ", "bgi", "mgi_filter_adapter", "2:adapter1", "manual=manual.pdf", S["mgi1_head"] + S["mgi1_tail"], "oligo:1-23")
    add("MGI_2", "BGI/DNBSEQ", "bgi", "mgi_filter_adapter", "2:adapter2", "manual=manual.pdf", S["mgi2"], "oligo:1-30")
    add("BC1001", "PacBio", "pacbio_barcode", "pacbio_barcode", "BC1001", pacbio, S["b1"], "barcode:1-16=BC1001", construct_class="barcode")
    add("BC1002", "PacBio", "pacbio_barcode", "pacbio_barcode", "BC1002", pacbio, S["b2"], "barcode:1-16=BC1002",
        verification="verified_secondary_source", construct_class="barcode")
    add("BC1001_hairpin", "PacBio", "pacbio_barcode_adapter", "pacbio_barcode_hairpin", "BC1001+NGB90001.1",
        f"{pacbio};hairpin=catalogue.fa;layout=layout.xlsx;reference=reference.pdf", S["b1"] + S["hairpin"],
        "barcode:1-16=BC1001;hairpin:17-46=NGB90001.1", verification="search_representation", construct_class="composite")
    add("PRIMER", "PacBio", "pacbio_smrtbell_adapter", "catalogue_record", "NGB90002.1", "catalogue=catalogue.fa", S["primer"], "oligo:1-20",
        construct_class="primer")
    return rows, expected


class PanelBuilderTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.sources, self.study = self.root / "remote", self.root / "study"
        write_sources(self.sources)
        self.study.mkdir()
        self.rows, self.expected = definition_rows()
        self.definition = self.write_definition(self.rows)
        self.pinned = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(self.sources.iterdir())}
        manifest = [dict(name=name, url=(self.sources / name).resolve().as_uri(), retrieved="2026-11-09T00:00:00+00:00", sha256=digest,
                         bytes=(self.sources / name).stat().st_size, purpose="synthetic fixture", redistribution="synthetic test data")
                    for name, digest in self.pinned.items()]
        (self.study / "sources.json").write_text(json.dumps(dict(kind="panel_sources", sources=manifest)), encoding="utf-8")

    def write_definition(self, rows):
        path = self.study / "panel_definition.tsv"
        path.unlink(missing_ok=True)
        write_table(path, rows, DEFINITION_SCHEMA)
        return path

    def invoke(self, arguments):
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            status = main(list(map(str, arguments)))
        return status, output.getvalue(), errors.getvalue()

    def test_every_rule_reproduces_the_definition(self):
        report = build_panel(self.definition, self.sources, self.root / "panel")
        fasta = self.root / "panel" / "queries.fasta"
        self.assertEqual(fasta.read_text(), "".join(f">{name}\n{sequence}\n" for name, sequence in self.expected.items()))
        self.assertEqual(read_fasta(fasta), self.expected)
        self.assertEqual(set(report["by_rule"]), set(RULES))
        self.assertEqual(report["by_verification"], {"verified": 21, "verified_secondary_source": 2, "search_representation": 1})
        self.assertEqual(report["canonical_parents"], len({canonical_id(sequence) for sequence in self.expected.values()}))
        self.assertEqual(set(report["text_extracts"]), {"vendor.pdf", "manual.pdf#page=2"})
        with (self.root / "panel" / "queries.tsv").open(encoding="utf-8") as handle:
            header = next(csv.reader(handle, delimiter="\t"))
        self.assertEqual(header, list(metadata_schema()["properties"]))
        self.assertNotIn("role", header)
        metadata = {row["query_id"]: row for row in read_table(self.root / "panel" / "queries.tsv", metadata_schema())}
        self.assertEqual(metadata["K40_NB01"]["family_id"], canonical_id(S["nb1"]))
        self.assertEqual(metadata["UDP0001_i7_adapter"]["family_id"], metadata["UDP0001_i7"]["family_id"])
        self.assertEqual((metadata["BC1001"]["source_id"], metadata["BC1002"]["source_id"]), ("official.zip", "mirror.fasta"))
        self.assertEqual([metadata[name]["architecture_status"] for name in ("K46_NB02", "BC1001_hairpin", "NB01_bare")],
                         ["verified", "verified", "not_applicable"])
        notes = {item["query_id"]: item["notes"] for item in report["derivations"]}
        self.assertEqual((notes["RAPID_RB01"]["ont_barcode_table"], notes["RAPID_RB02"]["ont_barcode_table"]), ("identical", "different"))
        self.assertEqual((notes["LIG_top"]["ont_document"], notes["LIG_bottom"]["ont_document"]), ("present", "absent"))
        self.assertEqual(notes["BC1001_hairpin"]["documented_full_adapter"], "yes")
        self.assertEqual(report["files"]["queries.fasta"], hashlib.sha256(fasta.read_bytes()).hexdigest())
        self.assertNotIn(S["nb1"], (self.root / "panel" / "panel_build.json").read_text())

    def test_output_feeds_controls_and_prepare_panel(self):
        build = build_panel(self.definition, self.sources, self.root / "panel")
        queries, metadata = self.root / "panel" / "queries.fasta", self.root / "panel" / "queries.tsv"
        controls = generate_controls(queries, self.root / "controls", 2)
        self.assertEqual(controls["canonical_parents"], build["canonical_parents"])
        panel = prepare_panel(queries, metadata, self.root / "controls", self.root / "prepared", 2)
        self.assertEqual((panel["physical_real_queries"], panel["canonical_parents"]), (len(self.expected), build["canonical_parents"]))

    def test_changed_or_missing_source_stops_before_writing(self):
        page = self.sources / "nextera.html"
        page.write_text(page.read_text(encoding="utf-8").replace("Read 2", "Read  2"), encoding="utf-8")
        with self.assertRaisesRegex(ContractError, "nextera.html: SHA-256 [0-9a-f]{64} differs from the pinned"):
            build_panel(self.definition, self.sources, self.root / "panel")
        (self.sources / "mirror.fasta").unlink()
        with self.assertRaisesRegex(ContractError, "mirror.fasta: missing"):
            build_panel(self.definition, self.sources, self.root / "panel")
        self.assertFalse((self.root / "panel").exists())

    def test_definition_disagreements_are_reported_per_query(self):
        rows = {row["query_id"]: dict(row) for row in self.rows}
        rows["TN5_A"]["sequence_sha256"] = sha256(S["tn5b"])
        rows["LIG_bottom"]["verification"] = "verified"
        rows["K46_NB02"]["segments"] = "front:1-14;barcode:15-38=NB01;rear:39-46"
        rows["MGI_2"]["length"] = 31
        self.write_definition(list(rows.values()))
        with self.assertRaises(ContractError) as caught:
            build_panel(self.definition, self.sources, self.root / "panel")
        message = str(caught.exception)
        self.assertIn(f"4 of {len(rows)} queries failed", message)
        for fragment in ("TN5_A: derived sequence_sha256", "LIG_bottom: derived verification verified_secondary_source != verified",
                         "K46_NB02: derived segments", "MGI_2: derived length 30 != 31"):
            self.assertIn(fragment, message)
        self.assertNotIn(S["tn5a"], message)
        self.assertFalse((self.root / "panel").exists())

    def test_definition_and_rule_contracts(self):
        cases = [("EHNS501_i5_adapter", dict(sources="table=vendor.pdf"), "needs exactly the source roles"),
                 ("TN5_A", dict(sources="document=unknown.html"), "missing from the manifest"),
                 ("NB01_bare", dict(transformation="U_to_T"), "transformations apply only to documented oligos"),
                 ("PRIMER", dict(index_name="NGB90002.1"), "index_name/index_side"),
                 ("TN5_A", dict(transformation="U_to_T"), "U_to_T declared for a DNA oligo"),
                 ("RA5", dict(transformation="none"), "RNA oligo requires the U_to_T transformation"),
                 ("PCR2", dict(element="Legacy Kit Z (Obsolete) > PCR Primers #3"), "No documented oligo follows"),
                 ("RAPID_RB02", dict(element="RB03"), "not a position of the Dorado SQK-RBK114-96 barcode list"),
                 ("BC1001", dict(element="BC1009"), "absent from the official barcode set and the mirror"),
                 ("LIG_top", dict(element="Kit_A_middle"), "Porechop has no adapter named Kit_A_middle"),
                 ("MGI_1", dict(element="2:adapter3"), "no sequence parameter adapter3")]
        for query_id, change, message in cases:
            with self.subTest(query_id=query_id, change=change):
                self.write_definition([dict(row, **change) if row["query_id"] == query_id else row for row in self.rows])
                with self.assertRaisesRegex(ContractError, message):
                    build_panel(self.definition, self.sources, self.root / "panel")
                self.assertFalse((self.root / "panel").exists())
        self.write_definition(self.rows)
        (self.root / "panel").mkdir()
        with self.assertRaisesRegex(ContractError, "Output directory already exists"):
            build_panel(self.definition, self.sources, self.root / "panel")

    def test_unreadable_pinned_source_fails_its_queries(self):
        (self.sources / "official.zip").write_bytes(b"not a zip archive\n")
        manifest = json.loads((self.study / "sources.json").read_text(encoding="utf-8"))
        for item in manifest["sources"]:
            if item["name"] == "official.zip":
                item.update(sha256=hashlib.sha256(b"not a zip archive\n").hexdigest(), bytes=18)
        (self.study / "sources.json").write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ContractError, r"3 of 24 queries failed.*BC1001: official.zip: unreadable as text \(BadZipFile"):
            build_panel(self.definition, self.sources, self.root / "panel")
        self.assertFalse((self.root / "panel").exists())

    def test_parsers_reject_ambiguous_or_unknown_layouts(self):
        page = f"x\n   【{S['mgi1_head']}\n【adapter1】  label\n   {S['mgi1_tail']}】\n【adapter2】 【{S['mgi2']}】\n【n】 【2】\n"
        self.assertEqual(parse_bracket_parameters(page), dict(adapter1=S["mgi1_head"] + S["mgi1_tail"], adapter2=S["mgi2"]))
        text = f"Kit Q PCR Primers 5′ {S['pcr1']} 5′ {S['pcr2']} Read 1 x"
        self.assertEqual(locate_oligo(text, "Kit Q > PCR Primers #2"), S["pcr2"])
        for locator, message in (("Kit R > PCR Primers", "anchor not found: Kit R"), ("Kit Q > Read 1", "No documented oligo"),
                                 ("Kit Q > PCR Primers #b", "Invalid locator ordinal")):
            with self.subTest(locator=locator), self.assertRaisesRegex(ContractError, message):
                locate_oligo(text, locator)
        with self.assertRaisesRegex(ContractError, "unrecognised column layout"):
            parse_illumina_indexes(f"UDP0009 {S['i7']} {S['i7']}")
        with self.assertRaisesRegex(ContractError, "lists UDP0009/i7 with different sequences"):
            parse_illumina_indexes(f"UDP0009 {S['i7']} {reverse_complement(S['i7'])} UDP0009 {S['i7b']} {reverse_complement(S['i7b'])}")
        dorado = (self.sources / "dorado.cpp").read_text(encoding="utf-8")
        self.assertEqual(dorado_kit(parse_dorado(dorado), "SQK-RBK114-96")["barcodes"], ["BC01", "RBK02"])
        with self.assertRaisesRegex(ContractError, "unexpected KitInfo layout"):
            dorado_kit(parse_dorado(dorado.replace('"NB96", true,', '"NB96",')), "SQK-NBD114-96")
        with self.assertRaisesRegex(ContractError, "Dorado kit not found"):
            dorado_kit(parse_dorado(dorado), "SQK-XYZ000")
        library = Library(self.sources, {name: {} for name in self.pinned})
        entry = dict(query_id="Q", element="Nextera Transposase Adapters > Read 2", transformation="A_tail_base")
        self.assertEqual(rule_documented_oligo(library, entry, dict(document="nextera.html")).segments, "a_tail:1-1;oligo:2-22")
        self.assertEqual(library.parsed("illumina_formulas", "nextera.html")["i5"], (S["nx_i5_front"], S["nx_i5_rear"]))
        with self.assertRaisesRegex(ContractError, "Several i7 adapter formulas"):
            parse_illumina_formulas(f"{S['ud_i7_front']}[i7]{S['ud_i7_rear']} {S['nx_i7_front']}[i7]{S['nx_i7_rear']}")

    def test_fetch_sources_accepts_pinned_files_and_never_replaces_conflicts(self):
        local = self.root / "local"
        local.mkdir()
        (local / "ont.html").write_bytes((self.sources / "ont.html").read_bytes())
        report = fetch_sources(self.study / "sources.json", local)
        self.assertEqual(report["status_counts"], dict(present=1, downloaded=len(self.pinned) - 1))
        self.assertEqual({path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in local.iterdir()}, self.pinned)
        (local / "dorado.cpp").write_text("locally edited\n", encoding="utf-8")
        with self.assertRaisesRegex(ContractError, "dorado.cpp: existing local file differs from the pinned SHA-256 and was left unchanged"):
            fetch_sources(self.study / "sources.json", local)
        self.assertEqual((local / "dorado.cpp").read_text(encoding="utf-8"), "locally edited\n")
        (self.sources / "porechop.py").write_text("changed upstream\n", encoding="utf-8")
        (self.sources / "reference.pdf").unlink()
        fresh = self.root / "fresh"
        with self.assertRaises(ContractError) as caught:
            fetch_sources(self.study / "sources.json", fresh)
        self.assertRegex(str(caught.exception), "porechop.py: downloaded content has SHA-256 [0-9a-f]{64}, not the pinned value, and was not stored")
        self.assertIn("reference.pdf: Download failed", str(caught.exception))
        self.assertEqual(sorted(path.name for path in fresh.iterdir()), sorted(set(self.pinned) - {"porechop.py", "reference.pdf"}))

    def test_cli_fetches_and_builds(self):
        local = self.root / "local"
        status, output, error = self.invoke(["fetch-sources", self.study / "sources.json", "--directory", local])
        self.assertEqual(status, 0, error)
        self.assertEqual(json.loads(output)["status_counts"], dict(downloaded=len(self.pinned)))
        status, output, error = self.invoke(["build-panel", self.definition, "--sources", local, "--output", self.root / "panel"])
        self.assertEqual(status, 0, error)
        result = json.loads(output)
        self.assertEqual((result["queries"], result["by_group"]), (len(self.expected), {"ONT": 7, "Illumina": 11, "BGI/DNBSEQ": 2, "PacBio": 4}))
        status, output, error = self.invoke(["build-panel", self.definition, "--sources", self.root / "empty", "--output", self.root / "other"])
        self.assertEqual((status, output), (1, ""))
        self.assertIn("Error: Pinned sources missing or changed", error)
        self.assertNotIn("Traceback", error)


if __name__ == "__main__":
    unittest.main()
