"""Derive a study query panel from pinned vendor and catalogue documents.

No vendor document and no derived sequence is shipped. `fetch_sources` downloads
every pinned source and checks its SHA-256; `build_panel` parses the local copies
with deterministic per-source parsers, assembles each query with its family rule
and accepts it only if length, segments, verification status and SHA-256 equal
the panel definition.
"""
import ast
import hashlib
import html
import http.client
import re
import shutil
import subprocess
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

from .contracts import ContractError, read_json, read_table, schema, write_table
from .design import digest
from .example import save_json
from .sequences import canonical_id, reverse_complement

NATIVE_KIT, RAPID_KIT = "SQK-NBD114-96", "SQK-RBK114-96"
USER_AGENTS = ("Mozilla/5.0", "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0", "curl/8.5.0")
VERIFICATIONS = ("verified", "verified_secondary_source", "search_representation")
MANIFEST_SCHEMA = dict(type="object", additionalProperties=False, required=["kind", "sources"], properties={
    "kind": {"const": "panel_sources"}, "description": {"type": "string"},
    "sources": {"type": "array", "minItems": 1, "items": dict(
        type="object", additionalProperties=False,
        required=["name", "url", "retrieved", "sha256", "bytes", "purpose", "redistribution"], properties={
            "name": {"type": "string", "pattern": r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$"},
            "url": {"type": "string", "pattern": r"^[a-z]+://\S+$"}, "retrieved": {"type": "string", "minLength": 1},
            "sha256": {"type": "string", "pattern": "^[a-f0-9]{64}$"}, "bytes": {"type": "integer", "minimum": 1},
            "purpose": {"type": "string", "minLength": 1}, "redistribution": {"type": "string", "minLength": 1}})}})


def load_manifest(path: Path) -> dict[str, dict]:
    entries = read_json(path, MANIFEST_SCHEMA)["sources"]
    manifest = {entry["name"]: entry for entry in entries}
    if len(manifest) != len(entries):
        raise ContractError(f"{path.name}: duplicate source names")
    return manifest


def download(url: str, timeout: float = 180) -> bytes:
    for agent in USER_AGENTS:
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": agent}), timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:  # some vendor firewalls reject individual user agents
            last = error
        except (OSError, ValueError, http.client.HTTPException) as error:
            raise ContractError(f"Download failed for {url}: {type(error).__name__}: {error}") from error
    raise ContractError(f"Download failed for {url}: {type(last).__name__}: {last}")


def collect_sources(manifest: Path, directory: Path, destination: Path, *, fetch: bool = True, allow_drift: bool = False) -> dict:
    """Snapshot every pinned source of `manifest` into `destination/files` for one study run.

    A local copy in `directory` is used when present; a missing file is downloaded (only with `fetch`) and kept in
    `directory` only when it matches its pin. A local or downloaded file whose SHA-256 differs from its pin is refused
    unless `allow_drift`: it is then snapshotted with its observed hash and listed as drift (a changed download is
    stored in the snapshot only, never in `directory`). `destination/sources.json` is the manifest of the snapshot
    (observed hashes): a panel built from it accepts a drifted file only if every derived sequence still matches its
    definition SHA-256. Writes `collection.json` with the per-source record and returns it.
    """
    entries = load_manifest(manifest)
    destination.mkdir(parents=True, exist_ok=False)
    (destination / "files").mkdir()
    records, problems, snapshot = [], [], []
    for name, entry in entries.items():
        local = directory / name
        record = dict(name=name, url=entry["url"], pinned_sha256=entry["sha256"], pinned_bytes=entry["bytes"])
        retrieved = entry["retrieved"]
        if local.is_file():
            body, origin = local.read_bytes(), "local"
        elif local.exists() or local.is_symlink():
            problems.append(f"{name}: local path is not a regular file")
            continue
        elif not fetch:
            problems.append(f"{name}: missing from {directory} and fetching is disabled")
            continue
        else:
            try:
                body, origin, retrieved = download(entry["url"]), "downloaded", datetime.now(timezone.utc).isoformat()
            except ContractError as error:
                problems.append(f"{name}: {error}")
                continue
        observed = hashlib.sha256(body).hexdigest()
        if observed == entry["sha256"]:
            status = "pinned"
            if origin == "downloaded":
                directory.mkdir(parents=True, exist_ok=True)
                partial = directory / (name + ".part")
                partial.write_bytes(body)
                partial.replace(local)
        elif allow_drift:
            status = "drift"
            if origin == "local":
                retrieved = "unknown_local_copy"
        else:
            problems.append(f"{name}: {origin} file has SHA-256 {observed}, not the pinned {entry['sha256']}"
                            + (" (local file left unchanged)" if origin == "local" else " (download not stored)"))
            continue
        (destination / "files" / name).write_bytes(body)
        records.append(dict(record, origin=origin, status=status, sha256=observed, bytes=len(body), retrieved=retrieved))
        snapshot.append(dict(entry, sha256=observed, bytes=len(body), retrieved=retrieved))
    if problems:
        raise ContractError("Pinned sources unavailable or changed: " + "; ".join(problems))
    save_json(destination / "sources.json", dict(kind="panel_sources", sources=snapshot,
                                                 description="Snapshot of the pinned sources of one study run (observed hashes)"))
    report = dict(manifest=str(manifest), manifest_sha256=digest(manifest), directory=str(directory), fetch=fetch,
                  allow_drift=allow_drift, sources=records, drift=[record for record in records if record["status"] == "drift"],
                  status_counts=dict(Counter(record["status"] for record in records)),
                  origin_counts=dict(Counter(record["origin"] for record in records)),
                  drift_rule="A drifted file is accepted only if every derived panel sequence still matches its definition SHA-256")
    save_json(destination / "collection.json", report)
    return report


def fetch_sources(manifest: Path, directory: Path) -> dict:
    """Download every pinned source into `directory` and verify its SHA-256.

    Existing files are kept only when their hash matches; a mismatching local file
    is never replaced, and a changed download is never stored.
    """
    entries = load_manifest(manifest)
    directory.mkdir(parents=True, exist_ok=True)
    records, problems = [], []
    for name, entry in entries.items():
        target = directory / name
        record = dict(name=name, url=entry["url"], sha256=entry["sha256"])
        if target.exists() or target.is_symlink():
            observed = digest(target) if target.is_file() else ""
            if observed == entry["sha256"]:
                records.append(dict(record, status="present"))
            else:
                problems.append(f"{name}: existing local file differs from the pinned SHA-256 and was left unchanged")
                records.append(dict(record, status="local_conflict", observed_sha256=observed))
            continue
        try:
            body = download(entry["url"])
        except ContractError as error:
            problems.append(f"{name}: {error}")
            records.append(dict(record, status="failed", error=str(error)))
            continue
        observed = hashlib.sha256(body).hexdigest()
        if observed != entry["sha256"]:
            problems.append(f"{name}: downloaded content has SHA-256 {observed}, not the pinned value, and was not stored")
            records.append(dict(record, status="changed", observed_sha256=observed, bytes=len(body)))
            continue
        partial = directory / (name + ".part")
        partial.write_bytes(body)
        partial.replace(target)
        records.append(dict(record, status="downloaded", bytes=len(body), retrieved_utc=datetime.now(timezone.utc).isoformat()))
    if problems:
        raise ContractError("Pinned sources unavailable or changed: " + "; ".join(problems))
    return dict(manifest=str(manifest), manifest_sha256=digest(manifest), directory=str(directory), sources=records,
                status_counts=dict(Counter(record["status"] for record in records)))


def verify_sources(directory: Path, entries: list[dict]) -> None:
    problems = []
    for entry in entries:
        path = directory / entry["name"]
        observed = digest(path) if path.is_file() else None
        if observed != entry["sha256"]:
            problems.append(f"{entry['name']}: " + ("missing" if observed is None else
                                                     f"SHA-256 {observed} differs from the pinned {entry['sha256']}"))
    if problems:
        raise ContractError("Pinned sources missing or changed: " + "; ".join(problems))


def normalize(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def html_text(raw: str) -> str:
    return normalize(html.unescape(re.sub(r"<[^>]+>", " ", raw)))


def pdf_text(path: Path, page: int | None = None, executable: str | None = None) -> str:
    """Layout text exactly as `pdftotext -layout` (poppler) writes it, optionally for one page."""
    executable = executable or shutil.which("pdftotext")
    if executable is None:
        raise ContractError("pdftotext (poppler-utils) is required to extract PDF sources")
    pages = [] if page is None else ["-f", str(page), "-l", str(page)]
    result = subprocess.run([executable, "-layout", "-enc", "UTF-8", *pages, str(path), "-"], capture_output=True, check=False)
    if result.returncode != 0:
        raise ContractError(f"pdftotext failed for {path.name}: {result.stderr.decode('utf-8', 'replace').strip()}")
    return result.stdout.decode("utf-8")


def pdftotext_version(executable: str | None = None) -> str:
    executable = executable or shutil.which("pdftotext")
    if executable is None:
        return ""
    result = subprocess.run([executable, "-v"], capture_output=True, check=False)
    return (result.stderr or result.stdout).decode("utf-8", "replace").splitlines()[0].strip()


def unique(pairs, what: str) -> dict:
    values = {}
    for key, value in pairs:
        if values.setdefault(key, value) != value:
            raise ContractError(f"{what} lists {'/'.join(key) if isinstance(key, tuple) else key} with different sequences")
    return values


def flank_statement(text: str, heading: str, label: str) -> tuple[str, str]:
    start = text.find(heading)
    match = re.compile(re.escape(label) + r" 5' - ([ACGT]+) - barcode - ([ACGT]+) - 3'").search(text, max(start, 0))
    if start < 0 or match is None:
        raise ContractError(f"ONT document: no '{label}' flank statement after '{heading}'")
    return match.group(1), match.group(2)


@dataclass
class OntDocument:
    text: str
    native: dict
    native_forward: tuple
    native_reverse: tuple
    rapid: tuple
    rapid_table: dict


def parse_ont_document(text: str) -> OntDocument:
    """ONT Chemistry Technical Document: native table (forward/reverse columns), barcode tables and flank statements."""
    text = normalize(text)
    rows = re.findall(r"(?<!\S)(NB\d{2}) ([ACGT]{12,}) ([ACGT]{12,})(?!\S)", text)
    if any(reverse != reverse_complement(forward) for _, forward, reverse in rows):
        raise ContractError("ONT native table: a reverse column is not the reverse complement of its forward column")
    labelled = re.findall(r"(?<!\S)((?:[A-Z0-9]+ / )*RB\d{2}(?: / [A-Z0-9]+)*) ([ACGT]{12,})(?!\S)", text)
    return OntDocument(text=text, native=unique(((name, forward) for name, forward, _ in rows), "ONT native table"),
                       native_forward=flank_statement(text, "Native barcode sequences", "Forward sequence:"),
                       native_reverse=flank_statement(text, "Native barcode sequences", "Reverse sequence:"),
                       rapid=flank_statement(text, "Rapid barcoding kits", "Barcode flanking sequence:"),
                       rapid_table=unique(((name, sequence) for label, sequence in labelled for name in label.split(" / ")),
                                          "ONT barcode table"))


@dataclass
class DoradoKits:
    text: str
    barcodes: dict
    constants: dict
    lists: dict


def parse_dorado(text: str) -> DoradoKits:
    """Dorado `barcode_kits.cpp` without line comments: barcode entries, flank constants and barcode lists."""
    code = re.sub(r"//[^\n]*", "", text)
    return DoradoKits(text=code,
                      barcodes=unique(re.findall(r'\{"([A-Za-z0-9_]+)",\s*"([ACGT]+)"\}', code), "Dorado barcode map"),
                      constants=unique(re.findall(r'const std::string (\w+) = "([ACGT]*)";', code), "Dorado constants"),
                      lists=unique(((name, re.findall(r'"([A-Za-z0-9_]+)"', body)) for name, body in
                                    re.findall(r"const std::vector<std::string> (\w+) = \{(.*?)\};", code, re.S)), "Dorado barcode lists"))


def braced(text: str, start: int) -> str:
    if start >= len(text) or text[start] != "{":
        raise ContractError("Dorado kit definition does not start with a brace")
    depth = 0
    for position in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[position], 0)
        if depth == 0:
            return text[start + 1:position]
    raise ContractError("Unbalanced braces in the Dorado kit definition")


def dorado_kit(dorado: DoradoKits, kit: str) -> dict:
    """Resolve a kit's first front/rear flank constants and barcode list (KitInfo fields 7, 10 and 15)."""
    text = dorado.text
    match = re.search(r'\{"' + re.escape(kit) + r'",\s*', text)
    if match is None:
        raise ContractError(f"Dorado kit not found: {kit}")
    start = match.end()
    if not text.startswith("{", start):
        variable = re.match(r"\w+", text[start:])
        definition = variable and re.search(r"const KitInfo " + re.escape(variable.group(0)) + r" = ", text)
        if not definition:
            raise ContractError(f"Dorado kit {kit}: KitInfo definition not found")
        start = definition.end()
    fields, depth, current = [], 0, ""
    for character in braced(text, start):
        depth += {"{": 1, "}": -1}.get(character, 0)
        if character == "," and depth == 0:
            fields.append(current.strip())
            current = ""
        else:
            current += character
    fields = [value for value in fields + [current.strip()] if value]
    if (len(fields) < 15 or not re.fullmatch(r'"[^"]*"', fields[0]) or any(value not in ("true", "false") for value in fields[1:6])
            or not fields[6].endswith("FRONT") or not fields[9].endswith("REAR") or fields[6] not in dorado.constants
            or fields[9] not in dorado.constants or fields[14] not in dorado.lists):
        raise ContractError(f"Dorado kit {kit}: unexpected KitInfo layout")
    return dict(front=dorado.constants[fields[6]], rear=dorado.constants[fields[9]], barcodes=dorado.lists[fields[14]],
                fields=f"{fields[6]},{fields[9]},{fields[14]}")


def parse_porechop(text: str) -> dict[str, str]:
    """Named constant sequences of Porechop `adapters.py`, read without executing the file."""
    pairs = []
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "Adapter":
            for keyword in node.keywords:
                value = keyword.value
                if (keyword.arg in {"start_sequence", "end_sequence", "both_ends_sequence"} and isinstance(value, ast.Tuple)
                        and len(value.elts) == 2 and all(isinstance(item, ast.Constant) and isinstance(item.value, str) for item in value.elts)):
                    pairs.append((value.elts[0].value, value.elts[1].value))
    return unique(pairs, "Porechop adapters")


def parse_illumina_indexes(text: str) -> dict[tuple[str, str], str]:
    """Index tables with bases in adapter; sample-sheet columns are used only to confirm the column layout.

    UD rows: i7 adapter, i7 sheet (= rc), i5 adapter, i5 forward sheet (= i5), i5 reverse sheet (= rc).
    Nextera rows: [H/N]7xx adapter, sheet (= rc); [E/H/N/S]5xx adapter, forward sheet, reverse sheet.
    """
    pairs = []
    for match in re.finditer(r"(?<!\S)(UDP\d{4}(?:V\d+)?|\[(?:[A-Z]/)*[A-Z]\]\d{3})((?: [ACGT]{6,})+)(?!\S)", normalize(text)):
        name, columns = re.sub(r"[\[\]/]", "", match.group(1)), match.group(2).split()
        if len(columns) == 5 and columns[1] == reverse_complement(columns[0]) and columns[3] == columns[2] == reverse_complement(columns[4]):
            pairs += [((name, "i7"), columns[0]), ((name, "i5"), columns[2])]
        elif len(columns) == 2 and columns[1] == reverse_complement(columns[0]):
            pairs.append(((name, "i7"), columns[0]))
        elif len(columns) == 3 and columns[1] == columns[0] == reverse_complement(columns[2]):
            pairs.append(((name, "i5"), columns[0]))
        else:
            raise ContractError(f"Illumina index row {match.group(1)} has an unrecognised column layout")
    return unique(pairs, "Illumina index tables")


def parse_illumina_formulas(text: str) -> dict[str, tuple[str, str]]:
    formulas = {}
    for side in ("i7", "i5"):
        found = set(re.findall(rf"(?<![A-Za-z])([ACGT]+)\[{side}\]([ACGT]+)(?![A-Za-z])", normalize(text)))
        if len(found) > 1:
            raise ContractError(f"Several {side} adapter formulas in one formula document")
        formulas.update({side: found.pop()} if found else {})
    return formulas


OLIGO = re.compile(r" ?(?:5['′’] ?)?(?:P-)?([ACGTU]{6,})(?:-3['′’])?(?= |$)")


def locate_oligo(text: str, locator: str) -> str:
    """Oligo following the last of the ' > '-separated anchors; '#n' selects the n-th oligo of a listed series.

    Anchors are literal, whitespace-normalised document text found in order; a 5′ marker and a 5′ phosphate
    notation (P-) precede the bases and are not bases.
    """
    text = normalize(text)
    path, _, ordinal = locator.partition(" #")
    if ordinal and not ordinal.isdigit():
        raise ContractError(f"Invalid locator ordinal: {locator}")
    position = 0
    for anchor in path.split(" > "):
        found = text.find(anchor, position)
        if found < 0:
            raise ContractError(f"Locator anchor not found: {anchor}")
        position = found + len(anchor)
    for _ in range(int(ordinal or 1)):
        match = OLIGO.match(text, position)
        if match is None:
            raise ContractError(f"No documented oligo follows the locator: {locator}")
        position = match.end()
    return match.group(1)


def parse_bracket_parameters(text: str) -> dict[str, str]:
    """Sequence-valued 【name】 【value】 parameters of a layout page; a value may wrap around its name line."""
    tokens = [token.strip() for token in re.findall(r"【[^【】\n]*】|【[^【】\n]*|[^【】\n]*】", text)]
    values = {}
    for index, token in enumerate(tokens):
        name = re.fullmatch(r"【(\w+)】", token)
        before = tokens[index - 1] if index else ""
        after = tokens[index + 1] if index + 1 < len(tokens) else ""
        head, tail, whole = re.fullmatch(r"【([ACGT]+)", before), re.fullmatch(r"([ACGT]+)】", after), re.fullmatch(r"【([ACGT]+)】", after)
        if name and head and tail:
            values[name.group(1)] = head.group(1) + tail.group(1)
        elif name and whole:
            values[name.group(1)] = whole.group(1)
    return values


def parse_fasta(text: str) -> dict[str, str]:
    """FASTA records keyed by the lower-case last '|' field of the identifier (e.g. bc1001, ngb00972.1:1-45)."""
    pairs, name, chunks = [], None, []
    for line in text.splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None:
                pairs.append((name, "".join(chunks).upper()))
            name, chunks = (line[1:].split() or [""])[0].split("|")[-1].lower(), []
        elif line.strip():
            if name is None:
                raise ContractError("FASTA content before the first header")
            chunks.append(line.strip())
    return unique(pairs, "FASTA source")


def fasta_record(records: dict[str, str], element: str, required: bool = True) -> str | None:
    key = element.lower()
    found = [sequence for name, sequence in records.items() if name == key or name.startswith(key + ":")]
    if len(found) > 1 or (required and not found):
        raise ContractError(f"FASTA source has {len(found)} records for {element}")
    return found[0] if found else None


def xlsx_rows(path: Path) -> list[dict[str, str]]:
    """Cell text of the first worksheet of an .xlsx workbook, keyed by column letter."""
    main = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    relation = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    with zipfile.ZipFile(path) as book:
        names = set(book.namelist())
        strings = (["".join(text.text or "" for text in item.iter(main + "t"))
                    for item in ElementTree.fromstring(book.read("xl/sharedStrings.xml"))] if "xl/sharedStrings.xml" in names else [])
        sheet = ElementTree.fromstring(book.read("xl/workbook.xml")).find(f"{main}sheets/{main}sheet")
        targets = {item.get("Id"): item.get("Target") for item in ElementTree.fromstring(book.read("xl/_rels/workbook.xml.rels"))}
        target = targets[sheet.get(relation)].lstrip("/")
        worksheet = ElementTree.fromstring(book.read(target if target.startswith("xl/") else "xl/" + target))
    rows = []
    for row in worksheet.iter(main + "row"):
        cells = {}
        for cell in row.iter(main + "c"):
            value = cell.find(main + "v")
            text = ("".join(item.text or "" for item in cell.iter(main + "t")) if cell.get("t") == "inlineStr" else
                    strings[int(value.text)] if cell.get("t") == "s" and value is not None else value.text if value is not None else "")
            cells[re.match(r"[A-Z]+", cell.get("r")).group(0)] = text.strip()
        rows.append(cells)
    return rows


def barcoded_adapters(rows: list[dict[str, str]]) -> list[tuple[str, str, str, str]]:
    """(right barcode, adapter, left barcode, full oligo without /modification/ notation) of a barcoded-adapter table."""
    labels = ("Right Barcode", "Adapter Sequence", "Left Barcode", "Barcoded Adapter Sequence")
    for index, row in enumerate(rows):
        columns = {value: column for column, value in row.items()}
        if all(label in columns for label in labels):
            return [tuple(re.sub(r"/[^/]*/", "", item.get(columns[label], "")).upper() for label in labels)
                    for item in rows[index + 1:] if all(item.get(columns[label]) for label in labels)]
    raise ContractError("Barcoded-adapter table header not found")


PARSERS = dict(ont=parse_ont_document, dorado=parse_dorado, porechop=parse_porechop, illumina_indexes=parse_illumina_indexes,
               illumina_formulas=parse_illumina_formulas, bracket_parameters=parse_bracket_parameters, fasta=parse_fasta)


class Library:
    """Lazily extracted and parsed text of the hash-verified local sources of one build."""

    def __init__(self, directory: Path, manifest: dict[str, dict], pdftotext: str | None = None):
        self.directory, self.manifest, self.cache, self.extracts = directory, manifest, {}, {}
        self.pdftotext = pdftotext

    def path(self, name: str) -> Path:
        if name not in self.manifest:
            raise ContractError(f"Source is not declared in the manifest: {name}")
        return self.directory / name

    def cached(self, key: tuple, factory):
        """Result of `factory`, computed once per source; a failure is also kept so dependent queries fail fast."""
        if key not in self.cache:
            try:
                self.cache[key] = factory()
            except ContractError as error:
                self.cache[key] = error
            except (ValueError, OSError, SyntaxError, KeyError, IndexError, zipfile.BadZipFile) as error:
                self.cache[key] = ContractError(f"{key[1]}: unreadable as {key[0]} ({type(error).__name__}: {error})")
        if isinstance(self.cache[key], ContractError):
            raise self.cache[key]
        return self.cache[key]

    def text(self, name: str, page: int | None = None) -> str:
        def extract():
            path = self.path(name)
            suffix = path.suffix.lower()
            if suffix == ".pdf":
                value = pdf_text(path, page, self.pdftotext)
                self.extracts[name + ("" if page is None else f"#page={page}")] = hashlib.sha256(value.encode("utf-8")).hexdigest()
                return value
            if suffix == ".zip":
                with zipfile.ZipFile(path) as archive:
                    members = [item for item in archive.namelist() if not item.endswith("/") and not item.startswith("__MACOSX")]
                    if len(members) != 1:
                        raise ContractError(f"{name}: expected exactly one archive member, found {len(members)}")
                    return archive.read(members[0]).decode("utf-8")
            value = path.read_text(encoding="utf-8")
            return html_text(value) if suffix in {".html", ".htm"} else value
        return self.cached(("text", name, page), extract)

    def parsed(self, kind: str, name: str, page: int | None = None):
        if kind == "barcoded_adapters":
            return self.cached((kind, name, page), lambda: barcoded_adapters(xlsx_rows(self.path(name))))
        return self.cached((kind, name, page), lambda: PARSERS[kind](self.text(name, page)))


@dataclass
class Derivation:
    sequence: str
    parts: list  # (segment kind, length, element name or "")
    verification: str
    core: str  # barcode/index of the query, or the whole oligo; keys the downstream family
    origin: str  # source role the sequence (or its barcode) was taken from
    notes: dict = field(default_factory=dict)

    @property
    def segments(self) -> str:
        items, start = [], 1
        for kind, length, name in self.parts:
            items.append(f"{kind}:{start}-{start + length - 1}" + (f"={name}" if name else ""))
            start += length
        return ";".join(items)


def native_barcode(library: Library, element: str, sources: dict) -> tuple[str, dict, OntDocument]:
    ont, dorado = library.parsed("ont", sources["ont"]), library.parsed("dorado", sources["dorado"])
    kit = dorado_kit(dorado, NATIVE_KIT)
    if element not in ont.native:
        raise ContractError(f"ONT native table has no {element}")
    if element not in kit["barcodes"] or dorado.barcodes.get(element) != ont.native[element]:
        raise ContractError(f"{element}: ONT native table and Dorado {NATIVE_KIT} disagree")
    return ont.native[element], kit, ont


def joined(longer: str, shorter: str, front: bool) -> str:
    longer, shorter = sorted((longer, shorter), key=len, reverse=True)
    if not (longer.endswith(shorter) if front else longer.startswith(shorter)):
        raise ContractError("ONT forward and reverse native contexts are inconsistent")
    return longer


def rule_ont_native_barcode(library, entry, sources):
    barcode, _, _ = native_barcode(library, entry["element"], sources)
    return Derivation(barcode, [("barcode", len(barcode), entry["element"])], "verified", barcode, "ont")


def rule_ont_native_context_k40(library, entry, sources):
    """ONT forward native context 5'-front-barcode-rear-3'; it must lie within the Dorado kit flanks."""
    barcode, kit, ont = native_barcode(library, entry["element"], sources)
    front, rear = ont.native_forward
    if not kit["front"].endswith(front) or kit["rear"] != rear:
        raise ContractError(f"ONT forward native context disagrees with the Dorado {NATIVE_KIT} flanks")
    return Derivation(front + barcode + rear, [("front", len(front), ""), ("barcode", len(barcode), entry["element"]),
                                              ("rear", len(rear), "")], "verified", barcode, "ont")


def rule_ont_native_context_k46(library, entry, sources):
    """Duplex union of the ONT forward and (reverse-complemented) reverse native contexts; must equal the Dorado kit flanks."""
    barcode, kit, ont = native_barcode(library, entry["element"], sources)
    (forward_front, forward_rear), (reverse_front, reverse_rear) = ont.native_forward, ont.native_reverse
    front = joined(forward_front, reverse_complement(reverse_rear), True)
    rear = joined(forward_rear, reverse_complement(reverse_front), False)
    if (front, rear) != (kit["front"], kit["rear"]):
        raise ContractError(f"Duplex native context differs from the Dorado {NATIVE_KIT} flanks")
    return Derivation(front + barcode + rear, [("front", len(front), ""), ("barcode", len(barcode), entry["element"]),
                                              ("rear", len(rear), "")], "verified", barcode, "ont")


def rule_ont_rapid_context_90(library, entry, sources):
    """Rapid context: ONT flank statement around the barcode bound to position RBnn of the Dorado kit list."""
    dorado, ont = library.parsed("dorado", sources["dorado"]), library.parsed("ont", sources["ont"])
    kit = dorado_kit(dorado, RAPID_KIT)
    position = re.fullmatch(r"RB(\d+)", entry["element"])
    if position is None or not 1 <= int(position.group(1)) <= len(kit["barcodes"]):
        raise ContractError(f"{entry['element']}: not a position of the Dorado {RAPID_KIT} barcode list")
    name = kit["barcodes"][int(position.group(1)) - 1]
    if name not in dorado.barcodes:
        raise ContractError(f"Dorado has no sequence for {name}")
    if ont.rapid != (kit["front"], kit["rear"]):
        raise ContractError(f"ONT Rapid flank statement disagrees with the Dorado {RAPID_KIT} flanks")
    barcode, (front, rear) = dorado.barcodes[name], ont.rapid
    table = ont.rapid_table.get(entry["element"])
    return Derivation(front + barcode + rear, [("front", len(front), ""), ("barcode", len(barcode), name), ("rear", len(rear), "")],
                      "verified", barcode, "dorado", dict(dorado_kit=f"{RAPID_KIT}:{kit['fields']}", ont_barcode_table=(
                          "absent" if table is None else "identical" if table == barcode else "different")))


def rule_porechop_adapter(library, entry, sources):
    """Named Porechop definition; verified by the vendor only if it occurs in the current ONT document."""
    sequence = library.parsed("porechop", sources["porechop"]).get(entry["element"])
    if sequence is None:
        raise ContractError(f"Porechop has no adapter named {entry['element']}")
    present = sequence in library.parsed("ont", sources["ont"]).text
    return Derivation(sequence, [("oligo", len(sequence), "")], "verified" if present else "verified_secondary_source", sequence,
                      "porechop", dict(ont_document="present" if present else "absent"))


def illumina_index(library, entry, sources) -> str:
    key = entry["element"], entry["index_side"]
    if entry["index_name"] != entry["element"] or entry["index_side"] not in ("i7", "i5"):
        raise ContractError("Illumina index rows need index_name equal to element and index_side i7 or i5")
    table = library.parsed("illumina_indexes", sources["table"])
    if key not in table:
        raise ContractError(f"Illumina index table has no {key[1]} index {key[0]}")
    return table[key]


def illumina_formula(library, name: str, side: str) -> tuple[str, str]:
    formulas = library.parsed("illumina_formulas", name)
    if side not in formulas:
        raise ContractError(f"{name}: no {side} adapter formula")
    return formulas[side]


def rule_illumina_index(library, entry, sources):
    index = illumina_index(library, entry, sources)
    return Derivation(index, [("index", len(index), entry["element"])], "verified", index, "table")


def rule_illumina_index_adapter(library, entry, sources):
    """Index bases in adapter embedded in the documented 5'-front[i7|i5]rear-3' adapter formula."""
    index = illumina_index(library, entry, sources)
    front, rear = illumina_formula(library, sources["formula"], entry["index_side"])
    return Derivation(front + index + rear, [("front", len(front), ""), ("index", len(index), entry["element"]), ("rear", len(rear), "")],
                      "verified", index, "table")


def rule_illumina_formula_segment(library, entry, sources):
    """Constant 5' (front) or 3' (rear) part of a documented index adapter formula, element 'i7:front'."""
    side, _, part = entry["element"].partition(":")
    if side not in ("i7", "i5") or part not in ("front", "rear"):
        raise ContractError(f"Invalid formula segment: {entry['element']}")
    sequence = dict(zip(("front", "rear"), illumina_formula(library, sources["formula"], side)))[part]
    return Derivation(sequence, [("oligo", len(sequence), "")], "verified", sequence, "formula")


def rule_documented_oligo(library, entry, sources):
    """Labelled oligo of a vendor document (see `locate_oligo`) with its declared transformation."""
    raw = locate_oligo(library.text(sources["document"]), entry["element"])
    transformation = entry["transformation"]
    if transformation == "U_to_T" and "U" not in raw:
        raise ContractError(f"{entry['query_id']}: U_to_T declared for a DNA oligo")
    sequence = raw.replace("U", "T") if transformation == "U_to_T" else "A" + raw if transformation == "A_tail_base" else raw
    if set(sequence) - set("ACGT"):
        raise ContractError(f"{entry['query_id']}: RNA oligo requires the U_to_T transformation")
    parts = [("a_tail", 1, ""), ("oligo", len(raw), "")] if transformation == "A_tail_base" else [("oligo", len(sequence), "")]
    return Derivation(sequence, parts, "verified", sequence, "document")


def rule_mgi_filter_adapter(library, entry, sources):
    """Default value of a 【adapterN】 parameter on one manual page, element 'page:parameter'."""
    page, _, name = entry["element"].partition(":")
    if not page.isdigit():
        raise ContractError(f"Invalid manual parameter element: {entry['element']}")
    sequence = library.parsed("bracket_parameters", sources["manual"], int(page)).get(name)
    if sequence is None:
        raise ContractError(f"Page {page} of {sources['manual']} has no sequence parameter {name}")
    return Derivation(sequence, [("oligo", len(sequence), "")], "verified", sequence, "manual")


def pacbio_barcode(library, name: str, sources: dict) -> tuple[str, str, str, dict]:
    official = fasta_record(library.parsed("fasta", sources["official"]), name, False)
    mirror = fasta_record(library.parsed("fasta", sources["mirror"]), name, False)
    if official and mirror and official != mirror:
        raise ContractError(f"{name}: official barcode set and mirror disagree")
    if official:
        return official, "verified", "official", dict(official_set="present", mirror="identical" if mirror else "absent")
    if mirror:
        return mirror, "verified_secondary_source", "mirror", dict(official_set="absent", mirror="present")
    raise ContractError(f"{name}: absent from the official barcode set and the mirror")


def rule_pacbio_barcode(library, entry, sources):
    barcode, status, origin, notes = pacbio_barcode(library, entry["element"], sources)
    return Derivation(barcode, [("barcode", len(barcode), entry["element"])], status, barcode, origin, notes)


def rule_pacbio_barcode_hairpin(library, entry, sources):
    """Search representation barcode + hairpin (element 'BCxxxx+accession'): the first part of the documented
    barcoded adapter 5'-barcode-hairpin-rc(barcode)-3', whose layout is checked on every documented example."""
    name, _, accession = entry["element"].partition("+")
    barcode, _, origin, notes = pacbio_barcode(library, name, sources)
    hairpin = fasta_record(library.parsed("fasta", sources["hairpin"]), accession)
    documented = library.parsed("barcoded_adapters", sources["layout"])
    if not documented or any(full != right + adapter + left or left != reverse_complement(right) or adapter != hairpin
                             for right, adapter, left, full in documented):
        raise ContractError("Documented barcoded adapters do not follow barcode-hairpin-reverse complement with this hairpin")
    notes["documented_full_adapter"] = "yes" if any(right == barcode for right, *_ in documented) else "no"
    return Derivation(barcode + hairpin, [("barcode", len(barcode), name), ("hairpin", len(hairpin), accession)],
                      "search_representation", barcode, origin, notes)


def rule_catalogue_record(library, entry, sources):
    sequence = fasta_record(library.parsed("fasta", sources["catalogue"]), entry["element"])
    return Derivation(sequence, [("oligo", len(sequence), "")], "verified", sequence, "catalogue")


RULES = {
    "ont_native_barcode": (rule_ont_native_barcode, ("ont", "dorado")),
    "ont_native_context_k40": (rule_ont_native_context_k40, ("ont", "dorado")),
    "ont_native_context_k46": (rule_ont_native_context_k46, ("ont", "dorado")),
    "ont_rapid_context_90": (rule_ont_rapid_context_90, ("dorado", "ont")),
    "porechop_adapter": (rule_porechop_adapter, ("porechop", "ont")),
    "illumina_index": (rule_illumina_index, ("table",)),
    "illumina_index_adapter": (rule_illumina_index_adapter, ("table", "formula")),
    "illumina_formula_segment": (rule_illumina_formula_segment, ("formula",)),
    "documented_oligo": (rule_documented_oligo, ("document",)),
    "mgi_filter_adapter": (rule_mgi_filter_adapter, ("manual",)),
    "pacbio_barcode": (rule_pacbio_barcode, ("official", "mirror")),
    "pacbio_barcode_hairpin": (rule_pacbio_barcode_hairpin, ("official", "mirror", "hairpin", "layout")),
    "catalogue_record": (rule_catalogue_record, ("catalogue",)),
}


SEGMENT = r"[a-z_]+:\d+-\d+(?:=[^;=]+)?"
SOURCE_REFERENCE = r"[a-z_]+=[A-Za-z0-9][A-Za-z0-9_.+-]*"
DEFINITION_SCHEMA = dict(type="object", additionalProperties=False, required=[
    "query_id", "group", "condition", "family", "construct_class", "length", "sequence_sha256", "kit_or_set", "orientation",
    "segments", "transformation", "rule", "element", "sources", "verification", "index_name", "index_side"], properties={
    "query_id": {"type": "string", "pattern": r"^[A-Za-z0-9_.:+-]+$"},
    **{name: {"type": "string", "minLength": 1} for name in ("group", "condition", "family", "kit_or_set", "orientation", "element")},
    "construct_class": {"enum": ["barcode", "adapter", "primer", "composite"]},
    "length": {"type": "integer", "minimum": 1},
    "sequence_sha256": {"type": "string", "pattern": "^[a-f0-9]{64}$"},
    "segments": {"type": "string", "pattern": f"^{SEGMENT}(?:;{SEGMENT})*$"},
    "transformation": {"enum": ["none", "U_to_T", "A_tail_base"]},
    "rule": {"enum": sorted(RULES)},
    "sources": {"type": "string", "pattern": f"^{SOURCE_REFERENCE}(?:;{SOURCE_REFERENCE})*$"},
    "verification": {"enum": list(VERIFICATIONS)},
    "index_name": {"type": "string"}, "index_side": {"enum": ["", "i7", "i5"]}})


def metadata_schema() -> dict:
    return schema("queries")


def metadata_row(entry: dict, derived: Derivation, sources: dict) -> dict:
    """Downstream query-metadata row (schemas/queries.json); adapt here to contract changes.

    family_id follows the panel importer: canonical ID of the barcode/index a query carries, or of the oligo itself.
    Every composite rule checks the architecture of its query against the documented layout, so derived composites
    carry a verified architecture; the source verification level stays in the definition and panel_build.json.
    """
    composite = entry["construct_class"] == "composite"
    return dict(query_id=entry["query_id"], group=entry["group"], condition=entry["condition"],
                construct_class=entry["construct_class"], family_id=canonical_id(derived.core), source_id=sources[derived.origin],
                provenance_status="verified", architecture_status="verified" if composite else "not_applicable")


def source_roles(entry: dict) -> tuple[dict, list]:
    roles, references = {}, []
    for item in entry["sources"].split(";"):
        role, _, name = item.partition("=")
        if role == "reference":
            references.append(name)
        elif role in roles:
            raise ContractError(f"{entry['query_id']}: source role {role} is declared twice")
        else:
            roles[role] = name
    required = RULES[entry["rule"]][1]
    if set(roles) != set(required):
        raise ContractError(f"{entry['query_id']}: rule {entry['rule']} needs exactly the source roles {', '.join(required)} "
                            "(plus optional reference=)")
    return roles, references


def read_definition(path: Path) -> list[dict]:
    entries = read_table(path, DEFINITION_SCHEMA)
    if not entries or len({entry["query_id"] for entry in entries}) != len(entries):
        raise ContractError(f"{path.name}: empty definition or duplicate query IDs")
    for entry in entries:
        if entry["transformation"] != "none" and entry["rule"] != "documented_oligo":
            raise ContractError(f"{entry['query_id']}: transformations apply only to documented oligos")
        index_rule = entry["rule"] in ("illumina_index", "illumina_index_adapter")
        if bool(entry["index_name"] and entry["index_side"]) != index_rule or (not index_rule and (entry["index_name"] or entry["index_side"])):
            raise ContractError(f"{entry['query_id']}: index_name/index_side are required for, and only for, Illumina index rules")
    return entries


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("x", encoding="ascii", newline="\n") as handle:
        handle.writelines(f">{name}\n{sequence}\n" for name, sequence in records.items())


def build_panel(definition: Path, sources_dir: Path, output_dir: Path, manifest: Path | None = None,
                pdftotext: str | None = None) -> dict:
    """Derive every query of `definition` from the pinned sources (manifest default: sources.json beside it).

    Writes queries.fasta, queries.tsv (downstream query metadata) and panel_build.json. Nothing is written
    if a source is missing or changed or any derived query disagrees with the definition. `pdftotext` names the
    executable used for PDF sources (default: found on PATH).
    """
    manifest = manifest or definition.with_name("sources.json")
    if output_dir.exists():
        raise ContractError(f"Output directory already exists: {output_dir}")
    pinned = load_manifest(manifest)
    entries = read_definition(definition)
    roles = {entry["query_id"]: source_roles(entry) for entry in entries}
    used = sorted({name for named, references in roles.values() for name in [*named.values(), *references]})
    undeclared = [name for name in used if name not in pinned]
    if undeclared:
        raise ContractError("Definition references sources missing from the manifest: " + ", ".join(undeclared))
    verify_sources(sources_dir, [pinned[name] for name in used])
    library = Library(sources_dir, pinned, pdftotext)
    sequences, rows, details, problems = {}, [], [], []
    for entry in entries:
        named = roles[entry["query_id"]][0]
        try:
            derived = RULES[entry["rule"]][0](library, entry, named)
        except ContractError as error:
            problems.append(f"{entry['query_id']}: {error}")
            continue
        observed = hashlib.sha256(derived.sequence.encode("ascii")).hexdigest()
        differences = [f"{label} {value} != {expected}" for label, value, expected in (
            ("length", len(derived.sequence), entry["length"]), ("sequence_sha256", observed, entry["sequence_sha256"]),
            ("segments", derived.segments, entry["segments"]), ("verification", derived.verification, entry["verification"]))
            if value != expected]
        if differences:
            problems.append(f"{entry['query_id']}: derived " + "; ".join(differences))
            continue
        sequences[entry["query_id"]] = derived.sequence
        rows.append(metadata_row(entry, derived, named))
        details.append(dict(query_id=entry["query_id"], rule=entry["rule"], element=entry["element"], sequence_sha256=observed,
                            canonical_id=canonical_id(derived.sequence), source=named[derived.origin], notes=derived.notes))
    if problems:
        raise ContractError(f"{len(problems)} of {len(entries)} queries failed derivation or verification: "
                            + " | ".join(problems[:25]) + (" | ..." if len(problems) > 25 else ""))
    output_dir.mkdir(parents=True)
    write_fasta(output_dir / "queries.fasta", sequences)
    write_table(output_dir / "queries.tsv", rows, metadata_schema())

    def counts(column):
        return dict(Counter(entry[column] for entry in entries))

    report = dict(
        kind="panel_build", research_freeze=False,
        scope="Reference sequences derived from pinned public documents; not a physical authentication of reagents",
        definition=str(definition), definition_sha256=digest(definition), manifest=str(manifest), manifest_sha256=digest(manifest),
        sources={name: {key: pinned[name][key] for key in ("sha256", "bytes", "url", "retrieved")} for name in used},
        text_extracts=library.extracts, pdftotext=pdftotext_version(pdftotext) if library.extracts else "",
        queries=len(sequences), distinct_sequences=len(set(sequences.values())),
        canonical_parents=len({canonical_id(sequence) for sequence in sequences.values()}),
        by_group=counts("group"), by_family=counts("family"), by_condition=counts("condition"), by_rule=counts("rule"),
        by_verification=counts("verification"),
        files={name: digest(output_dir / name) for name in ("queries.fasta", "queries.tsv")}, derivations=details)
    save_json(output_dir / "panel_build.json", report)
    return report
