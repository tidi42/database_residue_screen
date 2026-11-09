"""Code identity: per-stage import closures, package files, git state, interpreter, platform, dependencies and tools.

Every stage request binds the digest of the stage's code closure — the module defining its run function plus every
module of the package that it imports, transitively (static analysis of absolute and relative imports) — together with
the identities of the tools the stage declares, the Python version and implementation, the dependency versions and the
coarse platform. The package-wide digests and the git commit and tree are recorded, not bound.
"""
import ast
import functools
import hashlib
import importlib.metadata as metadata
import importlib.util
import os
import platform
import re
import shutil
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .common import StudyError, file_sha256, json_sha256, listing_sha256

PACKAGE_DIRECTORY = Path(__file__).resolve().parent.parent
MODES = ("production", "working")
EXCLUDED_EXTRAS = frozenset({"test", "tests", "testing", "dev", "develop", "docs", "doc", "lint"})
RESOURCE_MODULES = ("importlib.resources", "pkgutil")
_COMPILED = (".pyc", ".pyo")
_PARSED: dict[str, tuple] = {}


def run_module(function) -> str:
    """Name of the module defining a stage's run function (the start of its code closure)."""
    while isinstance(function, functools.partial):
        function = function.func
    module = getattr(function, "__module__", None)
    if not isinstance(module, str) or not module:
        raise StudyError(f"Cannot determine the module that defines {function!r}")
    return module


def package_files(package_dir: Path) -> tuple[dict, dict]:
    """SHA-256 of every module (.py) and every other package file, keyed by POSIX path relative to the package."""
    modules, data = {}, {}
    for root, directories, files in os.walk(package_dir):
        directories[:] = sorted(name for name in directories if name != "__pycache__" and not name.startswith("."))
        for name in sorted(files):
            if name.startswith(".") or name.endswith(_COMPILED):
                continue
            path = Path(root, name)
            (modules if name.endswith(".py") else data)[path.relative_to(package_dir).as_posix()] = file_sha256(path)
    return modules, data


class CodeIndex:
    """Static import graph of one package.

    ``closure(module)`` maps the module and every package module it imports, transitively, to their SHA-256, keyed by
    POSIX path relative to the import root (``database_residue_screen/search.py``). Imports inside functions count;
    parent-package ``__init__`` modules count only when imported explicitly; dynamic imports are not seen. If a module
    of the closure reads package resources (importlib.resources or pkgutil), every package data file is included."""

    def __init__(self, package_dir: Path):
        self.package_dir = Path(package_dir).resolve()
        self.package = self.package_dir.name
        self._modules: dict[str, tuple | None] = {}
        self._located: dict[str, Path | None] = {}
        self._data: dict[str, str] | None = None

    def inside(self, name: str) -> bool:
        return name == self.package or name.startswith(self.package + ".")

    def locate(self, name: str) -> Path | None:
        if name not in self._located:
            self._located[name] = self._find(name)
        return self._located[name]

    def _find(self, name: str) -> Path | None:
        if self.inside(name):
            base = self.package_dir.parent.joinpath(*name.split("."))
            for candidate in (base / "__init__.py", base.parent / (base.name + ".py")):
                if candidate.is_file():
                    return candidate
            return None
        location = getattr(sys.modules.get(name), "__file__", None)
        if location is None:
            try:
                spec = importlib.util.find_spec(name)
            except (ImportError, ValueError):
                spec = None
            location = spec.origin if spec is not None and spec.has_location else None
        return Path(location).resolve() if location and Path(location).is_file() else None

    def module(self, name: str) -> tuple | None:
        if name not in self._modules:
            path = self.locate(name)
            if path is None:
                self._modules[name] = None
            else:
                source = path.read_bytes()
                digest = hashlib.sha256(source).hexdigest()
                key = "/".join(name.split(".")) + ("/__init__.py" if path.name == "__init__.py" else ".py")
                self._modules[name] = (key, digest, *self._imports(name, path, _import_statements(source, digest, path)))
        return self._modules[name]

    def _imports(self, name: str, path: Path, statements: tuple) -> tuple[tuple[str, ...], bool]:
        package = name.split(".") if path.name == "__init__.py" else name.split(".")[:-1]
        found, resources = set(), False
        for level, module, names in statements:
            if module is None and not level:
                for alias in names:
                    resources = resources or alias.startswith(RESOURCE_MODULES)
                    if self.inside(alias):
                        found.add(alias)
                continue
            if level:
                if level - 1 > len(package):
                    continue
                base = ".".join(package[:len(package) - level + 1] + (module.split(".") if module else []))
            else:
                base = module
                resources = resources or base.startswith(RESOURCE_MODULES) or (base == "importlib" and "resources" in names)
            if not base or not self.inside(base):
                continue
            for alias in names:
                candidate = f"{base}.{alias}"
                found.add(candidate if alias != "*" and self.locate(candidate) is not None else base)
        return tuple(sorted(found)), resources

    def data_files(self) -> dict[str, str]:
        if self._data is None:
            self._data = {f"{self.package}/{name}": digest for name, digest in package_files(self.package_dir)[1].items()}
        return self._data

    def closure(self, start: str) -> dict[str, str]:
        result, pending, seen, resources = {}, [start], set(), False
        while pending:
            name = pending.pop()
            if name in seen:
                continue
            seen.add(name)
            entry = self.module(name)
            if entry is None:
                if name == start:
                    raise StudyError(f"Cannot locate the source of module {start!r}")
                continue
            key, digest, imports, reads = entry
            result[key] = digest
            resources = resources or reads
            pending.extend(imports)
        if resources:
            result.update(self.data_files())
        return dict(sorted(result.items()))


def _import_statements(source: bytes, digest: str, path: Path) -> tuple:
    """(level, module, names) of every import statement; ``module`` is None for plain ``import`` statements. Cached by
    content digest, so re-checking an unchanged closure only re-reads and hashes its files."""
    if digest not in _PARSED:
        try:
            tree = ast.parse(source, filename=str(path))
        except SyntaxError as error:
            raise StudyError(f"Cannot parse {path}: {error}") from error
        _PARSED[digest] = tuple(
            (0, None, tuple(alias.name for alias in node.names)) if isinstance(node, ast.Import)
            else (node.level, node.module or "", tuple(alias.name for alias in node.names))
            for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom)))
    return _PARSED[digest]


def stage_closure(index: CodeIndex, run) -> dict:
    module = run_module(run)
    closure = index.closure(module)
    return dict(module=module, closure=closure, closure_sha256=listing_sha256(closure.items()))


def package_version(package_dir: Path) -> str | None:
    try:
        tree = ast.parse((package_dir / "__init__.py").read_text(encoding="utf-8"))
    except (OSError, SyntaxError, UnicodeDecodeError):
        return None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                return node.value.value
    return None


def _normalise(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def _requirement(text: str, extras: bool) -> str | None:
    requirement, _, marker = text.partition(";")
    extra = re.search(r"""extra\s*==\s*["']([^"']+)["']""", marker)
    if extra and (not extras or _normalise(extra.group(1)) in EXCLUDED_EXTRAS):
        return None
    match = re.match(r"\s*([A-Za-z0-9][A-Za-z0-9._-]*)", requirement)
    return _normalise(match.group(1)) if match else None


def declared_requirements(package_dir: Path) -> list[str]:
    """Runtime and optional (non-test) requirements from the source tree's pyproject.toml, else installed metadata."""
    candidates = [package_dir.parent / "pyproject.toml"]
    if package_dir.parent.name == "src":
        candidates.append(package_dir.parent.parent / "pyproject.toml")
    for candidate in candidates:
        if candidate.is_file():
            with candidate.open("rb") as handle:
                project = tomllib.load(handle).get("project")
            if isinstance(project, dict):
                texts = list(project.get("dependencies", []))
                for group, items in sorted(project.get("optional-dependencies", {}).items()):
                    if _normalise(group) not in EXCLUDED_EXTRAS:
                        texts.extend(items)
                return sorted({name for name in (_requirement(text, False) for text in texts) if name})
    names = set()
    for distribution in metadata.packages_distributions().get(package_dir.name, []):
        try:
            requirements = metadata.requires(distribution) or []
        except metadata.PackageNotFoundError:
            continue
        names.update(name for name in (_requirement(text, True) for text in requirements) if name)
    return sorted(names)


def dependency_versions(names) -> dict:
    """Installed versions of the named distributions and, transitively, of what they require (None: not installed)."""
    declared, versions, pending = set(names), {}, sorted(names, reverse=True)
    while pending:
        name = pending.pop()
        if name in versions:
            continue
        try:
            distribution = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            if name in declared:
                versions[name] = None
            continue
        versions[name] = distribution.version
        pending.extend(child for child in (_requirement(text, False) for text in distribution.requires or []) if child and child not in versions)
    return dict(sorted(versions.items()))


def environment_identity(package_dir: Path) -> dict:
    """Interpreter, coarse platform and dependency versions (bound into every stage request)."""
    return dict(python=dict(implementation=platform.python_implementation(), version=platform.python_version()),
                platform=dict(system=platform.system(), machine=platform.machine(), libc=" ".join(part for part in platform.libc_ver() if part)),
                dependencies=dependency_versions(declared_requirements(package_dir)))


def _git(arguments: list[str], cwd: Path, environ: Mapping[str, str]):
    executable = shutil.which("git", path=environ.get("PATH", os.defpath))
    if executable is None:
        return None
    environment = {key: value for key, value in environ.items() if not key.startswith("GIT_") or key == "GIT_CEILING_DIRECTORIES"}
    environment.update(GIT_OPTIONAL_LOCKS="0", LC_ALL="C")
    try:
        return subprocess.run([executable, *arguments], cwd=cwd, env=environment, capture_output=True, timeout=300, check=False,
                              stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None


def _porcelain_paths(output: bytes) -> list[str]:
    tokens, paths, index = output.decode("utf-8", "replace").split("\0"), [], 0
    while index < len(tokens):
        token = tokens[index]
        index += 1
        if len(token) > 3:
            paths.append(token[3:])
            if "R" in token[:2] or "C" in token[:2]:
                index += 1
    return paths


def git_state(package_dir: Path, environ: Mapping[str, str], files=()) -> dict:
    """Commit, tree, cleanliness and file tracking of the git working tree containing the package."""
    top = _git(["rev-parse", "--show-toplevel"], package_dir, environ)
    if top is None:
        return dict(available=False, reason="git executable not found")
    if top.returncode != 0:
        message = top.stderr.decode("utf-8", "replace").strip().splitlines()
        return dict(available=False, reason=message[0] if message else "not inside a git working tree")
    root = Path(top.stdout.decode("utf-8").strip()).resolve()
    commit = _git(["rev-parse", "--verify", "--quiet", "HEAD^{commit}"], package_dir, environ)
    if commit is None or commit.returncode != 0:
        return dict(available=False, reason="git repository has no commits", root=str(root))
    tree = _git(["rev-parse", "--verify", "--quiet", "HEAD^{tree}"], package_dir, environ)
    status = _git(["status", "--porcelain=v1", "-z", "--untracked-files=all"], root, environ)
    listed = _git(["ls-files", "-z", "--cached"], package_dir, environ)
    if tree is None or status is None or listed is None or tree.returncode or status.returncode or listed.returncode:
        raise StudyError(f"Cannot read the git state of {root}")
    dirty = _porcelain_paths(status.stdout)
    return dict(available=True, root=str(root), commit=commit.stdout.decode().strip(), tree=tree.stdout.decode().strip(),
                clean=not dirty, dirty_paths=dirty, untracked_files=sorted(set(files) - set(listed.stdout.decode("utf-8", "replace").split("\0"))))


def git_ignored(root: Path, path: Path, environ: Mapping[str, str]) -> bool:
    result = _git(["check-ignore", "--quiet", "--no-index", path.relative_to(root).as_posix()], root, environ)
    return result is not None and result.returncode == 0


def tool_identity(name: str, specification: Mapping, environ: Mapping[str, str]) -> dict:
    """Resolved path, SHA-256 and version output of one configured external program."""
    command = specification["command"]
    if os.sep in command or (os.altsep and os.altsep in command):
        found = command if os.path.isfile(command) and os.access(command, os.X_OK) else None
    else:
        found = shutil.which(command, path=environ.get("PATH", os.defpath))
    if found is None:
        raise StudyError(f"Configured tool {name!r} was not found or is not executable: {command}")
    path = os.path.abspath(found)
    try:
        completed = subprocess.run([path, *specification["version_args"]], capture_output=True, timeout=300, check=False,
                                   env=dict(environ), stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError) as error:
        raise StudyError(f"Cannot run the version command of tool {name!r}: {error}") from error
    output = "\n".join(text for text in (completed.stdout.decode("utf-8", "replace").strip(),
                                         completed.stderr.decode("utf-8", "replace").strip()) if text)
    if completed.returncode != 0:
        raise StudyError(f"Version command of tool {name!r} failed with exit code {completed.returncode}: {output[:300]}")
    resolved = Path(path).resolve()
    return dict(command=command, path=path, resolved_path=str(resolved), sha256=file_sha256(resolved),
                version_args=list(specification["version_args"]), version=output)


@dataclass(frozen=True, eq=False)
class CodeIdentity:
    """``recorded``: mode, package-wide digests, git state and pipeline (not bound). ``environment`` and ``tools`` feed
    the per-stage code block of each request, together with the stage's closure digest."""
    mode: str
    package_dir: Path
    recorded: dict
    environment: dict
    tools: dict
    files: dict
    host: dict
    git_root: Path | None
    index: CodeIndex

    @property
    def sha256(self) -> str:
        return json_sha256(self.recorded)

    def stage_code(self, run, tools: Iterable[str] = (), fresh: bool = False) -> tuple[dict, dict]:
        """(bound code block of a stage request, recorded closure details). ``fresh`` re-reads the source files."""
        code = stage_closure(CodeIndex(self.package_dir) if fresh else self.index, run)
        bound = dict(closure_sha256=code["closure_sha256"], tools={name: self.tools[name] for name in sorted(tools)}, **self.environment)
        return bound, dict(module=code["module"], closure=code["closure"])


def collect_code_identity(mode: str, tools: Mapping[str, Mapping] = (), *, package_dir: Path | None = None, pipeline=None,
                          environ: Mapping[str, str] | None = None) -> CodeIdentity:
    """Describe the code that will run. Production mode refuses a dirty, untracked or foreign code tree and stage code
    defined outside the package."""
    if mode not in MODES:
        raise StudyError(f"Unknown mode {mode!r}; expected production or working")
    environ = os.environ if environ is None else environ
    package_dir = Path(package_dir or PACKAGE_DIRECTORY).resolve()
    if not (package_dir / "__init__.py").is_file():
        raise StudyError(f"Not a Python package directory: {package_dir}")
    index = CodeIndex(package_dir)
    modules, data = package_files(package_dir)
    git = git_state(package_dir, environ, [*modules, *data])
    foreign = sorted(name for name, module in list(sys.modules.items()) if index.inside(name)
                     and getattr(module, "__file__", None) and not Path(module.__file__).resolve().is_relative_to(package_dir))
    pipeline_record, outside = None, []
    if pipeline is not None:
        location = getattr(sys.modules.get(pipeline.module), "__file__", None)
        source = Path(location).resolve() if location else None
        pipeline_record = dict(name=pipeline.name, module=pipeline.module, source_sha256=file_sha256(source) if source else None,
                               in_package=bool(source and source.is_relative_to(package_dir)))
        outside = sorted({name for name in (run_module(stage.run) for stage in pipeline.ordered()) if not index.inside(name)})
    recorded = dict(
        mode=mode,
        package=dict(name=package_dir.name, version=package_version(package_dir), modules_sha256=listing_sha256(modules.items()),
                     module_count=len(modules), data_sha256=listing_sha256(data.items()), data_file_count=len(data)),
        git=dict(available=True, commit=git["commit"], tree=git["tree"], clean=git["clean"]) if git["available"]
        else dict(available=False, reason=git["reason"]),
        pipeline=pipeline_record)
    if mode == "production":
        problems = []
        if not git["available"]:
            problems.append(f"the package directory {package_dir} is not in a usable git working tree ({git['reason']})")
        else:
            if git["untracked_files"]:
                problems.append("the package is not importable from the git tree: package files not tracked by git: "
                                + ", ".join(git["untracked_files"][:10]))
            if not git["clean"]:
                problems.append(f"the git working tree {git['root']} is dirty ({len(git['dirty_paths'])} changed or untracked "
                                f"paths, e.g. {', '.join(git['dirty_paths'][:5])})")
        if foreign:
            problems.append("package modules were imported from outside the package directory: " + ", ".join(foreign[:10]))
        if pipeline_record is not None and not pipeline_record["in_package"]:
            problems.append(f"pipeline {pipeline_record['name']!r} is defined in {pipeline_record['module']}, outside the package")
        if outside:
            problems.append("stage functions are defined outside the package in " + ", ".join(outside))
        if problems:
            raise StudyError("Production mode refused: " + "; ".join(problems)
                             + ". Commit the code and run from that clean checkout, or use working mode (--working).")
    tool_records = {name: tool_identity(name, specification, environ) for name, specification in sorted(dict(tools).items())}
    host = dict(platform=platform.platform(), kernel_release=platform.release(), python_executable=sys.executable,
                python_build=sys.version, package_directory=str(package_dir), git_root=git.get("root"),
                git_dirty_paths=git.get("dirty_paths", [])[:200], foreign_package_modules=foreign)
    return CodeIdentity(mode, package_dir, recorded, environment_identity(package_dir), tool_records, dict(modules=modules, data=data),
                        host, Path(git["root"]) if git.get("root") else None, index)
