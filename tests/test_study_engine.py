import contextlib
import fcntl
import hashlib
import io
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
import uuid
from pathlib import Path
from unittest import mock

from database_residue_screen.cli import main
from database_residue_screen.study import (REDACTED, CodeIndex, Pipeline, Stage, StudyError, collect_code_identity, load_config,
                                           load_pipeline, register_pipeline, register_pipeline_reference, run_study, study_status)
from database_residue_screen.study.common import json_sha256
from database_residue_screen.study.demo import PIPELINE as DEMO, write_demo_study

CONTROL = {"failures": 0}
SEEN = {}
UTC_LINE = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z (INFO|WARNING|ERROR) ")


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(directory: Path, *arguments: str) -> str:
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    result = subprocess.run(["git", "-c", "user.name=Study Test", "-c", "user.email=study@example.org", "-c", "commit.gpgsign=false",
                             *arguments], cwd=directory, env=environment, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def new_pipeline(label: str, **options) -> Pipeline:
    return register_pipeline(Pipeline(f"test-{label}-{uuid.uuid4().hex[:8]}", **options))


def write_package(root: Path, name: str, files: dict) -> Path:
    package = root / name
    for relative, text in files.items():
        path = package / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(text))
    return package


def forget_modules(prefix: str) -> None:
    for name in [name for name in sys.modules if name.split(".")[0] == prefix]:
        del sys.modules[name]


def make_writable(root: Path) -> None:
    for directory, _, files in os.walk(root):
        os.chmod(directory, stat.S_IMODE(os.lstat(directory).st_mode) | 0o700)
        for name in files:
            path = os.path.join(directory, name)
            if not os.path.islink(path):
                os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) | 0o600)


class StudyCase(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.root = Path(temporary.name).resolve()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(make_writable, self.root)
        ceilings = [value for value in (os.environ.get("GIT_CEILING_DIRECTORIES"), str(self.root)) if value]
        self.environ = dict(os.environ, GIT_CEILING_DIRECTORIES=os.pathsep.join(ceilings))
        CONTROL["failures"] = 0
        SEEN.clear()

    def run_study(self, config, run, **options):
        options.setdefault("mode", "working")
        return run_study(config, run, environ=options.pop("environ", self.environ), **options)

    def demo(self, name="input", **options) -> Path:
        return write_demo_study(self.root / name, **options)

    def write_config(self, pipeline: Pipeline, body: str, *, name="study.toml", study_id="test") -> Path:
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f'[study]\nid = "{study_id}"\ntitle = "Test study"\npipeline = "{pipeline.name}"\n'
        path.write_text(header + textwrap.dedent(body), encoding="utf-8")
        return path

    def snapshot(self, run: Path) -> dict:
        files = {}
        for path in sorted(run.rglob("*")):
            if path.is_file() and path.parent != run:
                files[path.relative_to(run).as_posix()] = (path.stat().st_mtime_ns, sha256(path))
        return files

    def actions(self, result) -> dict:
        return {row["stage"]: row["action"] for row in result["stages"]}


class ConfigTests(StudyCase):
    def load(self, text, sections=None, environ=None, name="study.toml", secrets=()):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(text), encoding="utf-8")
        return load_config(path, DEMO.sections() if sections is None else sections, environ=environ or {}, secrets=secrets)

    def test_defaults_path_resolution_and_canonical_hash(self):
        config = self.load("""
            [study]
            id = "demo"
            title = "Demo"
            pipeline = "demo"
            [numbers]
            source = "data/numbers.txt"
        """)
        self.assertEqual(config.values["numbers"]["source"], str(self.root / "data/numbers.txt"))
        self.assertEqual(config.values["scaling"], {"factor": 1})
        self.assertEqual(config.values["report"], {"heading": "Demonstration report"})
        self.assertEqual(config.values["tools"], {})
        canonical = json.dumps(config.recorded, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
        self.assertEqual(config.sha256, hashlib.sha256(canonical).hexdigest())
        reordered = self.load("""
            # same content, other order, explicit defaults and comments
            [scaling]
            factor = 1
            [numbers]
            source = "./data/../data/numbers.txt"
            [study]
            pipeline = "demo"
            title = "Demo"
            id = "demo"
            [report]
            heading = "Demonstration report"
        """, name="other.toml")
        self.assertEqual(reordered.sha256, config.sha256)
        self.assertNotEqual(reordered.file_sha256, config.file_sha256)

    def test_relative_paths_resolve_against_the_config_directory(self):
        absolute = self.root / "elsewhere" / "n.txt"
        text = """
            [study]
            id = "demo"
            title = "Demo"
            pipeline = "demo"
            [numbers]
            source = "../data/numbers.txt"
        """
        with contextlib.chdir(self.root):
            first = self.load(text, name="nested/deeper/study.toml")
        second = self.load(text, name="nested/deeper/study.toml")
        self.assertEqual(first.values["numbers"]["source"], str(self.root / "nested/data/numbers.txt"))
        self.assertEqual(second.sha256, first.sha256)
        kept = self.load(text.replace("../data/numbers.txt", str(absolute)))
        self.assertEqual(kept.values["numbers"]["source"], str(absolute))

    def test_unknown_sections_and_keys_are_errors(self):
        base = '[study]\nid = "demo"\ntitle = "Demo"\npipeline = "demo"\n[numbers]\nsource = "n.txt"\n'
        for extra, expected in (('[serach]\nthreads = 4\n', "unknown section(s) [serach]"),
                                ('[scaling]\nfactr = 2\n', "unknown key(s) in [scaling]: factr"),
                                ('[report]\nheading = "x"\ncolour = "red"\n', "unknown key(s) in [report]: colour")):
            with self.subTest(extra=extra), self.assertRaises(StudyError) as caught:
                self.load(base + extra)
            self.assertIn(expected, str(caught.exception))
        with self.assertRaises(StudyError) as caught:
            self.load(base.replace('pipeline = "demo"', 'pipeline = "demo"\nowner = "me"'))
        self.assertIn("unknown key(s) in [study]: owner", str(caught.exception))

    def test_missing_and_invalid_values_are_errors(self):
        cases = (
            ('[study]\nid = "demo"\npipeline = "demo"\n[numbers]\nsource = "n"\n', "'title' is a required property"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n', "'numbers' is a required property"),
            ('[study]\nid = "has space"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = "n"\n', "study.id"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = "n"\n[scaling]\nfactor = 0\n', "scaling.factor"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = 5\n', "numbers.source"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = ""\n', "empty path"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = "n"\nwhen = 2026-01-01\n', "dates"),
            ('[study]\nid = "demo"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = "n"\n[scaling]\nfactor = nan\n', "non-finite"),
            ('[study\nid = "demo"\n', "invalid TOML"),
        )
        for text, expected in cases:
            with self.subTest(expected=expected), self.assertRaises(StudyError) as caught:
                self.load(text)
            self.assertIn(expected, str(caught.exception))

    def test_environment_substitution(self):
        text = """
            [study]
            id = "demo"
            title = "Costs $$5 for ${WHO}"
            pipeline = "demo"
            [numbers]
            source = "${DATA}/numbers.txt"
        """
        config = self.load(text, environ=dict(DATA=str(self.root / "d"), WHO="alice"))
        self.assertEqual(config.values["numbers"]["source"], str(self.root / "d/numbers.txt"))
        self.assertEqual(config.values["study"]["title"], "Costs $5 for alice")
        self.assertEqual(config.recorded, config.values)
        self.assertEqual(config.environment, dict(DATA=str(self.root / "d"), WHO="alice"))
        for environ, expected in ((dict(WHO="alice"), "numbers.source: environment variable DATA is not set"),
                                  (dict(WHO="alice", DATA=""), "environment variable DATA is empty")):
            with self.subTest(expected=expected), self.assertRaises(StudyError) as caught:
                self.load(text, environ=environ)
            self.assertIn(expected, str(caught.exception))
        for broken, expected in (("${1BAD}", "invalid environment variable reference"), ("${DATA", "unterminated")):
            with self.subTest(broken=broken), self.assertRaises(StudyError) as caught:
                self.load(text.replace("${DATA}", broken), environ=dict(DATA="x", WHO="y"))
            self.assertIn(expected, str(caught.exception))

    def test_secrets_are_substituted_but_recorded_redacted(self):
        sections = {"remote": {"type": "object", "properties": {"url": {"type": "string"}, "key": {"type": "string", "pattern": "^[a-z0-9]+$"},
                                                                "code": {"type": "string"}, "plain": {"type": "string"}, "who": {"type": "string"}}}}
        text = """
            [study]
            id = "demo"
            title = "Demo"
            pipeline = "demo"
            secrets = ["PLAIN_VALUE"]
            [remote]
            url = "https://example.org/x?api_key=${NCBI_API_KEY}"
            key = "${NCBI_API_KEY}"
            code = "${REMOTE_CODE}"
            plain = "${PLAIN_VALUE}"
            who = "${WHO}"
        """
        environ = dict(NCBI_API_KEY="abc123secret", REMOTE_CODE="code-77", PLAIN_VALUE="plain-99", WHO="alice")
        config = self.load(text, sections, environ, secrets=("REMOTE_CODE",))
        self.assertEqual(config.values["remote"], dict(url="https://example.org/x?api_key=abc123secret", key="abc123secret",
                                                       code="code-77", plain="plain-99", who="alice"))
        self.assertEqual(config.recorded["remote"], dict(url=f"https://example.org/x?api_key={REDACTED}", key=REDACTED, code=REDACTED,
                                                         plain=REDACTED, who="alice"))
        self.assertEqual(config.environment, dict(NCBI_API_KEY=REDACTED, PLAIN_VALUE=REDACTED, REMOTE_CODE=REDACTED, WHO="alice"))
        rotated = self.load(text, sections, dict(environ, NCBI_API_KEY="zzz999other"), secrets=("REMOTE_CODE",))
        self.assertEqual(rotated.sha256, config.sha256)
        renamed = self.load(text, sections, dict(environ, WHO="bob"), secrets=("REMOTE_CODE",))
        self.assertNotEqual(renamed.sha256, config.sha256)
        with self.assertRaises(StudyError) as caught:
            self.load(text, sections, dict(environ, NCBI_API_KEY="UPPER123secret"), secrets=("REMOTE_CODE",))
        self.assertIn(REDACTED, str(caught.exception))
        self.assertNotIn("UPPER123secret", str(caught.exception))

    def test_tools_are_normalised(self):
        config = self.load("""
            [study]
            id = "demo"
            title = "Demo"
            pipeline = "demo"
            [numbers]
            source = "n.txt"
            [tools]
            blastn = "blastn"
            local = "bin/tool.sh"
            Rscript = { command = "Rscript", version_args = ["--version"] }
        """)
        self.assertEqual(config.values["tools"], {
            "Rscript": {"command": "Rscript", "version_args": ["--version"]},
            "blastn": {"command": "blastn", "version_args": ["-version"]},
            "local": {"command": str(self.root / "bin/tool.sh"), "version_args": ["--version"]}})
        with self.assertRaises(StudyError):
            self.load('[study]\nid = "d"\ntitle = "D"\npipeline = "demo"\n[numbers]\nsource = "n"\n[tools]\nx = { comand = "y" }\n')


class PipelineTests(unittest.TestCase):
    def test_topological_order_uses_registration_order_for_ties(self):
        self.assertEqual([stage.name for stage in DEMO.ordered()], ["ingest", "scale", "describe", "report"])
        pipeline = Pipeline("order-check")
        for name, dependencies in (("c", ("a",)), ("a", ()), ("b", ()), ("d", ("c", "b"))):
            pipeline.register(Stage(name, lambda context: None, dependencies))
        self.assertEqual([stage.name for stage in pipeline.ordered()], ["a", "c", "b", "d"])
        self.assertEqual(pipeline.ancestors("d"), {"a", "b", "c"})

    def test_invalid_pipelines_are_rejected(self):
        cyclic = Pipeline("cyclic")
        cyclic.register(Stage("a", lambda context: None, ("b",)))
        cyclic.register(Stage("b", lambda context: None, ("a",)))
        with self.assertRaisesRegex(StudyError, "dependency cycle"):
            cyclic.ordered()
        dangling = Pipeline("dangling")
        dangling.register(Stage("a", lambda context: None, ("missing",)))
        with self.assertRaisesRegex(StudyError, "unknown stage"):
            dangling.ordered()
        undeclared = Pipeline("undeclared")
        undeclared.register(Stage("a", lambda context: None, config_sections=("search",)))
        with self.assertRaisesRegex(StudyError, r"undeclared section \[search\]"):
            undeclared.ordered()
        with self.assertRaisesRegex(StudyError, "no stages"):
            Pipeline("empty").ordered()
        for invalid in (lambda: Stage("Bad-Name", lambda context: None), lambda: Stage("a", lambda context: None, ("a",)),
                        lambda: Stage("a", "not callable"), lambda: Pipeline("x").section("study", {}),
                        lambda: Pipeline("x").section("s", {"type": "no-such-type"}), lambda: Pipeline("bad name")):
            with self.assertRaises(StudyError):
                invalid()
        duplicate = Pipeline("duplicate")
        duplicate.register(Stage("a", lambda context: None))
        with self.assertRaisesRegex(StudyError, "already registered"):
            duplicate.register(Stage("a", lambda context: None))

    def test_selection(self):
        self.assertEqual(DEMO.select(), ({"ingest", "scale", "describe", "report"}, set()))
        self.assertEqual(DEMO.select(until="scale"), ({"ingest", "scale"}, set()))
        self.assertEqual(DEMO.select(only=["report"]), ({"report"}, {"ingest", "scale", "describe"}))
        self.assertEqual(DEMO.select(only=["describe", "report"]), ({"describe", "report"}, {"ingest", "scale"}))
        with self.assertRaisesRegex(StudyError, "mutually exclusive"):
            DEMO.select(until="scale", only=["report"])
        with self.assertRaisesRegex(StudyError, "Unknown stage"):
            DEMO.select(until="nope")

    def test_pipeline_lookup_and_registration(self):
        self.assertIs(load_pipeline("demo"), DEMO)
        self.assertIs(load_pipeline("database_residue_screen.study.demo:PIPELINE"), DEMO)
        with self.assertRaisesRegex(StudyError, "Unknown pipeline 'nope'.*demo"):
            load_pipeline("nope")
        with self.assertRaisesRegex(StudyError, "Cannot import"):
            load_pipeline("no_such_module_xyz:PIPELINE")
        name = f"alias-{uuid.uuid4().hex[:8]}"
        register_pipeline_reference(name, "database_residue_screen.study.demo:PIPELINE")
        with self.assertRaisesRegex(StudyError, "defines pipeline 'demo'"):
            load_pipeline(name)
        pipeline = new_pipeline("lookup")
        self.assertIs(load_pipeline(pipeline.name), pipeline)
        self.assertEqual(pipeline.module, __name__)
        with self.assertRaisesRegex(StudyError, "already registered"):
            register_pipeline(Pipeline(pipeline.name))


CLOSURE_PACKAGE = {
    "__init__.py": '__version__ = "1.0"\n',
    "pipeline.py": """\
        from database_residue_screen.study import Pipeline

        from . import stage_one, stage_three, stage_two

        PIPELINE = Pipeline("PKG")
        PIPELINE.stage(name="one")(stage_one.run)
        PIPELINE.stage(name="two")(stage_two.run)
        PIPELINE.stage(name="three", depends_on=("two",))(stage_three.run)
        """,
    "stage_one.py": """\
        from .helper_one import double


        def run(context):
            (context.output / "one.txt").write_text(str(double(1)))
        """,
    "helper_one.py": "def double(value):\n    return 2 * value\n",
    "stage_two.py": """\
        import PKG.helper_two
        from .sub import VALUE


        def run(context):
            from . import late
            (context.output / "two.txt").write_text(str(PKG.helper_two.triple(VALUE) + late.OFFSET))
        """,
    "helper_two.py": "def triple(value):\n    return 3 * value\n",
    "sub/__init__.py": "from .deep import VALUE\n",
    "sub/deep.py": "VALUE = 7\n",
    "late.py": "OFFSET = 1\n",
    "stage_three.py": """\
        def run(context):
            (context.output / "three.txt").write_text(context.upstream["two"].path("two.txt").read_text())
        """,
    "unrelated.py": "X = 1\n",
    "reader.py": "from importlib.resources import files\n",
    "data/table.json": "{}\n",
}


class CodeClosureTests(StudyCase):
    def setUp(self):
        super().setUp()
        self.name = f"closurefixture_{uuid.uuid4().hex[:10]}"
        self.package = write_package(self.root / "src", self.name,
                                     {path: text.replace("PKG", self.name) for path, text in CLOSURE_PACKAGE.items()})
        sys.path.insert(0, str(self.root / "src"))
        self.addCleanup(sys.path.remove, str(self.root / "src"))
        self.addCleanup(forget_modules, self.name)
        self.config = self.root / "study.toml"
        self.config.write_text(f'[study]\nid = "closure"\ntitle = "Closure"\npipeline = "{self.name}.pipeline:PIPELINE"\n')
        self.run_directory = self.root / "run"

    def study(self, **options):
        return self.run_study(self.config, self.run_directory, package_dir=self.package, **options)

    def edit(self, relative):
        path = self.package / relative
        path.write_text(path.read_text() + "# edited\n")

    def current(self, stage):
        document = json.loads((self.run_directory / "run.json").read_text())["stages"][stage]
        return next(entry for entry in document["attempts"] if entry["attempt"] == document["current_attempt"])

    def test_closure_follows_relative_absolute_package_and_function_level_imports(self):
        import database_residue_screen
        index, name = CodeIndex(self.package), self.name
        self.assertEqual(set(index.closure(f"{name}.stage_one")), {f"{name}/stage_one.py", f"{name}/helper_one.py"})
        self.assertEqual(set(index.closure(f"{name}.stage_two")), {f"{name}/stage_two.py", f"{name}/helper_two.py",
                                                                    f"{name}/sub/__init__.py", f"{name}/sub/deep.py", f"{name}/late.py"})
        self.assertEqual(set(index.closure(f"{name}.stage_three")), {f"{name}/stage_three.py"})
        self.assertEqual(index.closure(f"{name}.stage_one")[f"{name}/helper_one.py"], sha256(self.package / "helper_one.py"))
        self.assertEqual(set(index.closure(f"{name}.reader")), {f"{name}/reader.py", f"{name}/data/table.json"})
        self.assertLessEqual({f"{name}/pipeline.py", f"{name}/stage_one.py", f"{name}/stage_two.py"}, set(index.closure(f"{name}.pipeline")))
        with self.assertRaisesRegex(StudyError, "Cannot locate the source of module"):
            index.closure(f"{name}.missing")
        package = CodeIndex(Path(database_residue_screen.__file__).parent)
        self.assertEqual(set(package.closure("database_residue_screen.study.demo_stages")), {"database_residue_screen/study/demo_stages.py"})
        self.assertIn("database_residue_screen/study/engine.py", package.closure("database_residue_screen.study.demo"))

    def test_changes_outside_a_stage_closure_keep_its_request_and_changes_inside_change_it(self):
        self.assertEqual(self.actions(self.study()), dict(one="executed", two="executed", three="executed"))
        before, digests = self.snapshot(self.run_directory), {stage: self.current(stage)["request_digest"] for stage in ("one", "two", "three")}
        for relative in ("unrelated.py", "pipeline.py", "__init__.py"):
            self.edit(relative)
        self.assertEqual(self.actions(self.study()), dict(one="skipped", two="skipped", three="skipped"))
        self.assertEqual(self.snapshot(self.run_directory), before)
        self.edit("helper_two.py")
        with self.assertRaises(StudyError) as caught:
            self.study()
        self.assertEqual(str(caught.exception),
                         f"Stage 'two': the request differs from its current attempt 02_two/attempt_1 (completed); changed fields: "
                         f"code.closure_sha256 (changed modules: {self.name}/helper_two.py). Use --rerun-changed to run a new attempt "
                         "(earlier attempts are kept), or use a new run directory.")
        self.assertEqual(self.actions(self.study(rerun_changed=True)), dict(one="skipped", two="executed", three="executed"))
        self.assertEqual(self.current("one")["request_digest"], digests["one"])
        self.assertNotEqual(self.current("two")["request_digest"], digests["two"])
        record = json.loads((self.run_directory / "02_two/attempt_2/stage.json").read_text())
        self.assertEqual(record["code"]["closure"][f"{self.name}/helper_two.py"], sha256(self.package / "helper_two.py"))
        self.assertEqual(record["request"]["code"]["closure_sha256"], json.loads((self.run_directory / "run.json").read_text())
                         ["stages"]["two"]["attempts"][1]["code_closure_sha256"])
        self.assertEqual(json.loads((self.run_directory / "03_three/attempt_2/stage.json").read_text())["request"]["upstream"]["two"]["attempt"],
                         "02_two/attempt_2")
        self.assertEqual({key: value for key, value in self.snapshot(self.run_directory).items() if key in before}, before)
        self.edit("sub/deep.py")
        self.assertEqual(self.actions(self.study(rerun_changed=True)), dict(one="skipped", two="executed", three="executed"))
        self.edit("helper_one.py")
        self.assertEqual(self.actions(self.study(rerun_changed=True)), dict(one="executed", two="skipped", three="skipped"))

    def test_verify_code_checks_every_current_attempt(self):
        from database_residue_screen.study import format_status
        from database_residue_screen.study import status as status_module
        self.study()
        check = dict(verify_code=True, package_dir=self.package, environ=self.environ)
        self.assertTrue(study_status(self.run_directory, **check)["code"]["ok"])
        self.edit("helper_two.py")
        report = study_status(self.run_directory, **check)
        rows = {row["stage"]: row for row in report["code"]["stages"]}
        self.assertFalse(report["code"]["ok"])
        self.assertEqual([(rows[stage]["code"], rows[stage]["upstream"], rows[stage]["ok"]) for stage in ("one", "two", "three")],
                         [("ok", "-", True), ("CHANGED", "-", False), ("ok", "ok", True)])
        self.assertEqual(rows["two"]["details"], [f"changed code: {self.name}/helper_two.py"])
        self.assertRegex(format_status(report), rf"(?m)^two\s+1\s+working\s+CHANGED\s+-\s+ok\s+ok\s+ok\s+-$")
        self.study(rerun_changed=True)
        self.assertTrue(study_status(self.run_directory, **check)["code"]["ok"])
        real = status_module.environment_identity

        def drifted(package_dir):
            environment = real(package_dir)
            environment["dependencies"] = dict(environment["dependencies"], numpy="0.0.1")
            environment["python"] = dict(environment["python"], version="0.0")
            return environment

        with mock.patch.object(status_module, "environment_identity", drifted):
            report = study_status(self.run_directory, **check)
        row = report["code"]["stages"][0]
        self.assertEqual((row["code"], row["dependencies"], row["python"], row["platform"], row["ok"]), ("ok", "CHANGED", "CHANGED", "ok", False))
        self.assertTrue(any(detail.startswith("changed dependencies: numpy ") and detail.endswith(" -> 0.0.1") for detail in row["details"]))


class RunTests(StudyCase):
    def flaky_pipeline(self):
        pipeline = new_pipeline("flaky")
        pipeline.section("data", {"type": "object", "required": ["source"], "properties": {"source": {"type": "string", "format": "path"}}})
        pipeline.section("flaky", {"type": "object", "properties": {"label": {"type": "string", "default": "first"}}})

        @pipeline.stage(config_sections=("data",), input_files=lambda config: [config["data"]["source"]])
        def load(context):
            (context.output / "copy.txt").write_text(Path(context.config["data"]["source"]).read_text())

        @pipeline.stage(depends_on=("load",), config_sections=("flaky",))
        def flaky(context):
            (context.output / "partial.txt").write_text(context.config["flaky"]["label"])
            if CONTROL["failures"] > 0:
                CONTROL["failures"] -= 1
                raise RuntimeError("transient failure")
            (context.output / "done.txt").write_text(context.upstream["load"].path("copy.txt").read_text())

        @pipeline.stage(depends_on=("flaky",))
        def final(context):
            (context.output / "final.txt").write_text(context.upstream["flaky"].relative)

        (self.root / "source.txt").write_text("payload\n")
        return self.write_config(pipeline, '[data]\nsource = "source.txt"\n')

    def test_first_run_executes_and_identical_rerun_skips_without_touching_outputs(self):
        config, run = self.demo(), self.root / "run"
        first = self.run_study(config, run)
        self.assertEqual(self.actions(first), dict(ingest="executed", scale="executed", describe="executed", report="executed"))
        self.assertEqual({(row["status"], row["current_attempt"]) for row in first["stages"]}, {("completed", 1)})
        before = self.snapshot(run)
        second = self.run_study(config, run)
        self.assertEqual(set(self.actions(second).values()), {"skipped"})
        self.assertEqual(self.snapshot(run), before)
        self.assertEqual(sorted(path.name for path in run.iterdir()),
                         ["01_ingest", "02_scale", "03_describe", "04_report", "run.json", "run.lock", "run.log"])
        document = json.loads((run / "run.json").read_text())
        self.assertEqual([invocation["status"] for invocation in document["invocations"]], ["completed", "completed"])
        self.assertEqual(document["invocations"][1]["skipped"], ["ingest", "scale", "describe", "report"])
        self.assertIn("Built from `02_scale/attempt_1`, `03_describe/attempt_1`.", (run / "04_report/attempt_1/output/report.md").read_text())

    def test_stage_and_run_records(self):
        import platform
        import database_residue_screen
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        document = json.loads((run / "run.json").read_text())
        ingest_record = run / "01_ingest/attempt_1/stage.json"
        ingest, scale = json.loads(ingest_record.read_text()), json.loads((run / "02_scale/attempt_1/stage.json").read_text())
        self.assertEqual({key: scale[key] for key in ("schema_version", "study_id", "pipeline", "stage", "attempt", "directory", "mode",
                                                       "status", "error", "summary")},
                         dict(schema_version="2.0.0", study_id="demo", pipeline="demo", stage="scale", attempt=1,
                              directory="02_scale/attempt_1", mode="working", status="completed", error=None, summary=dict(factor=2, rows=8)))
        package = Path(database_residue_screen.__file__).parent
        closure = {"database_residue_screen/study/demo_stages.py": sha256(package / "study/demo_stages.py")}
        identity = document["code_identity"]
        environment = {key: identity[key] for key in ("python", "platform", "dependencies")}
        code = dict(closure_sha256=hashlib.sha256("".join(f"{digest}  {path}\n" for path, digest in sorted(closure.items())).encode()).hexdigest(),
                    tools={}, **environment)
        request, numbers = scale["request"], self.root / "input/numbers.txt"
        self.assertEqual(request, dict(stage="scale", config={"scaling": {"factor": 2}}, inputs={}, code=code,
                                       upstream={"ingest": {"attempt": "01_ingest/attempt_1", "stage_json_sha256": sha256(ingest_record)}}))
        self.assertEqual(scale["code"], dict(module="database_residue_screen.study.demo_stages", closure=closure))
        self.assertEqual(ingest["request"]["config"], {"numbers": {"source": str(numbers)}})
        self.assertEqual(ingest["request"]["inputs"], {str(numbers): dict(size=numbers.stat().st_size, sha256=sha256(numbers))})
        self.assertEqual(scale["request_digest"], json_sha256(request))
        self.assertEqual(environment["python"], dict(implementation=platform.python_implementation(), version=platform.python_version()))
        self.assertEqual(set(environment["platform"]), {"system", "machine", "libc"})
        self.assertIn("jsonschema", environment["dependencies"])
        self.assertEqual(scale["code_identity"], {key: identity[key] for key in ("mode", "package", "git", "pipeline")})
        self.assertEqual((identity["mode"], identity["package"]["name"], identity["package"]["version"]),
                         ("working", "database_residue_screen", database_residue_screen.__version__))
        self.assertEqual(identity["git"]["available"], document["host"]["git_root"] is not None)
        self.assertEqual((identity["pipeline"]["name"], identity["pipeline"]["in_package"]), ("demo", True))
        self.assertEqual(document["code_files"]["modules"]["study/engine.py"], sha256(package / "study/engine.py"))
        self.assertEqual(identity["package"]["module_count"], len(document["code_files"]["modules"]))
        self.assertIn("schemas/study.json", document["code_files"]["data"])
        output, log = run / "02_scale/attempt_1/output/scaled.tsv", run / "02_scale/attempt_1/stage.log"
        started_record = run / "02_scale/attempt_1/request.json"
        self.assertEqual(scale["outputs"], [dict(path="scaled.tsv", size=output.stat().st_size, sha256=sha256(output))])
        self.assertEqual(scale["log"], dict(path="stage.log", size=log.stat().st_size, sha256=sha256(log)))
        self.assertEqual(scale["request_json"], dict(path="request.json", size=started_record.stat().st_size, sha256=sha256(started_record)))
        self.assertTrue(scale["started_utc"].endswith("Z") and scale["ended_utc"].endswith("Z"))
        self.assertLessEqual(scale["started_utc"], scale["ended_utc"])
        self.assertGreaterEqual(scale["duration_seconds"], 0)
        started = json.loads(started_record.read_text())
        self.assertEqual((started["request"], started["request_digest"], started["code"]), (request, scale["request_digest"], scale["code"]))
        self.assertEqual({key: document[key] for key in ("study_id", "title", "mode", "config_file", "config_file_sha256")},
                         dict(study_id="demo", title="Demonstration study", mode="working", config_file=str(config),
                              config_file_sha256=sha256(config)))
        self.assertEqual(document["config"]["numbers"]["source"], str(numbers))
        self.assertEqual(document["config_sha256"], json_sha256(document["config"]))
        self.assertEqual([(stage["name"], stage["directory"], stage["module"]) for stage in document["pipeline"]["stages"]],
                         [(name, f"0{index}_{name}", "database_residue_screen.study.demo_stages")
                          for index, name in enumerate(("ingest", "scale", "describe", "report"), 1)])
        self.assertEqual(document["stages"]["scale"], dict(
            directory="02_scale", status="completed", current_attempt=1, latest_attempt=1, attempts=[dict(
                attempt=1, status="completed", mode="working", request_digest=scale["request_digest"],
                code_closure_sha256=code["closure_sha256"], stage_json_sha256=sha256(run / "02_scale/attempt_1/stage.json"),
                started_utc=scale["started_utc"], ended_utc=scale["ended_utc"], duration_seconds=scale["duration_seconds"])]))
        invocation = document["invocations"][0]
        self.assertEqual((invocation["selection"], invocation["git_commit"], invocation["package_modules_sha256"]),
                         (dict(until=None, only=[], rerun_changed=False), identity["git"].get("commit"), identity["package"]["modules_sha256"]))
        for line in (run / "run.log").read_text().splitlines() + (run / "01_ingest/attempt_1/stage.log").read_text().splitlines():
            self.assertRegex(line, UTC_LINE)
        self.assertIn("study.stage.ingest: ingested 8 integers", (run / "01_ingest/attempt_1/stage.log").read_text())

    def test_until_and_only(self):
        config, run = self.demo(), self.root / "run"
        result = self.run_study(config, run, until="scale")
        self.assertEqual(self.actions(result), dict(ingest="executed", scale="executed", describe="not_selected", report="not_selected"))
        self.assertFalse((run / "03_describe").exists())
        with self.assertRaisesRegex(StudyError, "--only report needs completed dependencies; not complete: describe"):
            self.run_study(config, run, only=["report"])
        self.assertFalse((run / "04_report").exists())
        result = self.run_study(config, run, only=["describe", "report"])
        self.assertEqual(self.actions(result), dict(ingest="verified", scale="verified", describe="executed", report="executed"))
        self.assertEqual(self.actions(self.run_study(config, run, only="report"))["report"], "skipped")

    def test_failed_attempt_is_kept_and_rerun_creates_a_new_attempt(self):
        config, run = self.flaky_pipeline(), self.root / "run"
        CONTROL["failures"] = 1
        with self.assertRaisesRegex(StudyError, "Stage 'flaky' failed in attempt 1: transient failure"):
            self.run_study(config, run)
        attempt = run / "02_flaky/attempt_1"
        failed = json.loads((attempt / "stage.json").read_text())
        self.assertEqual((failed["status"], failed["error"]["type"], failed["error"]["message"]), ("failed", "RuntimeError", "transient failure"))
        self.assertIn("Traceback", failed["error"]["traceback"])
        self.assertEqual([entry["path"] for entry in failed["outputs"]], ["partial.txt"])
        self.assertIn("transient failure", (attempt / "stage.log").read_text())
        self.assertEqual(stat.S_IMODE(attempt.stat().st_mode) & 0o222, 0)
        self.assertFalse((run / "03_final").exists())
        document = json.loads((run / "run.json").read_text())
        self.assertEqual((document["stages"]["flaky"]["status"], document["stages"]["flaky"]["current_attempt"]), ("failed", None))
        self.assertEqual((document["invocations"][-1]["status"], document["invocations"][-1]["failed"]), ("failed", "flaky"))
        before = self.snapshot(run)
        result = self.run_study(config, run)
        self.assertEqual(self.actions(result), dict(load="skipped", flaky="executed", final="executed"))
        self.assertEqual((run / "03_final/attempt_1/output/final.txt").read_text(), "02_flaky/attempt_2")
        self.assertEqual(json.loads((run / "03_final/attempt_1/stage.json").read_text())["request"]["upstream"]["flaky"]["attempt"],
                         "02_flaky/attempt_2")
        self.assertEqual({key: value for key, value in self.snapshot(run).items() if key in before}, before)
        flaky = json.loads((run / "run.json").read_text())["stages"]["flaky"]
        self.assertEqual((flaky["status"], flaky["current_attempt"], [entry["status"] for entry in flaky["attempts"]]),
                         ("completed", 2, ["failed", "completed"]))

    def test_changed_request_after_failure(self):
        config, run = self.flaky_pipeline(), self.root / "run"
        CONTROL["failures"] = 1
        with self.assertRaises(StudyError):
            self.run_study(config, run)
        config.write_text(config.read_text() + '[flaky]\nlabel = "second"\n')
        with self.assertRaisesRegex(StudyError, r"the request differs from its latest attempt 02_flaky/attempt_1 \(failed\); "
                                                r"changed fields: config\.flaky\.label\. Use --rerun-changed"):
            self.run_study(config, run)
        self.assertEqual([path.name for path in (run / "02_flaky").iterdir()], ["attempt_1"])
        CONTROL["failures"] = 1
        with self.assertRaisesRegex(StudyError, "failed in attempt 2"):
            self.run_study(config, run, rerun_changed=True)
        result = self.run_study(config, run)
        self.assertEqual(self.actions(result), dict(load="skipped", flaky="executed", final="executed"))
        self.assertEqual((run / "02_flaky/attempt_3/output/partial.txt").read_text(), "second")

    def test_incomplete_attempts_get_a_new_attempt(self):
        from database_residue_screen.study import runner
        config, run = self.flaky_pipeline(), self.root / "run"
        original = runner.write_json_atomic

        def killed(path, value):
            if path.name == "stage.json" and value["stage"] == "flaky":
                raise SystemExit("killed")
            original(path, value)

        with mock.patch.object(runner, "write_json_atomic", killed), self.assertRaises(SystemExit):
            self.run_study(config, run)
        self.assertTrue((run / "02_flaky/attempt_1/request.json").is_file())
        self.assertFalse((run / "02_flaky/attempt_1/stage.json").exists())
        self.assertEqual(study_status(run)["stages"][1]["status"], "incomplete")
        self.assertEqual(json.loads((run / "run.json").read_text())["invocations"][-1]["status"], "failed")
        (run / "03_final/attempt_1/output").mkdir(parents=True)
        result = self.run_study(config, run)
        self.assertEqual(self.actions(result), dict(load="skipped", flaky="executed", final="executed"))
        self.assertEqual(sorted(path.name for path in (run / "02_flaky").iterdir()), ["attempt_1", "attempt_2"])
        self.assertEqual(sorted(path.name for path in (run / "03_final").iterdir()), ["attempt_1", "attempt_2"])
        self.assertEqual(json.loads((run / "03_final/attempt_2/stage.json").read_text())["request"]["upstream"]["flaky"]["attempt"],
                         "02_flaky/attempt_2")

    def test_changed_config_is_refused_by_default_naming_the_fields(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        record, before = (run / "run.json").read_bytes(), self.snapshot(run)
        config.write_text(config.read_text().replace("factor = 2", "factor = 3"))
        for options, subject in (({}, "Stage 'scale'"), ({"until": "ingest"}, "Stage 'scale'"), ({"only": ["report"]}, "Stage 'scale' (required by --only)")):
            with self.subTest(options=options), self.assertRaises(StudyError) as caught:
                self.run_study(config, run, **options)
            self.assertIn(f"{subject}: the request differs from its current attempt 02_scale/attempt_1 (completed); changed fields: "
                          "config.scaling.factor.", str(caught.exception))
        self.assertEqual((run / "run.json").read_bytes(), record)
        self.assertEqual(self.snapshot(run), before)
        self.assertIn("ERROR study.engine: invocation stopped: Stage 'scale'", (run / "run.log").read_text())

    def test_rerun_changed_runs_changed_stages_and_their_descendants(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        before = self.snapshot(run)
        config.write_text(config.read_text().replace("factor = 2", "factor = 3"))
        result = self.run_study(config, run, rerun_changed=True)
        self.assertEqual(self.actions(result), dict(ingest="skipped", scale="executed", describe="skipped", report="executed"))
        self.assertEqual({row["stage"]: row["current_attempt"] for row in result["stages"]}, dict(ingest=1, scale=2, describe=1, report=2))
        self.assertEqual({key: value for key, value in self.snapshot(run).items() if key in before}, before)
        report = json.loads((run / "04_report/attempt_2/stage.json").read_text())["request"]["upstream"]
        self.assertEqual({name: value["attempt"] for name, value in report.items()}, dict(scale="02_scale/attempt_2", describe="03_describe/attempt_1"))
        self.assertIn("scaled total 93", (run / "04_report/attempt_2/output/report.md").read_text())
        scale = json.loads((run / "run.json").read_text())["stages"]["scale"]
        self.assertEqual((scale["current_attempt"], [entry["attempt"] for entry in scale["attempts"]]), (2, [1, 2]))
        self.assertNotEqual(scale["attempts"][0]["request_digest"], scale["attempts"][1]["request_digest"])
        self.assertEqual(set(self.actions(self.run_study(config, run)).values()), {"skipped"})
        rows = {row["stage"]: row for row in study_status(run)["stages"]}
        self.assertEqual((rows["scale"]["status"], rows["scale"]["current_attempt"], len(rows["scale"]["attempts"])), ("completed", 2, 2))

    def test_unselected_descendants_stay_stale_until_rerun(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        config.write_text(config.read_text().replace("factor = 2", "factor = 3"))
        result = self.run_study(config, run, until="scale", rerun_changed=True)
        self.assertEqual(self.actions(result), dict(ingest="skipped", scale="executed", describe="not_selected", report="not_selected"))
        with self.assertRaisesRegex(StudyError, r"Stage 'report': the request differs from its current attempt 04_report/attempt_1 "
                                                r"\(completed\); changed fields: upstream\.scale\.attempt, upstream\.scale\.stage_json_sha256\."):
            self.run_study(config, run)
        code = study_status(run, verify_code=True)["code"]
        self.assertEqual([(row["stage"], row["upstream"], row["ok"]) for row in code["stages"]],
                         [("ingest", "-", True), ("scale", "ok", True), ("describe", "ok", True), ("report", "STALE", False)])
        self.assertIn("not computed from the current attempt of: scale", code["stages"][3]["details"])
        self.assertEqual(self.actions(self.run_study(config, run, rerun_changed=True))["report"], "executed")
        self.assertTrue(study_status(run, verify_code=True)["code"]["ok"])

    def test_config_of_stages_without_attempts_may_still_change(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run, until="scale")
        config.write_text(config.read_text().replace('heading = "Demonstration report"', 'heading = "Revised heading"'))
        result = self.run_study(config, run)
        self.assertEqual(self.actions(result), dict(ingest="skipped", scale="skipped", describe="executed", report="executed"))
        self.assertTrue((run / "04_report/attempt_1/output/report.md").read_text().startswith("# Revised heading"))
        document = json.loads((run / "run.json").read_text())
        self.assertEqual(document["config"]["report"]["heading"], "Revised heading")
        self.assertNotEqual(document["invocations"][0]["config_sha256"], document["invocations"][1]["config_sha256"])

    def test_changed_input_file(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        numbers = self.root / "input/numbers.txt"
        numbers.write_text(numbers.read_text() + "7\n")
        with self.assertRaises(StudyError) as caught:
            self.run_study(config, run)
        self.assertIn("Stage 'ingest': the request differs from its current attempt 01_ingest/attempt_1 (completed)", str(caught.exception))
        self.assertIn(f'inputs["{numbers}"].sha256', str(caught.exception))
        result = self.run_study(config, run, rerun_changed=True)
        self.assertEqual(set(self.actions(result).values()), {"executed"})
        numbers.unlink()
        with self.assertRaisesRegex(StudyError, "Stage 'ingest': Input does not exist"):
            self.run_study(config, run, rerun_changed=True)

    def test_modified_upstream_records_are_refused_even_with_rerun_changed(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run)
        record = run / "01_ingest/attempt_1/stage.json"
        make_writable(run / "01_ingest")
        value = json.loads(record.read_text())
        value["summary"]["note"] = "edited"
        record.write_text(json.dumps(value, indent=2))
        for options in ({}, {"rerun_changed": True}):
            with self.subTest(options=options), self.assertRaisesRegex(
                    StudyError, r"Stage 'scale': attempt_1 used upstream 01_ingest/attempt_1, whose stage.json has changed; "
                                r"changed fields: upstream\.ingest\.stage_json_sha256\. The run directory was modified"):
                self.run_study(config, run, **options)
        shutil.rmtree(run / "01_ingest")
        with self.assertRaisesRegex(StudyError, "attempt_1 used upstream 01_ingest/attempt_1, whose stage.json is missing"):
            self.run_study(config, run, rerun_changed=True)
        self.assertFalse((run / "01_ingest").exists())

    def test_outputs_are_read_only_and_modifications_are_detected(self):
        config, run = self.demo(), self.root / "run"
        self.run_study(config, run, until="ingest")
        attempt = run / "01_ingest/attempt_1"
        for path in [attempt, *attempt.rglob("*")]:
            self.assertEqual(stat.S_IMODE(path.lstat().st_mode) & 0o222, 0, path)
        numbers = attempt / "output/numbers.json"
        original = numbers.read_bytes()
        make_writable(attempt)
        numbers.write_bytes(original.replace(b"3", b"7", 1))
        self.assertEqual(self.actions(self.run_study(config, run, until="ingest"))["ingest"], "skipped")
        with self.assertRaisesRegex(StudyError, r"Completed attempt 01_ingest/attempt_1 was modified after completion "
                                                r"\(content of numbers.json changed\)"):
            self.run_study(config, run)
        report = study_status(run, verify=True)
        self.assertEqual((report["stages"][0]["verified"], report["stages"][0]["problems"]), (False, ["content of numbers.json changed"]))
        numbers.write_bytes(original + b"\n")
        with self.assertRaisesRegex(StudyError, "size of numbers.json changed"):
            self.run_study(config, run, until="ingest")
        numbers.write_bytes(original)
        (attempt / "output/extra.txt").write_text("x")
        with self.assertRaisesRegex(StudyError, "unexpected output extra.txt"):
            self.run_study(config, run, until="ingest")
        (attempt / "output/extra.txt").unlink()
        self.assertEqual(self.actions(self.run_study(config, run))["scale"], "executed")
        self.assertTrue(study_status(run, verify=True)["stages"][0]["verified"])


class RunDirectoryTests(StudyCase):
    def test_misbehaving_stages_fail_cleanly(self):
        pipeline = new_pipeline("misbehaving")
        pipeline.section("data", {"type": "object", "required": ["source"], "properties": {"source": {"type": "string", "format": "path"}}})

        @pipeline.stage()
        def link(context):
            (context.output / "target.txt").write_text("x")
            os.symlink("target.txt", context.output / "alias.txt")

        @pipeline.stage()
        def listing(context):
            return ["not", "a", "mapping"]

        @pipeline.stage(config_sections=("data",), input_files=lambda config: [config["data"]["source"]])
        def mutate(context):
            with open(context.config["data"]["source"], "a", encoding="utf-8") as handle:
                handle.write("more\n")

        @pipeline.stage()
        def undeclared(context):
            context.tool("blastn")

        (self.root / "source.txt").write_text("data\n")
        config = self.write_config(pipeline, '[data]\nsource = "source.txt"\n')
        for stage, expected in (("link", "alias.txt is a symbolic link"), ("listing", "returns None or a JSON-serialisable mapping"),
                                ("mutate", "input changed while the stage was running"), ("undeclared", "did not declare tool 'blastn'")):
            run = self.root / f"run-{stage}"
            with self.subTest(stage=stage), self.assertRaisesRegex(StudyError, expected):
                self.run_study(config, run, only=[stage])
            self.assertEqual(json.loads(next(run.glob(f"*_{stage}/attempt_1/stage.json")).read_text())["status"], "failed")
        link_record = json.loads(next((self.root / "run-link").glob("*_link/attempt_1/stage.json")).read_text())
        self.assertIn(dict(path="alias.txt", size=len("target.txt"), kind="symbolic link"), link_record["outputs"])

    def test_secrets_are_never_written(self):
        pipeline = new_pipeline("secrets", secrets=("REMOTE_CODE",))
        pipeline.section("remote", {"type": "object", "required": ["url", "key", "code", "who"],
                                    "properties": {name: {"type": "string"} for name in ("url", "key", "code", "who")}})

        @pipeline.stage(config_sections=("remote",))
        def fetch(context):
            remote = context.config["remote"]
            SEEN.update(remote)
            context.log.info("using key %s and code %s for %s", remote["key"], remote["code"], remote["who"])
            if CONTROL["failures"] > 0:
                CONTROL["failures"] -= 1
                raise RuntimeError(f"request to {remote['url']} failed")
            (context.output / "fetched.txt").write_text("ok\n")
            return dict(echo=remote["code"])

        config = self.write_config(pipeline, """
            [remote]
            url = "https://example.org/x?api_key=${NCBI_API_KEY}"
            key = "${NCBI_API_KEY}"
            code = "${REMOTE_CODE}"
            who = "${WHO}"
        """)
        keys = ["key-0123456789abcdef", "key-fedcba9876543210", "key-third-value-000"]
        environ = dict(self.environ, NCBI_API_KEY=keys[0], REMOTE_CODE="code-secret-77", WHO="alice")
        run = self.root / "run"
        CONTROL["failures"] = 1
        with self.assertRaises(StudyError) as caught:
            self.run_study(config, run, environ=environ)
        self.assertIn(f"api_key={REDACTED}", str(caught.exception))
        self.assertNotIn(keys[0], str(caught.exception))
        self.assertEqual((SEEN["key"], SEEN["code"]), (keys[0], "code-secret-77"))
        self.assertEqual(self.actions(self.run_study(config, run, environ=dict(environ, NCBI_API_KEY=keys[1])))["fetch"], "executed")
        self.assertEqual(self.actions(self.run_study(config, run, environ=dict(environ, NCBI_API_KEY=keys[2])))["fetch"], "skipped")
        written = b"".join(path.read_bytes() for path in run.rglob("*") if path.is_file())
        for secret in (*keys, "code-secret-77"):
            self.assertNotIn(secret.encode(), written)
        self.assertIn(f"using key {REDACTED} and code {REDACTED} for alice", (run / "01_fetch/attempt_1/stage.log").read_text())
        document = json.loads((run / "run.json").read_text())
        self.assertEqual(document["environment"], dict(NCBI_API_KEY=REDACTED, REMOTE_CODE=REDACTED, WHO="alice"))
        self.assertEqual(document["config"]["remote"], dict(url=f"https://example.org/x?api_key={REDACTED}", key=REDACTED,
                                                            code=REDACTED, who="alice"))
        self.assertEqual(json.loads((run / "01_fetch/attempt_2/stage.json").read_text())["summary"], dict(echo=REDACTED))

    def test_tools_are_bound_only_into_the_stages_that_declare_them(self):
        tool = self.root / "bin/fake-tool"
        tool.parent.mkdir()
        tool.write_text('#!/bin/sh\nif [ "$1" = "--version" ]; then echo "fake-tool 1.2.3"; exit 0; fi\necho "processed $1"\n')
        tool.chmod(0o755)
        pipeline = new_pipeline("tools")

        @pipeline.stage(tools=("fake",))
        def use(context):
            completed = subprocess.run([context.tool("fake"), "item"], capture_output=True, text=True, check=True)
            (context.output / "tool.txt").write_text(completed.stdout)

        @pipeline.stage()
        def plain(context):
            (context.output / "plain.txt").write_text("no tool")

        config, run = self.write_config(pipeline, '[tools]\nfake = "bin/fake-tool"\n'), self.root / "run"
        self.run_study(config, run)
        identity = dict(command=str(tool), path=str(tool), resolved_path=str(tool), sha256=sha256(tool), version_args=["--version"],
                        version="fake-tool 1.2.3")
        self.assertEqual(json.loads((run / "run.json").read_text())["code_identity"]["tools"]["fake"], identity)
        self.assertEqual(json.loads((run / "01_use/attempt_1/stage.json").read_text())["request"]["code"]["tools"], {"fake": identity})
        self.assertEqual(json.loads((run / "02_plain/attempt_1/stage.json").read_text())["request"]["code"]["tools"], {})
        self.assertEqual((run / "01_use/attempt_1/output/tool.txt").read_text(), "processed item\n")
        tool.write_text(tool.read_text().replace("1.2.3", "1.2.4"))
        with self.assertRaises(StudyError) as caught:
            self.run_study(config, run)
        self.assertIn("Stage 'use': the request differs", str(caught.exception))
        self.assertIn("changed fields: code.tools.fake.sha256, code.tools.fake.version.", str(caught.exception))
        code = study_status(run, verify_code=True, environ=self.environ)["code"]["stages"]
        self.assertEqual([(row["stage"], row["tools"], row["ok"]) for row in code], [("use", "CHANGED", False), ("plain", "-", True)])
        self.assertIn("changed tools: fake (sha256, version)", code[0]["details"])
        self.assertEqual(self.actions(self.run_study(config, run, rerun_changed=True)), dict(use="executed", plain="skipped"))
        with self.assertRaisesRegex(StudyError, r"uses tool\(s\) fake; configure them in \[tools\]"):
            self.run_study(self.write_config(pipeline, "", name="untooled.toml"), self.root / "run-untooled")
        with self.assertRaisesRegex(StudyError, "Configured tool 'fake' was not found"):
            self.run_study(self.write_config(pipeline, '[tools]\nfake = "bin/missing"\n', name="missing.toml"), self.root / "run-missing")

    def test_run_directory_rules(self):
        config, run = self.demo(), self.root / "run"
        foreign = self.root / "foreign"
        foreign.mkdir()
        (foreign / "notes.txt").write_text("x")
        with self.assertRaisesRegex(StudyError, "is not empty and has no run.json"):
            self.run_study(config, foreign)
        self.assertEqual([path.name for path in foreign.iterdir()], ["notes.txt"])
        self.run_study(config, run, until="ingest")
        with self.assertRaisesRegex(StudyError, "belongs to study 'demo', not 'other'"):
            self.run_study(self.demo("other-input", study_id="other"), run)
        pipeline = new_pipeline("layout")
        pipeline.register(Stage("single", lambda context: None))
        with self.assertRaisesRegex(StudyError, "belongs to pipeline 'demo', not 'test-layout-"):
            self.run_study(self.write_config(pipeline, "", study_id="demo"), run)
        (run / "stray").mkdir()
        with self.assertRaisesRegex(StudyError, "Unexpected directory .*stray in the run directory"):
            self.run_study(config, run)
        (run / "stray").rmdir()
        with open(run / "run.lock", encoding="utf-8") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.assertTrue(study_status(run)["running"])
            with self.assertRaisesRegex(StudyError, "in use by another run-study process"):
                self.run_study(config, run)
        self.assertFalse(study_status(run)["running"])
        self.assertEqual(self.actions(self.run_study(config, run))["ingest"], "skipped")

    def test_stage_directories_persist_when_the_pipeline_changes(self):
        name = f"test-evolving-{uuid.uuid4().hex[:8]}"

        def write(context):
            (context.output / "out.txt").write_text(context.name)

        def version(stages):
            pipeline = Pipeline(name)
            for stage, dependencies in stages:
                pipeline.register(Stage(stage, write, dependencies))
            register_pipeline(pipeline, replace=True)

        version([("alpha", ()), ("gamma", ("alpha",))])
        config, run = self.write_config(load_pipeline(name), ""), self.root / "run"
        self.run_study(config, run)
        version([("alpha", ()), ("beta", ("alpha",)), ("gamma", ("alpha",))])
        self.assertEqual(self.actions(self.run_study(config, run)), dict(alpha="skipped", beta="executed", gamma="skipped"))
        self.assertEqual(sorted(path.name for path in run.iterdir() if path.is_dir()), ["01_alpha", "02_gamma", "03_beta"])
        version([("alpha", ()), ("beta", ("alpha",))])
        self.assertEqual(self.actions(self.run_study(config, run)), dict(alpha="skipped", beta="skipped"))
        document = json.loads((run / "run.json").read_text())
        self.assertEqual((document["pipeline"]["orphan_directories"], list(document["stages"])), (["02_gamma"], ["alpha", "beta"]))
        self.assertIn("WARNING study.engine: directories of stages no longer in pipeline", (run / "run.log").read_text())
        self.assertTrue((run / "02_gamma/attempt_1/stage.json").is_file())


class ProductionModeTests(StudyCase):
    def setUp(self):
        super().setUp()
        self.package_name = f"studyfixture_{uuid.uuid4().hex[:10]}"
        self.repository = self.root / "repository"
        self.package = write_package(self.repository / "src", self.package_name, {
            "__init__.py": '__version__ = "2.0"\n',
            "stages.py": """\
                from pathlib import Path


                def copy(context):
                    (context.output / "copy.txt").write_text(Path(context.config["data"]["source"]).read_text())
                """,
            "pipeline.py": f"""\
                from database_residue_screen.study import Pipeline

                from . import stages

                PIPELINE = Pipeline("{self.package_name}")
                PIPELINE.section("data", {{"type": "object", "required": ["source"],
                                          "properties": {{"source": {{"type": "string", "format": "path"}}}}}})
                PIPELINE.stage(config_sections=("data",), input_files=lambda config: [config["data"]["source"]])(stages.copy)
                """})
        (self.repository / ".gitignore").write_text("__pycache__/\nruns/\n")
        git(self.repository, "init", "-q")
        git(self.repository, "add", "-A")
        git(self.repository, "commit", "-q", "-m", "fixture")
        sys.path.insert(0, str(self.repository / "src"))
        self.addCleanup(sys.path.remove, str(self.repository / "src"))
        self.addCleanup(forget_modules, self.package_name)
        (self.root / "source.txt").write_text("payload\n")
        self.config = self.root / "study.toml"
        self.config.write_text(f'[study]\nid = "prod"\ntitle = "Production"\npipeline = "{self.package_name}.pipeline:PIPELINE"\n'
                               '[data]\nsource = "source.txt"\n')

    def production(self, run=None, config=None, package=None, **options):
        return self.run_study(config or self.config, run or self.root / "run", mode="production", package_dir=package or self.package,
                              **options)

    def working(self, **options):
        return self.run_study(self.config, self.root / "run", mode="working", package_dir=self.package, **options)

    def test_clean_tree_runs_and_records_the_commit(self):
        self.assertEqual(self.actions(self.production()), {"copy": "executed"})
        record = json.loads((self.root / "run/01_copy/attempt_1/stage.json").read_text())
        self.assertEqual(record["mode"], "production")
        self.assertEqual(record["code_identity"]["git"], dict(available=True, commit=git(self.repository, "rev-parse", "HEAD"),
                                                             tree=git(self.repository, "rev-parse", "HEAD^{tree}"), clean=True))
        self.assertEqual(record["code_identity"]["package"]["version"], "2.0")
        self.assertEqual(record["code"], dict(module=f"{self.package_name}.stages",
                                              closure={f"{self.package_name}/stages.py": sha256(self.package / "stages.py")}))
        document = json.loads((self.root / "run/run.json").read_text())
        self.assertEqual((document["mode"], document["invocations"][0]["mode"], document["invocations"][0]["git_clean"]),
                         ("production", "production", True))
        self.assertEqual(document["code_files"]["modules"], {name: sha256(self.package / name) for name in ("__init__.py", "pipeline.py", "stages.py")})
        self.assertEqual(self.actions(self.production()), {"copy": "skipped"})

    def test_production_reuses_only_production_attempts(self):
        self.production()
        self.assertEqual(self.actions(self.working()), {"copy": "skipped"})
        source = self.package / "stages.py"
        source.write_text(source.read_text() + "\n# reviewed\n")
        self.assertEqual(self.actions(self.working(rerun_changed=True)), {"copy": "executed"})
        git(self.repository, "commit", "-q", "-am", "review")
        with self.assertRaisesRegex(StudyError, r"Stage 'copy': the request differs from its current attempt 01_copy/attempt_2 \(completed\); "
                                                r"changed fields: mode \(production mode does not reuse a working-mode attempt\)\."):
            self.production()
        self.assertEqual(self.actions(self.production(rerun_changed=True)), {"copy": "executed"})
        modes = [entry["mode"] for entry in json.loads((self.root / "run/run.json").read_text())["stages"]["copy"]["attempts"]]
        self.assertEqual(modes, ["production", "working", "production"])

    def test_dirty_tree_is_refused(self):
        source = self.package / "stages.py"
        source.write_text(source.read_text() + "\n# local edit\n")
        with self.assertRaises(StudyError) as caught:
            self.production()
        self.assertIn("Production mode refused", str(caught.exception))
        self.assertIn(f"is dirty (1 changed or untracked paths, e.g. src/{self.package_name}/stages.py)", str(caught.exception))
        self.assertFalse((self.root / "run").exists())
        git(self.repository, "checkout", "--", ".")
        (self.repository / "notes.txt").write_text("untracked\n")
        with self.assertRaisesRegex(StudyError, "is dirty .*notes.txt"):
            self.production()
        (self.repository / "notes.txt").unlink()
        self.assertEqual(self.actions(self.production()), {"copy": "executed"})

    def test_working_mode_records_a_dirty_tree(self):
        source = self.package / "stages.py"
        source.write_text(source.read_text() + "\n# local edit\n")
        self.working()
        record = json.loads((self.root / "run/01_copy/attempt_1/stage.json").read_text())
        self.assertEqual((record["mode"], record["code_identity"]["mode"], record["code_identity"]["git"]["clean"]), ("working", "working", False))
        self.assertEqual(record["host"]["git_dirty_paths"], [f"src/{self.package_name}/stages.py"])
        self.assertEqual(json.loads((self.root / "run/run.json").read_text())["mode"], "working")

    def test_package_outside_git_untracked_package_and_foreign_stage_code_are_refused(self):
        outside = self.root / "outside" / self.package_name
        shutil.copytree(self.package, outside, ignore=shutil.ignore_patterns("__pycache__"))
        with self.assertRaisesRegex(StudyError, "not in a usable git working tree"):
            self.production(package=outside)
        ignoring = self.root / "ignoring"
        shutil.copytree(self.package, ignoring / "src" / self.package_name, ignore=shutil.ignore_patterns("__pycache__"))
        (ignoring / ".gitignore").write_text("src/\n")
        git(ignoring, "init", "-q")
        git(ignoring, "add", "-A")
        git(ignoring, "commit", "-q", "-m", "only the ignore file")
        with self.assertRaisesRegex(StudyError, "package is not importable from the git tree: package files not tracked by git: __init__.py"):
            self.production(package=ignoring / "src" / self.package_name)
        with self.assertRaises(StudyError) as caught:
            self.production(config=self.demo())
        self.assertIn("pipeline 'demo' is defined in database_residue_screen.study.demo, outside the package", str(caught.exception))
        self.assertIn("stage functions are defined outside the package in database_residue_screen.study.demo_stages", str(caught.exception))
        self.assertEqual(self.actions(self.run_study(self.demo("again"), self.root / "run", package_dir=self.package)),
                         dict(ingest="executed", scale="executed", describe="executed", report="executed"))

    def test_run_directory_inside_the_repository_must_be_ignored(self):
        with self.assertRaisesRegex(StudyError, "inside the git working tree .* but not ignored by git"):
            self.production(run=self.repository / "output")
        self.assertFalse((self.repository / "output").exists())
        self.assertEqual(self.actions(self.production(run=self.repository / "runs/first")), {"copy": "executed"})
        self.assertEqual(git(self.repository, "status", "--porcelain"), "")
        self.assertEqual(self.actions(self.production(run=self.repository / "runs/first")), {"copy": "skipped"})

    def test_identity_collection_without_a_run(self):
        identity = collect_code_identity("production", package_dir=self.package, environ=self.environ)
        self.assertEqual(identity.recorded["git"]["commit"], git(self.repository, "rev-parse", "HEAD"))
        self.assertEqual(identity.sha256, json_sha256(identity.recorded))
        self.assertEqual((identity.git_root, identity.mode), (self.repository, "production"))
        with self.assertRaisesRegex(StudyError, "Unknown mode"):
            collect_code_identity("draft", package_dir=self.package, environ=self.environ)


class StudyCliTests(StudyCase):
    def invoke(self, *arguments):
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            status = main([str(argument) for argument in arguments])
        return status, output.getvalue(), errors.getvalue()

    def test_run_study_and_study_status(self):
        config, run = self.demo(), self.root / "run"
        status, output, errors = self.invoke("run-study", config, "--output", run, "--until", "scale", "--working")
        self.assertEqual(status, 0, errors)
        self.assertEqual([(row["stage"], row["action"]) for row in json.loads(output)["stages"]],
                         [("ingest", "executed"), ("scale", "executed"), ("describe", "not_selected"), ("report", "not_selected")])
        self.assertIn("stage scale: attempt 1 completed", errors)
        status, table, _ = self.invoke("study-status", run, "--verify-code")
        self.assertEqual(status, 1)
        self.assertRegex(table, r"(?m)^describe\s+-\s+-\s+NOT RUN\s+-\s+-\s+-\s+-\s+-$")
        status, output, errors = self.invoke("run-study", config, "--output", run, "--only", "describe", "report", "--working")
        self.assertEqual(status, 0, errors)
        self.assertEqual(json.loads(output)["mode"], "working")
        status, table, _ = self.invoke("study-status", run)
        self.assertEqual(status, 0)
        document = json.loads((run / "run.json").read_text())
        lines = table.splitlines()
        self.assertEqual(lines[0], "Study demo: Demonstration study")
        self.assertIn(f"pipeline demo | last invocation in working mode | config {document['config_sha256'][:12]}", lines[1])
        self.assertRegex(lines[2], r"^STAGE\s+STATUS\s+CURRENT\s+ATTEMPTS\s+DURATION\s+REQUEST$")
        for line, name in zip(lines[3:], ("ingest", "scale", "describe", "report"), strict=True):
            cells = line.split()
            self.assertEqual(cells[:4], [name, "completed", "1", "1"])
            self.assertRegex(cells[4], r"^\d+\.\ds$")
            self.assertEqual(cells[5], document["stages"][name]["attempts"][0]["request_digest"][:12])
        status, table, _ = self.invoke("study-status", run, "--verify", "--verify-code")
        self.assertEqual(status, 0, table)
        self.assertRegex(table.splitlines()[3], r"^ingest\s+completed\s+1\s+1\s+\S+\s+[0-9a-f]{12}\s+ok$")
        self.assertIn(": all identical", table)
        self.assertRegex(table, r"(?m)^report\s+1\s+working\s+ok\s+-\s+ok\s+ok\s+ok\s+ok$")
        config.write_text(config.read_text().replace("factor = 2", "factor = 4"))
        status, _, errors = self.invoke("run-study", config, "--output", run, "--working")
        self.assertEqual(status, 1)
        self.assertIn("changed fields: config.scaling.factor. Use --rerun-changed", errors)
        status, output, errors = self.invoke("run-study", config, "--output", run, "--working", "--rerun-changed")
        self.assertEqual(status, 0, errors)
        self.assertEqual({row["stage"]: row["current_attempt"] for row in json.loads(output)["stages"]}, dict(ingest=1, scale=2, describe=1, report=2))
        target = run / "03_describe/attempt_1/output/description.json"
        make_writable(target.parent)
        target.write_text(target.read_text().replace("31", "32"))
        status, table, _ = self.invoke("study-status", run, "--verify")
        self.assertEqual(status, 1)
        self.assertIn("MODIFIED", table)
        self.assertIn("describe: content of description.json changed", table)

    def test_status_of_failed_and_unstarted_stages_and_durations(self):
        from database_residue_screen.study.status import format_duration
        self.assertEqual([format_duration(value) for value in (None, 4.25, 125, 7300, 200000)], ["-", "4.2s", "2m05s", "2h01m", "2d07h"])
        config, run = self.demo(numbers=()), self.root / "run"
        status, _, errors = self.invoke("run-study", config, "--output", run, "--working")
        self.assertEqual(status, 1)
        self.assertIn("Error: Stage 'describe' failed in attempt 1: min() ", errors)
        self.assertNotIn("Traceback", errors)
        status, table, _ = self.invoke("study-status", run)
        rows = {line.split()[0]: line.split() for line in table.splitlines()[3:]}
        self.assertEqual(rows["describe"][1:4], ["failed", "-", "1"])
        self.assertEqual(rows["report"][1:], ["not_started", "-", "0", "-", "-"])

    def test_errors_are_reported_without_traceback(self):
        status, _, errors = self.invoke("run-study", self.root / "missing.toml", "--output", self.root / "run", "--working")
        self.assertEqual(status, 1)
        self.assertIn("Error: Cannot read study configuration", errors)
        self.assertNotIn("Traceback", errors)
        status, _, errors = self.invoke("study-status", self.root)
        self.assertEqual(status, 1)
        self.assertIn("is not a study run directory", errors)
        status, _, errors = self.invoke("run-study", self.demo(), "--output", self.root / "run", "--only", "report", "--working")
        self.assertEqual(status, 1)
        self.assertIn("not complete: ingest, scale, describe", errors)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(["run-study", "c.toml", "--output", "r", "--until", "scale", "--only", "report"])
        from database_residue_screen.cli import parser
        args = parser().parse_args(["run-study", "c.toml", "--output", "r", "--only", "a", "b", "--only", "c", "--rerun-changed"])
        self.assertEqual((args.only, args.until, args.working, args.rerun_changed), (["a", "b", "c"], None, False, True))
        args = parser().parse_args(["study-status", "r", "--verify-code"])
        self.assertEqual((args.verify, args.verify_code), (False, True))


if __name__ == "__main__":
    unittest.main()
