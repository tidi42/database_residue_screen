"""The reporting stages registered into a study pipeline and run through the study engine on synthetic upstream outputs."""
import hashlib
import json
import os
import shutil
import tempfile
import unittest
import uuid
from pathlib import Path

import reporting_fixture
from database_residue_screen.pipeline import reporting
from database_residue_screen.pipeline.stages.report_inputs import LAYOUT
from database_residue_screen.study import Pipeline, StudyError, register_pipeline, run_study

UPSTREAM = ("sources", "panel", "controls", "design", "database", "search", "analysis", "metadata", "events", "taxonomy", "insdc",
            "positional", "record_cases", "census", "exact_scan", "scan_validation", "cap_review", "target_limit", "benchmark", "inference",
            "validation")
SCRIPTS = Path(__file__).resolve().parents[1] / "studies" / "paper3" / "figures"
RSCRIPT = shutil.which("Rscript")
OUTPUTS: dict[str, Path] = {}


def upstream(context):
    """Synthetic upstream stage: places the fixture output of the stage of the same name in the stage's output layout."""
    name, source, output = context.name, OUTPUTS[context.name], context.output
    if name == "panel":
        definition = source / "panel_definition.tsv"
        output.joinpath("panel.json").write_text(json.dumps(dict(definition=str(definition), definition_sha256=hashlib.sha256(
            definition.read_bytes()).hexdigest())), encoding="utf-8")
    elif name == "benchmark":
        sizes = {}
        for run in sorted(source.glob("benchmark_*")):
            label = run.name.removeprefix("benchmark_")
            shutil.copytree(run, output / "runs" / label)
            sizes[label] = {}
        output.joinpath("benchmark.json").write_text(json.dumps(dict(sizes=sizes)), encoding="utf-8")
    elif name == "inference":
        shutil.copytree(source, output / "precision")
        (output / "specification").mkdir()
        shutil.copyfile(source / "specification.json", output / "specification" / "precision_specification.json")
    elif name in LAYOUT:
        shutil.copytree(source, output / LAYOUT[name][0], dirs_exist_ok=True)
    else:
        shutil.copytree(source, output, dirs_exist_ok=True)


def screen_pipeline() -> Pipeline:
    pipeline = Pipeline(f"screen-test-{uuid.uuid4().hex[:8]}")
    for name in UPSTREAM:
        pipeline.stage(name)(upstream)
    return register_pipeline(reporting.register(pipeline))


class RegistrationTests(unittest.TestCase):
    def test_sections_stages_and_dependencies(self):
        pipeline = screen_pipeline()
        stages = {stage.name: stage for stage in pipeline.ordered()}
        self.assertEqual([stage.name for stage in pipeline.ordered()][-4:], ["metrics", "reports", "figure_sources", "figures"])
        self.assertEqual(set(stages["metrics"].depends_on), {"design", "search", "taxonomy"})
        self.assertTrue(set(stages["reports"].depends_on) <= set(UPSTREAM) | {"metrics"})
        self.assertEqual(stages["figures"].tools, ("Rscript",))
        self.assertEqual(stages["figures"].config_sections, ("figures",))
        sections = pipeline.sections()
        self.assertEqual(sorted(sections), ["figures", "metrics", "reports"])
        self.assertEqual(sections["reports"]["properties"]["draw_prefixes"]["default"], [5, 10, 15])
        with self.assertRaises(StudyError):
            reporting.register(pipeline)


class EngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name).resolve()
        study = reporting_fixture.build(cls.root / "study")
        OUTPUTS.update({name: study.get(name, study["taxonomy"]) for name in UPSTREAM})
        for name in ("sources", "database"):
            directory = cls.root / "study" / name
            directory.mkdir()
            (directory / f"{name}.json").write_text("{}\n")
            OUTPUTS[name] = directory

    @classmethod
    def tearDownClass(cls):
        for directory, _, files in os.walk(cls.root):
            os.chmod(directory, 0o700)
            for name in files:
                os.chmod(os.path.join(directory, name), 0o600)
        cls.temporary.cleanup()

    def test_reporting_stages_run_through_the_engine(self):
        rscript = RSCRIPT
        until = None
        if rscript is None:
            # Without R the figures stage cannot run; a stand-in executable satisfies the tool declaration.
            rscript = str(self.root / "Rscript")
            Path(rscript).write_text("#!/bin/sh\necho 'R stand-in'\n")
            os.chmod(rscript, 0o755)
            until = "figure_sources"
        pipeline = screen_pipeline()
        config = self.root / "study.toml"
        config.write_text(f'[study]\nid = "reporting-test"\ntitle = "Reporting stages"\npipeline = "{pipeline.name}"\n\n'
                          f'[tools]\nRscript = "{rscript}"\n\n[figures]\nscripts = "{SCRIPTS}"\n', encoding="utf-8")
        environ = dict(os.environ, GIT_CEILING_DIRECTORIES=str(self.root))
        run = self.root / "run"
        result = run_study(config, run, until=until, mode="working", environ=environ)
        actions = {row["stage"]: row["action"] for row in result["stages"]}
        expected = ["metrics", "reports", "figure_sources"] + ([] if until else ["figures"])
        self.assertTrue(all(actions[name] == "executed" for name in expected), actions)
        record = json.loads(next(run.glob("*_metrics/attempt_1/stage.json")).read_text())
        self.assertEqual(record["status"], "completed")
        closure = set(record["code"]["closure"])
        self.assertIn("database_residue_screen/reporting/metrics.py", closure)
        self.assertIn("database_residue_screen/pipeline/stages/metrics.py", closure)
        for module in ("database_residue_screen/reporting/figure_sources.py", "database_residue_screen/reporting/reports.py",
                       "database_residue_screen/pipeline/reporting.py", "database_residue_screen/study/engine.py"):
            self.assertNotIn(module, closure)
        self.assertEqual(record["summary"]["counts"]["canonical_parents"], 16)
        reports = json.loads(next(run.glob("*_reports/attempt_1/output/reports.json")).read_text())
        self.assertEqual(reports["skipped"], [])
        sources = json.loads(next(run.glob("*_figure_sources/attempt_1/output/manifest.json")).read_text())
        self.assertEqual(sources["contract"], "publication_figure_sources_v1")
        if until is None:
            figures = json.loads(next(run.glob("*_figures/attempt_1/stage.json")).read_text())
            self.assertIn(str(SCRIPTS), figures["request"]["inputs"])
            self.assertEqual(figures["summary"]["figures"]["figureS4"], "Figure S2")
        again = run_study(config, run, until=until, mode="working", environ=environ)
        self.assertTrue(all(row["action"] == "skipped" for row in again["stages"]))


if __name__ == "__main__":
    unittest.main()
