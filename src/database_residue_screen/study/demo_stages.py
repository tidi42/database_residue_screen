"""Stage functions of the demonstration pipeline.

Stage code lives apart from pipeline registration so that each stage's code closure holds only what the stage
imports (here nothing else from the package), not the engine imported by the registering module.
"""
import json
from pathlib import Path


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def ingest(context):
    """Parse the configured integers."""
    values = [int(token) for token in Path(context.config["numbers"]["source"]).read_text(encoding="utf-8").split()]
    (context.output / "numbers.json").write_text(json.dumps(values) + "\n", encoding="utf-8")
    context.log.info("ingested %d integers", len(values))
    return dict(values=len(values))


def scale(context):
    """Multiply every integer by the configured factor."""
    values = _read_json(context.upstream["ingest"].path("numbers.json"))
    factor = context.config["scaling"]["factor"]
    rows = "".join(f"{value}\t{value * factor}\n" for value in values)
    (context.output / "scaled.tsv").write_text("value\tscaled\n" + rows, encoding="utf-8")
    return dict(rows=len(values), factor=factor)


def describe(context):
    """Count, minimum, maximum and total of the integers."""
    values = _read_json(context.upstream["ingest"].path("numbers.json"))
    description = dict(count=len(values), minimum=min(values), maximum=max(values), total=sum(values))
    (context.output / "description.json").write_text(json.dumps(description, indent=2) + "\n", encoding="utf-8")
    return description


def report(context):
    """Markdown report naming the exact upstream attempts it was built from."""
    description = _read_json(context.upstream["describe"].path("description.json"))
    scaled = context.upstream["scale"].path("scaled.tsv").read_text(encoding="utf-8").splitlines()[1:]
    lines = [f"# {context.config['report']['heading']}", "", f"Study `{context.study_id}`: {context.title}.", "",
             f"{description['count']} integers from {description['minimum']} to {description['maximum']} "
             f"(total {description['total']}); scaled total {sum(int(row.split()[1]) for row in scaled)}.", "",
             "Built from " + ", ".join(f"`{attempt.relative}`" for attempt in context.upstream.values()) + "."]
    (context.output / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
