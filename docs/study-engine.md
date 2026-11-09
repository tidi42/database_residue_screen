# Study engine

One command runs a whole study — days of BLAST searches, metadata retrieval,
scans, benchmarks, reports and figures — from **one TOML configuration file**,
into a run directory whose attempts are immutable and carry full provenance.
Code: [src/database_residue_screen/study](../src/database_residue_screen/study).

```bash
residue-screen run-study study.toml --output runs/paper3                  # production (clean git tree required)
residue-screen run-study study.toml --output runs/dev --working           # working mode (dirty tree allowed, recorded)
residue-screen run-study study.toml --output runs/paper3 --until search
residue-screen run-study study.toml --output runs/paper3 --only figures report
residue-screen run-study study.toml --output runs/paper3 --rerun-changed  # accept changed requests: new attempts
residue-screen study-status runs/paper3 [--verify] [--verify-code]
```

`--until S` runs S and everything it depends on; `--only S ...` runs exactly
those stages, whose other dependencies must already be complete and unchanged.
The status table shows, per stage, the status of the latest attempt
(`completed`, `failed`, `incomplete`, `running`, `not_started`), the current
attempt, the number of attempts, and the duration and request digest prefix of
the latest attempt. `--verify` re-hashes the outputs of every current attempt;
`--verify-code` checks every current attempt against the code now installed
(see [Release check](#release-check)). Both exit with code 1 on any mismatch.

## Concepts

**Pipelines, stage modules and code closures.** A pipeline declares
configuration-section schemas and stages. A stage has a name, `depends_on`,
the `config_sections` it reads, optional `input_files(config)` (external files
or directories hashed into the request on every invocation — name evidence or
manifest files rather than multi-GB databases), optional `tools` and a run
function. The **code closure** of a stage is the module defining its run
function (`run.__module__`) plus every module of the package that it imports,
transitively: absolute and relative imports, including imports inside
functions, found by static analysis. Parent-package `__init__` modules count
only when imported explicitly, dynamic imports are not seen, and if a closure
module reads package resources (`importlib.resources`, `pkgutil`) all package
data files (the JSON schemas) join the closure. Keep stage functions in
modules that import only what they need, and register them from a separate
pipeline module: the registering module and the engine it imports are then
outside every closure, and a fix in a report module leaves the closure — and
the completed attempts — of a BLAST stage untouched.

```python
# database_residue_screen/paper3_stages/search.py  (stage code: imports only what the stage uses)
from ..search import run_search

def search(context):
    blastn = context.tool("blastn")                    # the exact binary that was hashed
    panel = context.upstream["panel"].output           # exact upstream attempt_K/output
    threads = context.config["search"]["threads"]      # only declared sections are visible
    context.log.info("searching with %d threads", threads)
    ...                                                # write only below context.output
    return {"queries": 123}                            # optional JSON summary for stage.json

# database_residue_screen/paper3.py  (registration only)
from .paper3_stages import search
from .study import Pipeline

PIPELINE = Pipeline("paper3", secrets=("NCBI_API_KEY",))
PIPELINE.section("search", {"type": "object", "required": ["queries"], "properties": {
    "queries": {"type": "string", "format": "path"},              # resolved against the config directory
    "threads": {"type": "integer", "minimum": 1, "default": 8}}})  # defaults are filled in and recorded
PIPELINE.stage(depends_on=("panel",), config_sections=("search",), tools=("blastn",),
               input_files=lambda config: [config["search"]["queries"]])(search.search)
```

The context also offers `attempt`, `mode`, `study_id`, `title` and, per
dependency, `directory`, `output`, `relative`, `record` and `path(name)`.
Stages run in topological order (ties: registration order). `[study] pipeline`
names the stage set: a registered name (`register_pipeline(p)` in-process, an
entry of `PIPELINE_REFERENCES` or `register_pipeline_reference("paper3",
"database_residue_screen.paper3:PIPELINE")`, imported on first use) or a direct
`"package.module:ATTRIBUTE"`. The built-in `demo` pipeline
([demo.py](../src/database_residue_screen/study/demo.py) registers the stages
of [demo_stages.py](../src/database_residue_screen/study/demo_stages.py);
`write_demo_study(directory)` writes inputs) needs no BLAST, network or R.

**Configuration.** `[study]` has `id`, `title`, `pipeline` and optional
`secrets` (environment variable names). `[tools]` maps tool names to a command
(`blastn = "blastn"`, looked up on `PATH`; values with `/` are paths) or to
`{ command = "...", version_args = ["..."] }` (default `-version` for NCBI
BLAST+ programs, `--version` otherwise). Every other table must be a section
read by at least one stage; unknown sections and keys are errors. Strings may
use `${VAR}` (unset or empty: error; `$$` is a literal `$`). Values from secret
variables — those listed in `[study] secrets` or `Pipeline(secrets=...)`, and
any variable with a `_`-separated word KEY, APIKEY, TOKEN, SECRET, PASSWORD,
PASSWD, PASSPHRASE or CREDENTIAL(S), such as `NCBI_API_KEY` — reach the stages
but are recorded, hashed and logged as `<redacted>` (also inside error messages
and summaries). TOML dates and non-finite numbers are rejected. The config hash
is the SHA-256 of the canonical JSON (sorted keys, no whitespace, UTF-8) of the
resolved, redacted configuration.

**Modes.** *Production* refuses to start unless the imported package lies in a
git working tree, every package file is tracked, the whole tree is clean (no
modified, staged or untracked non-ignored paths), and the pipeline and every
stage function are defined inside the package. *Working* mode always runs and
is recorded as `working` in every `run.json`, `request.json` and `stage.json`.
Working invocations may reuse production attempts; production invocations never
reuse working-mode attempts (they count as changed, field `mode`). A run
directory inside the repository must be git-ignored (e.g. `runs/`).

## Files

```
RUN_DIR/run.json                          run record, rewritten atomically after every stage
RUN_DIR/run.log                           run-level log, UTC timestamps
RUN_DIR/run.lock                          advisory lock of the running invocation (one at a time)
RUN_DIR/NN_stage/attempt_K/request.json   request written when the attempt starts
RUN_DIR/NN_stage/attempt_K/stage.log      attempt log (engine and stage messages), UTC
RUN_DIR/NN_stage/attempt_K/output/        the only directory the stage writes
RUN_DIR/NN_stage/attempt_K/stage.json     final record (completed or failed)
```

`NN` is the stage's topological position when the run directory is created; a
stage keeps its directory, stages added later get the next free number, and
directories of stages removed from the pipeline are kept untouched (listed as
`orphan_directories`). A finished attempt (completed or failed) is made
read-only and never modified; to delete a run directory use
`chmod -R u+w RUN_DIR` first. Outputs must be regular files.

## Requests and rerun semantics

Each stage's **request** binds: the configuration sections it reads (redacted
form), the `stage.json` SHA-256 of the exact upstream attempts, the SHA-256 and
size of its input files, and a code block — the SHA-256 of its code closure
(over the sorted `"<sha256>  <module path>\n"` lines), the identities of the
tools it declares, the Python implementation and version, the dependency
versions and the coarse platform (system, machine, libc). The git commit, tree
and cleanliness, the package-wide digests and the mode are recorded but not
bound. A stage's **current attempt** is its latest completed attempt; it is the
one downstream stages use.

| Stage situation | Default | `--rerun-changed` |
| --- | --- | --- |
| no attempts | run `attempt_1` | same |
| current attempt has an identical request | skip | skip |
| latest attempt failed or incomplete with an identical request | new attempt | new attempt |
| request differs from the current attempt (or, if none completed, from the latest) | stop before anything executes, naming the changed fields | new attempt; earlier attempts stay for audit; every selected descendant reruns because its upstream digest changes |

Without `--rerun-changed`, all checks — for every stage that has attempts,
selected or not, and for stages whose upstream will produce a new attempt — run
before the first stage executes, so a refusal cannot waste days of upstream
work. With `--rerun-changed`, only the selected stages and their dependencies
are examined; unselected descendants stay stale until a later invocation
reruns them (a default invocation then stops on them, and `--verify-code`
reports them as `STALE`). Always refused: a different study id or pipeline
name, an attempt record or referenced upstream attempt that changed or
vanished, and completed outputs that were modified (sizes are checked on every
invocation; contents are re-hashed before a dependent stage runs and by
`study-status --verify`). Inputs, the stage's closure files and its tools are
re-checked after the stage runs and the closure also before it starts; a change
fails the attempt or stops the invocation. Each invocation should run in a fresh
process (the CLI does), so that the code executed is the code hashed.

## Release check

Before tagging a release, from the checkout to be tagged:

```bash
residue-screen study-status runs/paper3 --verify --verify-code   # every stage: ok
residue-screen run-study study.toml --output runs/paper3           # every stage: skipped
```

`--verify-code` recomputes each stage's code closure from the current package
source and reports per stage whether the current attempt was produced by
identical stage code (`CODE`), whether the declared tools (path, SHA-256,
version output), the dependency versions, Python and platform still match, and
whether the attempt was computed from the current attempts of its dependencies
(`UPSTREAM`). Stages without a completed attempt are reported as `NOT RUN`.
The second command confirms that configuration and inputs are unchanged too.

## Provenance fields

`stage.json`: `study_id`, `pipeline`, `stage`, `attempt`, `directory`, `mode`,
`status`, `started_utc`, `ended_utc`, `duration_seconds`, `request` (`stage`,
`config`, `upstream` {`attempt`, `stage_json_sha256`}, `inputs` {path: `size`,
`sha256` [, `files` for directories]}, `code` {`closure_sha256`, `tools`,
`python`, `platform`, `dependencies`}), `request_digest` (SHA-256 of the
canonical request), `code` (`module` and the `closure` listing: module path →
SHA-256), `code_identity` (recorded: `mode`, `package` name, version and
package-wide `modules_sha256`/`data_sha256`, `git` commit, tree and clean flag,
`pipeline` module and source digest), `outputs` (every file: `path`, `size`,
`sha256`), `outputs_sha256`, `log` and `request_json` (`path`, `size`,
`sha256`), `summary`, `error` (`type`, `message`, `traceback`; secrets
redacted) and `host` (full platform string, kernel release, interpreter path
and build, git dirty paths).

`run.json`: `study_id`, `title`, `created_utc`, `updated_utc`, `mode` of the
latest invocation, `pipeline` (name, module, stages with directories, modules
and declarations, `orphan_directories`), `config_file`, `config_file_sha256`,
`config_sha256`, `config` (resolved, redacted), `environment` (variables used;
secrets `<redacted>`), `code_identity` (recorded identity plus Python,
platform, dependencies and all configured tools of the latest invocation),
`code_files` (SHA-256 of every package module and data file), `host`, `stages`
(per stage: status of the latest attempt, `current_attempt`, `latest_attempt`
and every attempt with status, mode, request and closure digests,
`stage_json_sha256`, times and duration) and `invocations` (selection
including `rerun_changed`, mode, config digest, git commit and cleanliness,
package digest, executed/skipped/failed stages, outcome). Refused invocations
change nothing and are only logged in `run.log`.
