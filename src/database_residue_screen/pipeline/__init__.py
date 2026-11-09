"""Study pipelines of the framework.

`screen` (screen.py) runs a complete reference-panel screen from one study configuration (``[study] pipeline =
"screen"``). It only declares configuration sections and registers stages; every stage function lives in its own
module under `stages`, which imports only what the stage uses, so that each stage's code closure stays small.
"""
