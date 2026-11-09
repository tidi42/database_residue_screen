"""Stage ``figures``: render the publication figures with Rscript (reporting.figures).

Upstream mapping: ``figure_sources`` (contract publication_figure_sources_v1). The R scripts and build script are the
directory named by ``[figures] scripts`` (default: ``figures`` beside the study configuration); it is hashed as a stage input.
"""
from ...contracts import ContractError
from ...reporting.figures import build_figures


def figures(context):
    """Publication figures (R)"""
    source = context.upstream["figure_sources"].output
    if not (source / "manifest.json").is_file():
        raise ContractError("Stage figure_sources: manifest.json missing")
    return build_figures(source, context.output, scripts=context.config["figures"]["scripts"], rscript=context.tool("Rscript"),
                         log=context.log.info)
