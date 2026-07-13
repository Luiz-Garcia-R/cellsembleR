#' cellsembleR: Cell Signature Analysis for Bulk Transcriptomics
#'
#' cellsembleR provides a unified framework for estimating, analyzing,
#' and visualizing biological signatures from bulk transcriptomic data.
#'
#' Signature estimation can be performed using reference-based methods,
#' including EPIC, xCell, and MCP-counter, or through user-defined gene
#' signatures. Downstream functions support statistical comparisons,
#' association analyses, dimensionality reduction, and visualization.
#'
#' The package is designed to transform inferred cellular signatures into
#' interpretable biological features that can be integrated into
#' downstream analyses.
#'
#' ## Getting started
#'
#' The typical workflow consists of:
#'
#' 1. Estimate signatures with one of the `*_run()` functions.
#' 2. Compare experimental groups with `cellsemble_compare()`.
#' 3. Visualize results using the plotting functions.
#'
#' See the package README and function documentation for complete workflows
#' and examples.
#'
#' @docType package
#' @name cellsembleR
#'
"_PACKAGE"
