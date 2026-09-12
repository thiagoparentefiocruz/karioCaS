#' karioCaS: Kraken Confidence Scores for Reliable Domain-Specific
#' Microbiota Inference and Discovery
#'
#' karioCaS provides a framework to explore the effect of the Kraken2
#' confidence score (CS) on metagenomic classification. Reports produced at
#' several CS values are imported into a single
#' \code{TreeSummarizedExperiment}, and taxonomic stability is then evaluated
#' across stringency levels, separately for Bacteria, Archaea, Eukaryota and
#' Viruses, to identify an objective CS threshold per domain.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item \code{\link{import_karioCaS}()}: import the per-CS reports into a
#'     \code{TreeSummarizedExperiment}.
#'   \item \code{\link{taxa_retention}()}: taxa retention curves across CS and
#'     the Stability Index used to select the optimal CS per domain.
#'   \item \code{\link{reads_per_taxa}()}: saturation analysis and optimal
#'     minimum reads per domain.
#'   \item \code{\link{retrieve_selected_taxa}()}: build the final,
#'     domain-specific taxonomic mosaic.
#'   \item \code{\link{upset_kariocas}()}, \code{\link{group_upset}()},
#'     \code{\link{heatmaps_karioCaS}()}, \code{\link{taxa_resolution}()}:
#'     exploratory and comparative visualisations.
#' }
#'
#' @section Acknowledgement of AI assistance:
#' The scientific concept, the method design and its validation are the
#' author's. The package code was written with the assistance of large
#' language models, which also contributed to the formulation of some of the
#' calculations. All code was reviewed and tested by the author.
#'
#' @keywords internal
"_PACKAGE"
