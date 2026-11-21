#' FSNA: Fuzzy Social Network Analysis with Triangular Fuzzy Ties
#'
#' The package implements two families of methods for directed social networks
#' whose ties are triangular fuzzy numbers.
#'
#' \strong{Indices of Porreca, Maturo and Ventre (2025).} Fuzzy degree as an
#' FTOWA of the ties (\code{\link{ijar_degree}}), best paths selected by a
#' ranking score (\code{\link{ijar_best_paths}}), fuzzy closeness
#' (\code{\link{ijar_closeness}}) and betweenness on best paths
#' (\code{\link{ijar_betweenness}}). \code{variant = "definition"} follows the
#' definitions of the paper; \code{variant = "alternative"} changes the path
#' strength, the treatment of equal scores and the normalization, and each
#' choice can also be set individually.
#' Any of the ranking scores of \code{\link{tfn_score}} can replace the centre
#' of gravity.
#'
#' \strong{Exact fuzzy centralities.} For every index that is continuous and
#' non-decreasing in the tie weights, \code{\link{fuzzy_centrality}} computes the
#' \eqn{\alpha}-cuts of the fuzzy centrality defined by the extension principle
#' from two crisp evaluations per level. Actors are compared through possibility
#' and necessity degrees (\code{\link{possibility_matrix}}), certified orderings
#' (\code{\link{certified_pairs}}), rank intervals (\code{\link{rank_intervals}})
#' and the fuzzy set of the most central actor (\code{\link{top_set}}), and
#' ranked by nine methods (\code{\link{rank_fuzzy}}).
#'
#' \strong{Data.} Readers of the questionnaire data of two departments, which
#' are available on request and not distributed with the package
#' (\code{\link{fsna_questionnaire}}, \code{\link{fsn_department}}), a
#' simulated survey network for the examples (\code{\link{fsn_simulated}}), importers
#' (\code{\link{read_fuzzy_questionnaire}}, \code{\link{fsn_from_ratings}},
#' \code{\link{read_ratings}}, \code{\link{fuzzify_responses}}), simulators
#' (\code{\link{simulate_fsn}}, \code{\link{simulate_ratings}}), a simulation
#' study (\code{\link{simulation_study}}) and a complete analysis with report
#' (\code{\link{fsna_analyze}}).
#'
#' @section Main entry points:
#' \tabular{lll}{
#'   \strong{Task} \tab \strong{Function} \tab \strong{Returns} \cr
#'   Build a network from ratings \tab \code{\link{fsn_from_ratings}} \tab \code{fsn} \cr
#'   Use the departmental data (on request) \tab \code{\link{fsn_department}} \tab \code{fsn} \cr
#'   Simulated survey network \tab \code{\link{fsn_simulated}} \tab \code{fsn} \cr
#'   Simulate a network \tab \code{\link{simulate_fsn}} \tab \code{fsn} with latent network \cr
#'   Indices of the 2025 paper \tab \code{\link{ijar_degree}}, \code{\link{ijar_closeness}}, \code{\link{ijar_betweenness}} \tab TFN per actor \cr
#'   Exact fuzzy centrality \tab \code{\link{fuzzy_centrality}} \tab \code{fuzzy_centrality} (alpha-cuts) \cr
#'   Compare actors \tab \code{\link{possibility_matrix}}, \code{\link{certified_pairs}}, \code{\link{top_set}} \tab degrees, relations \cr
#'   Rank actors \tab \code{\link{rank_fuzzy}}, \code{\link{compare_rankings}} \tab ranks \cr
#'   Draw \tab \code{\link{plot.fuzzy_centrality}}, \code{\link{plot_fsn_glyphs}}, \code{\link{plot_fuzzy_scatter}} \tab figures \cr
#'   Everything at once \tab \code{\link{fsna_analyze}} \tab figures, tables, HTML report \cr
#'   Collect new data \tab \code{\link{make_collector}}, \code{\link{read_ratings}}, \code{\link{fuzzify_responses}} \tab web form, data frame
#' }
#'
#' @section What is exact and what is approximate:
#' \tabular{lll}{
#'   \strong{Quantity} \tab \strong{Status} \tab \strong{Conditions} \cr
#'   \eqn{\alpha}-cuts of \code{\link{fuzzy_centrality}} \tab exact at the grid levels \tab isotone index (Theorem 1 of Porreca and Maturo, 2026) \cr
#'   Triangle \eqn{(c(L), c(M), c(R))} \tab exact at support and core; everywhere for strength \tab cores of the ties single values (no unknown ties) \cr
#'   Marginal \eqn{\Pi}, \eqn{N}, \eqn{\mu_{top}} \tab exact for strength (interpolation); numerical estimate with bracket otherwise (bisection) \tab see \code{possibility_matrix()} \cr
#'   Marginal vs joint degrees \tab equal for in/out-strength and OWA of in/out ties \tab directed networks; bounds otherwise \cr
#'   Certified ordering (\code{\link{certified_pairs}}) \tab exact, from the cuts \tab holds for every compatible network \cr
#'   Rank intervals \tab attained for in/out-strength; bounds otherwise \tab positions after tie-breaking \cr
#'   Monte Carlo cuts (\code{\link{mc_centrality}}) \tab inner approximation \tab contained in the exact cuts \cr
#'   IJAR indices \tab defined by ranking TFNs with a score \tab not the extension principle (except degree with equal weights)
#' }
#' A ranking by a score (e.g. the centre of gravity) is a total order; only the
#' certified pairs are orderings that hold for every network compatible with the
#' ratings.
#'
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed
#'   fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings
#'   of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#'
#'   Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social
#'   network analysis: Theory and application in a university department
#'   collaboration network. \emph{International Journal of Approximate Reasoning}
#'   176, 109319. \doi{10.1016/j.ijar.2024.109319}
#'
#' Start with \code{vignette("FSNA")}.
#' @keywords internal
"_PACKAGE"
