# ---------------------------------------------------------------------------
# Questionnaire data of the two departments (not distributed with the package)
# ---------------------------------------------------------------------------

# identifiers of all the members of the two departments in the original
# processing (department A: 1-53; department B: 54-119 except 115)
.department_members <- list(A = 1:53, B = setdiff(54:119, 115))
.fsna_data_cache <- new.env(parent = emptyenv())

#' Questionnaire data of the two departments
#'
#' The collaboration ratings of Porreca, Maturo and Ventre (2025) and Porreca
#' and Maturo (2026) are \strong{not distributed} with the package: the
#' structure of the networks could allow the identification of the respondents
#' even with numerical codes. They are available from the corresponding author
#' upon reasonable request, within the terms of the informed consent. Once the
#' raw export of the questionnaire (a semicolon-separated file) is available,
#' set its path with \code{options(FSNA.questionnaire = "path/to/file.csv")} or
#' the environment variable \code{FSNA_QUESTIONNAIRE}; then
#' \code{fsna_questionnaire()} reads it and \code{\link{fsn_department}} builds
#' the networks, and the examples, tests and vignettes that use the real data
#' are run.
#'
#' The ratings were collected with a web questionnaire among the professors of
#' two departments of the University "G. d'Annunzio" of Chieti-Pescara: the
#' Department of Economic Studies (\code{"A"}), whose network is analysed in
#' Porreca, Maturo and Ventre (2025), and the Department of Medical, Oral and
#' Biotechnological Sciences (\code{"B"}). Following the Dynamic Fuzzy Rating
#' Tracker of Calcagni and Lombardi (2014), each respondent rated the
#' collaboration with colleagues by moving the mouse on a semicircular
#' five-point scale (Never, Seldom, Occasionally, Often, Constantly, coded 1 to
#' 5), and the trajectory of the pointer was recorded until the final click.
#' Each rating is the TFN (minimum, final position, maximum) of the trajectory
#' (\code{\link{tfn_from_trajectory}}). Ratings whose trajectory could not be
#' reconstructed (52 of 275 in department A, 61 of 316 in B) are dropped.
#'
#' @param file Path of the raw export; by default \code{fsna_data_path()}.
#' @param refresh Logical: read the file again even if a copy with the same
#'   content is cached in the session. The cache is keyed by the path and the
#'   MD5 checksum of the file, so a modified file is always read again.
#' @return \code{fsna_questionnaire()}: a data frame with one row per rating and
#'   columns \code{department} (\code{"A"} or \code{"B"}), \code{order} (order
#'   of submission), \code{rater}, \code{rated} (numerical codes), \code{l},
#'   \code{m}, \code{r} (TFN divided by 5), \code{l_raw}, \code{m_raw},
#'   \code{r_raw} (TFN on the scale 1-5), \code{button} and \code{n_points}.
#'   \code{fsna_data_path()}: the path set by the user, or \code{""}.
#'   \code{fsna_data_available()}: \code{TRUE} if the file exists.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures
#'   in social network analysis: Theory and application in a university
#'   department collaboration network. \emph{International Journal of Approximate
#'   Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#'
#'   Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#'
#'   Calcagni A., Lombardi L. (2014). Dynamic Fuzzy Rating Tracker (DYFRAT): a
#'   novel methodology for modeling real-time dynamic cognitive processes in
#'   rating scales. \emph{Applied Soft Computing} 24, 948-961.
#' @examples
#' fsna_data_available()
#' if (fsna_data_available()) {
#'   d <- fsna_questionnaire()
#'   table(d$department)
#'   with(d, tapply(r - l, department, mean))   # mean width of the ratings
#' }
#' @section Methodology:
#' The data illustrate the methods; they are not a sample from a population of departments. About half of the members responded, only respondents send ties, and ratings were given to a subset of colleagues: the in-centralities describe the ratings given by the respondents, and treating unrated pairs as null ties is an assumption (see \code{\link{set_unknown}}). See Section 4 of Porreca and Maturo (2026).
#' @export
fsna_questionnaire <- function(file = fsna_data_path(), refresh = FALSE) {
  if (!nzchar(file) || !file.exists(file))
    stop("the questionnaire data are not distributed with FSNA; they are available from the ",
         "corresponding author upon reasonable request. Set their path with ",
         "options(FSNA.questionnaire = \"path/to/file.csv\"); see ?fsna_questionnaire")
  # the cache is keyed by the content of the file (MD5), so that a new version
  # of the file at the same path is read again
  key <- paste(normalizePath(file), unname(tools::md5sum(file)), sep = "|")
  if (refresh || is.null(.fsna_data_cache[[key]])) {
    raw <- read_fuzzy_questionnaire(file)
    # department codes of the export: 1 = Economic Studies (A), 0 = Medical,
    # Oral and Biotechnological Sciences (B)
    raw$department <- ifelse(raw$department == 1, "A", "B")
    d <- raw[, c("department", "order", "rater", "rated", "l", "m", "r",
                 "l_raw", "m_raw", "r_raw", "button", "n_points")]
    rownames(d) <- NULL
    .fsna_data_cache[[key]] <- d
  }
  .fsna_data_cache[[key]]
}

#' @rdname fsna_questionnaire
#' @export
fsna_data_path <- function() {
  p <- getOption("FSNA.questionnaire", default = "")
  if (is.null(p) || !nzchar(p)) p <- Sys.getenv("FSNA_QUESTIONNAIRE", unset = "")
  p
}

#' @rdname fsna_questionnaire
#' @export
fsna_data_available <- function() {
  p <- fsna_data_path(); nzchar(p) && file.exists(p)
}

#' The collaboration networks of the two departments
#'
#' Builds the fuzzy collaboration network of a department from the
#' questionnaire data, which are not distributed with the package (see
#' \code{\link{fsna_questionnaire}}).
#' @param department \code{"A"} (Economic Studies, the network of Porreca, Maturo
#'   and Ventre, 2025) or \code{"B"} (Medical, Oral and Biotechnological
#'   Sciences).
#' @param members \code{"rated"} (default): the actors with at least one tie of
#'   positive support (51 and 55 actors with \code{scale = "unit"}, 51 and 56
#'   with \code{scale = "ijar"}, since with the latter a rating Never is a tie); \code{"all"}: all the members of the department (53 and 65);
#'   the choice affects the indices that are averages over \eqn{n - 1} actors.
#' @param scale Mapping of the ratings from the scale \eqn{[1, 5]} to
#'   \eqn{[0, 1]}: \code{"unit"}, \eqn{x \mapsto (x - 1)/4}, so that the answer
#'   Never is the absence of collaboration, as in Porreca and Maturo (2026);
#'   \code{"ijar"}, \eqn{x \mapsto x/5}, as in Porreca, Maturo and Ventre (2025),
#'   which gives strength 0.2 to Never. By default \code{"ijar"} when
#'   \code{members = "all"} and
#'   \code{"unit"} otherwise.
#' @param keep Rule for pairs rated twice, see \code{\link{fsn_from_ratings}}.
#' @param unknown \code{"none"} (default: pairs without a rating are null ties,
#'   as in the original study), \code{"nonrespondents"} or \code{"unrated"}:
#'   see \code{\link{set_unknown}}.
#' @param file Path of the raw export of the questionnaire, see
#'   \code{\link{fsna_questionnaire}}.
#' @return An \code{\link{fsn}} object.
#' @examples
#' if (fsna_data_available()) {
#'   A <- fsn_department("A")
#'   A
#'   summary(A)
#'   fsn_department("B", members = "all")
#' }
#' @export
#' @section Methodology:
#' The default (\code{members = "rated"}) excludes members without ties, which do not carry information on the network; \code{members = "all"} keeps every member of the department. Ratings whose trajectory could not be reconstructed (52 of 275 in department A, 61 of 316 in B) are dropped and the pairs are null ties; later ratings of the same colleague replace earlier ones. With \code{scale = "unit"} and the default arguments the networks are those analysed in Section 4 of Porreca and Maturo (2026): 51 actors (24 respondents) and 55 actors (28 respondents).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' if (fsna_data_available()) {
#'   A <- fsn_department("A"); B <- fsn_department("B")
#'   print(rbind(A = c(n = n_actors(A), respondents = sum(A$respondents)),
#'               B = c(n = n_actors(B), respondents = sum(B$respondents))))
#'   fsn_department("A", keep = "mean")   # average repeated ratings
#' }
fsn_department <- function(department = c("A", "B"), members = c("rated", "all"), scale = NULL,
                           keep = "last", unknown = c("none", "nonrespondents", "unrated"),
                           file = fsna_data_path()) {
  unknown <- match.arg(unknown)
  department <- match.arg(department); members <- match.arg(members)
  if (is.null(scale)) scale <- if (members == "all") "ijar" else "unit"
  scale <- match.arg(scale, c("unit", "ijar"))
  d <- fsna_questionnaire(file)
  d <- d[d$department == department, ]
  if (scale == "unit") {
    d$l <- (d$l_raw - 1) / 4; d$m <- (d$m_raw - 1) / 4; d$r <- (d$r_raw - 1) / 4
  }
  act <- if (members == "all") .department_members[[department]] else NULL
  desc <- c(A = "Department of Economic Studies (A)", B = "Department of Medical, Oral and Biotechnological Sciences (B)")
  x <- fsn_from_ratings(d, from = "rater", to = "rated", order = "order", keep = keep, actors = act,
                        drop_isolated = members == "rated",
                        info = list(description = paste0(desc[[department]], if (scale == "unit") " [ratings mapped by (x - 1)/4]" else " [ratings divided by 5]"),
                                    department = department, scale = scale))
  if (unknown != "none") x <- set_unknown(x, unknown)
  x
}
