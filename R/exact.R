# ---------------------------------------------------------------------------
# Exact fuzzy centralities by the extension principle
# ---------------------------------------------------------------------------

#' Exact fuzzy centrality by the extension principle
#'
#' For a centrality index \eqn{c} that is continuous and non-decreasing in every
#' tie weight on a box containing the support of the network, the fuzzy
#' centrality defined by the extension principle,
#' \eqn{\mu_{\tilde c_i}(v) = \sup\{\min_{kl}\mu_{kl}(W_{kl}) : c_i(W) = v\}},
#' has \eqn{\alpha}-cuts
#' \deqn{[\underline c_i(\alpha), \overline c_i(\alpha)] = [c_i(\underline W_\alpha),\, c_i(\overline W_\alpha)],}
#' where \eqn{\underline W_\alpha} and \eqn{\overline W_\alpha} collect the lower
#' and upper endpoints of the \eqn{\alpha}-cuts of the ties. Two crisp
#' evaluations per level are therefore enough. The class includes in- and
#' out-strength, ordered weighted averages of the ties, closeness based on any
#' continuous t-norm and Katz centrality; betweenness and eigenvector centrality
#' are excluded.
#'
#' The result is computed on a grid of levels. When the cores of the ties are
#' single values, the triangle \eqn{(c(L), c(M), c(R))} coincides with the fuzzy
#' centrality at \eqn{\alpha \in \{0, 1\}}; for affine indices (strength) the
#' endpoints are linear in \eqn{\alpha}, and the fuzzy centrality is triangular
#' and equal to that triangle. With unknown ties (\code{\link{set_unknown}}) the
#' 1-cut is an interval: the fuzzy centrality of an affine index is then
#' trapezoidal, and the triangle (built on the midpoint of the core) is an
#' approximation also at the core.
#'
#' @param x An \code{\link{fsn}} object.
#' @param index Index name or function, see \code{\link{centrality_index}}.
#' @param alpha Increasing grid of levels from 0 to 1.
#' @param ... Parameters of the index (e.g. \code{tnorm}, \code{max_steps},
#'   \code{w}, \code{beta}, \code{normalize}). Strength and Katz centrality
#'   are divided by \eqn{n - 1} unless \code{normalize = "none"}. For Katz
#'   centrality, \code{beta} defaults to \eqn{0.9/\rho(A)}, with \eqn{A} the
#'   binary adjacency matrix of the ties with positive support, as in Porreca
#'   and Maturo (2026).
#' @param check Logical, check that the cuts are nested (a violation indicates an
#'   index that is not isotone).
#' @return An object of class \code{"fuzzy_centrality"}: a list with matrices
#'   \code{lower} and \code{upper} (actors by levels), the grid \code{alpha}, the
#'   triangle \code{triangle} (a \code{\link{tfn}}), the network \code{network}
#'   (used to compute exact cuts at further levels, see
#'   \code{\link{possibility_matrix}}), the index specification and
#'   the labels.
#' @seealso \code{\link{triangular_approx}}, \code{\link{possibility_matrix}},
#'   \code{\link{rank_fuzzy}}, \code{\link{plot.fuzzy_centrality}}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' fc <- fuzzy_centrality(A, "in_closeness_prod")
#' fc
#' summary(fc)
#' plot(fc, actors = 6)
#' fuzzy_centrality(A, "in_strength", normalize = "n-1")
#' @export
#' @section Methodology:
#' Theorem (Section 4.1 of the article). Let \eqn{c} be continuous and non-decreasing in every tie on a box containing the support of the network. For every \eqn{\alpha \in [0, 1]}, the \eqn{\alpha}-cut of the fuzzy centrality defined by the extension principle is \eqn{[c_i(\underline W_\alpha), c_i(\overline W_\alpha)]}. The proof combines the commutation of the extension principle with \eqn{\alpha}-cuts for continuous functions on compact sets (Nguyen, 1978) with monotonicity, which places the minimum and the maximum of \eqn{c_i} over the box at its lowest and highest vertices.
#'
#' Consequences. (i) The fuzzy centrality is computed exactly with two crisp evaluations per level, without optimization or sampling. (ii) When the cores of the ties are single values, the triangle \eqn{(c(L), c(M), c(R))} coincides with it at the support and at the core; for nonlinear indices it is only approximate at intermediate levels (\code{\link{triangular_approx}}). With unknown ties the 1-cut is an interval and, for affine indices, the fuzzy centrality is trapezoidal. (iii) If the latent crisp network lies in the \eqn{\alpha}-cut of the fuzzy network, its centrality lies in the \eqn{\alpha}-cut of the fuzzy centrality (Proposition 2 of Porreca and Maturo, 2026). (iv) In undirected networks the box contains only symmetric networks; its lowest and highest vertices are symmetric, so the cuts remain exact, but the centralities of two actors share ties and their comparisons are bounds. The grid of levels affects only the resolution, not the exactness of the computed cuts. See Theorem 1 and Corollary 1 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- fsn_example()
#' # indices of the registry
#' sapply(c("in_strength", "in_closeness_min", "in_closeness_prod", "in_katz"),
#'        function(ix) round(fc_score(fuzzy_centrality(x, ix)), 3))
#' # parameters: step limit, OWA weights, Katz attenuation
#' fuzzy_centrality(x, "in_closeness", tnorm = "lukasiewicz", max_steps = 2)
#' fuzzy_centrality(x, "in_owa", w = c(0.4, 0.3, 0.2, 0.1, 0, 0, 0))
#' fuzzy_centrality(x, "out_katz", beta = 0.2)
#' # finer grid
#' fc <- fuzzy_centrality(x, "in_closeness_prod", alpha = seq(0, 1, by = 0.01))
#' plot(fc, actors = c("B", "D", "E"))
#' head(as.data.frame(fc))
fuzzy_centrality <- function(x, index = "in_strength", alpha = seq(0, 1, by = 0.05), ..., check = TRUE) {
  if (!is.fsn(x)) stop("'x' must be an fsn object")
  if (abs(alpha[1]) > 1e-12 || abs(alpha[length(alpha)] - 1) > 1e-12 || is.unsorted(alpha))
    stop("'alpha' must be an increasing grid from 0 to 1")
  spec <- centrality_index(index, ...)
  if (isTRUE(spec$needs_beta) && is.null(spec$params$beta)) spec <- .katz_beta(spec, x)
  if (isTRUE(spec$needs_beta)) spec$label <- sprintf("%s (beta = %.4g)", spec$label, spec$params$beta)
  n <- n_actors(x); G <- length(alpha)
  lo <- up <- matrix(NA_real_, n, G, dimnames = list(x$labels, format(alpha)))
  for (g in seq_len(G)) {
    ct <- fsn_cut(x, alpha[g])
    lo[, g] <- spec$fun(ct$lower); up[, g] <- spec$fun(ct$upper)
  }
  if (check) {
    tol <- 1e-8 * max(1, max(abs(up)))
    if (any(lo > up + tol)) stop("lower endpoint above upper endpoint: the index is not isotone")
    if (G > 1 && (any(diff(t(lo)) < -tol) || any(diff(t(up)) > tol)))
      stop("the cuts are not nested: the index is not isotone")
  }
  c1 <- fsn_cut(x, 1)
  tri <- tfn(lo[, 1], pmin(pmax(spec$fun((c1$lower + c1$upper) / 2), lo[, G]), up[, G]), up[, 1], names = x$labels)
  # comparisons of two actors are exact only if their indices depend on disjoint
  # sets of non-interactive ties: never in undirected networks, where a tie is
  # shared by its two endpoints
  if (!x$directed) spec$exact_pairs <- FALSE
  interval_core <- any(up[, G] - lo[, G] > 1e-10 * max(1, max(abs(up))))
  structure(list(lower = lo, upper = up, alpha = alpha, triangle = tri, index = spec,
                 labels = x$labels, respondents = x$respondents, directed = x$directed,
                 interval_core = interval_core, network = x,
                 cache = new.env(parent = emptyenv())), class = "fuzzy_centrality")
}

#' @export
print.fuzzy_centrality <- function(x, n = 10, digits = 4, ...) {
  cat("Exact fuzzy centrality:", x$index$label, "\n")
  cat(length(x$labels), "actors,", length(x$alpha), "levels;",
      if (isTRUE(x$interval_core)) "cores are intervals (unknown ties): the triangle is an approximation also at the core"
      else if (x$index$affine) "affine index (triangular, exact on the endpoints)" else "nonlinear index", "\n")
  d <- data.frame(actor = x$labels, support_lo = x$lower[, 1], core_lo = x$lower[, length(x$alpha)],
                  core_up = x$upper[, length(x$alpha)], support_up = x$upper[, 1],
                  cog = fc_score(x, "cog"))
  d <- d[order(-d$cog), ]; rownames(d) <- NULL
  print(utils::head(d, n), digits = digits)
  invisible(x)
}

#' @export
`[.fuzzy_centrality` <- function(x, i) {
  n <- length(x$labels)
  if (is.character(i)) {
    j <- match(i, x$labels)
    if (anyNA(j)) stop("unknown actor(s): ", paste(i[is.na(j)], collapse = ", "))
    i <- j
  } else if (is.logical(i)) {
    if (length(i) != n) stop("a logical selector must have one value per actor")
    i <- which(i)
  } else {
    i <- as.integer(i)
    if (anyNA(i) || any(i == 0)) stop("invalid numerical selector")
    i <- seq_len(n)[i]
    if (anyNA(i)) stop("selector out of range")
  }
  map <- if (is.null(x$actor_map)) seq_len(n) else x$actor_map
  x$actor_map <- map[i]
  x$lower <- x$lower[i, , drop = FALSE]; x$upper <- x$upper[i, , drop = FALSE]
  x$triangle <- x$triangle[i]; x$labels <- x$labels[i]
  if (!is.null(x$respondents)) x$respondents <- x$respondents[i]
  G <- ncol(x$lower)
  x$interval_core <- any(x$upper[, G] - x$lower[, G] > 1e-10 * max(1, max(abs(x$upper))))
  # 'network' and the cache of cuts refer to the whole network: the centralities
  # of the selected actors are those computed on it, not on the induced subnetwork
  x
}

#' @export
as.data.frame.fuzzy_centrality <- function(x, row.names = NULL, optional = FALSE, ...) {
  G <- length(x$alpha); n <- length(x$labels)
  data.frame(actor = rep(x$labels, G), alpha = rep(x$alpha, each = n),
             lower = as.vector(x$lower), upper = as.vector(x$upper), stringsAsFactors = FALSE)
}

# trapezoidal weights on the alpha grid
.trapz_w <- function(a) { G <- length(a); if (G == 1) return(1); d <- diff(a); (c(d, 0) + c(0, d)) / 2 }

#' Ranking scores of exact fuzzy centralities
#'
#' Computes the scores of \code{\link{tfn_score}} for fuzzy numbers given by
#' their \eqn{\alpha}-cuts, by trapezoidal integration over the grid:
#' centre of gravity \eqn{\int (a^{+2}_\alpha - a^{-2}_\alpha)/2\, d\alpha \big/
#' \int (a^+_\alpha - a^-_\alpha) d\alpha}, core (midpoint of the 1-cut), Yager's
#' index \eqn{\int (a^-_\alpha + a^+_\alpha)/2\, d\alpha}, possibilistic mean
#' \eqn{\int \alpha (a^-_\alpha + a^+_\alpha) d\alpha} (\code{"graded_mean"}) and
#' the Liou-Wang total integral value.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param method See \code{\link{fuzzy_rank_methods}}.
#' @param lambda Optimism index for \code{"liou_wang"}.
#' @param triangle Logical: compute the score of the triangular approximation.
#' @return Named numeric vector.
#' @examples
#' fc <- fuzzy_centrality(fsn_example(), "in_closeness_min")
#' cbind(exact = fc_score(fc, "cog"), triangle = fc_score(fc, "cog", triangle = TRUE))
#' @export
#' @section Methodology:
#' Scores are integrals over the \eqn{\alpha}-cuts, approximated by the trapezoidal rule on the grid; for strength they coincide with the closed forms of \code{\link{tfn_score}} up to the integration error. Comparing the score of the exact fuzzy centrality with that of its triangle shows how much the triangular approximation changes a ranking. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_min")
#' d <- fc_score(fc, "cog") - fc_score(fc, "cog", triangle = TRUE)
#' summary(d)
#' cor(fc_score(fc, "cog"), fc_score(fc, "cog", triangle = TRUE), method = "kendall")
fc_score <- function(x, method = "cog", lambda = 0.5, triangle = FALSE) {
  method <- match.arg(method, fuzzy_rank_methods())
  .check_lambda(lambda)
  if (triangle) return(tfn_score(x$triangle, method, lambda))
  a <- x$alpha; w <- .trapz_w(a); lo <- x$lower; up <- x$upper; G <- length(a)
  out <- switch(method,
                cog = {
                  area <- as.vector((up - lo) %*% w)
                  mom <- as.vector(((up^2 - lo^2) / 2) %*% w)
                  ifelse(area > 1e-14, mom / pmax(area, 1e-300), (lo[, G] + up[, G]) / 2)
                },
                core = (lo[, G] + up[, G]) / 2,
                yager = as.vector(((lo + up) / 2) %*% w),
                graded_mean = as.vector(((lo + up) * matrix(a, nrow(lo), G, byrow = TRUE)) %*% w),
                liou_wang = as.vector(((1 - lambda) * lo + lambda * up) %*% w))
  names(out) <- x$labels; out
}

#' Triangular approximation of an exact fuzzy centrality
#'
#' The triangle \eqn{(c(L), c(M), c(R))} built on the endpoints coincides with the
#' exact fuzzy centrality at the support and, when the cores of the ties are
#' single values, at the core. With unknown ties (\code{\link{set_unknown}}) the
#' 1-cut of the fuzzy centrality is an interval: the triangle, whose core is the
#' centrality of the midpoint network, is then an approximation also at
#' \eqn{\alpha = 1}, and for affine indices the exact fuzzy centrality is
#' trapezoidal. This function returns
#' the triangle and the error of its \eqn{\alpha}-cuts at every level, relative to
#' the width of the exact support.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @return A list with the triangle (\code{tfn}), the matrices of its cut
#'   endpoints, the matrix \code{error} (largest absolute endpoint error divided
#'   by the support width, actors by levels) and \code{max_error} by actor.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' ta <- triangular_approx(fc)
#' summary(ta$max_error)
#' @export
#' @section Methodology:
#' Earlier computations of fuzzy centralities return a TFN. When the endpoints are computed correctly, the triangle \eqn{(c(L), c(M), c(R))} is exact at the support and at the core but not in between for nonlinear indices, whose exact membership has curved sides and may be markedly asymmetric. The error is reported relative to the width of the support.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' ta <- triangular_approx(fc)
#' round(apply(ta$error, 2, max), 3)       # largest error by level
#' plot(fc, actors = names(sort(ta$max_error, decreasing = TRUE))[1:3])
triangular_approx <- function(x) {
  ct <- alpha_cut(x$triangle, x$alpha)
  width <- pmax(x$upper[, 1] - x$lower[, 1], 1e-14)
  err <- pmax(abs(ct$lower - x$lower), abs(ct$upper - x$upper)) / width
  err[x$upper[, 1] - x$lower[, 1] <= 1e-14, ] <- 0
  dimnames(err) <- dimnames(x$lower)
  list(triangle = x$triangle, lower = ct$lower, upper = ct$upper, error = err,
       max_error = apply(err, 1, max))
}

#' Monte Carlo inner approximation of a fuzzy centrality
#'
#' Draws \code{B} networks uniformly in the box \eqn{\mathcal W_\alpha} and
#' returns, for each actor, the minimum and maximum of the centrality over the
#' draws. The result is an inner approximation of the exact \eqn{\alpha}-cut:
#' extreme values are attained at the vertices of the box and are rarely
#' sampled.
#' @param x An \code{\link{fsn}} object.
#' @param index Index name or function.
#' @param B Number of draws.
#' @param alpha Level.
#' @param ... Parameters of the index.
#' @param keep_draws Logical, return the draws.
#' @return A list with vectors \code{lower} and \code{upper} and, optionally, the
#'   matrix of draws.
#' @examples
#' x <- fsn_simulated(n = 20, respondents = 10)
#' mc <- mc_centrality(x, "in_strength", B = 200)
#' fc <- fuzzy_centrality(x, "in_strength")
#' summary((mc$upper - mc$lower) / (fc$upper[, 1] - fc$lower[, 1]))
#' @export
#' @section Methodology:
#' Sampling networks in the box and taking the range of the centralities gives an inner approximation of the \eqn{\alpha}-cut: the extremes are attained only at the vertices of a box of dimension equal to the number of ties, which uniform sampling almost never reaches. The approximation therefore underestimates the vagueness of the centralities, increasingly so as the network grows.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- fsn_simulated(n = 20, respondents = 10)
#' fc <- fuzzy_centrality(x, "in_strength")
#' w <- fc$upper[, 1] - fc$lower[, 1]
#' sapply(c(10, 100, 1000), function(B) {
#'   mc <- mc_centrality(x, "in_strength", B = B)
#'   mean((mc$upper - mc$lower)[w > 0] / w[w > 0])
#' })
mc_centrality <- function(x, index = "in_strength", B = 500, alpha = 0, ..., keep_draws = FALSE) {
  spec <- centrality_index(index, ...)
  if (isTRUE(spec$needs_beta) && is.null(spec$params$beta)) spec <- .katz_beta(spec, x)
  ct <- fsn_cut(x, alpha); n <- n_actors(x)
  D <- matrix(NA_real_, B, n, dimnames = list(NULL, x$labels))
  for (b in seq_len(B)) {
    U <- matrix(stats::runif(n * n), n, n)
    if (!x$directed) U[lower.tri(U)] <- t(U)[lower.tri(U)]   # one draw per undirected tie
    W <- ct$lower + (ct$upper - ct$lower) * U
    D[b, ] <- spec$fun(W)
  }
  out <- list(lower = apply(D, 2, min), upper = apply(D, 2, max), B = B, alpha = alpha)
  if (keep_draws) out$draws <- D
  out
}

#' Crisp centrality of the latent network of a simulated fuzzy network
#'
#' @param x A network from \code{\link{simulate_fsn}}.
#' @param index Index name or function.
#' @param ... Parameters of the index (the same used for the fuzzy centrality;
#'   for Katz centrality pass the \code{beta} of the fuzzy centrality).
#' @return Numeric vector.
#' @examples
#' x <- simulate_fsn(30, 0.2, seed = 1)
#' fc <- fuzzy_centrality(x, "in_strength")
#' mean(coverage(fc, latent_centrality(x, "in_strength")))
#' @export
#' @section Methodology:
#' The centrality of the latent crisp network of a simulated fuzzy network is the reference against which coverage and ranking accuracy are evaluated in \code{\link{simulation_study}}.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- simulate_fsn(40, 0.25, kappa = 1, seed = 6)
#' fc <- fuzzy_centrality(x, "in_strength")
#' truth <- latent_centrality(x, "in_strength")
#' discordance(fc_score(fc, "cog"), truth)
latent_centrality <- function(x, index = "in_strength", ...) {
  if (is.null(x$latent)) stop("'x' has no latent network")
  spec <- centrality_index(index, ...)
  if (isTRUE(spec$needs_beta) && is.null(spec$params$beta)) spec <- .katz_beta(spec, x)
  stats::setNames(spec$fun(x$latent), x$labels)
}

#' Coverage of crisp values by the alpha-cuts of a fuzzy centrality
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param value Numeric vector, one value per actor.
#' @param alpha Level (must belong to the grid).
#' @return Logical vector.
#' @examples
#' x <- simulate_fsn(30, 0.2, kappa = 1, seed = 1)
#' fc <- fuzzy_centrality(x, "in_closeness_min")
#' mean(coverage(fc, latent_centrality(x, "in_closeness_min")))
#' @export
#' @section Methodology:
#' Coverage at level \eqn{\alpha} is guaranteed when the latent network lies in \eqn{\mathcal W_\alpha}; otherwise it is an empirical property, usually high for indices that aggregate many ties because the errors on individual ties partly cancel.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- simulate_fsn(40, 0.3, kappa = 1, seed = 7)
#' fc <- fuzzy_centrality(x, "in_closeness_min")
#' truth <- latent_centrality(x, "in_closeness_min")
#' sapply(c(0, 0.5, 1), function(a) mean(coverage(fc, truth, a)))
coverage <- function(x, value, alpha = 0) {
  g <- .alpha_index(x, alpha)
  tol <- 1e-10
  value >= x$lower[, g] - tol & value <= x$upper[, g] + tol
}

.alpha_index <- function(x, alpha) {
  g <- which(abs(x$alpha - alpha) < 1e-9)
  if (!length(g)) stop("alpha = ", alpha, " is not in the grid of the fuzzy centrality")
  g[1]
}

#' @export
summary.fuzzy_centrality <- function(object, alpha = 0.5, method = "cog", ...) {
  x <- object; G <- length(x$alpha); g <- .alpha_index(x, alpha)
  ri0 <- rank_intervals(x, 0); ria <- rank_intervals(x, alpha)
  d <- data.frame(actor = x$labels,
                  support_lo = x$lower[, 1], support_up = x$upper[, 1],
                  cut_lo = x$lower[, g], cut_up = x$upper[, g],
                  core_lo = x$lower[, G], core_up = x$upper[, G],
                  core_midpoint = (x$lower[, G] + x$upper[, G]) / 2,
                  score = fc_score(x, method), stringsAsFactors = FALSE)
  d$rank <- rank(-d$score, ties.method = "min")
  d$best_rank_0 <- ri0$best; d$worst_rank_0 <- ri0$worst
  d$best_rank_a <- ria$best; d$worst_rank_a <- ria$worst
  d$top_membership <- top_set(x)
  d <- d[order(d$rank), ]; rownames(d) <- NULL
  attr(d, "alpha") <- alpha; attr(d, "method") <- method; attr(d, "index") <- x$index$label
  class(d) <- c("summary.fuzzy_centrality", "data.frame"); d
}

#' @export
print.summary.fuzzy_centrality <- function(x, n = 15, digits = 3, ...) {
  cat("Exact fuzzy centrality:", attr(x, "index"), "\n")
  cat("Ranked by", attr(x, "method"), "; cut at alpha =", attr(x, "alpha"),
      "; rank intervals at alpha = 0 and at alpha =", attr(x, "alpha"), "\n")
  print(utils::head(as.data.frame(unclass(x)), n), digits = digits)
  invisible(x)
}
