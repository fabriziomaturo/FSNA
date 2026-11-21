# ---------------------------------------------------------------------------
# Possibilistic comparison of actors, certified orderings, rankings
# ---------------------------------------------------------------------------

# supremum of the levels at which D >= 0 (or D > 0), D non-increasing in alpha;
# D: matrix (items x levels). rule: "interpolate", "upper", "grid"
.sup_level <- function(D, alpha, strict = FALSE, rule = "upper") {
  G <- length(alpha)
  ok <- if (strict) D > 0 else D >= 0
  k <- if (G == 1) as.integer(ok[, 1]) else apply(ok, 1, function(z) if (any(z)) max(which(z)) else 0L)
  out <- numeric(nrow(D))
  out[k == G] <- 1
  mid <- which(k > 0 & k < G)
  if (length(mid)) {
    a0 <- alpha[k[mid]]; a1 <- alpha[k[mid] + 1]
    out[mid] <- switch(rule,
                       grid = a0, upper = a1,
                       interpolate = {
                         d0 <- D[cbind(mid, k[mid])]; d1 <- D[cbind(mid, k[mid] + 1)]
                         a0 + (a1 - a0) * d0 / pmax(d0 - d1, 1e-300)
                       })
  }
  out
}

.rule <- function(x, bracket) {
  bracket <- match.arg(bracket, c("auto", "interpolate", "bisection", "upper", "grid"))
  if (bracket == "interpolate" && !isTRUE(x$index$affine))
    stop("linear interpolation is exact only for affine indices (strength); ",
         "use bracket = \"bisection\" (default for this index) or \"upper\"")
  if (bracket == "bisection" && is.null(x$network))
    stop("bisection needs the network: recompute the object with fuzzy_centrality()")
  if (bracket == "auto")
    if (isTRUE(x$index$affine)) "interpolate" else if (!is.null(x$network)) "bisection" else "upper"
  else bracket
}

# exact cuts of the index at an arbitrary level, for all the actors of the
# network on which the centrality was computed; cached by level in the object's
# environment, which is shared by its subsets (the cuts do not depend on the
# selection of actors)
.cut_fun <- function(x) {
  cache <- if (is.environment(x$cache)) x$cache else new.env(parent = emptyenv())
  function(a) {
    key <- sprintf("cut_%.17g", a)
    if (is.null(cache[[key]])) {
      ct <- fsn_cut(x$network, a)
      cache[[key]] <- list(lower = x$index$fun(ct$lower), upper = x$index$fun(ct$upper))
    }
    cache[[key]]
  }
}

# position of the actors of x in the network on which the centrality was computed
.actor_map <- function(x) if (is.null(x$actor_map)) seq_along(x$labels) else x$actor_map

.check_tol <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0)
    stop("'tol' must be a single finite positive number")
  tol
}

# Suprema sup{alpha : D_ij(alpha) >= 0}, D_ij = upper_i - lower_j, for the pairs
# (ii, jj) of rows of x. Returns the estimate and a bracket [lo, hi] that
# contains the supremum of the marginal degree: lo is a level at which the
# condition holds (or 0 when it fails at 0), hi a level at which it fails (or
# the supremum itself when it is attained exactly).
.sup_bounds <- function(x, ii, jj, rule, tol = 1e-4) {
  a <- x$alpha; G <- length(a)
  D <- x$upper[ii, , drop = FALSE] - x$lower[jj, , drop = FALSE]
  ok <- D >= 0
  k <- if (G == 1) as.integer(ok[, 1]) else apply(ok, 1, function(z) if (any(z)) max(which(z)) else 0L)
  est <- lo <- hi <- numeric(length(ii))
  est[k == G] <- lo[k == G] <- hi[k == G] <- 1
  mid <- which(k > 0 & k < G)
  if (!length(mid)) return(list(estimate = est, lower = lo, upper = hi))
  a0 <- a[k[mid]]; a1 <- a[k[mid] + 1]
  if (rule == "interpolate") {
    d0 <- D[cbind(mid, k[mid])]; d1 <- D[cbind(mid, k[mid] + 1)]
    v <- a0 + (a1 - a0) * d0 / pmax(d0 - d1, 1e-300)
    est[mid] <- lo[mid] <- hi[mid] <- v
  } else if (rule %in% c("grid", "upper")) {
    lo[mid] <- a0; hi[mid] <- a1; est[mid] <- if (rule == "grid") a0 else a1
  } else {
    blo <- a0; bhi <- a1
    map <- .actor_map(x); cut_at <- .cut_fun(x)
    pi_all <- map[ii[mid]]; pj_all <- map[jj[mid]]
    max_iter <- ceiling(log2(max(a1 - a0) / tol)) + 2
    it <- 0
    while (any(act <- (bhi - blo) > tol) && it < max_iter) {
      it <- it + 1
      m <- (blo[act] + bhi[act]) / 2
      stag <- m <= blo[act] | m >= bhi[act]          # no representable midpoint
      if (all(stag)) break
      w <- which(act)[!stag]; m <- m[!stag]
      um <- unique(m); cuts <- lapply(um, cut_at); idx <- match(m, um)
      UP <- matrix(vapply(cuts, `[[`, numeric(length(cuts[[1]]$upper)), "upper"), ncol = length(um))
      LO <- matrix(vapply(cuts, `[[`, numeric(length(cuts[[1]]$lower)), "lower"), ncol = length(um))
      d <- UP[cbind(pi_all[w], idx)] - LO[cbind(pj_all[w], idx)]
      blo[w[d >= 0]] <- m[d >= 0]; bhi[w[d < 0]] <- m[d < 0]
    }
    lo[mid] <- blo; hi[mid] <- bhi; est[mid] <- blo
  }
  list(estimate = est, lower = lo, upper = hi)
}

#' Possibility and necessity degrees of the orderings of actors
#'
#' For two actors \eqn{i, j} with fuzzy centralities \eqn{\tilde c_i, \tilde c_j}
#' (Dubois and Prade, 1983):
#' \deqn{\Pi(i \succeq j) = \sup\{\alpha : \overline c_i(\alpha) \ge \underline c_j(\alpha)\},\qquad
#'       N(i \succ j) = 1 - \Pi(j \succeq i)}
#' (Definition 2 and Eq. 7 of Porreca and Maturo, 2026).
#' These are \emph{marginal} degrees, computed from the fuzzy centralities of the
#' two actors separately. They coincide with the degrees of the joint
#' statement on the network when the two indices depend on disjoint sets of
#' ties (in- or out-strength, OWA of the ties); otherwise \eqn{\Pi} is an upper
#' bound and \eqn{N} a lower bound of the joint degrees.
#'
#' The supremum lies between the last level of the grid at which the condition
#' holds and the next one (Algorithm 1 of Porreca and Maturo, 2026). With
#' \code{bracket = "auto"} it is found by linear interpolation for affine indices
#' (strength), which is exact, and by bisection on \eqn{\alpha} inside the
#' bracket for the other indices, computing the exact cuts at each trial level
#' with two evaluations of the crisp index; the bisection stops when the bracket
#' is shorter than \code{tol} and returns the last level at which the condition
#' holds, so that a grid level, in particular 1, is returned only if the
#' condition holds there. Trial levels are shared across pairs, and the exact
#' cuts computed at them are cached in the object.
#'
#' The bisection gives a numerical approximation of the marginal degree, not a
#' bound: the returned level can lie below the supremum by up to \code{tol}. The
#' result therefore also contains, for every pair, a bracket
#' \eqn{[\Pi_{lo}, \Pi_{hi}]} that contains the marginal possibility (the level
#' at which the condition was last found to hold and the level at which it was
#' found to fail), and the corresponding necessity bracket
#' \eqn{[1 - \Pi_{hi}(j \succeq i), 1 - \Pi_{lo}(j \succeq i)]}. Statements that
#' must hold for certain, such as a dominance with \eqn{N = 1}, should use the
#' conservative ends (\code{necessity_lower}, \code{possibility_upper}); a
#' necessity estimate equal to one does not imply \code{necessity_lower = 1}.
#' A guaranteed marginal necessity equal to one (\code{necessity_lower}
#' exactly 1) means that the \eqn{\alpha}-cuts of the two actors are strictly
#' separated at every \eqn{\alpha > 0}; the closed supports may still touch at
#' an endpoint, where the membership is zero. It is therefore implied by, but
#' weaker than, certification at \eqn{\alpha = 0} (\code{\link{certified_pairs}}).
#' Comparisons with one are exact: values within \eqn{10^{-12}} of one are
#' reported separately as numerically equal to one, not as guaranteed.
#' Linear interpolation (\code{bracket = "interpolate"}) is refused for
#' indices that are not affine, for which it is not exact.
#' \code{bracket = "upper"} takes the upper end of the grid bracket, a quick
#' alternative that keeps \eqn{\Pi} an upper bound and \eqn{N} a lower bound of
#' the marginal degrees within one grid step. \code{bracket = "grid"} takes the
#' largest grid level at which the condition holds; it underestimates \eqn{\Pi}
#' and overestimates \eqn{N}, and is provided only for comparison.
#'
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param bracket \code{"auto"}, \code{"interpolate"}, \code{"bisection"},
#'   \code{"upper"} or \code{"grid"}.
#' @param tol Tolerance of the bisection on \eqn{\alpha} (default \eqn{10^{-4}}, as in
#'   Porreca and Maturo, 2026).
#' @return An object of class \code{"fuzzy_comparison"}: a list with the matrices
#'   \code{possibility} (estimate of \eqn{\Pi(i \succeq j)}, row \eqn{i},
#'   column \eqn{j}), \code{possibility_lower} and \code{possibility_upper}
#'   (numerical bracket of the marginal possibility; equal to the estimate when it
#'   is computed exactly, by interpolation or at the ends of the grid),
#'   \code{necessity} (estimate of \eqn{N(i \succ j)}, zero on the diagonal),
#'   \code{necessity_lower} and \code{necessity_upper} (bracket of the marginal
#'   necessity), the rule used, the tolerance of the bisection and a flag
#'   \code{exact}: \code{TRUE} when the marginal degrees coincide with the joint
#'   ones (in- and out-strength and OWA of the incoming or outgoing ties, in
#'   directed networks) and are computed by interpolation or bisection; the
#'   numerical brackets still apply.
#' @references Dubois D., Prade H. (1983). Ranking fuzzy numbers in the setting
#'   of possibility theory. \emph{Information Sciences} 30, 183-224.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' cmp <- possibility_matrix(fc)
#' cmp
#' round(cmp$necessity[1:5, 1:5], 2)
#' @export
#' @section Methodology:
#' For fuzzy quantities \eqn{\tilde c_i} and \eqn{\tilde c_j}, the possibility that \eqn{i} is at least as central as \eqn{j} is \eqn{\Pi(i \succeq j) = \sup_{u \ge v}\min(\mu_i(u), \mu_j(v))}, which for fuzzy numbers equals the largest \eqn{\alpha} at which the upper endpoint of \eqn{i} reaches the lower endpoint of \eqn{j}. The necessity \eqn{N(i \succ j) = 1 - \Pi(j \succeq i)} is the degree to which the ordering is implied by the information. \eqn{\Pi = 1} and \eqn{N = 0} mean that the ordering is plausible but not supported; \eqn{N = 1} means that it holds at every level above zero.
#'
#' These are marginal degrees, computed from each actor's fuzzy centrality separately. When the centralities of two actors depend on common ties, the joint possibility of the ordering may be smaller than the marginal one; marginal degrees are then bounds. The suprema lie between two grid levels: interpolation is exact for strength; for the other indices the supremum is located by bisection, since the exact cuts can be computed at any level with two evaluations of the crisp index (Algorithm 1 of Porreca and Maturo, 2026). Taking the suprema on the grid overestimates necessity: when the endpoints cross between \eqn{\alpha = 0} and the first positive level it gives \eqn{N = 1} to an ordering that is not certified. See Sections 2.4 and 2.5 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' cmp <- possibility_matrix(fc)
#' top3 <- names(sort(fc$lower[, ncol(fc$lower)], decreasing = TRUE))[1:3]   # highest cores
#' round(cmp$possibility[top3, top3], 3)
#' round(cmp$necessity[top3, top3], 3)
#' # grid suprema overestimate necessity
#' sum(possibility_matrix(fc, "grid")$necessity > cmp$necessity + 1e-12)
#' plot(cmp, "possibility")
possibility_matrix <- function(x, bracket = "auto", tol = 1e-4) {
  rule <- .rule(x, bracket); n <- length(x$labels)
  if (rule == "bisection") .check_tol(tol)
  key <- paste("P", rule, if (rule == "bisection") format(tol, digits = 17) else "",
               paste(.actor_map(x), collapse = ","), sep = "|")
  res <- if (is.environment(x$cache)) x$cache[[key]] else NULL
  if (is.null(res)) {
    ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
    off <- ii != jj
    b <- .sup_bounds(x, ii[off], jj[off], rule, tol)
    mk <- function(v) { m <- matrix(1, n, n); m[off] <- v; m }
    res <- list(P = mk(b$estimate), Plo = mk(b$lower), Phi = mk(b$upper))
    if (is.environment(x$cache)) x$cache[[key]] <- res
  }
  P <- res$P; Plo <- res$Plo; Phi <- res$Phi
  N <- 1 - t(P); Nlo <- 1 - t(Phi); Nhi <- 1 - t(Plo)
  diag(N) <- diag(Nlo) <- diag(Nhi) <- 0
  dn <- list(x$labels, x$labels)
  dimnames(P) <- dimnames(Plo) <- dimnames(Phi) <- dimnames(N) <- dimnames(Nlo) <- dimnames(Nhi) <- dn
  structure(list(possibility = P, possibility_lower = Plo, possibility_upper = Phi,
                 necessity = N, necessity_lower = Nlo, necessity_upper = Nhi, rule = rule,
                 exact = isTRUE(x$index$exact_pairs) && !isFALSE(x$directed) && rule %in% c("interpolate", "bisection"),
                 tol = if (rule == "bisection") tol else NA_real_,
                 index = x$index$label, labels = x$labels), class = "fuzzy_comparison")
}

#' @export
print.fuzzy_comparison <- function(x, ...) {
  n <- length(x$labels); off <- row(x$necessity) != col(x$necessity)
  cat("Possibilistic comparison of", n, "actors:", x$index, "\n")
  cat("Supremum rule:", x$rule,
      if (x$rule == "bisection") paste0("(tol = ", format(x$tol), ")") else "", "-",
      if (x$exact) "marginal degrees equal to the joint ones" else "marginal degrees (bounds of the joint ones)", "\n")
  w <- max(x$possibility_upper[off] - x$possibility_lower[off], 0)
  cat("Maximum width of the numerical brackets:", format(signif(w, 3)), "\n")
  cat("Ordered pairs with N = 1: guaranteed (lower end of the bracket exactly 1)",
      sum(x$necessity_lower[off] == 1), "; estimate numerically 1 (>= 1 - 1e-12)",
      sum(x$necessity[off] >= 1 - 1e-12), "\n")
  cat("Unordered pairs with Pi = 1 in both directions:",
      sum(x$possibility_lower[off] == 1 & t(x$possibility_lower)[off] == 1) / 2, "\n")
  invisible(x)
}

#' Certified orderings
#'
#' The ordering \eqn{i \succ j} is \emph{certified} at level \eqn{\alpha} when the
#' \eqn{\alpha}-cuts are strictly separated,
#' \eqn{\underline c_i(\alpha) > \overline c_j(\alpha)}: then \eqn{c_i(W) > c_j(W)}
#' for every network \eqn{W} in the box \eqn{\mathcal W_\alpha}. At
#' \eqn{\alpha = 0} the ordering holds for every network compatible with the
#' supports of the ties. A necessity degree equal to one does not guarantee a
#' strict ordering, because the supports may touch.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param alpha Level (in the grid).
#' @return \code{certified_pairs}: logical matrix, \code{TRUE} in row \eqn{i},
#'   column \eqn{j} when \eqn{i \succ j} is certified. \code{certification}: a
#'   list with the proportion of unordered pairs certified at each level and the
#'   number of ordered pairs with \eqn{N = 1} whose supports are not strictly
#'   separated: with the estimate of \code{\link{possibility_matrix}}
#'   numerically equal to one (\code{N1_not_separated}), with its conservative
#'   lower end exactly equal to one (\code{N1_not_separated_guaranteed}: pairs
#'   whose cuts are separated at every \eqn{\alpha > 0} but whose supports touch
#'   at an endpoint; usually zero, not zero by construction) and with the
#'   suprema taken as the largest grid level (\code{N1_not_separated_grid}).
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' C0 <- certified_pairs(fc, 0)
#' sum(C0)
#' certification(fc)
#' @export
#' @section Methodology:
#' An ordering \eqn{i \succ j} is certified at level \eqn{\alpha} when \eqn{\underline c_i(\alpha) > \overline c_j(\alpha)}: by the theorem of \code{\link{fuzzy_centrality}}, \eqn{c_i(W) > c_j(W)} for every network \eqn{W} in \eqn{\mathcal W_\alpha}, and in particular for the latent network whenever it lies in the box. A necessity degree equal to one does not imply certification, because supports may touch without overlapping. An ordering certified at \eqn{\alpha = 0} agrees with the ordering of the centres of gravity of the triangles, because every parameter of one triangle exceeds every parameter of the other; an ordering certified only at a positive level is conditional on the corresponding box and may disagree with it (e.g. \eqn{(0, 0.1, 1)} and \eqn{(0.2, 0.21, 0.22)} at \eqn{\alpha = 0.95}). Certification adds the separation of reliable orderings from undecided ones. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' C0 <- certified_pairs(fc, 0); C5 <- certified_pairs(fc, 0.5)
#' c(alpha0 = sum(C0), alpha05 = sum(C5))
#' # actors certified below the actor with the highest core at alpha = 0.5
#' best <- names(which.max(fc$lower[, ncol(fc$lower)]))
#' names(which(C5[best, ]))
#' certification(fc, levels = seq(0, 1, by = 0.25))$proportion_certified
certified_pairs <- function(x, alpha = 0) {
  g <- .alpha_index(x, alpha)
  C <- outer(x$lower[, g], x$upper[, g], ">")
  diag(C) <- FALSE; dimnames(C) <- list(x$labels, x$labels); C
}

#' @rdname certified_pairs
#' @param levels Levels at which the proportion of certified pairs is computed.
#' @export
certification <- function(x, levels = c(0, 0.5, 1)) {
  n <- length(x$labels); np <- n * (n - 1) / 2
  prop <- vapply(levels, function(a) sum(certified_pairs(x, a)) / np, numeric(1))
  names(prop) <- paste0("alpha=", levels)
  C0 <- certified_pairs(x, 0); off <- row(C0) != col(C0)
  pm <- possibility_matrix(x); Ng <- possibility_matrix(x, "grid")$necessity
  list(proportion_certified = prop,
       N1_not_separated = sum(pm$necessity[off] >= 1 - 1e-12 & !C0[off]),
       N1_not_separated_guaranteed = sum(pm$necessity_lower[off] == 1 & !C0[off]),
       N1_not_separated_grid = sum(Ng[off] >= 1 - 1e-12 & !C0[off]))
}

#' Fuzzy set of the most central actor
#'
#' Membership of actor \eqn{i} in the fuzzy set of the most central actor:
#' \eqn{\mu_{\mathrm{top}}(i) = \sup\{\alpha : \overline c_i(\alpha) \ge \max_{j \ne i}\underline c_j(\alpha)\}
#' = \min_{j \ne i} \Pi(i \succeq j)}, the possibility that \eqn{i} is the most
#' central actor. The marginal degree coincides with the joint one for in- and
#' out-strength and OWA of the incoming or outgoing ties in directed networks,
#' and is an upper bound of the joint possibility for the other indices. When the
#' possibility degrees are computed by bisection the membership is a numerical
#' estimate with a bracket, returned as attributes.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param bracket,tol See \code{\link{possibility_matrix}}.
#' @param alpha Optional level: if given, the actors in the \eqn{\alpha}-cut of
#'   the fuzzy set are returned.
#' @param which With \code{alpha}: membership used for the cut,
#'   \code{"estimate"} (default), \code{"guaranteed"} (lower end of the
#'   numerical bracket: actors certainly in the cut of the marginal fuzzy set) or
#'   \code{"possible"} (upper end: actors that may be in it).
#' @return Without \code{alpha}, the named vector of estimated membership
#'   degrees, with attributes \code{lower} and \code{upper} (numerical
#'   bracket); with \code{alpha}, a character vector of actors.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' sort(top_set(fc), decreasing = TRUE)[1:5]
#' top_set(fc, alpha = 0.5)
#' @export
#' @section Methodology:
#' The fuzzy set of the most central actor gives, for each actor, the possibility of being the most central. Its \eqn{\alpha}-cut is the set of actors that cannot be excluded from the top position at level \eqn{\alpha}, a set-valued answer that replaces the single winner of a crisp ranking. For in- and out-strength and OWA of the incoming or outgoing ties in directed networks the marginal degrees coincide with the joint ones; otherwise (other indices, undirected networks) they are upper bounds and the cuts may contain too many actors. When the degrees are computed by bisection, they are numerical estimates within the tolerance: \code{which = "guaranteed"} and \code{"possible"} use the two ends of the numerical bracket. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' mu <- top_set(fc)
#' sort(mu[mu > 0], decreasing = TRUE)
#' sapply(c(0, 0.25, 0.5, 0.75, 1), function(a) length(top_set(fc, alpha = a)))
top_set <- function(x, bracket = "auto", alpha = NULL, tol = 1e-4,
                    which = c("estimate", "guaranteed", "possible")) {
  which <- match.arg(which)
  pm <- possibility_matrix(x, bracket, tol)
  mu <- apply(pm$possibility, 1, min)
  mlo <- apply(pm$possibility_lower, 1, min); mhi <- apply(pm$possibility_upper, 1, min)
  if (is.null(alpha)) return(structure(mu, lower = mlo, upper = mhi))
  v <- switch(which, estimate = mu, guaranteed = mlo, possible = mhi)
  names(v)[v >= alpha - 1e-12]
}

#' Possible ranks of the actors
#'
#' At level \eqn{\alpha}, the best possible position of actor \eqn{i} is one plus
#' the number of actors whose lower endpoint exceeds the upper endpoint of
#' \eqn{i}, and the worst possible position is \eqn{n} minus the number of actors
#' whose upper endpoint is below the lower endpoint of \eqn{i}. These are
#' \emph{positions in a complete ordering}, i.e. after breaking ties: two actors
#' with the same crisp centrality share the positions 1 and 2 and both get the
#' interval \eqn{[1, 2]}, whereas \code{\link{rank_fuzzy}} gives tied scores the
#' same rank (\code{ties.method = "min"}). The interval is attained by a
#' compatible network for in- and out-strength in directed networks and is a bound otherwise.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param alpha Level in the grid.
#' @return A data frame with \code{actor}, \code{best} and \code{worst}.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' head(rank_intervals(fc, 0.5))
#' @export
#' @section Methodology:
#' The best and worst ranks compatible with the \eqn{\alpha}-cuts summarise how far the position of an actor is determined by the information. For in- and out-strength in directed networks the extremes are attained by a compatible network, since the centralities of different actors depend on disjoint ties; for the other indices and in undirected networks, where a tie is shared by its two endpoints, they are bounds. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength")
#' ri <- rank_intervals(fc, 0.5)
#' ri$width <- ri$worst - ri$best
#' head(ri[order(ri$best), ], 10)
rank_intervals <- function(x, alpha = 0) {
  g <- .alpha_index(x, alpha); lo <- x$lower[, g]; up <- x$upper[, g]; n <- length(lo)
  best <- 1 + vapply(seq_len(n), function(i) sum(lo[-i] > up[i]), numeric(1))
  worst <- n - vapply(seq_len(n), function(i) sum(up[-i] < lo[i]), numeric(1))
  data.frame(actor = x$labels, best = best, worst = worst, stringsAsFactors = FALSE)
}

#' Rank actors by their fuzzy centralities
#'
#' Ranks the actors with one of several methods:
#' \describe{
#'   \item{score-based}{\code{"cog"}, \code{"core"}, \code{"yager"},
#'     \code{"graded_mean"}, \code{"liou_wang"} (see \code{\link{tfn_score}} and
#'     \code{\link{fc_score}});}
#'   \item{\code{"possibility"}}{net possibility flow
#'     \eqn{\sum_j [\Pi(i\succeq j) - \Pi(j \succeq i)]};}
#'   \item{\code{"necessity"}}{number of actors \eqn{j} with \eqn{N(i \succ j) = 1}
#'     guaranteed (lower end of the numerical bracket exactly equal to one) minus the
#'     number with \eqn{N(j \succ i) = 1} guaranteed, ties broken by the sum of
#'     the estimated necessity degrees;}
#'   \item{\code{"certified"}}{number of actors certified below \eqn{i} at
#'     \eqn{\alpha = 0} minus the number certified above: a heuristic net
#'     dominance score, not the certified partial order (see Methodology);}
#'   \item{\code{"top"}}{membership in the fuzzy set of the most central actor
#'     (\code{\link{top_set}}).}
#' }
#' @param x A \code{\link{fuzzy_centrality}} object, a \code{\link{tfn}} vector or
#'   an \code{"ijar_index"} object (score-based methods only for the last two).
#' @param method Ranking method.
#' @param lambda Optimism index for \code{"liou_wang"}.
#' @param triangle Logical: for a fuzzy centrality, rank its triangular
#'   approximation (score-based methods).
#' @return A data frame with \code{actor}, \code{score} and \code{rank} (1 = most
#'   central), ordered by rank.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' head(rank_fuzzy(fc, "cog"))
#' head(rank_fuzzy(fc, "necessity"))
#' compare_rankings(fc)
#' @export
#' @section Methodology:
#' Nine methods are available. Score-based methods order all actors; methods based on comparisons use the possibility and necessity degrees and leave ties when the information does not discriminate. \code{"certified"} returns the net count of certified dominances at \eqn{\alpha = 0}, a heuristic summary: two actors whose ordering is not certified may receive different counts because of their comparisons with third actors, so the resulting ranking is not the certified partial order. The partial order itself is the relation returned by \code{\link{certified_pairs}}, in which uncertified pairs are incomparable. Agreement between methods indicates that the ranking does not depend on the ranking criterion, which Porreca, Maturo and Ventre (2025) set to the centre of gravity. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' meth <- c("cog", "graded_mean", "necessity", "certified", "top")
#' r <- lapply(meth, function(m) head(rank_fuzzy(fc, m), 5)$actor)
#' names(r) <- meth; r
#' # ranking a vector of TFNs
#' rank_fuzzy(tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32), names = c("a", "b")), "yager")
rank_fuzzy <- function(x, method = "cog", lambda = 0.5, triangle = FALSE) {
  methods <- c(fuzzy_rank_methods(), "possibility", "necessity", "certified", "top")
  method <- match.arg(method, methods)
  if (inherits(x, "ijar_index")) x <- as_tfn(x)
  if (is.tfn(x)) {
    if (!method %in% fuzzy_rank_methods()) stop("method '", method, "' requires a fuzzy_centrality object")
    s <- tfn_score(x, method, lambda); lab <- names(x); if (is.null(lab)) lab <- seq_along(s)
  } else {
    lab <- x$labels
    s <- switch(method,
                possibility = { P <- possibility_matrix(x)$possibility; rowSums(P) - colSums(P) },
                necessity = { pm <- possibility_matrix(x); N <- pm$necessity
                              N1 <- pm$necessity_lower == 1   # guaranteed N = 1 (exact)
                              rowSums(N1) - colSums(N1) + 1e-6 * (rowSums(N) - colSums(N)) },
                certified = { C <- certified_pairs(x, 0); rowSums(C) - colSums(C) },
                top = as.numeric(top_set(x)),
                if (triangle) tfn_score(x$triangle, method, lambda) else fc_score(x, method, lambda))
  }
  d <- data.frame(actor = lab, score = as.numeric(s), rank = rank(-s, ties.method = "min"),
                  stringsAsFactors = FALSE)
  d <- d[order(d$rank), ]; rownames(d) <- NULL; attr(d, "method") <- method; d
}

#' Compare rankings obtained with different methods
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param methods Ranking methods (see \code{\link{rank_fuzzy}}).
#' @param cor_method \code{"kendall"} or \code{"spearman"}.
#' @return A list with the data frame of ranks (actors by methods) and the
#'   correlation matrix between methods; correlations with a method that gives
#'   all actors the same rank are \code{NA}.
#' @examples
#' compare_rankings(fuzzy_centrality(fsn_example(), "in_closeness_min"))
#' @export
#' @section Methodology:
#' Kendall's \eqn{\tau} between the rankings of different methods measures the sensitivity of conclusions to the ranking criterion.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x2 <- fsn_simulated(n = 20, respondents = 10, seed = 2)
#' cr <- compare_rankings(fuzzy_centrality(x2, "in_strength"))
#' round(cr$correlation, 2)
#' head(cr$ranks)
compare_rankings <- function(x, methods = c("cog", "core", "yager", "graded_mean", "possibility", "necessity", "top"),
                             cor_method = c("kendall", "spearman")) {
  cor_method <- match.arg(cor_method)
  R <- vapply(methods, function(m) { r <- rank_fuzzy(x, m); r$rank[match(x$labels, r$actor)] }, numeric(length(x$labels)))
  R <- matrix(R, ncol = length(methods), dimnames = list(x$labels, methods))
  const <- apply(R, 2, function(v) length(unique(v[!is.na(v)])) < 2)
  C <- matrix(NA_real_, ncol(R), ncol(R), dimnames = list(methods, methods))
  if (any(!const)) C[!const, !const] <- stats::cor(R[, !const, drop = FALSE], method = cor_method)
  diag(C) <- 1
  list(ranks = R, correlation = C)
}

#' Proportion of discordant pairs between two rankings
#'
#' Over the pairs of actors whose reference values differ, the proportion that
#' the scores order the other way; pairs with equal scores count one half.
#' @param score Numeric vector of scores.
#' @param truth Numeric vector of reference values (e.g. latent centralities).
#' @param pairs Optional logical matrix restricting the pairs considered (upper
#'   or lower triangle is irrelevant: \code{pairs[i, j] | pairs[j, i]} is used).
#' @return A proportion, with the number of pairs as attribute.
#' @examples
#' discordance(c(3, 2, 1), c(3, 1, 2))
#' @export
#' @section Methodology:
#' The proportion of discordant pairs is \eqn{(1 - \tau)/2} in the absence of ties, where \eqn{\tau} is Kendall's coefficient. Restricting it to the pairs that are not certified shows where the errors of a crisp ranking concentrate.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- simulate_fsn(40, 0.3, kappa = 2, seed = 8)
#' fc <- fuzzy_centrality(x, "in_closeness_prod")
#' truth <- latent_centrality(x, "in_closeness_prod")
#' C <- certified_pairs(fc, 0)
#' c(all = discordance(fc_score(fc), truth),
#'   not_certified = discordance(fc_score(fc), truth, pairs = !(C | t(C))),
#'   certified = discordance(fc_score(fc), truth, pairs = C))
discordance <- function(score, truth, pairs = NULL) {
  n <- length(score); ut <- upper.tri(matrix(0, n, n))
  ds <- outer(score, score, "-"); dt <- outer(truth, truth, "-")
  keep <- ut & dt != 0
  if (!is.null(pairs)) keep <- keep & (pairs | t(pairs))
  if (!any(keep)) return(structure(NA_real_, n_pairs = 0L))
  v <- ifelse(ds[keep] == 0, 0.5, as.numeric(sign(ds[keep]) != sign(dt[keep])))
  structure(mean(v), n_pairs = sum(keep))
}
