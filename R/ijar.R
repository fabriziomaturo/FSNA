# ---------------------------------------------------------------------------
# Fuzzy centrality indices of Porreca, Maturo and Ventre (2025, IJAR 176, 109319)
# ---------------------------------------------------------------------------

#' Fuzzy Triangular Ordered Weighted Average (FTOWA)
#'
#' Orders the TFNs with a ranking score and combines them with TFN arithmetic:
#' \eqn{\mathrm{FTOWA}(A_1, \dots, A_k) = \bigoplus_h w_h \odot A_{(h)}}
#' (Porreca, Maturo and Ventre, 2025, Eq. 8).
#'
#' The same permutation is applied to the three parameters. The extension
#' principle would instead order the left endpoints, the cores and the right
#' endpoints separately; the two computations coincide when the weights are
#' equal or the order of the TFNs does not depend on the parameter, and differ
#' in general (see \code{\link{owa_extension}}).
#'
#' @param x A \code{\link{tfn}} vector.
#' @param w Weights: a numeric vector of the same length as \code{x}, or one of
#'   \code{"mean"} (equal weights \eqn{1/k}), \code{"max"} (first ordered TFN),
#'   \code{"min"} (last ordered TFN).
#' @param rank_method Ranking score used to order the TFNs, see
#'   \code{\link{tfn_score}}.
#' @param decreasing Logical: order from the largest (default, as in Section 4.1
#'   of the paper) or from the smallest (Section 3.3).
#' @param lambda Optimism index for \code{rank_method = "liou_wang"}.
#' @return A TFN of length one.
#' @examples
#' a <- tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32))
#' ftowa(a, "max")                 # chosen by the CoG: (0, 0.2, 1)
#' owa_extension(a, c(1, 0))       # extension principle: (0.3, 0.31, 1)
#' ftowa(a, "max", rank_method = "core")
#' @export
#' @section Methodology:
#' The FTOWA (Porreca, Maturo and Ventre, 2025, Eq. 8) first orders the TFNs by a ranking score and then applies the ordered weights to the three parameters with the same permutation. The extension of the crisp OWA (Yager, 1988) instead orders the left endpoints, the cores and the right endpoints separately. The two coincide when the weights are equal (the order is irrelevant) or when the ordering by score agrees with the orderings of all parameters, and may differ otherwise, even at the support and at the core. With weights \eqn{(1, 0, \dots, 0)} the FTOWA selects the TFN with the largest score, whereas the extension returns the componentwise maximum. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @examples
#' a <- tfn(c(0.1, 0.3, 0.0), c(0.5, 0.45, 0.2), c(0.9, 0.5, 0.6))
#' w <- c(0.5, 0.3, 0.2)
#' ftowa(a, w)                       # ordered by CoG
#' ftowa(a, w, rank_method = "core") # ordered by core
#' owa_extension(a, w)               # extension principle
#' ftowa(a, "mean"); tfn_sum(a, 1/3) # equal weights: identical
ftowa <- function(x, w = "mean", rank_method = "cog", decreasing = TRUE, lambda = 0.5) {
  k <- length(x)
  w <- .owa_weights(w, k)
  o <- order(tfn_score(x, rank_method, lambda), decreasing = decreasing)
  tfn_sum(x[o], w)
}

#' OWA of TFNs by the extension principle
#'
#' Ordered weighted average computed separately on the \eqn{\alpha}-cut
#' endpoints, which is the extension of the crisp OWA (a continuous and
#' non-decreasing function of its arguments). For TFN arguments the result at
#' the support and at the core is obtained by ordering the left endpoints, the
#' cores and the right endpoints separately; at intermediate levels the result
#' is piecewise linear and is returned by \code{\link{fuzzy_centrality}} with an
#' OWA index.
#' @param x A \code{\link{tfn}} vector.
#' @param w Weights (applied to the values in decreasing order).
#' @return A TFN of length one: the triangle (support and core) of the extension.
#' @examples
#' owa_extension(tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32)), c(0.7, 0.3))
#' @export
#' @section Methodology:
#' For TFN arguments, the crisp OWA is a continuous non-decreasing function, so the endpoints of the \eqn{\alpha}-cuts of the extension are the OWA of the endpoints of the \eqn{\alpha}-cuts. Between the support and the core the endpoints are piecewise linear in \eqn{\alpha}, because the order of the arguments may change with the level.
#' @examples
#' a <- tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32))
#' owa_extension(a, c(1, 0))   # componentwise maximum
#' owa_extension(a, c(0, 1))   # componentwise minimum
owa_extension <- function(x, w) {
  y <- unclass(x); w <- .owa_weights(w, length(y$l))
  tfn(sum(w * sort(y$l, decreasing = TRUE)), sum(w * sort(y$m, decreasing = TRUE)),
      sum(w * sort(y$r, decreasing = TRUE)))
}

.owa_weights <- function(w, k) {
  if (is.character(w)) {
    w <- match.arg(w, c("mean", "max", "min"))
    w <- switch(w, mean = rep(1 / k, k), max = c(1, rep(0, k - 1)), min = c(rep(0, k - 1), 1))
  }
  if (length(w) != k) stop("the weights must have the same length as the TFNs (", k, ")")
  if (any(w < 0)) stop("weights must be non-negative")
  w
}

.ijar_variant <- function(variant, ...) {
  variant <- match.arg(variant, c("definition", "alternative"))
  def <- if (variant == "alternative")
    list(path_tfn = "mean", normalize_closeness = "rescaled", tie_break = "indirect",
         direct_rule = "transposed", normalize_betweenness = "relative")
  else
    list(path_tfn = "weakest", normalize_closeness = "mean", tie_break = "shorter",
         direct_rule = "same_pair", normalize_betweenness = "relative")
  a <- list(...)
  for (nm in names(a)) if (!is.null(a[[nm]])) def[[nm]] <- a[[nm]]
  def
}

# object returned by the IJAR indices
.ijar_result <- function(x, lab, index, rank_method, lambda, extra = list()) {
  y <- unclass(x)
  d <- data.frame(actor = lab, l = y$l, m = y$m, r = y$r, stringsAsFactors = FALSE)
  d$score <- tfn_score(tfn(d$l, d$m, d$r), rank_method, lambda)
  d$rank <- rank(-d$score, ties.method = "min")
  attr(d, "index") <- index; attr(d, "rank_method") <- rank_method
  for (nm in names(extra)) attr(d, nm) <- extra[[nm]]
  class(d) <- c("ijar_index", "data.frame"); d
}

#' @export
print.ijar_index <- function(x, n = 10, digits = 4, ...) {
  cat("Fuzzy centrality (Porreca, Maturo and Ventre, 2025):", attr(x, "index"), "\n")
  if (!is.null(attr(x, "settings"))) cat("Settings:", attr(x, "settings"), "\n")
  cat("Ranking by", attr(x, "rank_method"), "; first", min(n, nrow(x)), "of", nrow(x), "actors\n")
  d <- as.data.frame(unclass(x))[order(x$rank), ]; rownames(d) <- NULL
  print(utils::head(d, n), digits = digits)
  invisible(x)
}

#' @rdname ijar_degree
#' @param object An \code{"ijar_index"} object.
#' @export
as_tfn <- function(object) {
  tfn(object$l, object$m, object$r, names = object$actor)
}

#' Fuzzy degree centrality (FTOWA of the ties)
#'
#' Fuzzy in-degree (Eq. 16), out-degree (Eq. 19) and total degree (Eq. 21) of
#' Porreca, Maturo and Ventre (2025): the FTOWA of the \eqn{n - 1} incoming or
#' outgoing ties of each actor, null ties included, ordered by
#' \code{rank_method}.
#'
#' With equal weights the FTOWA is the mean of the ties and does not depend on
#' their order: it is the in- or out-strength divided by \eqn{n - 1}, and it
#' coincides with the extension principle, so that the result is exact.
#' \code{weights = "mean"} gives equal weights \eqn{1/(n-1)}, so that the index
#' is the mean of the ties; \code{weights = "n"} gives equal weights \eqn{1/n}.
#' With unequal weights the TFNs are ordered
#' by their score and the result differs in general from the extension; see
#' \code{\link{fuzzy_centrality}} with \code{index = "in_owa"} for the exact
#' computation.
#'
#' @param x An \code{\link{fsn}} object.
#' @param mode \code{"in"}, \code{"out"} or \code{"total"}.
#' @param weights \code{"mean"} (\eqn{1/(n-1)}), \code{"n"} (\eqn{1/n}),
#'   \code{"max"}, \code{"min"} or a numeric vector of length \eqn{n - 1}.
#' @param rank_method Score used to order the ties and to rank the actors, see
#'   \code{\link{tfn_score}}.
#' @param lambda Optimism index for \code{"liou_wang"}.
#' @return An object of class \code{"ijar_index"}: a data frame with the actor,
#'   the TFN \code{(l, m, r)}, its score and the rank. \code{as_tfn()} converts it
#'   to a \code{\link{tfn}} vector.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures
#'   in social network analysis: Theory and application in a university
#'   department collaboration network. \emph{International Journal of Approximate
#'   Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @seealso \code{\link{ijar_closeness}}, \code{\link{ijar_betweenness}},
#'   \code{\link{fuzzy_centrality}}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' ijar_degree(A, "in")
#' print(ijar_degree(A, "out"), n = 5)
#' # ordered weighted average with decreasing weights
#' w <- rev(seq_len(n_actors(A) - 1)); w <- w / sum(w)
#' ijar_degree(A, "in", weights = w)
#' @export
#' @section Methodology:
#' Definition (Porreca, Maturo and Ventre, 2025, Eqs. 16 and 19): the fuzzy in-degree of actor \eqn{i} is \eqn{\bigoplus_j w_j \odot \tilde O_j}, where \eqn{\tilde O_j} are the \eqn{n-1} incoming ties ordered decreasingly by their score. With equal weights the index is the in-strength divided by \eqn{n-1}; its \eqn{\alpha}-cuts are exact and coincide with \code{fuzzy_centrality(x, "in_strength", normalize = "n-1")}. With unequal weights it differs from the extension principle, because the order of the ties is fixed by their score.
#'
#' Equal weights \eqn{1/(n-1)} (\code{"mean"}) and \eqn{1/n} (\code{"n"}) give indices that differ by the constant factor \eqn{(n-1)/n} and the same ranking. Out-degree uses the outgoing ties and measures the activity of an actor as a rater, which in questionnaire data reflects response behaviour. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' d <- ijar_degree(A, "in")
#' plot(d, n = 5)
#' # exact counterpart
#' fc <- fuzzy_centrality(A, "in_strength", normalize = "n-1")
#' all.equal(unname(fc$lower[, 1]), d$l)
#' # ranking by another score
#' dg <- as.data.frame(unclass(ijar_degree(A, "in", rank_method = "graded_mean")))
#' head(dg[order(dg$rank), ], 3)
#' ijar_degree(A, "total")
ijar_degree <- function(x, mode = c("in", "out", "total"), weights = "mean", rank_method = "cog",
                        lambda = 0.5) {
  mode <- match.arg(mode); n <- n_actors(x)
  if (mode == "total") {
    a <- ijar_degree(x, "in", weights, rank_method, lambda); b <- ijar_degree(x, "out", weights, rank_method, lambda)
    tot <- as_tfn(a) + as_tfn(b)
    return(.ijar_result(tot, x$labels, "fuzzy total degree", rank_method, lambda))
  }
  if (is.character(weights) && weights[1] == "n") weights <- rep(1 / n, n - 1)
  if (is.character(weights) && weights[1] == "mean") weights <- rep(1 / (n - 1), n - 1)
  res <- lapply(seq_len(n), function(i) {
    idx <- setdiff(seq_len(n), i)
    t <- if (mode == "in") tfn(x$L[idx, i], x$M[idx, i], x$R[idx, i]) else tfn(x$L[i, idx], x$M[i, idx], x$R[i, idx])
    ftowa(t, weights, rank_method, TRUE, lambda)
  })
  .ijar_result(do.call(c, res), x$labels, paste0("fuzzy ", mode, "-degree"), rank_method, lambda)
}

#' Best paths between actors (Porreca, Maturo and Ventre, 2025, Def. 11)
#'
#' For every ordered pair of actors, the best path with at most
#' \code{max_steps} steps is the path whose weakest tie, as ranked by
#' \code{rank_method}, is the strongest (the minimum is selected with
#' \eqn{w = (0, \dots, 0, 1)} in Eq. 24). The fuzzy strength assigned to the
#' path depends on \code{path_tfn}:
#' \describe{
#'   \item{\code{"weakest"}}{the TFN of the weakest tie, as stated in the paper
#'     (FTOWA with \eqn{w = (0,\dots,0,1)}, Definition 5);}
#'   \item{\code{"mean"}}{the componentwise mean of the ties of the path;}
#'   \item{\code{"extension"}}{the componentwise minimum of the ties, which is the
#'     extension of the minimum at the support and at the core for the selected
#'     path (the selection itself is still made by \code{rank_method}).}
#' }
#' A best path with null strength (no path with positive score) gives a null
#' TFN.
#'
#' Paths are found by dynamic programming over the number of steps. With
#' \code{tie_break = "shorter"} (default) a longer path replaces a shorter one
#' only when its score is strictly greater, so that among paths with equal score
#' the shorter is kept; with \code{"indirect"} an indirect path replaces a
#' shorter one when its score is greater or equal, and among intermediate nodes
#' with equal score the first is taken. \code{"indirect"} is defined for at most
#' two steps.
#'
#' @param x An \code{\link{fsn}} object.
#' @param max_steps Maximum number of steps \eqn{s}. \code{Inf} means \eqn{n - 1}.
#' @param rank_method Score used to rank the ties, see \code{\link{tfn_score}}.
#' @param path_tfn \code{"weakest"}, \code{"mean"} or \code{"extension"}.
#' @param tie_break \code{"shorter"} or \code{"indirect"}.
#' @param lambda Optimism index for \code{"liou_wang"}.
#' @return A list with the matrices \code{L}, \code{M}, \code{R} of the path
#'   strengths (rows: origin, columns: destination), \code{score} (score of the
#'   weakest tie of the best path), \code{length} (number of steps) and the
#'   information needed by \code{\link{best_path_route}}.
#' @examples
#' x <- fsn_example()
#' bp <- ijar_best_paths(x, max_steps = 3)
#' round(bp$score, 2)
#' best_path_route(bp, "A", "D")
#' @export
#' @section Methodology:
#' In a fuzzy network every pair is connected, and the shortest path is replaced by the \emph{best} path: the one whose weakest tie, as ranked by a score, is strongest (Definition 11 of Porreca, Maturo and Ventre, 2025). The selection is a widest (bottleneck) path problem on the scores of the ties and is solved by dynamic programming over the number of steps.
#'
#' Two features distinguish this computation from the extension principle. First, the path is selected by a score; under the extension the best path may change across the realizations of the ties, and the strength of the best path is \eqn{\max_\pi \min_{(k,l) \in \pi} W_{kl}} computed separately at each endpoint. Second, the fuzzy strength assigned to a path is the TFN of its weakest tie, as in the definition of the paper, or, with \code{path_tfn = "mean"}, the componentwise mean of its ties, rather than the componentwise minimum. For crisp ties the selected path is the classical maximin path; its strength is the bottleneck value for \code{"weakest"} and \code{"extension"}, but the mean of the ties for \code{"mean"} (e.g. 0.75 rather than 0.6 for ties 0.9 and 0.6). The exact counterpart is the closeness with the minimum t-norm in \code{\link{fuzzy_centrality}}. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' x <- fsn_example()
#' p_w <- ijar_best_paths(x, 3, path_tfn = "weakest")
#' p_e <- ijar_best_paths(x, 3, path_tfn = "extension")
#' best_path_route(p_w, "A", "D")
#' rbind(weakest = c(p_w$L["A", "D"], p_w$M["A", "D"], p_w$R["A", "D"]),
#'       extension = c(p_e$L["A", "D"], p_e$M["A", "D"], p_e$R["A", "D"]))
#' # number of steps of the best paths
#' table(p_w$length[row(p_w$length) != col(p_w$length)])
ijar_best_paths <- function(x, max_steps = 2, rank_method = "cog", path_tfn = c("weakest", "mean", "extension"),
                            tie_break = c("shorter", "indirect"), lambda = 0.5) {
  path_tfn <- match.arg(path_tfn); tie_break <- match.arg(tie_break)
  n <- n_actors(x); s <- if (is.infinite(max_steps)) n - 1 else as.integer(max_steps)
  if (!is.numeric(max_steps) || length(max_steps) != 1 || is.na(max_steps) || max_steps < 1 ||
      (is.finite(max_steps) && max_steps != round(max_steps)))
    stop("max_steps must be a positive integer or Inf")
  if (tie_break == "indirect" && s > 2)
    stop("tie_break = 'indirect' is defined for max_steps <= 2; ",
         "with more steps it could select walks with repeated actors")
  sc <- .score_matrix(x$L, x$M, x$R, rank_method, lambda); diag(sc) <- 0
  L <- x$L; M <- x$M; R <- x$R
  best <- sc; diag(best) <- -Inf
  len <- matrix(1L, n, n)
  # aggregates of the current best path
  if (path_tfn == "mean") { aL <- L; aM <- M; aR <- R; cnt <- matrix(1, n, n) }
  else { aL <- L; aM <- M; aR <- R }
  wS <- sc                                   # score of the weakest tie (for "weakest")
  pred <- vector("list", s); changed <- vector("list", s)
  ind2 <- NULL
  for (k in seq_len(s)[-1]) {
    nb <- best; nlen <- len; nL <- aL; nM <- aM; nR <- aR; nS <- wS
    if (path_tfn == "mean") ncnt <- cnt
    pk <- matrix(NA_integer_, n, n); ch <- matrix(FALSE, n, n)
    if (k == 2) ind2 <- list(value = matrix(0, n, n), via = matrix(NA_integer_, n, n))
    for (j in seq_len(n)) {
      cand <- pmin(best, matrix(sc[, j], n, n, byrow = TRUE))   # cand[i, h]: i -> ... -> h -> j
      cand[, j] <- -Inf
      cand[!is.finite(cand)] <- -Inf
      h <- max.col(cand, ties.method = "first")
      v <- cand[cbind(seq_len(n), h)]
      if (k == 2) { ind2$value[, j] <- pmax(v, 0); ind2$via[, j] <- h }
      cur <- best[, j]
      upd <- if (tie_break == "indirect") v >= cur else v > cur
      upd[j] <- FALSE; upd[!is.finite(v)] <- FALSE
      if (!any(upd)) next
      i <- which(upd); hh <- h[upd]
      nb[i, j] <- v[upd]; nlen[i, j] <- len[cbind(i, hh)] + 1L; pk[i, j] <- hh; ch[i, j] <- TRUE
      if (path_tfn == "mean") {
        nL[i, j] <- aL[cbind(i, hh)] + L[hh, j]; nM[i, j] <- aM[cbind(i, hh)] + M[hh, j]
        nR[i, j] <- aR[cbind(i, hh)] + R[hh, j]; ncnt[i, j] <- cnt[cbind(i, hh)] + 1
      } else if (path_tfn == "extension") {
        nL[i, j] <- pmin(aL[cbind(i, hh)], L[hh, j]); nM[i, j] <- pmin(aM[cbind(i, hh)], M[hh, j])
        nR[i, j] <- pmin(aR[cbind(i, hh)], R[hh, j])
      } else {
        prevS <- wS[cbind(i, hh)]; tie_weaker <- sc[hh, j] < prevS
        nL[i, j] <- ifelse(tie_weaker, L[hh, j], aL[cbind(i, hh)])
        nM[i, j] <- ifelse(tie_weaker, M[hh, j], aM[cbind(i, hh)])
        nR[i, j] <- ifelse(tie_weaker, R[hh, j], aR[cbind(i, hh)])
        nS[i, j] <- pmin(prevS, sc[hh, j])
      }
    }
    best <- nb; len <- nlen; aL <- nL; aM <- nM; aR <- nR; wS <- nS
    if (path_tfn == "mean") cnt <- ncnt
    pred[[k]] <- pk; changed[[k]] <- ch
  }
  if (path_tfn == "mean") { aL <- aL / cnt; aM <- aM / cnt; aR <- aR / cnt }
  null <- !(best > 0); diag(null) <- TRUE
  aL[null] <- 0; aM[null] <- 0; aR[null] <- 0
  best[null] <- 0; len[null] <- 0L
  dimnames(aL) <- dimnames(aM) <- dimnames(aR) <- dimnames(best) <- dimnames(len) <- list(x$labels, x$labels)
  structure(list(L = aL, M = aM, R = aR, score = best, length = len, pred = pred, changed = changed,
                 indirect2 = ind2, direct_score = sc, labels = x$labels, max_steps = s,
                 settings = list(rank_method = rank_method, path_tfn = path_tfn, tie_break = tie_break)),
            class = "ijar_paths")
}

#' Route of a best path
#'
#' @param paths An object returned by \code{\link{ijar_best_paths}}.
#' @param from,to Actor labels or indices.
#' @return Character vector with the actors on the best path (empty if the best
#'   path has null strength).
#' @examples
#' best_path_route(ijar_best_paths(fsn_example(), 3), "A", "D")
#' @export
#' @section Methodology:
#' Reconstructs the sequence of actors of a best path from the dynamic programming tables of \code{\link{ijar_best_paths}}. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @examples
#' bp <- ijar_best_paths(fsn_simulated(n = 20, respondents = 10), 2)
#' best_path_route(bp, "a01", "a02")
best_path_route <- function(paths, from, to) {
  lab <- paths$labels
  i <- if (is.character(from)) match(from, lab) else from
  j <- if (is.character(to)) match(to, lab) else to
  if (paths$score[i, j] <= 0) return(character(0))
  lab[.route(paths, i, j, paths$max_steps)]
}

.route <- function(paths, i, j, k) {
  while (k > 1 && !paths$changed[[k]][i, j]) k <- k - 1
  if (k == 1) return(c(i, j))
  h <- paths$pred[[k]][i, j]
  c(.route(paths, i, h, k - 1), j)
}

#' Fuzzy closeness centrality (Porreca, Maturo and Ventre, 2025)
#'
#' The fuzzy in-closeness of actor \eqn{i} is the FTOWA (with equal weights, the
#' mean) of the strengths of the best paths from the other actors to \eqn{i};
#' the out-closeness uses the best paths from \eqn{i}; the total closeness is the
#' sum of the two. Best paths are computed by \code{\link{ijar_best_paths}}.
#'
#' Two sets of computational choices are available. \code{variant = "definition"}
#' (default) follows the definitions of Section 3 of the paper: path strength
#' equal to the weakest tie, shorter paths kept when scores are equal, and plain
#' mean over the other actors. \code{variant = "alternative"} is defined for at
#' most two steps and uses the componentwise mean of the ties of a best path as
#' its strength, prefers indirect paths when scores are equal
#' (\code{tie_break = "indirect"}) and divides the mean once more by \eqn{n - 1}
#' (\code{normalize = "rescaled"}); the rescaling does not change the ranking.
#' Individual settings override the variant, so that the effect of each choice
#' can be examined separately.
#'
#' @param x An \code{\link{fsn}} object.
#' @param mode \code{"in"}, \code{"out"} or \code{"total"}.
#' @param max_steps Maximum number of steps of the best paths.
#' @param variant \code{"definition"} or \code{"alternative"}.
#' @param rank_method Score used to rank the ties, the paths and the actors.
#' @param path_tfn,tie_break Override the variant, see \code{\link{ijar_best_paths}}.
#' @param normalize \code{"mean"} or \code{"rescaled"} (mean divided by \eqn{n - 1}); overrides the variant.
#' @param lambda Optimism index for \code{"liou_wang"}.
#' @param paths Optional precomputed \code{\link{ijar_best_paths}} object.
#' @return An \code{"ijar_index"} data frame.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' print(ijar_closeness(A, "in"), n = 5)
#' print(ijar_closeness(A, "out"), n = 5)
#' # alternative computational choices
#' print(ijar_closeness(A, "in", variant = "alternative"), n = 5)
#' @export
#' @section Methodology:
#' Definition (Porreca, Maturo and Ventre, 2025, Eqs. 27-31): the fuzzy closeness of an actor is the FTOWA, with equal weights the mean, of the strengths of its best paths to (out) or from (in) the other actors. It measures how strongly an actor is connected to the others through their best routes, not a distance; higher values mean higher centrality.
#'
#' The index inherits the two departures from the extension principle described in \code{\link{ijar_best_paths}}, so its support and core may differ from those of the exact fuzzy closeness. \code{variant = "alternative"} changes the path strength, the treatment of equal scores and the normalization, as described above. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' alt <- ijar_closeness(A, "in", variant = "alternative")
#' def <- ijar_closeness(A, "in", variant = "definition")
#' cor(alt$score, def$score, method = "kendall")
#' # one setting at a time: mean path strength, plain mean
#' ijar_closeness(A, "in", path_tfn = "mean", normalize = "mean")
#' # exact counterpart (minimum t-norm, two steps)
#' fc <- fuzzy_centrality(A, "in_closeness_min", max_steps = 2)
#' plot(fc, actors = 5)
ijar_closeness <- function(x, mode = c("in", "out", "total"), max_steps = 2, variant = c("definition", "alternative"),
                           rank_method = "cog", path_tfn = NULL, tie_break = NULL, normalize = NULL,
                           lambda = 0.5, paths = NULL) {
  mode <- match.arg(mode); variant <- match.arg(variant)
  st <- .ijar_variant(variant, path_tfn = path_tfn, tie_break = tie_break, normalize_closeness = normalize)
  n <- n_actors(x)
  if (is.null(paths)) paths <- ijar_best_paths(x, max_steps, rank_method, st$path_tfn, st$tie_break, lambda)
  if (mode == "total") {
    a <- ijar_closeness(x, "in", max_steps, variant, rank_method, path_tfn, tie_break, normalize, lambda, paths)
    b <- ijar_closeness(x, "out", max_steps, variant, rank_method, path_tfn, tie_break, normalize, lambda, paths)
    return(.ijar_result(as_tfn(a) + as_tfn(b), x$labels, "fuzzy total closeness", rank_method, lambda,
                        list(settings = attr(a, "settings"))))
  }
  if (!st$normalize_closeness %in% c("mean", "rescaled")) stop("normalize must be 'mean' or 'rescaled'")
  use_in <- mode == "in"
  f <- function(A, i) if (use_in) mean(A[-i, i]) else mean(A[i, -i])
  den <- if (st$normalize_closeness == "rescaled") n - 1 else 1
  res <- tfn(vapply(seq_len(n), function(i) f(paths$L, i), 0) / den,
             vapply(seq_len(n), function(i) f(paths$M, i), 0) / den,
             vapply(seq_len(n), function(i) f(paths$R, i), 0) / den)
  set <- sprintf("s = %d, path strength = %s, tie break = %s, normalization = %s", paths$max_steps,
                 paths$settings$path_tfn, paths$settings$tie_break, st$normalize_closeness)
  .ijar_result(res, x$labels, paste0("fuzzy ", mode, "-closeness"), rank_method, lambda, list(settings = set))
}

#' Betweenness on best paths (Porreca, Maturo and Ventre, 2025, Eq. 25)
#'
#' Counts how often each actor is an intermediate node of the best path between
#' two other actors (\code{\link{ijar_best_paths}}). The index is crisp: it is a
#' relative frequency computed on paths selected by the ranking score, and it
#' depends on that criterion. Betweenness is not monotone in the tie weights, so
#' the exact methods of \code{\link{fuzzy_centrality}} do not apply to it.
#'
#' \code{normalize = "relative"} divides the counts by their total, so that the
#' values sum to one; \code{"pairs"} divides by the
#' number \eqn{(n-1)(n-2)} of ordered pairs not involving the actor (Eq. 23);
#' \code{"none"} returns the counts.
#'
#' \strong{Ties between best paths.} Several paths can share the same best
#' score and number of steps. With \code{ties = "split"} (default of
#' \code{variant = "definition"}) every pair contributes \eqn{\sigma_{st}(i)/\sigma_{st}},
#' the share of its optimal paths through \eqn{i} (Freeman, 1977), where the
#' optimal paths are the paths with the largest score of the weakest tie and,
#' among them, the smallest number of steps (at most \code{max_steps}). They are
#' exactly the shortest paths from \eqn{s} to \eqn{t} in the graph of the ties
#' whose score is at least the best score of the pair, and they are counted by
#' breadth-first search; the result does not depend on the order of the actors.
#' With \code{ties = "first"} a single representative path is used; the result
#' then depends on the labelling of the actors when optimal paths tie.
#'
#' \code{direct_rule} decides when the best two-step path of a pair is preferred
#' to its direct tie. With \code{"same_pair"} (default) the intermediate node of
#' the best two-step path from \eqn{i} to \eqn{j} is counted unless the direct tie
#' from \eqn{i} to \eqn{j} has a higher score, as in the definition of the paper.
#' With \code{"transposed"} (two steps only) the comparison is made with the
#' direct tie from \eqn{j} to \eqn{i}, which makes the count depend on the
#' strength of the ties in both directions; it is provided as a sensitivity
#' option for directed networks.
#'
#' @inheritParams ijar_closeness
#' @param normalize \code{"relative"}, \code{"pairs"} or \code{"none"}.
#' @param direct_rule \code{"same_pair"} or \code{"transposed"}; overrides the variant.
#' @param tie_break \code{"shorter"} or \code{"indirect"}; overrides the variant, see \code{\link{ijar_best_paths}}.
#' @param ties \code{"split"} (share each pair among all its optimal paths) or
#'   \code{"first"} (one representative path). Default: \code{"split"} for
#'   \code{variant = "definition"}, \code{"first"} for \code{"alternative"}.
#' @return A data frame with the actor, the (fractional) number of best paths through it, the
#'   betweenness and its rank, of class \code{"ijar_betweenness"}.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' print(ijar_betweenness(A), n = 5)
#' print(ijar_betweenness(A, variant = "alternative"), n = 5)
#' @export
#' @section Methodology:
#' Betweenness on best paths (Porreca, Maturo and Ventre, 2025, Eq. 25) counts how often an actor is an intermediate node of the best path between two other actors. It is a crisp number, because the paths are selected by a ranking score; a different score may select different paths and change the index. Betweenness is not monotone in the tie weights (strengthening a tie can divert best paths away from an actor), so it has no exact fuzzy counterpart in the class of isotone indices. When several best paths tie, the default divides the contribution of a pair among them (Freeman, 1977), so that the index does not depend on the labelling of the actors.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @references Freeman L.C. (1977). A set of measures of centrality based on betweenness. \emph{Sociometry} 40, 35-41. \doi{10.2307/3033543}
#' @examples
#' # two equally good routes from A to D: the pair is shared by B and C
#' m <- matrix(0, 4, 4); m[1, 2:3] <- 0.8; m[2:3, 4] <- 0.8
#' dm <- fsn(M = m, labels = c("A", "B", "C", "D"))
#' ijar_betweenness(dm, normalize = "none")
#' ijar_betweenness(dm, normalize = "none", ties = "first")   # one representative
#' A <- fsn_simulated(n = 20, respondents = 10)
#' b_def <- ijar_betweenness(A)
#' b_tr <- ijar_betweenness(A, direct_rule = "transposed", tie_break = "shorter", ties = "first")
#' cor(b_def$betweenness, b_tr$betweenness, method = "spearman")
#' plot(b_def, n = 20)
#' # sensitivity to the ranking score
#' b_core <- ijar_betweenness(A, rank_method = "core")
#' cor(b_def$betweenness, b_core$betweenness, method = "spearman")
ijar_betweenness <- function(x, max_steps = 2, variant = c("definition", "alternative"), rank_method = "cog",
                             normalize = c("relative", "pairs", "none"), direct_rule = NULL,
                             tie_break = NULL, ties = NULL, lambda = 0.5, paths = NULL) {
  variant <- match.arg(variant); normalize <- match.arg(normalize)
  st <- .ijar_variant(variant, direct_rule = direct_rule, tie_break = tie_break)
  if (is.null(ties)) ties <- if (st$tie_break == "shorter" && st$direct_rule == "same_pair") "split" else "first"
  ties <- match.arg(ties, c("split", "first"))
  if (!st$direct_rule %in% c("same_pair", "transposed")) stop("direct_rule must be 'same_pair' or 'transposed'")
  if (ties == "split" && (st$tie_break != "shorter" || st$direct_rule != "same_pair"))
    stop("ties = 'split' requires tie_break = 'shorter' and direct_rule = 'same_pair'")
  n <- n_actors(x)
  if (is.null(paths)) paths <- ijar_best_paths(x, max_steps, rank_method, "weakest", st$tie_break, lambda)
  counts <- numeric(n)
  if (ties == "split") {
    sc <- paths$direct_score; diag(sc) <- 0; s_max <- paths$max_steps
    for (i in seq_len(n)) for (j in seq_len(n)) {
      v <- paths$score[i, j]
      if (i == j || v <= 0) next
      A <- sc >= v; diag(A) <- FALSE
      fw <- .bfs_count(A, i, s_max); bw <- .bfs_count(t(A), j, s_max)
      d <- fw$dist[j]
      if (!is.finite(d) || d < 2) next
      on <- is.finite(fw$dist) & is.finite(bw$dist) & fw$dist + bw$dist == d
      on[c(i, j)] <- FALSE
      if (any(on)) counts[on] <- counts[on] + fw$sigma[on] * bw$sigma[on] / fw$sigma[j]
    }
  } else if (st$direct_rule == "transposed") {
    if (paths$max_steps != 2) stop("direct_rule = 'transposed' is defined for max_steps = 2 only")
    v <- paths$indirect2$value; z <- paths$indirect2$via
    direct_better <- paths$direct_score > v
    excl <- t(direct_better)
    ok <- v > 0 & !excl; diag(ok) <- FALSE
    counts <- tabulate(z[ok], n)
  } else {
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i == j || paths$score[i, j] <= 0 || paths$length[i, j] < 2) next
      r <- .route(paths, i, j, paths$max_steps)
      mid <- unique(r[-c(1, length(r))]); mid <- setdiff(mid, c(i, j))
      counts[mid] <- counts[mid] + 1
    }
  }
  b <- switch(normalize, relative = if (sum(counts)) counts / sum(counts) else counts,
              pairs = if (n > 2) counts / ((n - 1) * (n - 2)) else rep(NA_real_, n), none = counts)
  d <- data.frame(actor = x$labels, paths_through = counts, betweenness = b,
                  rank = rank(-b, ties.method = "min", na.last = "keep"), stringsAsFactors = FALSE)
  attr(d, "settings") <- sprintf("s = %d, direct rule = %s, ties = %s, normalization = %s", paths$max_steps,
                                 st$direct_rule, ties, normalize)
  class(d) <- c("ijar_betweenness", "data.frame"); d
}

# breadth-first search with path counting on a logical adjacency matrix,
# up to max_d steps: distances (Inf if unreachable) and numbers of shortest paths
.bfs_count <- function(A, s, max_d) {
  n <- nrow(A); dist <- rep(Inf, n); sigma <- numeric(n)
  dist[s] <- 0; sigma[s] <- 1; front <- s; d <- 0
  while (length(front) && d < max_d) {
    d <- d + 1
    reach <- colSums(A[front, , drop = FALSE] * sigma[front])
    new <- which(reach > 0 & is.infinite(dist))
    if (!length(new)) break
    dist[new] <- d; sigma[new] <- reach[new]; front <- new
  }
  list(dist = dist, sigma = sigma)
}

#' @export
print.ijar_betweenness <- function(x, n = 10, digits = 4, ...) {
  cat("Betweenness on best paths (crisp):", attr(x, "settings"), "\n")
  d <- as.data.frame(unclass(x))[order(x$rank), ]; rownames(d) <- NULL
  print(utils::head(d, n), digits = digits); invisible(x)
}

#' Strength of a path from its fuzzy ties
#'
#' Compares the strength of a path computed as in Porreca, Maturo and Ventre
#' (2025), i.e. the TFN of the weakest tie identified by a ranking score, with
#' the extension of the minimum, obtained separately on the left endpoints, the
#' cores and the right endpoints, and with the componentwise mean.
#' @param ties A \code{\link{tfn}} vector with the ties of the path.
#' @param rule \code{"weakest"}, \code{"extension"} or \code{"mean"}.
#' @param rank_method Ranking score for \code{"weakest"}.
#' @return A TFN of length one.
#' @examples
#' p <- tfn(c(0.1, 0.3), c(0.5, 0.45), c(0.9, 0.5))
#' path_strength(p, "weakest")    # (0.3, 0.45, 0.5)
#' path_strength(p, "extension")  # (0.1, 0.45, 0.5)
#' @export
#' @section Methodology:
#' For a path with ties \eqn{\tilde W_1, \dots, \tilde W_k}, the extension of the minimum has \eqn{\alpha}-cuts \eqn{[\min_h \underline w_{h,\alpha}, \min_h \overline w_{h,\alpha}]}; for TFNs its support and core are the componentwise minimum. Selecting the weakest tie by a score returns one of the ties, which overstates the left endpoint whenever another tie has a lower left endpoint but a higher score. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @examples
#' p <- tfn(c(0.1, 0.3), c(0.5, 0.45), c(0.9, 0.5))
#' rbind(as.matrix(path_strength(p, "weakest")), as.matrix(path_strength(p, "extension")),
#'       as.matrix(path_strength(p, "mean")))
path_strength <- function(ties, rule = c("weakest", "extension", "mean"), rank_method = "cog") {
  rule <- match.arg(rule); y <- unclass(ties)
  switch(rule,
         weakest = ties[which.min(tfn_score(ties, rank_method))],
         extension = tfn(min(y$l), min(y$m), min(y$r)),
         mean = tfn(mean(y$l), mean(y$m), mean(y$r)))
}
