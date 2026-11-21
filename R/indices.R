# ---------------------------------------------------------------------------
# Crisp centrality indices that are continuous and non-decreasing in the ties
# ---------------------------------------------------------------------------

#' Crisp centrality indices for weighted networks
#'
#' Centrality indices of a crisp weighted directed network \eqn{W} (rows:
#' senders, columns: receivers) that are continuous and non-decreasing in every
#' tie weight (isotone). They are the building blocks of the exact fuzzy
#' centralities of \code{\link{fuzzy_centrality}}.
#'
#' \describe{
#'   \item{\code{centrality_strength}}{in-strength \eqn{\sum_k W_{ki}},
#'     out-strength \eqn{\sum_k W_{ik}} or their sum.}
#'   \item{\code{centrality_owa}}{ordered weighted average (Yager, 1988) of the
#'     \eqn{n-1} incoming or outgoing ties, ordered decreasingly, with fixed
#'     non-negative weights; the crisp counterpart of the FTOWA degree.}
#'   \item{\code{centrality_closeness}}{mean over the other actors of the
#'     strength of the best path toward (in) or from (out) the actor, where the
#'     strength of a path is the t-norm \eqn{T} of its ties and the best path
#'     maximizes it, with at most \code{max_steps} steps. With \eqn{T} the
#'     minimum this is the crisp counterpart of the best-path closeness of
#'     Porreca, Maturo and Ventre (2025).}
#'   \item{\code{centrality_katz}}{Katz centrality
#'     \eqn{\sum_{k \ge 1} \beta^k (W^\top)^k \mathbf 1} for in-centrality
#'     (walks ending at the actor) and \eqn{\sum_{k\ge1}\beta^k W^k \mathbf 1}
#'     for out-centrality; isotone on the box of networks with spectral radius
#'     below \eqn{1/\beta}.}
#' }
#' Betweenness and eigenvector centrality are not monotone in the ties and are
#' not included.
#'
#' @param W Non-negative square matrix.
#' @param mode \code{"in"}, \code{"out"} or (strength only) \code{"total"}.
#' @param normalize For strength: \code{"none"}, \code{"n-1"} (divide by
#'   \eqn{n - 1}) or \code{"n"}; for Katz:
#'   \code{"none"} or \code{"n-1"}. The crisp functions default to
#'   \code{"none"}; \code{\link{centrality_index}}, hence
#'   \code{\link{fuzzy_centrality}}, defaults to \code{"n-1"}, as in
#'   Proposition 1 of Porreca and Maturo (2026). The normalization does not
#'   affect orderings, certification or possibility degrees.
#' @return Numeric vector, one value per actor.
#' @references Yager R.R. (1988). On ordered weighted averaging aggregation
#'   operators in multicriteria decisionmaking. \emph{IEEE Transactions on
#'   Systems, Man, and Cybernetics} 18, 183-190.
#'
#'   Katz L. (1953). A new status index derived from sociometric analysis.
#'   \emph{Psychometrika} 18, 39-43.
#' @examples
#' W <- fsn_example()$M
#' centrality_strength(W, "in")
#' centrality_owa(W, w = c(0.5, 0.3, 0.2, 0, 0, 0, 0))
#' centrality_closeness(W, "in", tnorm = "product")
#' centrality_katz(W, beta = 0.3)
#' @name crisp_indices
#' @section Methodology:
#' These crisp indices are continuous and non-decreasing in every tie weight (isotone), which is the condition under which the \eqn{\alpha}-cuts of the fuzzy centrality are obtained from two crisp evaluations (\code{\link{fuzzy_centrality}}). Strength is linear; the OWA is piecewise linear; the t-norm closeness is the mean of best-path strengths in the \eqn{(\max, T)} semiring, isotone because max and every t-norm are non-decreasing; Katz centrality is a power series with non-negative coefficients, isotone wherever it converges. Betweenness and eigenvector centrality are not isotone. The choice of the t-norm determines how strength decays along a path: the minimum retains the weakest tie, the product multiplies the ties, the Lukasiewicz t-norm subtracts their deficits from one. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' W <- fsn_simulated(n = 20, respondents = 10)$M
#' cl <- sapply(c("min", "product", "lukasiewicz"), function(t) centrality_closeness(W, "in", t))
#' round(cor(cl, method = "kendall"), 2)
#' # step limit
#' one_step <- centrality_closeness(W, "in", "min", max_steps = 1)
#' summary(one_step - centrality_strength(W, "in") / (nrow(W) - 1))
#' # Katz: beta must be below the inverse spectral radius
#' rho <- max(Mod(eigen(W, only.values = TRUE)$values))
#' head(centrality_katz(W, beta = 0.5 / rho))
#' # default of fuzzy_centrality: 0.9 / spectral radius of the binary adjacency
#' A01 <- (W > 0) * 1
#' head(centrality_katz(W, beta = 0.9 / max(Mod(eigen(A01, only.values = TRUE)$values))))
NULL

#' @rdname crisp_indices
#' @export
centrality_strength <- function(W, mode = c("in", "out", "total"), normalize = c("none", "n-1", "n")) {
  mode <- match.arg(mode); normalize <- match.arg(normalize); n <- nrow(W)
  diag(W) <- 0
  v <- switch(mode, `in` = colSums(W), out = rowSums(W), total = colSums(W) + rowSums(W))
  v / switch(normalize, none = 1, `n-1` = n - 1, n = n)
}

#' @rdname crisp_indices
#' @param w Non-negative OWA weights of length \eqn{n - 1}, applied to the ties in
#'   decreasing order.
#' @export
centrality_owa <- function(W, w, mode = c("in", "out")) {
  mode <- match.arg(mode); n <- nrow(W)
  if (length(w) != n - 1) stop("'w' must have length n - 1")
  if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0))
    stop("OWA weights must be finite and non-negative: negative weights make the index non-isotone")
  A <- if (mode == "in") W else t(W)            # column i = ties of actor i
  vapply(seq_len(n), function(i) sum(w * sort(A[-i, i], decreasing = TRUE)), numeric(1))
}

#' @rdname crisp_indices
#' @param tnorm \code{"min"} (Goedel), \code{"product"} or \code{"lukasiewicz"}.
#' @param max_steps Maximum number of steps of the paths (\code{Inf}: no limit).
#' @export
centrality_closeness <- function(W, mode = c("in", "out"), tnorm = c("min", "product", "lukasiewicz"),
                                 max_steps = Inf) {
  mode <- match.arg(mode); S <- best_path_strength(W, tnorm, max_steps); n <- nrow(W)
  if (mode == "in") (colSums(S) - diag(S)) / (n - 1) else (rowSums(S) - diag(S)) / (n - 1)
}

#' Strength of the best paths under a t-norm
#'
#' Computes \eqn{S_{kj} = \max_\pi T(W_{\pi})}, the largest t-norm composition of
#' the ties along a path from \eqn{k} to \eqn{j} with at most \code{max_steps}
#' steps. Since \eqn{T \le \min}, cycles never increase the strength, so the
#' maximum over walks equals the maximum over paths. Uses the Floyd-Warshall
#' recursion in the \eqn{(\max, T)} semiring when there is no step limit, and
#' \code{max_steps} rounds of the Bellman-Ford recursion otherwise.
#' @param W Non-negative matrix with values in \eqn{[0, 1]}.
#' @param tnorm \code{"min"}, \code{"product"} or \code{"lukasiewicz"}.
#' @param max_steps Maximum number of steps.
#' @return Matrix of best-path strengths, with zero diagonal.
#' @examples
#' round(best_path_strength(fsn_example()$M, "min"), 2)
#' @export
#' @section Methodology:
#' For a t-norm \eqn{T}, the strength of a path is the \eqn{T}-composition of its ties and the best-path strength is the maximum over paths. Since \eqn{T(a, b) \le \min(a, b)}, cycles never increase the strength and the maximum over walks equals the maximum over paths. The Floyd-Warshall recursion (Floyd, 1962) in the \eqn{(\max, T)} semiring computes it in \eqn{O(n^3)}. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' W <- fsn_example()$M
#' S1 <- best_path_strength(W, "min", max_steps = 1)
#' S3 <- best_path_strength(W, "min")
#' round(S3 - S1, 2)     # gain from indirect paths
best_path_strength <- function(W, tnorm = c("min", "product", "lukasiewicz"), max_steps = Inf) {
  tnorm <- match.arg(tnorm); n <- nrow(W); diag(W) <- 0
  if (!is.numeric(max_steps) || length(max_steps) != 1 || is.na(max_steps) || max_steps < 1 ||
      (is.finite(max_steps) && max_steps != round(max_steps)))
    stop("max_steps must be a positive integer or Inf")
  if (any(W > 1 + 1e-12)) stop("t-norm closeness requires tie values in [0, 1]")
  Tn <- switch(tnorm, min = pmin, product = function(a, b) a * b,
               lukasiewicz = function(a, b) pmax(0, a + b - 1))
  S <- W
  if (is.infinite(max_steps) || max_steps >= n - 1) {
    # Floyd recursion in the (max, T) semiring, restricted at each step to the
    # actors that reach h and are reached from h (sparse networks)
    for (h in seq_len(n)) {
      ri <- which(S[, h] > 0); cj <- which(S[h, ] > 0)
      if (!length(ri) || !length(cj)) next
      S[ri, cj] <- pmax(S[ri, cj], outer(S[ri, h], S[h, cj], Tn))
    }
  } else if (max_steps > 1) {
    for (k in 2:max_steps) {
      S1 <- S
      for (h in seq_len(n)) S1 <- pmax(S1, Tn(matrix(S[, h], n, n), matrix(W[h, ], n, n, byrow = TRUE)))
      S <- S1
    }
  }
  diag(S) <- 0; dimnames(S) <- dimnames(W); S
}

#' @rdname crisp_indices
#' @param beta Attenuation factor; must satisfy \eqn{\beta\rho(W) < 1}.
#'   In \code{\link{fuzzy_centrality}} the default is
#'   \eqn{\beta = 0.9/\rho(A)}, with \eqn{A} the binary adjacency matrix of the
#'   ties with positive support, as in Porreca and Maturo (2026): since every
#'   compatible network satisfies \eqn{W \le A}, the series converges on the
#'   whole support box.
#' @export
centrality_katz <- function(W, beta, mode = c("in", "out"), normalize = c("none", "n-1")) {
  mode <- match.arg(mode); normalize <- match.arg(normalize); n <- nrow(W); diag(W) <- 0
  if (!is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta < 0)
    stop("beta must be a single finite non-negative number")
  rho <- max(Mod(eigen(W, only.values = TRUE)$values))
  if (beta * rho >= 1) stop(sprintf("beta * spectral radius = %.3f >= 1: Katz centrality is not defined", beta * rho))
  A <- if (mode == "in") t(W) else W
  x <- solve(diag(n) - beta * A, rep(1, n)) - 1
  as.numeric(x) / if (normalize == "n-1") n - 1 else 1
}

#' Registry of isotone centrality indices
#'
#' Returns the specification of an index for \code{\link{fuzzy_centrality}}.
#' Available names: \code{"in_strength"}, \code{"out_strength"},
#' \code{"total_strength"}, \code{"in_owa"}, \code{"out_owa"},
#' \code{"in_closeness"}, \code{"out_closeness"} (t-norm chosen with
#' \code{tnorm}), the shortcuts \code{"in_closeness_min"},
#' \code{"in_closeness_prod"}, \code{"in_closeness_luk"} and the corresponding
#' \code{out_} versions, \code{"in_katz"}, \code{"out_katz"}. A user-defined
#' isotone index can be supplied as a function of the weight matrix.
#'
#' @param index Name of the index or a function \code{function(W)} returning a
#'   vector with one value per actor.
#' @param ... Parameters of the index: \code{normalize} (strength and Katz,
#'   default \code{"n-1"}), \code{w} (OWA), \code{tnorm} and \code{max_steps}
#'   (closeness), \code{beta} (Katz; default \eqn{0.9/\rho(A)}, see
#'   \code{\link{centrality_katz}}).
#' @return A list with the function \code{fun}, the label, and the flags
#'   \code{affine} (the index is affine in the ties, so its \eqn{\alpha}-cut
#'   endpoints are linear in \eqn{\alpha} and the triangle on the endpoints is
#'   exact) and \code{exact_pairs} (the indices of two different actors depend on
#'   disjoint sets of ties, so that marginal possibility and necessity degrees are
#'   exact).
#' @examples
#' str(centrality_index("in_closeness_prod", max_steps = 2))
#' my_index <- function(W) colSums(sqrt(W))
#' fuzzy_centrality(fsn_example(), my_index)
#' @export
#' @section Methodology:
#' The registry records two properties used by the other functions: \code{affine} (strength), for which the triangle on the endpoints is exact and possibility degrees can be located by linear interpolation; and \code{exact_pairs} (in- and out-strength, OWA), for which the centralities of two actors depend on disjoint sets of ties, so that marginal possibility and necessity degrees coincide with the joint ones. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' spec <- centrality_index("in_katz", beta = 0.05)
#' spec$label; spec$affine; spec$exact_pairs
#' spec$fun(fsn_example()$M)
centrality_index <- function(index, ...) {
  a <- list(...)
  if (is.function(index))
    return(list(fun = index, label = if (!is.null(a$label)) a$label else "user-defined index",
                name = "user", affine = isTRUE(a$affine), exact_pairs = isTRUE(a$exact_pairs), params = a))
  alias <- c(in_closeness_min = "min", in_closeness_prod = "product", in_closeness_luk = "lukasiewicz",
             out_closeness_min = "min", out_closeness_prod = "product", out_closeness_luk = "lukasiewicz")
  if (index %in% names(alias)) { a$tnorm <- alias[[index]]; index <- sub("_(min|prod|luk)$", "", index) }
  mode <- sub("_.*$", "", index); type <- sub("^[a-z]+_", "", index)
  if (!mode %in% c("in", "out", "total") || !type %in% c("strength", "owa", "closeness", "katz") ||
      (mode == "total" && type != "strength"))
    stop("unknown index '", index, "'")
  switch(type,
         strength = {
           nz <- if (is.null(a$normalize)) "n-1" else a$normalize
           list(fun = function(W) centrality_strength(W, mode, nz), name = index,
                label = paste0(mode, "-strength"), affine = TRUE, exact_pairs = mode != "total", params = a)
         },
         owa = {
           if (is.null(a$w)) stop("an OWA index requires the weights 'w'")
           if (!is.numeric(a$w) || any(!is.finite(a$w)) || any(a$w < 0))
             stop("OWA weights must be finite and non-negative: negative weights make the index non-isotone")
           list(fun = function(W) centrality_owa(W, a$w, mode), name = index, label = paste0(mode, "-OWA"),
                affine = FALSE, exact_pairs = TRUE, params = a)
         },
         closeness = {
           tn <- if (is.null(a$tnorm)) "min" else a$tnorm
           st <- if (is.null(a$max_steps)) Inf else a$max_steps
           list(fun = function(W) centrality_closeness(W, mode, tn, st), name = index,
                label = sprintf("%s-closeness (%s t-norm%s)", mode, tn,
                                if (is.finite(st)) paste0(", s = ", st) else ""),
                affine = FALSE, exact_pairs = FALSE, params = c(a, list(tnorm = tn, max_steps = st)))
         },
         katz = {
           nz <- if (is.null(a$normalize)) "n-1" else a$normalize
           list(fun = function(W) centrality_katz(W, a$beta, mode, nz), name = index,
                label = paste0(mode, "-Katz"), affine = FALSE, exact_pairs = FALSE, needs_beta = TRUE, params = a)
         })
}

#' Check numerically that an index is isotone
#'
#' Draws random networks in the support box of a fuzzy network, increases one
#' random tie at a time and checks that no centrality decreases.
#' @param x An \code{\link{fsn}} object.
#' @param index Index name or function (see \code{\link{centrality_index}}).
#' @param n_check Number of random checks.
#' @param tol Tolerance.
#' @param ... Parameters of the index.
#' @return Logical, with attributes \code{max_decrease} (largest decrease found) and
#'   \code{n_checked} (number of random pairs of networks compared; zero, with a
#'   note, when no tie can vary). A numerical check on random networks, not a
#'   proof of monotonicity.
#' @examples
#' check_isotone(fsn_example(), "in_closeness_prod", n_check = 20)
#' check_isotone(fsn_example(), function(W) -colSums(W), n_check = 20)
#' @export
#' @section Methodology:
#' The exact computation is valid only for isotone indices. This function provides a numerical check, not a proof: it can detect a violation but cannot establish monotonicity. \code{\link{fuzzy_centrality}} also checks that the computed cuts are nested. See Section 2.4 of Porreca and Maturo (2026).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- fsn_example()
#' check_isotone(x, function(W) colSums(W)^2)   # isotone on non-negative ties
#' check_isotone(x, function(W) colSums(W) - rowSums(W), n_check = 30)  # not isotone
check_isotone <- function(x, index, n_check = 100, tol = 1e-9, ...) {
  spec <- centrality_index(index, ...)
  if (isTRUE(spec$needs_beta) && is.null(spec$params$beta)) spec <- .katz_beta(spec, x)
  ct <- fsn_cut(x, 0); E <- which(ct$upper > ct$lower); worst <- 0
  if (!length(E)) {
    # no tie can vary: the box reduces to one network and there is nothing to check
    return(structure(TRUE, max_decrease = 0, n_checked = 0L,
                     note = "no tie with positive spread: monotonicity not checked"))
  }
  for (k in seq_len(n_check)) {
    W <- ct$lower + (ct$upper - ct$lower) * matrix(stats::runif(length(x$M)), nrow(x$M))
    e <- E[sample.int(length(E), 1)]
    W2 <- W; W2[e] <- W[e] + stats::runif(1) * (ct$upper[e] - W[e])
    d <- spec$fun(W2) - spec$fun(W); worst <- min(worst, min(d))
  }
  structure(worst >= -tol, max_decrease = -worst, n_checked = as.integer(n_check))
}

.katz_beta <- function(spec, x, share = 0.9) {
  A <- (fsn_cut(x, 0)$upper > 0) * 1; diag(A) <- 0
  rho <- max(Mod(eigen(A, only.values = TRUE)$values))
  beta <- if (rho > 0) share / rho else 0.5
  pr <- spec$params; pr$beta <- beta; pr$label <- NULL
  do.call(centrality_index, c(list(spec$name), pr))
}

#' Fuzzy spectral radius of a fuzzy network
#'
#' Exact \eqn{\alpha}-cuts of the spectral radius of the network, a
#' network-level index that is continuous and non-decreasing in every tie.
#' @param x An \code{\link{fsn}} object.
#' @param alpha Grid of levels.
#' @return A data frame with columns \code{alpha}, \code{lower}, \code{upper}.
#' @export
#' @section Methodology:
#' The spectral radius of a non-negative matrix is continuous and non-decreasing in each entry (Perron-Frobenius theory), so by Theorem 1 of Porreca and Maturo (2026) its \eqn{\alpha}-cut is \eqn{[\rho(W^L_\alpha), \rho(W^U_\alpha)]} (Proposition 1(e) of that paper). It bounds the attenuation factor admissible for Katz centrality on every compatible network.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' fuzzy_spectral_radius(fsn_example(), alpha = c(0, 0.5, 1))
#' sr <- fuzzy_spectral_radius(fsn_simulated(n = 20, respondents = 10))
#' plot(c(sr$lower, rev(sr$upper)), c(sr$alpha, rev(sr$alpha)), type = "l",
#'      xlab = "spectral radius", ylab = expression(alpha))
fuzzy_spectral_radius <- function(x, alpha = seq(0, 1, by = 0.05)) {
  rho <- function(W) { diag(W) <- 0; max(Mod(eigen(W, only.values = TRUE)$values)) }
  r <- t(vapply(alpha, function(a) { ct <- fsn_cut(x, a); c(rho(ct$lower), rho(ct$upper)) }, numeric(2)))
  data.frame(alpha = alpha, lower = r[, 1], upper = r[, 2])
}
