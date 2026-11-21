# ---------------------------------------------------------------------------
# Fuzzy social networks: class, constructors, accessors
# ---------------------------------------------------------------------------

#' Fuzzy social networks with triangular fuzzy ties
#'
#' Creates a directed (or undirected) fuzzy social network in which the tie from
#' actor \eqn{i} to actor \eqn{j} is the triangular fuzzy number
#' \eqn{(L_{ij}, M_{ij}, R_{ij})}. Absent ties are null fuzzy numbers
#' \eqn{(0, 0, 0)}, so that every pair of actors is connected to some degree, as
#' in Porreca, Maturo and Ventre (2025). Rows are senders (raters) and columns
#' receivers (rated actors): in-centralities are computed on the columns.
#'
#' Under the epistemic reading adopted by the exact methods of the package, each
#' TFN is a possibility distribution on the ill-known strength of a crisp tie,
#' and the network is a fuzzy set of crisp weighted networks.
#'
#' @param L,M,R Numeric square matrices with the left endpoints, the cores and
#'   the right endpoints of the ties. \code{L} and \code{R} default to \code{M}
#'   (crisp ties). Missing values are treated as absent ties.
#' @param labels Actor labels (default: row names of \code{M} or \code{1:n}).
#' @param directed Logical; if \code{FALSE} the three matrices must be symmetric.
#' @param respondents Optional logical vector or labels identifying the actors
#'   that sent ratings (respondents of a questionnaire).
#' @param info Optional list of metadata.
#' @param unknown Optional logical matrix flagging pairs whose tie is unknown
#'   (e.g. ties that non-respondents would have sent). Under the possibilistic
#'   reading an unknown tie is the vacuous possibility distribution on
#'   \code{unknown_range}: every value is fully possible, so its
#'   \eqn{\alpha}-cuts are the whole range at every level. The exact methods
#'   (\code{\link{fuzzy_centrality}}, \code{\link{mc_centrality}}) propagate
#'   this ignorance; the indices of Porreca, Maturo and Ventre (2025) treat
#'   unknown ties as null. See \code{\link{set_unknown}}.
#' @param unknown_range Range of the values of an unknown tie.
#' @param na Treatment of ties with missing values (in any of \code{L},
#'   \code{M}, \code{R}): \code{"zero"} (null ties, with a warning),
#'   \code{"unknown"} (unknown ties, see \code{\link{set_unknown}}) or
#'   \code{"error"}.
#' @return An object of class \code{"fsn"}: a list with components \code{L},
#'   \code{M}, \code{R}, \code{labels}, \code{directed}, \code{respondents},
#'   \code{info} and, if present, \code{unknown} and \code{unknown_range}.
#' @seealso \code{\link{fsn_from_ratings}}, \code{\link{simulate_fsn}},
#'   \code{\link{fsn_cut}}
#' @examples
#' M <- matrix(c(0, .6, .2,
#'               .5, 0, .8,
#'               0, .7, 0), 3, byrow = TRUE)
#' x <- fsn(L = pmax(M - 0.1, 0) * (M > 0), M = M, R = pmin(M + 0.1, 1) * (M > 0),
#'          labels = c("Ann", "Bob", "Cid"))
#' x
#' summary(x)
#' fsn_tie(x, "Ann", "Bob")
#' @export
#' @section Methodology:
#' A directed fuzzy social network (Porreca, Maturo and Ventre, 2025, Definition 4) is a set of actors and a matrix of TFN ties \eqn{\tilde W_{kl} = (L_{kl}, M_{kl}, R_{kl})}. Since absent ties are null fuzzy numbers, every pair of actors is connected to some degree and the network can be treated as complete.
#'
#' Under the possibilistic reading each tie constrains the latent crisp tie \eqn{W_{kl}}, and the network is a fuzzy set of crisp weighted networks with membership \eqn{\min_{kl}\mu_{kl}(W_{kl})}. Its \eqn{\alpha}-cut is the box \eqn{\mathcal W_\alpha = \{W : \underline W_\alpha \le W \le \overline W_\alpha\}} (\code{\link{fsn_cut}}), on which the exact centralities are computed.
#'
#' Pairs whose tie is not known, typically because the potential rater did not respond, can be declared unknown (\code{\link{set_unknown}}); they are then vacuous possibility distributions, whose \eqn{\alpha}-cuts are the whole range at every level.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' # a network with an unknown tie
#' M <- matrix(c(0, .6, 0, .5, 0, .8, 0, .7, 0), 3, byrow = TRUE)
#' U <- matrix(FALSE, 3, 3); U[3, 1] <- TRUE
#' y <- fsn(L = 0.8 * M, M = M, R = pmin(1.2 * M, 1), labels = c("a", "b", "c"), unknown = U)
#' y
#' fsn_cut(y, 0.5)
#' # undirected networks require symmetric matrices
#' S <- (M + t(M)) / 2
#' fsn(S, S, S, directed = FALSE)
fsn <- function(L = M, M, R = M, labels = NULL, directed = TRUE, respondents = NULL, info = list(),
                unknown = NULL, unknown_range = c(0, 1), na = c("zero", "unknown", "error")) {
  na <- match.arg(na)
  M <- as.matrix(M); L <- as.matrix(L); R <- as.matrix(R)
  n <- nrow(M)
  if (ncol(M) != n || any(dim(L) != dim(M)) || any(dim(R) != dim(M)))
    stop("L, M and R must be square matrices of the same size")
  if (any(is.infinite(c(L, M, R)))) stop("tie values must be finite")
  isna <- is.na(L) | is.na(M) | is.na(R); diag(isna) <- FALSE
  if (any(isna)) {
    if (na == "error") stop(sum(isna), " ties contain missing values")
    if (na == "zero") warning(sum(isna), " ties with missing values are treated as null ties (0, 0, 0); ",
                              "use na = 'unknown' to treat them as unknown", call. = FALSE)
    if (na == "unknown") unknown <- (if (is.null(unknown)) matrix(FALSE, n, n) else as.matrix(unknown)) | isna
  }
  for (nm in c("L", "M", "R")) { z <- get(nm); z[isna] <- 0; z[is.na(z)] <- 0; storage.mode(z) <- "double"; assign(nm, z) }
  absent <- M == 0 & L == 0 & R == 0
  if (any(L < 0)) stop("negative tie values are not allowed")
  tol <- 1e-10
  if (any(L > M + tol | M > R + tol)) stop("every tie must satisfy L <= M <= R")
  diag(L) <- 0; diag(M) <- 0; diag(R) <- 0
  if (is.null(labels)) labels <- if (!is.null(rownames(M))) rownames(M) else as.character(seq_len(n))
  labels <- as.character(labels)
  if (length(labels) != n || anyDuplicated(labels)) stop("'labels' must be n distinct values")
  dimnames(L) <- dimnames(M) <- dimnames(R) <- list(labels, labels)
  if (!directed && (!isSymmetric(unname(M)) || !isSymmetric(unname(L)) || !isSymmetric(unname(R))))
    stop("an undirected network requires symmetric matrices")
  if (!is.null(respondents)) {
    if (!is.logical(respondents)) respondents <- labels %in% as.character(respondents)
    names(respondents) <- labels
  }
  out <- list(L = L, M = M, R = R, labels = labels, directed = directed, respondents = respondents, info = info)
  if (!is.null(unknown) && any(unknown)) {
    unknown <- as.matrix(unknown) & TRUE; diag(unknown) <- FALSE
    if (any(dim(unknown) != dim(M))) stop("'unknown' must have the dimensions of M")
    if (any((M > 0 | R > 0) & unknown)) stop("a pair cannot have both a rating and an unknown tie")
    if (!directed && !isSymmetric(unname(unknown)))
      stop("an undirected network requires a symmetric 'unknown' matrix")
    if (!is.numeric(unknown_range) || length(unknown_range) != 2 || anyNA(unknown_range) ||
        unknown_range[1] < 0 || unknown_range[2] > 1 || unknown_range[1] > unknown_range[2])
      stop("'unknown_range' must be an interval [a, b] with 0 <= a <= b <= 1")
    dimnames(unknown) <- list(labels, labels)
    out$unknown <- unknown; out$unknown_range <- unknown_range
  }
  structure(out, class = "fsn")
}

#' Declare ties as unknown
#'
#' In questionnaire data, a pair without a rating may mean that there is no
#' collaboration (a null tie, the convention of Porreca, Maturo and Ventre, 2025)
#' or that the tie is unknown, because the potential rater did not respond. This
#' function flags unknown ties, which the exact methods treat as vacuous
#' possibility distributions on \code{range}: the resulting fuzzy centralities
#' contain all the values compatible with any answer of the non-respondents.
#' @param x An \code{\link{fsn}} object.
#' @param which \code{"nonrespondents"}: the ties that non-respondents would
#'   have sent; \code{"unrated"}: every pair without a rating; \code{"none"}:
#'   remove the flags.
#' @param range Range of the values of an unknown tie.
#' @return An \code{\link{fsn}} object.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' Au <- set_unknown(A, "nonrespondents")
#' Au
#' fc0 <- fuzzy_centrality(A, "in_strength")
#' fc1 <- fuzzy_centrality(Au, "in_strength")
#' summary(fc1$upper[, 1] - fc1$lower[, 1])
#' @export
#' @section Methodology:
#' An unknown tie is the vacuous possibility distribution on \code{range}: every value is fully possible, so its \eqn{\alpha}-cut is the whole range at every level. For an isotone index this is handled exactly, by setting unknown ties to the lower end of the range in \eqn{\underline W_\alpha} and to the upper end in \eqn{\overline W_\alpha}. The comparison with the default convention (null ties) shows how much a ranking depends on the treatment of non-response. Unknown ties are ignored (treated as null) by the indices of Porreca, Maturo and Ventre (2025).
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' Au <- set_unknown(A, "nonrespondents")
#' f0 <- fuzzy_centrality(A, "in_strength"); f1 <- fuzzy_centrality(Au, "in_strength")
#' c(null = certification(f0)$proportion_certified[1],
#'   unknown = certification(f1)$proportion_certified[1])
#' op <- par(mfrow = c(1, 2))
#' plot(f0, "intervals", actors = 10, main = "Unrated = null")
#' plot(f1, "intervals", actors = 10, main = "Non-respondents unknown")
#' par(op)
set_unknown <- function(x, which = c("nonrespondents", "unrated", "none"), range = c(0, 1)) {
  which <- match.arg(which)
  E <- x$M > 0 | x$R > 0
  U <- switch(which,
              none = matrix(FALSE, nrow(E), ncol(E)),
              unrated = !E,
              nonrespondents = {
                if (is.null(x$respondents)) stop("the network has no information on respondents")
                matrix(!x$respondents, nrow(E), ncol(E)) & !E
              })
  diag(U) <- FALSE
  fsn(x$L, x$M, x$R, labels = x$labels, directed = x$directed, respondents = x$respondents,
      info = x$info, unknown = U, unknown_range = range)
}

#' @rdname fsn
#' @param x An object.
#' @export
is.fsn <- function(x) inherits(x, "fsn")

#' Number of actors of a fuzzy social network
#' @param x An \code{\link{fsn}} object.
#' @return Integer.
#' @examples
#' n_actors(fsn_example())
#' @export
#' @section Methodology:
#' Number of actors \eqn{n}; indices averaged over the other actors divide by \eqn{n - 1}, so including or excluding isolated members changes their values. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @examples
#' n_actors(fsn_simulated(n = 20, respondents = 10))
#' n_actors(drop_isolated(fsn_simulated(n = 20, respondents = 10)))   # actors with at least one tie
n_actors <- function(x) nrow(x$M)

#' @export
print.fsn <- function(x, ...) {
  E <- x$M > 0 | x$R > 0; n <- n_actors(x)
  cat(if (x$directed) "Directed" else "Undirected", "fuzzy social network with", n, "actors and",
      sum(E) / if (x$directed) 1 else 2, "non-null ties\n")
  if (!is.null(x$respondents)) cat("Respondents:", sum(x$respondents), "of", n, "\n")
  if (!is.null(x$unknown)) cat("Unknown ties (vacuous on [", x$unknown_range[1], ", ", x$unknown_range[2], "]): ",
                               sum(x$unknown), "\n", sep = "")
  if (sum(E)) {
    w <- (x$R - x$L)[E]
    cat(sprintf("Density %.3f; cores in [%.3f, %.3f]; mean support width %.3f\n",
                mean(E[row(E) != col(E)]), min(x$M[E]), max(x$M[E]), mean(w)))
  }
  if (!is.null(x$info$description)) cat(x$info$description, "\n")
  invisible(x)
}

#' @export
summary.fsn <- function(object, ...) {
  x <- object; E <- x$M > 0 | x$R > 0
  ties <- data.frame(l = x$L[E], m = x$M[E], r = x$R[E])
  sp <- tfn_spread(tfn(ties$l, ties$m, ties$r))
  out <- list(n = n_actors(x), n_ties = sum(E), density = mean(E[row(E) != col(E)]),
              directed = x$directed,
              respondents = if (is.null(x$respondents)) NA else sum(x$respondents),
              isolated = x$labels[rowSums(E) + colSums(E) == 0],
              reciprocity = if (x$directed && sum(E)) sum(E & t(E)) / sum(E) else NA,
              cores = summary(ties$m), widths = summary(sp$width),
              share_asymmetric = mean(abs(sp$asymmetry) > 0.25 * pmax(sp$width, 1e-12)),
              in_strength_core = colSums(x$M), out_strength_core = rowSums(x$M))
  class(out) <- "summary.fsn"; out
}

#' @export
print.summary.fsn <- function(x, ...) {
  cat("Fuzzy social network:", x$n, "actors,", x$n_ties, "non-null ties, density", round(x$density, 3), "\n")
  if (!is.na(x$respondents)) cat("Respondents:", x$respondents, "\n")
  if (length(x$isolated)) cat("Isolated actors:", paste(x$isolated, collapse = ", "), "\n")
  if (!is.na(x$reciprocity)) cat("Share of reciprocated ties:", round(x$reciprocity, 3), "\n")
  cat("Cores of the ties:\n"); print(x$cores)
  cat("Support widths r - l:\n"); print(x$widths)
  cat("Share of ties with |(r-m)-(m-l)| > (r-l)/4:", round(x$share_asymmetric, 3), "\n")
  invisible(x)
}

#' Extract ties as triangular fuzzy numbers
#'
#' @param x An \code{\link{fsn}} object.
#' @param from,to Actor labels or indices.
#' @return \code{fsn_tie}: a \code{\link{tfn}} of length one. \code{fsn_ties}: a
#'   data frame with one row per non-null tie (\code{from}, \code{to}, \code{l},
#'   \code{m}, \code{r}).
#' @examples
#' x <- fsn_example()
#' fsn_tie(x, "A", "B")
#' head(fsn_ties(x))
#' @export
#' @section Methodology:
#' Ties are stored with rows as senders and columns as receivers; \code{fsn_ties} returns the non-null ties in long format, which is the input of \code{\link{fsn_from_ratings}}.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' head(fsn_ties(A))
#' plot(fsn_tie(A, "a05", "a01"))
#' nrow(fsn_ties(A, include_null = TRUE))   # n (n - 1) ordered pairs
fsn_tie <- function(x, from, to) {
  i <- .actor_index(x, from); j <- .actor_index(x, to)
  tfn(x$L[i, j], x$M[i, j], x$R[i, j], names = paste(x$labels[i], x$labels[j], sep = "->"))
}

#' @rdname fsn_tie
#' @param include_null Logical, include null ties.
#' @export
fsn_ties <- function(x, include_null = FALSE) {
  n <- n_actors(x); idx <- which(row(x$M) != col(x$M), arr.ind = FALSE)
  E <- (x$M > 0 | x$R > 0)
  keep <- if (include_null) idx else idx[E[idx]]
  rr <- row(x$M)[keep]; cc <- col(x$M)[keep]
  data.frame(from = x$labels[rr], to = x$labels[cc], l = x$L[keep], m = x$M[keep], r = x$R[keep],
             stringsAsFactors = FALSE)
}

.actor_index <- function(x, a) {
  if (is.character(a)) { i <- match(a, x$labels); if (anyNA(i)) stop("unknown actor(s): ", paste(a[is.na(i)], collapse = ", ")); i }
  else as.integer(a)
}

#' Subset the actors of a fuzzy social network
#'
#' @param x An \code{\link{fsn}} object.
#' @param i Actor labels, indices or a logical vector.
#' @param ... Unused.
#' @return An \code{\link{fsn}} with the ties among the selected actors.
#' @examples
#' fsn_example()[c("A", "B", "C")]
#' @export
#' @section Methodology:
#' Subsetting keeps the ties among the selected actors. Centralities of a subnetwork differ from those of the same actors in the whole network, since paths through excluded actors and the normalization by \eqn{n - 1} change.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' resp <- A[A$respondents]       # network among respondents only
#' resp
#' summary(resp)
`[.fsn` <- function(x, i, ...) {
  if (is.logical(i)) i <- which(i) else i <- .actor_index(x, i)
  y <- fsn(x$L[i, i, drop = FALSE], x$M[i, i, drop = FALSE], x$R[i, i, drop = FALSE],
           labels = x$labels[i], directed = x$directed,
           respondents = if (is.null(x$respondents)) NULL else x$respondents[i], info = x$info,
           unknown = if (is.null(x$unknown)) NULL else x$unknown[i, i, drop = FALSE],
           unknown_range = if (is.null(x$unknown_range)) c(0, 1) else x$unknown_range)
  # metadata of simulated networks, reduced to the selected actors
  if (!is.null(x$latent)) y$latent <- x$latent[i, i, drop = FALSE]
  if (is.matrix(x$sigma)) y$sigma <- x$sigma[i, i, drop = FALSE]
  if (!is.null(x$blocks)) y$blocks <- x$blocks[i]
  y
}

#' Remove isolated actors
#'
#' Removes the actors that neither send nor receive a non-null or unknown tie.
#' @param x An \code{\link{fsn}} object.
#' @param keep_respondents Logical: keep respondents without ties (a respondent
#'   who rated nobody declares null ties and is informative).
#' @return An \code{\link{fsn}} object.
#' @examples
#' summary(drop_isolated(fsn_example()))
#' @export
#' @section Methodology:
#' Actors without non-null ties do not affect sums but change averages over \eqn{n - 1} actors.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' setdiff(A$labels, drop_isolated(A)$labels)
drop_isolated <- function(x, keep_respondents = FALSE) {
  E <- x$M > 0 | x$R > 0
  if (!is.null(x$unknown)) E <- E | x$unknown
  keep <- rowSums(E) + colSums(E) > 0
  if (keep_respondents && !is.null(x$respondents)) keep <- keep | x$respondents
  x[keep]
}

#' Alpha-cuts of the ties of a fuzzy social network
#'
#' The \eqn{\alpha}-cut of the network is the box
#' \eqn{\mathcal W_\alpha = \{W: \underline W_\alpha \le W \le \overline W_\alpha\}}
#' of crisp weighted networks whose ties lie in the \eqn{\alpha}-cuts of the
#' fuzzy ties. Unknown ties (\code{\link{set_unknown}}) span their whole range
#' at every level.
#' @param x An \code{\link{fsn}} object.
#' @param alpha A level in \eqn{[0, 1]}.
#' @return A list with the matrices \code{lower} and \code{upper}.
#' @examples
#' fsn_cut(fsn_example(), 0.5)$lower[1:3, 1:3]
#' @export
#' @section Methodology:
#' Returns the lower and upper endpoint networks \eqn{\underline W_\alpha} and \eqn{\overline W_\alpha}, whose entries are the endpoints of the \eqn{\alpha}-cuts of the ties. For an isotone index \eqn{c} the \eqn{\alpha}-cut of the fuzzy centrality is \eqn{[c(\underline W_\alpha), c(\overline W_\alpha)]} (\code{\link{fuzzy_centrality}}). Unknown ties span their whole range at every level. See Section 2.4 of Porreca and Maturo (2026).
#' @examples
#' x <- fsn_example()
#' ct <- fsn_cut(x, 0.5)
#' centrality_strength(ct$lower, "in"); centrality_strength(ct$upper, "in")
#' # equals the 0.5-cut of the exact fuzzy in-strength
#' fc <- fuzzy_centrality(x, "in_strength")
#' cbind(fc$lower[, "0.50"], fc$upper[, "0.50"])
fsn_cut <- function(x, alpha) {
  if (length(alpha) != 1 || alpha < 0 || alpha > 1) stop("alpha must be a single value in [0, 1]")
  lo <- x$L + alpha * (x$M - x$L); up <- x$R - alpha * (x$R - x$M)
  if (!is.null(x$unknown)) { lo[x$unknown] <- x$unknown_range[1]; up[x$unknown] <- x$unknown_range[2] }
  list(lower = lo, upper = up)
}

#' Convert to and from other representations
#'
#' \code{as_long} returns the non-null ties in long format; \code{as_igraph}
#' converts the network into an \pkg{igraph} object with edge attributes
#' \code{l}, \code{m}, \code{r}, \code{cog} and \code{weight} (the core);
#' \code{fsn_from_igraph} builds a network from edge attributes.
#' @param x An \code{\link{fsn}} object.
#' @param g An \pkg{igraph} graph.
#' @param l,m,r Names of the edge attributes with the TFN parameters.
#' @return A data frame, an \pkg{igraph} object or an \code{\link{fsn}} object.
#' @examples
#' head(as_long(fsn_example()))
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- as_igraph(fsn_example())
#'   fsn_from_igraph(g)
#' }
#' @export
#' @section Methodology:
#' Conversions keep the three TFN parameters as edge attributes. In \pkg{igraph} the edge attribute \code{weight} is the core, so crisp algorithms applied to the converted graph analyse the network of the cores and ignore vagueness.
#' @examples
#' x <- fsn_example()
#' d <- as_long(x)
#' x2 <- fsn_from_ratings(d, actors = x$labels, drop_isolated = FALSE)
#' all.equal(x$M, x2$M)
as_long <- function(x) fsn_ties(x)

#' @rdname as_long
#' @export
as_igraph <- function(x) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("package 'igraph' is required")
  e <- fsn_ties(x)
  if (!x$directed && nrow(e)) {
    # one edge per unordered pair: the matrices of an undirected network are symmetric
    ip <- match(e$from, x$labels); jp <- match(e$to, x$labels)
    e <- e[ip < jp, , drop = FALSE]
  }
  g <- igraph::graph_from_data_frame(e[, c("from", "to")], directed = x$directed,
                                     vertices = data.frame(name = x$labels))
  igraph::E(g)$l <- e$l; igraph::E(g)$m <- e$m; igraph::E(g)$r <- e$r
  igraph::E(g)$cog <- (e$l + e$m + e$r) / 3; igraph::E(g)$weight <- e$m
  g
}

#' @rdname as_long
#' @export
fsn_from_igraph <- function(g, l = "l", m = "m", r = "r") {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("package 'igraph' is required")
  lab <- igraph::V(g)$name; if (is.null(lab)) lab <- as.character(seq_len(igraph::vcount(g)))
  ed <- igraph::as_edgelist(g, names = FALSE); n <- length(lab)
  L <- M <- R <- matrix(0, n, n)
  ea <- function(a) igraph::edge_attr(g, a)
  L[ed] <- ea(l); M[ed] <- ea(m); R[ed] <- ea(r)
  if (!igraph::is_directed(g)) { L[ed[, 2:1]] <- ea(l); M[ed[, 2:1]] <- ea(m); R[ed[, 2:1]] <- ea(r) }
  fsn(L, M, R, labels = lab, directed = igraph::is_directed(g))
}

#' Build a fuzzy social network from ratings in long format
#'
#' Each row of \code{data} is a rating given by a rater to a rated actor,
#' expressed as a TFN. Pairs without a rating are null ties; when the same pair
#' is rated more than once, \code{keep} decides which rating is used.
#'
#' @param data A data frame.
#' @param from,to Names of the columns with the rater and the rated actor.
#' @param l,m,r Names of the columns with the TFN parameters.
#' @param order Optional name of a column giving the temporal order of the
#'   ratings (default: row order).
#' @param keep \code{"last"} (later ratings replace earlier ones, as in the
#'   original study), \code{"first"} or \code{"mean"} (componentwise mean).
#' @param actors Optional vector with all the members of the network, including
#'   those who neither rated nor were rated; by default the actors appearing in
#'   \code{data}.
#' @param drop_isolated Logical, remove actors without non-null ties; respondents
#'   are kept. Rows with an empty or missing rated actor, or with missing TFN
#'   values, are participation records: the rater is a respondent but no tie is
#'   created (this is how \code{\link{make_collector}} records respondents who
#'   rate nobody).
#' @param scale Optional positive number by which the TFN parameters are divided
#'   (e.g. 5 for ratings on \eqn{[1, 5]}).
#' @param info Optional list of metadata.
#' @return An \code{\link{fsn}} object, with the respondents (actors that appear
#'   as raters) recorded.
#' @examples
#' r <- data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "a", "a"),
#'                 l = c(.2, .5, .6, .1), m = c(.4, .6, .8, .3), r = c(.5, .9, .9, .3))
#' fsn_from_ratings(r)
#' # ratings simulated as in a questionnaire survey
#' sr <- simulate_ratings(n = 30, seed = 1)
#' head(sr)
#' fsn_from_ratings(sr, from = "rater", to = "rated", order = "order")
#' @export
#' @section Methodology:
#' Builds the network from one row per rating. Pairs without a rating become null ties: this is the convention of the original study, and it is an assumption, since a missing rating may mean no collaboration or no response. Use \code{\link{set_unknown}} to treat the ties that non-respondents would have sent, or all unrated pairs, as unknown. When a pair is rated more than once the last rating is kept by default.
#' @examples
#' # three-point answers on a 0-10 scale
#' f <- system.file("extdata", "three_point_example.csv", package = "FSNA")
#' d <- fuzzify_responses(read_ratings(f), "three_point", scale = c(0, 10))
#' fsn_from_ratings(d, drop_isolated = FALSE)
#' # repeated ratings: last, first or mean
#' r <- data.frame(from = c("a", "a", "b"), to = c("b", "b", "a"),
#'                 l = c(.2, .4, .5), m = c(.3, .6, .6), r = c(.5, .7, .8))
#' sapply(c("last", "first", "mean"), function(k) fsn_from_ratings(r, keep = k)$M["a", "b"])
fsn_from_ratings <- function(data, from = "from", to = "to", l = "l", m = "m", r = "r",
                             order = NULL, keep = c("last", "first", "mean"), actors = NULL,
                             drop_isolated = TRUE, scale = NULL, info = list()) {
  keep <- match.arg(keep)
  d <- as.data.frame(data)
  miss <- setdiff(c(from, to, l, m, r), names(d))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
  if (!is.null(order)) d <- d[base::order(d[[order]]), , drop = FALSE]
  # participation records: a rater with no target or no rating (e.g. a respondent who
  # rated nobody) is a respondent, not a tie
  part <- is.na(d[[to]]) | trimws(as.character(d[[to]])) == "" | is.na(d[[l]]) | is.na(d[[m]]) | is.na(d[[r]])
  raters <- unique(d[[from]][!is.na(d[[from]])])
  d <- d[!part, , drop = FALSE]
  d <- d[as.character(d[[from]]) != as.character(d[[to]]), , drop = FALSE]
  key <- paste(d[[from]], d[[to]], sep = "\r")
  if (keep == "mean") {
    agg <- stats::aggregate(d[, c(l, m, r)], by = list(key = key), FUN = mean)
    first <- d[!duplicated(key), c(from, to)]; first$key <- key[!duplicated(key)]
    d <- merge(first, agg, by = "key")
  } else {
    d <- d[if (keep == "last") !duplicated(key, fromLast = TRUE) else !duplicated(key), , drop = FALSE]
  }
  if (is.null(actors)) {
    actors <- unique(c(raters, d[[from]], d[[to]]))
    if (is.numeric(actors)) actors <- sort(actors)
  }
  lab <- as.character(actors); n <- length(lab)
  i <- match(as.character(d[[from]]), lab); j <- match(as.character(d[[to]]), lab)
  if (anyNA(i) || anyNA(j)) stop("some raters or rated actors are not in 'actors'")
  s <- if (is.null(scale)) 1 else scale
  L <- M <- R <- matrix(0, n, n)
  L[cbind(i, j)] <- d[[l]] / s; M[cbind(i, j)] <- d[[m]] / s; R[cbind(i, j)] <- d[[r]] / s
  x <- fsn(L, M, R, labels = lab, respondents = lab %in% as.character(raters), info = info)
  if (drop_isolated) drop_isolated(x, keep_respondents = TRUE) else x
}

#' A small example network
#'
#' A directed fuzzy network with eight actors used in the examples, in which the
#' direct tie from \code{A} to \code{D} is weaker than the indirect route through
#' \code{B} (the situation of Fig. 3 in Porreca, Maturo and Ventre, 2025).
#' @return An \code{\link{fsn}} object.
#' @examples
#' x <- fsn_example()
#' plot(x)
#' @export
#' @section Methodology:
#' A small network for examples, in which the direct tie from A to D is weak while the indirect route through B is strong, so that best paths (\code{\link{ijar_best_paths}}) and t-norm closeness differ from the direct ties.
#' @examples
#' x <- fsn_example()
#' round(x$M, 2)
#' best_path_route(ijar_best_paths(x, 3), "A", "D")
fsn_example <- function() {
  lab <- LETTERS[1:8]
  e <- data.frame(
    from = c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E", "F", "F", "G", "H", "H", "B", "D", "G"),
    to   = c("B", "D", "D", "C", "A", "E", "E", "B", "F", "C", "G", "A", "H", "A", "F", "A", "G", "E"),
    l    = c(.70, .00, .75, .30, .40, .50, .20, .55, .60, .10, .45, .05, .30, .15, .55, .50, .25, .35),
    m    = c(.80, .10, .85, .45, .55, .60, .40, .70, .70, .30, .60, .20, .45, .25, .70, .65, .35, .50),
    r    = c(.90, .30, .95, .60, .60, .80, .55, .75, .85, .40, .70, .45, .70, .40, .80, .70, .50, .60))
  fsn_from_ratings(e, actors = lab, info = list(description = "Example network (fsn_example)"))
}
