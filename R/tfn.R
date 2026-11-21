# ---------------------------------------------------------------------------
# Triangular fuzzy numbers: class, arithmetic, alpha-cuts, ranking scores
# ---------------------------------------------------------------------------

#' Triangular fuzzy numbers
#'
#' Creates a vector of triangular fuzzy numbers (TFNs) \eqn{(l, m, r)} with
#' \eqn{l \le m \le r}, where \eqn{[l, r]} is the support and \eqn{m} the core.
#' The membership function is \eqn{(x - l)/(m - l)} on \eqn{[l, m]},
#' \eqn{(r - x)/(r - m)} on \eqn{[m, r]} and zero elsewhere; degenerate sides
#' (\eqn{l = m} or \eqn{m = r}) and crisp numbers (\eqn{l = m = r}) are allowed.
#'
#' Arithmetic follows the pseudo-operations of Porreca, Maturo and Ventre (2025,
#' Eqs. 5-7): \code{+} and \code{*} act componentwise on the three parameters and
#' multiplication by a non-negative scalar rescales them. Addition and
#' multiplication by a non-negative scalar coincide with the extension principle;
#' componentwise multiplication coincides with it only for non-negative TFNs at
#' the support and at the core.
#'
#' @param l,m,r Numeric vectors (recycled) with left endpoint, core and right
#'   endpoint. Alternatively \code{l} can be a matrix or data frame with three
#'   columns.
#' @param names Optional names.
#' @return An object of class \code{"tfn"}: a list with numeric components
#'   \code{l}, \code{m}, \code{r} and optional \code{names}.
#' @seealso \code{\link{tfn_score}}, \code{\link{alpha_cut}}, \code{\link{ftowa}}
#' @examples
#' a <- tfn(c(0.1, 0.3), c(0.5, 0.45), c(0.9, 0.5), names = c("tie1", "tie2"))
#' a
#' a + a
#' 0.5 * a
#' tfn_score(a, "cog")
#' alpha_cut(a, 0.5)
#' plot(a)
#' @export
#' @section Methodology:
#' A triangular fuzzy number (TFN) \eqn{\tilde A = (l, m, r)} is a normal, convex fuzzy subset of the real line with piecewise linear membership: it is fully possible at the core \eqn{m} and impossible outside the support \eqn{[l, r]}. In the package a TFN tie is read epistemically, as a possibility distribution on the ill-known strength of a crisp tie (Dubois and Prade, 1988): \eqn{\mu(w)} is the degree to which the value \eqn{w} is compatible with the available information.
#'
#' The \eqn{\alpha}-cut \eqn{[l + \alpha(m-l),\, r - \alpha(r-m)]} is the set of values with possibility at least \eqn{\alpha}; every computation of the exact methods is carried out on \eqn{\alpha}-cuts. The pseudo-operations of Porreca, Maturo and Ventre (2025) act on the three parameters. For sums and non-negative scalar multiples they coincide with the extension principle; the componentwise product coincides with it only at the support and at the core of non-negative TFNs, since the true product has curved sides.
#'
#'
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' # degenerate TFNs: crisp numbers and one-sided triangles
#' tfn(0.4, 0.4, 0.4)
#' tfn(0, 0, 0.3)
#' # from a matrix with columns l, m, r
#' tfn(cbind(c(0.1, 0.2), c(0.3, 0.5), c(0.6, 0.7)), names = c("x", "y"))
#' # subsetting, combining, conversion
#' a <- tfn(c(0.1, 0.3, 0.5), c(0.2, 0.4, 0.6), c(0.4, 0.5, 0.9), names = c("p", "q", "s"))
#' a[2:3]
#' c(a, tfn(0, 0.1, 0.2))
#' as.data.frame(a)
#' # the sum of TFNs is exact (extension principle); the mean is a scaled sum
#' tfn_sum(a, w = rep(1/3, 3))
tfn <- function(l, m, r, names = NULL) {
  if (missing(m) && (is.matrix(l) || is.data.frame(l))) {
    x <- as.matrix(l)
    if (ncol(x) != 3) stop("a matrix argument must have three columns (l, m, r)")
    if (is.null(names)) names <- rownames(x)
    l <- x[, 1]; m <- x[, 2]; r <- x[, 3]
  }
  if (any(!is.finite(c(l, m, r)))) stop("TFN parameters must be finite (no NA, NaN or Inf)")
  k <- max(length(l), length(m), length(r))
  l <- rep_len(as.numeric(l), k); m <- rep_len(as.numeric(m), k); r <- rep_len(as.numeric(r), k)
  tol <- 1e-10 * pmax(1, abs(m))
  if (any(l > m + tol | m > r + tol, na.rm = TRUE))
    stop("each TFN must satisfy l <= m <= r")
  structure(list(l = l, m = m, r = r, names = if (is.null(names)) NULL else as.character(names)),
            class = "tfn")
}

#' @rdname tfn
#' @param x An object.
#' @export
is.tfn <- function(x) inherits(x, "tfn")

#' @export
length.tfn <- function(x) length(unclass(x)$l)

#' @export
names.tfn <- function(x) unclass(x)$names

#' @export
`names<-.tfn` <- function(x, value) {
  y <- unclass(x); y$names <- if (is.null(value)) NULL else as.character(value)
  structure(y, class = "tfn")
}

#' @export
`[.tfn` <- function(x, i) {
  y <- unclass(x)
  if (is.character(i)) i <- match(i, y$names)
  tfn(y$l[i], y$m[i], y$r[i], names = if (is.null(y$names)) NULL else y$names[i])
}

#' @export
c.tfn <- function(...) {
  a <- list(...)
  nm <- unlist(lapply(a, function(z) { n <- names(z); if (is.null(n)) rep("", length(z)) else n }))
  tfn(unlist(lapply(a, function(z) unclass(z)$l)), unlist(lapply(a, function(z) unclass(z)$m)),
      unlist(lapply(a, function(z) unclass(z)$r)), names = if (all(nm == "")) NULL else nm)
}

#' @export
as.data.frame.tfn <- function(x, row.names = NULL, optional = FALSE, ...) {
  y <- unclass(x)
  data.frame(l = y$l, m = y$m, r = y$r, row.names = if (is.null(row.names)) y$names else row.names)
}

#' @export
as.matrix.tfn <- function(x, ...) {
  y <- unclass(x); out <- cbind(l = y$l, m = y$m, r = y$r); rownames(out) <- y$names; out
}

#' @export
format.tfn <- function(x, digits = 3, ...) {
  y <- unclass(x)
  out <- sprintf("(%s, %s, %s)", formatC(y$l, digits = digits, format = "f"),
                 formatC(y$m, digits = digits, format = "f"), formatC(y$r, digits = digits, format = "f"))
  names(out) <- y$names; out
}

#' @export
print.tfn <- function(x, digits = 3, ...) {
  cat("Triangular fuzzy numbers (l, m, r):", length(x), "\n")
  if (length(x)) print(as.data.frame(x), digits = digits)
  invisible(x)
}

#' @export
Ops.tfn <- function(e1, e2) {
  op <- .Generic
  if (op %in% c("+", "*")) {
    if (is.tfn(e1) && is.tfn(e2)) {
      a <- unclass(e1); b <- unclass(e2)
      f <- match.fun(op)
      return(tfn(f(a$l, b$l), f(a$m, b$m), f(a$r, b$r), names = a$names))
    }
    if (op == "*") {
      s <- if (is.tfn(e1)) e2 else e1; a <- unclass(if (is.tfn(e1)) e1 else e2)
      if (any(s < 0)) stop("scalar multiplication is defined here for non-negative scalars")
      return(tfn(s * a$l, s * a$m, s * a$r, names = a$names))
    }
    s <- if (is.tfn(e1)) e2 else e1; a <- unclass(if (is.tfn(e1)) e1 else e2)
    return(tfn(a$l + s, a$m + s, a$r + s, names = a$names))
  }
  if (op == "/" && is.tfn(e1) && is.numeric(e2)) {
    if (any(e2 <= 0)) stop("division by a positive scalar only")
    a <- unclass(e1); return(tfn(a$l / e2, a$m / e2, a$r / e2, names = a$names))
  }
  stop("operation '", op, "' is not defined for TFNs")
}

#' Sum of triangular fuzzy numbers
#'
#' @param x A \code{\link{tfn}} vector.
#' @param w Optional non-negative weights.
#' @return A TFN of length one, \eqn{\bigoplus_k w_k \odot x_k}.
#' @examples
#' tfn_sum(tfn(c(0, 1), c(1, 2), c(2, 3)))
#' @export
#' @section Methodology:
#' Computes \eqn{\bigoplus_k w_k \odot \tilde A_k = (\sum_k w_k l_k, \sum_k w_k m_k, \sum_k w_k r_k)}. Because the weighted sum is continuous and increasing in each argument for non-negative weights, this coincides with the extension principle: the \eqn{\alpha}-cut of the result is the weighted sum of the \eqn{\alpha}-cuts. This is the operation behind the fuzzy degree and strength, which are therefore exact.
#' @examples
#' ties <- tfn(c(0.2, 0.5, 0), c(0.4, 0.6, 0), c(0.5, 0.9, 0))
#' tfn_sum(ties)                    # in-strength of an actor receiving these ties
#' tfn_sum(ties, w = 1 / 3)         # mean, as in the fuzzy degree with weights 1/k
tfn_sum <- function(x, w = 1) {
  y <- unclass(x); w <- rep_len(w, length(y$l))
  if (any(w < 0)) stop("weights must be non-negative")
  tfn(sum(w * y$l), sum(w * y$m), sum(w * y$r))
}

#' Membership degree of real values in triangular fuzzy numbers
#'
#' @param x A \code{\link{tfn}} vector.
#' @param v Numeric vector of values.
#' @return A matrix with one row per TFN and one column per value.
#' @examples
#' tfn_membership(tfn(0, 0.5, 1), c(0.25, 0.5, 0.9))
#' @export
#' @section Methodology:
#' Evaluates the membership function \eqn{\mu(v) = (v-l)/(m-l)} on \eqn{[l, m]}, \eqn{(r-v)/(r-m)} on \eqn{[m, r]} and 0 elsewhere, with the conventions \eqn{\mu(m) = 1} for degenerate sides. Under the possibilistic reading, \eqn{\mu(v)} is the possibility that the latent tie strength equals \eqn{v}.
#' @examples
#' x <- tfn(c(0.2, 0.5), c(0.4, 0.5), c(0.8, 0.5), names = c("vague", "crisp"))
#' round(tfn_membership(x, seq(0, 1, by = 0.1)), 2)
#' curve(tfn_membership(tfn(0.2, 0.4, 0.8), x)[1, ], 0, 1, ylab = "membership")
tfn_membership <- function(x, v) {
  y <- unclass(x)
  out <- vapply(v, function(u) {
    up <- ifelse(y$m > y$l, (u - y$l) / (y$m - y$l), as.numeric(u >= y$l))
    dn <- ifelse(y$r > y$m, (y$r - u) / (y$r - y$m), as.numeric(u <= y$r))
    mu <- pmin(up, dn); mu[u < y$l | u > y$r] <- 0; pmax(0, pmin(1, mu))
  }, numeric(length(y$l)))
  matrix(out, nrow = length(y$l), dimnames = list(y$names, NULL))
}

#' Alpha-cuts of triangular fuzzy numbers
#'
#' The \eqn{\alpha}-cut of \eqn{(l, m, r)} is
#' \eqn{[l + \alpha(m - l),\, r - \alpha(r - m)]}; for \eqn{\alpha = 0} it is the
#' closure of the support.
#'
#' @param x A \code{\link{tfn}} vector.
#' @param alpha Numeric vector of levels in \eqn{[0, 1]}.
#' @return A list with matrices \code{lower} and \code{upper} (one row per TFN,
#'   one column per level).
#' @examples
#' alpha_cut(tfn(0, 0.5, 1), c(0, 0.5, 1))
#' @export
#' @section Methodology:
#' The \eqn{\alpha}-cut \eqn{A_\alpha = \{v : \mu(v) \ge \alpha\}} of a TFN is the closed interval \eqn{[l + \alpha(m-l),\, r - \alpha(r-m)]}; for \eqn{\alpha = 0} the closure of the support is returned. The family of \eqn{\alpha}-cuts characterizes a fuzzy number (Negoita and Ralescu's representation), and the extension principle commutes with \eqn{\alpha}-cuts for continuous functions (Nguyen, 1978). This is the basis of \code{\link{fuzzy_centrality}}.
#' @examples
#' x <- tfn(0.2, 0.5, 0.6)
#' ct <- alpha_cut(x, seq(0, 1, by = 0.25))
#' cbind(alpha = ct$alpha, lower = ct$lower[1, ], upper = ct$upper[1, ])
#' # nested cuts: widths decrease with alpha
#' ct$upper - ct$lower
alpha_cut <- function(x, alpha) {
  if (any(alpha < 0 | alpha > 1)) stop("alpha must lie in [0, 1]")
  y <- unclass(x)
  lo <- outer(y$l, rep(1, length(alpha))) + outer(y$m - y$l, alpha)
  up <- outer(y$r, rep(1, length(alpha))) - outer(y$r - y$m, alpha)
  dimnames(lo) <- dimnames(up) <- list(y$names, format(alpha))
  list(lower = lo, upper = up, alpha = alpha)
}

#' Ranking scores of fuzzy numbers
#'
#' Real-valued indices used to rank fuzzy numbers. For a TFN \eqn{(l, m, r)}:
#' \describe{
#'   \item{\code{"cog"}}{centre of gravity (centroid) of the membership function,
#'     \eqn{(l + m + r)/3}; the criterion used by Porreca, Maturo and Ventre (2025).}
#'   \item{\code{"core"}}{the core \eqn{m}.}
#'   \item{\code{"yager"}}{Yager's index, the mean over \eqn{\alpha} of the
#'     midpoints of the \eqn{\alpha}-cuts, \eqn{(l + 2m + r)/4} (Yager, 1981).}
#'   \item{\code{"graded_mean"}}{graded mean integration representation
#'     \eqn{(l + 4m + r)/6} (Chen and Hsieh, 1999); it coincides with the
#'     possibilistic mean of Carlsson and Fuller (2001),
#'     \eqn{\int_0^1 \alpha (a^-_\alpha + a^+_\alpha) d\alpha}.}
#'   \item{\code{"liou_wang"}}{total integral value of Liou and Wang (1992) with
#'     optimism index \code{lambda}: \eqn{(1-\lambda)(l+m)/2 + \lambda(m+r)/2}.}
#' }
#' The same definitions are applied to general fuzzy numbers given by their
#' \eqn{\alpha}-cuts in \code{\link{fuzzy_centrality}} objects, through numerical
#' integration over \eqn{\alpha}.
#'
#' @param x A \code{\link{tfn}} vector.
#' @param method One of \code{"cog"}, \code{"core"}, \code{"yager"},
#'   \code{"graded_mean"}, \code{"liou_wang"}.
#' @param lambda Optimism index in \eqn{[0, 1]} for \code{"liou_wang"}.
#' @return Numeric vector of scores (larger = ranked higher).
#' @references
#' Yager R.R. (1981). A procedure for ordering fuzzy subsets of the unit interval.
#' \emph{Information Sciences} 24, 143-161.
#'
#' Liou T.-S., Wang M.-J.J. (1992). Ranking fuzzy numbers with integral value.
#' \emph{Fuzzy Sets and Systems} 50, 247-255.
#'
#' Carlsson C., Fuller R. (2001). On possibilistic mean value and variance of
#' fuzzy numbers. \emph{Fuzzy Sets and Systems} 122, 315-326.
#' @examples
#' x <- tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32))
#' sapply(fuzzy_rank_methods(), function(m) tfn_score(x, m))
#' @export
#' @section Methodology:
#' Ranking fuzzy numbers requires a criterion, and different criteria can order the same TFNs differently. Score-based methods map each fuzzy number to a real value. The centre of gravity weights each value by its membership; Yager's index averages the midpoints of the \eqn{\alpha}-cuts with uniform weight on \eqn{\alpha}; the graded mean (possibilistic mean) weights the cuts by \eqn{\alpha}, giving more importance to the most plausible values; the Liou-Wang index weights the right and left parts by an optimism index \eqn{\lambda}; the core ignores vagueness altogether. All are linear in \eqn{(l, m, r)} for TFNs, with weights \eqn{(1,1,1)/3}, \eqn{(1,2,1)/4}, \eqn{(1,4,1)/6}, \eqn{((1-\lambda), 1, \lambda)/2} and \eqn{(0,1,0)}.
#'
#' A score orders all pairs, including pairs whose fuzzy numbers overlap and whose ordering is therefore not supported by the data. The comparison-based methods of \code{\link{rank_fuzzy}} make this uncertainty explicit.
#' @examples
#' # two TFNs ordered differently by different scores
#' x <- tfn(c(0, 0.3), c(0.2, 0.31), c(1, 0.32), names = c("vague", "precise"))
#' s <- sapply(fuzzy_rank_methods(), function(m) tfn_score(x, m))
#' s
#' apply(s, 2, function(v) names(which.max(v)))
#' # optimism index of Liou and Wang
#' sapply(c(0, 0.5, 1), function(l) tfn_score(x, "liou_wang", lambda = l))
tfn_score <- function(x, method = "cog", lambda = 0.5) {
  method <- match.arg(method, fuzzy_rank_methods())
  .check_lambda(lambda)
  y <- unclass(x)
  out <- switch(method,
                cog = (y$l + y$m + y$r) / 3,
                core = y$m,
                yager = (y$l + 2 * y$m + y$r) / 4,
                graded_mean = (y$l + 4 * y$m + y$r) / 6,
                liou_wang = (1 - lambda) * (y$l + y$m) / 2 + lambda * (y$m + y$r) / 2)
  names(out) <- y$names; out
}

#' Names of the available score-based ranking methods
#' @return Character vector.
#' @examples
#' fuzzy_rank_methods()
#' @export
#' @section Methodology:
#' Lists the score-based methods available for TFNs and for exact fuzzy centralities. The comparison-based methods \code{"possibility"}, \code{"necessity"}, \code{"certified"} and \code{"top"}, which require a \code{\link{fuzzy_centrality}} object, are documented in \code{\link{rank_fuzzy}}. See Section 2.4 of Porreca and Maturo (2026).
#' @examples
#' x <- tfn(c(0.1, 0.2, 0.3), c(0.5, 0.4, 0.45), c(0.6, 0.9, 0.5))
#' sapply(fuzzy_rank_methods(), function(m) rank(-tfn_score(x, m)))
fuzzy_rank_methods <- function() c("cog", "core", "yager", "graded_mean", "liou_wang")

# score matrix of the ties of a network from three matrices L, M, R
.check_lambda <- function(lambda)
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) || lambda < 0 || lambda > 1)
    stop("lambda must be a single number in [0, 1]")

.score_matrix <- function(L, M, R, method = "cog", lambda = 0.5) {
  method <- match.arg(method, fuzzy_rank_methods())
  .check_lambda(lambda)
  switch(method,
         cog = (L + M + R) / 3, core = M, yager = (L + 2 * M + R) / 4,
         graded_mean = (L + 4 * M + R) / 6,
         liou_wang = (1 - lambda) * (L + M) / 2 + lambda * (M + R) / 2)
}

#' Spread and vagueness measures of triangular fuzzy numbers
#'
#' @param x A \code{\link{tfn}} vector.
#' @return A data frame with the support width \eqn{r - l}, the left and right
#'   spreads, the asymmetry \eqn{(r - m) - (m - l)} and the ambiguity
#'   \eqn{\int_0^1 \alpha (a^+_\alpha - a^-_\alpha) d\alpha = (r - l)/6}
#'   (Delgado, Vila and Voxman, 1998).
#' @examples
#' tfn_spread(tfn(c(0, 0.2), c(0.5, 0.3), c(1, 0.9)))
#' @export
#' @section Methodology:
#' The width \eqn{r - l} measures the overall vagueness of a rating; the asymmetry \eqn{(r-m)-(m-l)} is positive when the plausible values extend more to the right of the core. In mouse-tracking data the core is the final position of the pointer, which rarely lies at the midpoint of the explored range, so ratings are typically asymmetric. The ambiguity \eqn{\int_0^1 \alpha\,(a^+_\alpha - a^-_\alpha)\,d\alpha} weights the width of each cut by its level (Delgado, Vila and Voxman, 1998); for a TFN it equals one sixth of the width.
#' @examples
#' ties <- fsn_ties(fsn_simulated(n = 20, respondents = 10))
#' sp <- tfn_spread(tfn(ties$l, ties$m, ties$r))
#' summary(sp)
#' hist(sp$asymmetry, main = "Asymmetry of the ratings", xlab = "(r - m) - (m - l)")
tfn_spread <- function(x) {
  y <- unclass(x)
  data.frame(width = y$r - y$l, left_spread = y$m - y$l, right_spread = y$r - y$m,
             asymmetry = (y$r - y$m) - (y$m - y$l), ambiguity = (y$r - y$l) / 6,
             row.names = y$names)
}

#' Plot triangular fuzzy numbers
#'
#' Draws the membership functions.
#'
#' @param x A \code{\link{tfn}} vector.
#' @param col,lty Colours and line types (recycled).
#' @param labels Logical, label the cores with the names.
#' @param xlim,xlab,ylab,main Graphical parameters.
#' @param legend Logical, add a legend when names are available.
#' @param ... Further arguments passed to \code{\link[graphics]{plot.default}}.
#' @return \code{x}, invisibly.
#' @examples
#' plot(tfn(c(0.1, 0.3), c(0.5, 0.45), c(0.9, 0.5), names = c("a", "b")))
#' @export
#' @section Methodology:
#' Draws the piecewise linear membership functions. Overlapping triangles indicate pairs of values whose ordering is not determined by the information: see \code{\link{possibility_matrix}}.
#' @examples
#' x <- tfn(c(0.1, 0.3, 0.2), c(0.5, 0.45, 0.3), c(0.9, 0.5, 0.35), names = c("a", "b", "c"))
#' plot(x, labels = TRUE, main = "Three fuzzy ratings")
#' plot(x, col = "grey30", lty = 1:3, legend = FALSE)
plot.tfn <- function(x, col = NULL, lty = 1, labels = FALSE, xlim = NULL, xlab = "Value",
                     ylab = "Membership", main = "", legend = TRUE, ...) {
  y <- unclass(x); k <- length(y$l)
  if (is.null(col)) col <- fsna_palette(k)
  col <- rep_len(col, k); lty <- rep_len(lty, k)
  if (is.null(xlim)) { rg <- range(c(y$l, y$r)); xlim <- rg + c(-1, 1) * 0.05 * max(diff(rg), 1e-8) }
  graphics::plot(NA, xlim = xlim, ylim = c(0, 1.05), xlab = xlab, ylab = ylab, main = main, las = 1, ...)
  for (i in seq_len(k)) graphics::lines(c(y$l[i], y$m[i], y$r[i]), c(0, 1, 0), col = col[i], lty = lty[i], lwd = 1.6)
  if (labels && !is.null(y$names)) graphics::text(y$m, 1, y$names, pos = 3, cex = 0.7, col = col)
  if (legend && !is.null(y$names) && k <= 15)
    graphics::legend("topright", legend = y$names, col = col, lty = lty, lwd = 1.6, bty = "n", cex = 0.8)
  invisible(x)
}

#' Colour palette used by the package
#'
#' Colour-blind safe palette (Okabe-Ito), interpolated when more colours are
#' needed.
#' @param n Number of colours.
#' @return Character vector of colours.
#' @examples
#' fsna_palette(5)
#' @export
#' @section Methodology:
#' Okabe-Ito colour-blind safe palette used throughout the graphics of the package, interpolated when more than eight colours are requested.
#' @examples
#' pal <- fsna_palette(12)
#' barplot(rep(1, 12), col = pal, border = NA, axes = FALSE)
fsna_palette <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base[1:7])(n)
}
