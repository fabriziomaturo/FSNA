# ---------------------------------------------------------------------------
# Graphics (base R)
# ---------------------------------------------------------------------------

.alpha_col <- function(col, a) { k <- max(length(col), length(a)); mapply(function(cc, aa) grDevices::adjustcolor(cc, alpha.f = aa), rep_len(col, k), rep_len(a, k), USE.NAMES = FALSE) }

#' Plot a fuzzy social network
#'
#' Draws the actors on a circle (or with a user layout) and the non-null ties as
#' arrows whose width is proportional to the core and whose transparency
#' increases with the width of the support.
#' @param x An \code{\link{fsn}} object.
#' @param layout Optional two-column matrix of coordinates.
#' @param vertex_col Colours of the actors (default: respondents dark, others light).
#' @param vertex_size Vector of sizes (e.g. a centrality); rescaled.
#' @param min_core Ties with core below this value are not drawn.
#' @param main Title.
#' @param label_cex Size of the labels.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' plot(fsn_example())
#' A <- fsn_simulated(n = 20, respondents = 10)
#' plot(A, vertex_size = colSums(A$M), min_core = 0.6)
#' @export
#' @section Methodology:
#' The width of an arrow is proportional to the core of the tie and its transparency increases with the width of the support, so that vague ties are drawn lighter.
#' @examples
#' B <- fsn_simulated(n = 20, respondents = 10, seed = 2)
#' plot(B, vertex_size = fc_score(fuzzy_centrality(B, "in_strength")), min_core = 0.5,
#'      main = "Simulated network, ties with core above 0.5")
plot.fsn <- function(x, layout = NULL, vertex_col = NULL, vertex_size = NULL, min_core = 0,
                     main = "", label_cex = 0.7, ...) {
  n <- n_actors(x)
  if (is.null(layout)) { th <- 2 * pi * (seq_len(n) - 1) / n; layout <- cbind(cos(th), sin(th)) }
  if (is.null(vertex_col)) vertex_col <- if (is.null(x$respondents)) "#0072B2" else ifelse(x$respondents, "#0072B2", "#9ecae1")
  vs <- if (is.null(vertex_size)) rep(1.6, n) else 0.8 + 2.4 * (vertex_size - min(vertex_size)) / max(diff(range(vertex_size)), 1e-12)
  op <- graphics::par(mar = c(0.5, 0.5, if (nzchar(main)) 2 else 0.5, 0.5)); on.exit(graphics::par(op))
  graphics::plot(layout, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, main = main,
                 xlim = range(layout[, 1]) * 1.12, ylim = range(layout[, 2]) * 1.12)
  E <- which(x$M > min_core & row(x$M) != col(x$M), arr.ind = TRUE)
  if (nrow(E)) {
    w <- x$M[E]; spr <- (x$R - x$L)[E]; al <- 0.85 - 0.6 * spr / max(spr, 1e-12)
    p0 <- layout[E[, 1], , drop = FALSE]; p1 <- layout[E[, 2], , drop = FALSE]
    sh <- 0.06; d <- p1 - p0; len <- sqrt(rowSums(d^2)); d <- d / pmax(len, 1e-12)
    graphics::arrows(p0[, 1] + sh * d[, 1], p0[, 2] + sh * d[, 2], p1[, 1] - sh * d[, 1], p1[, 2] - sh * d[, 2],
                     length = 0.05, lwd = 0.4 + 2.5 * w, col = .alpha_col("grey30", al))
  }
  graphics::points(layout, pch = 21, bg = vertex_col, col = "white", cex = vs)
  graphics::text(layout * 1.07, labels = x$labels, cex = label_cex)
  invisible(x)
}

#' Plot exact fuzzy centralities
#'
#' \code{type = "membership"} draws the membership functions of the exact fuzzy
#' centralities of the selected actors (solid) and of their triangular
#' approximations (dashed). \code{type = "intervals"} draws, for all actors
#' ordered by a ranking score, the support (thin), the \eqn{\alpha}-cut at
#' \code{alpha} (thick) and the core (point, or short bar when the core is an
#' interval because of unknown ties), with the actors in the
#' \code{alpha}-cut of the fuzzy set of the most central actor highlighted.
#' \code{type = "profile"} draws the width of the \eqn{\alpha}-cuts against
#' \eqn{\alpha}.
#' @param x A \code{\link{fuzzy_centrality}} object.
#' @param type \code{"membership"}, \code{"intervals"} or \code{"profile"}.
#' @param actors Actors to show: labels, indices, or a number \eqn{k} (the
#'   \eqn{k} actors with the highest score). Default: 6 for membership plots, all
#'   otherwise.
#' @param triangle Logical, add the triangular approximation (membership plot).
#' @param alpha Level of the thick segments (interval plot).
#' @param method Ranking score used to order and select actors.
#' @param col Colours.
#' @param main,xlab Titles.
#' @param legend Logical.
#' @param ... Further graphical parameters.
#' @return \code{x}, invisibly.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_closeness_prod")
#' plot(fc)
#' plot(fc, type = "intervals")
#' plot(fc, type = "profile", actors = 5)
#' @export
#' @section Methodology:
#' Membership plots show the exact membership functions; dashed triangles show the triangular approximation, which coincides with the exact function at the support and at the core. The interval plot shows supports, \eqn{\alpha}-cuts and cores, and highlights the actors in the \eqn{\alpha}-cut of the fuzzy set of the most central actor. The profile plot shows how quickly vagueness decreases with the level.
#' @examples
#' fc <- fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_katz")
#' top3 <- names(sort(fc$lower[, ncol(fc$lower)], decreasing = TRUE))[1:3]
#' plot(fc, actors = top3, col = c("black", "grey40", "grey70"))
#' plot(fc, "intervals", alpha = 0.75, method = "graded_mean")
#' plot(fc, "profile", actors = 8)
plot.fuzzy_centrality <- function(x, type = c("membership", "intervals", "profile"), actors = NULL,
                                  triangle = TRUE, alpha = 0.5, method = "cog", col = NULL,
                                  main = NULL, xlab = NULL, legend = TRUE, ...) {
  type <- match.arg(type)
  sc <- fc_score(x, method)
  if (is.null(actors)) actors <- if (type == "intervals") length(x$labels) else 6
  sel <- if (is.numeric(actors) && length(actors) == 1) order(-sc)[seq_len(min(actors, length(sc)))]
         else if (is.character(actors)) match(actors, x$labels) else actors
  if (is.null(xlab)) xlab <- x$index$label
  if (is.null(main)) main <- ""
  k <- length(sel); if (is.null(col)) col <- fsna_palette(k)
  G <- length(x$alpha)
  if (type == "membership") {
    xl <- range(c(x$lower[sel, 1], x$upper[sel, 1])); xl <- xl + c(-1, 1) * 0.04 * max(diff(xl), 1e-8)
    graphics::plot(NA, xlim = xl, ylim = c(0, 1.04), xlab = xlab, ylab = "Membership", main = main, las = 1, ...)
    for (h in seq_len(k)) {
      i <- sel[h]
      graphics::lines(c(x$lower[i, ], rev(x$upper[i, ])), c(x$alpha, rev(x$alpha)), col = col[h], lwd = 1.8)
      if (triangle && (!x$index$affine || isTRUE(x$interval_core))) {
        t <- unclass(x$triangle)
        graphics::lines(c(t$l[i], t$m[i], t$r[i]), c(0, 1, 0), col = col[h], lty = 2)
      }
    }
    if (legend) graphics::legend("topright", legend = x$labels[sel], col = col, lwd = 1.8, bty = "n", cex = 0.8)
  } else if (type == "intervals") {
    g <- .alpha_index(x, alpha); o <- sel[order(sc[sel])]; y <- seq_along(o)
    top <- x$labels %in% top_set(x, alpha = alpha)
    graphics::plot(NA, xlim = range(c(x$lower[o, 1], x$upper[o, 1])), ylim = c(0.5, length(o) + 0.5),
                   yaxt = "n", xlab = xlab, ylab = "", main = main, ...)
    graphics::axis(2, at = y, labels = x$labels[o], las = 1, cex.axis = max(0.4, min(0.8, 25 / length(o))))
    cl <- ifelse(top[o], "#D55E00", "#0072B2")
    graphics::segments(x$lower[o, 1], y, x$upper[o, 1], y, col = .alpha_col(cl, 0.5), lwd = 1)
    graphics::segments(x$lower[o, g], y, x$upper[o, g], y, col = cl, lwd = 3)
    graphics::points((x$lower[o, G] + x$upper[o, G]) / 2, y, pch = 19, cex = 0.6)
    ic <- (x$upper[o, G] - x$lower[o, G]) > 1e-10
    if (any(ic)) graphics::segments(x$lower[o, G][ic], y[ic], x$upper[o, G][ic], y[ic], lwd = 5, col = "black")
    if (legend) graphics::legend("bottomright", bty = "n", cex = 0.75, lwd = c(1, 3, 3), col = c("grey50", "#0072B2", "#D55E00"),
                                 legend = c("support", paste0(alpha, "-cut"), paste0("in the ", alpha, "-cut of the top set")))
  } else {
    W <- x$upper[sel, , drop = FALSE] - x$lower[sel, , drop = FALSE]
    graphics::matplot(x$alpha, t(W), type = "l", lty = 1, lwd = 1.6, col = col, xlab = expression(alpha),
                      ylab = "Width of the alpha-cut", main = main, las = 1, ...)
    if (legend) graphics::legend("topright", legend = x$labels[sel], col = col, lwd = 1.6, bty = "n", cex = 0.8)
  }
  invisible(x)
}

#' Heat map of possibility or necessity degrees
#' @param x A \code{"fuzzy_comparison"} object (\code{\link{possibility_matrix}}).
#' @param which \code{"necessity"} or \code{"possibility"}.
#' @param order Optional ordering of the actors (default: by row sums).
#' @param main Title.
#' @param ... Further arguments to \code{\link[graphics]{image}}.
#' @return \code{x}, invisibly.
#' @examples
#' cmp <- possibility_matrix(fuzzy_centrality(fsn_simulated(n = 20, respondents = 10), "in_strength"))
#' plot(cmp)
#' @export
#' @section Methodology:
#' Rows and columns are ordered by the row sums of the degrees, so that a nearly triangular pattern indicates a well-determined ranking and grey blocks near the diagonal indicate groups of actors that cannot be ordered. See Section 2.4 of Porreca and Maturo (2026).
#' @examples
#' x2 <- fsn_simulated(n = 20, respondents = 10, seed = 2)
#' cmp <- possibility_matrix(fuzzy_centrality(x2, "in_closeness_prod"))
#' plot(cmp, "necessity")
#' plot(cmp, "possibility")
plot.fuzzy_comparison <- function(x, which = c("necessity", "possibility"), order = NULL, main = NULL, ...) {
  which <- match.arg(which); A <- x[[which]]
  if (is.null(order)) order <- base::order(-rowSums(A))
  A <- A[order, order]; n <- nrow(A)
  if (is.null(main)) main <- if (which == "necessity") "N(row > column)" else "Pi(row >= column)"
  op <- graphics::par(mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
  graphics::image(seq_len(n), seq_len(n), t(A[n:1, ]), col = grDevices::hcl.colors(50, "Blues 3", rev = TRUE),
                  zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = main, ...)
  cx <- max(0.35, min(0.8, 20 / n))
  graphics::axis(1, at = seq_len(n), labels = colnames(A), las = 2, cex.axis = cx)
  graphics::axis(2, at = seq_len(n), labels = rev(rownames(A)), las = 1, cex.axis = cx)
  graphics::box()
  invisible(x)
}

#' Fuzzy scatter plot of two centralities
#'
#' Each actor is drawn as an ellipse centred at the scores of two fuzzy
#' centralities, with semi-axes equal to half the widths of the supports (or of
#' the \eqn{\alpha}-cuts), as in Figs. 10 and 11 of Porreca, Maturo and Ventre
#' (2025). The colour may encode a third, crisp, variable (e.g. betweenness).
#' The ellipses are descriptive glyphs of vagueness: they are not confidence
#' regions, nor joint \eqn{\alpha}-cuts of the two centralities, and, being
#' centred at a score, they do not reproduce the endpoints of asymmetric cuts.
#' The two centralities are matched by actor label; \code{x} and \code{y} must
#' refer to the same actors, and a named \code{colour} is matched by label too.
#' @param x,y \code{\link{fuzzy_centrality}} objects, \code{"ijar_index"} objects
#'   or \code{\link{tfn}} vectors for the same actors.
#' @param colour Optional numeric vector for the colour scale.
#' @param alpha Level of the cuts used for the semi-axes (fuzzy centralities only).
#' @param method Score for the centres.
#' @param labels Logical, label the actors.
#' @param xlab,ylab,main Titles.
#' @param colour_label Legend title for the colour scale.
#' @param ... Further graphical parameters.
#' @return Invisibly, the data frame of centres and semi-axes.
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' fi <- fuzzy_centrality(A, "in_strength", normalize = "n-1")
#' fo <- fuzzy_centrality(A, "out_strength", normalize = "n-1")
#' b <- ijar_betweenness(A)
#' plot_fuzzy_scatter(fi, fo, colour = b$betweenness, colour_label = "Betweenness")
#' @export
#' @section Methodology:
#' Draws fuzzy scatter plots of the kind used by Porreca, Maturo and Ventre (2025, Figs. 10 and 11). The semi-axes of an ellipse are half the widths of the supports (or of the \eqn{\alpha}-cuts); large ellipses indicate actors whose centralities are vague in one or both dimensions.
#' @examples
#' x <- simulate_fsn(30, 0.25, seed = 9)
#' fi <- fuzzy_centrality(x, "in_closeness_prod"); fo <- fuzzy_centrality(x, "out_closeness_prod")
#' plot_fuzzy_scatter(fi, fo, alpha = 0.5, xlab = "in-closeness", ylab = "out-closeness")
#' # IJAR indices
#' A <- fsn_simulated(n = 20, respondents = 10)
#' plot_fuzzy_scatter(ijar_closeness(A, "in"), ijar_closeness(A, "out"),
#'                    colour = ijar_betweenness(A)$betweenness, labels = FALSE)
plot_fuzzy_scatter <- function(x, y, colour = NULL, alpha = 0, method = "cog", labels = TRUE,
                               xlab = "x", ylab = "y", main = "", colour_label = "", ...) {
  get <- function(z) {
    if (inherits(z, "ijar_index")) z <- as_tfn(z)
    if (is.tfn(z)) { u <- unclass(z); return(list(c = tfn_score(z, method), h = (u$r - u$l) / 2, lab = u$names)) }
    g <- .alpha_index(z, alpha)
    list(c = fc_score(z, method), h = (z$upper[, g] - z$lower[, g]) / 2, lab = z$labels)
  }
  a <- get(x); b <- get(y)
  # align the two axes by actor: both must identify the same set of actors
  if (!is.null(a$lab) && !is.null(b$lab)) {
    if (anyDuplicated(a$lab) || anyDuplicated(b$lab)) stop("duplicated actor identifiers in 'x' or 'y'")
    if (!setequal(a$lab, b$lab)) stop("'x' and 'y' do not refer to the same actors: ",
                                      paste(utils::head(c(setdiff(a$lab, b$lab), setdiff(b$lab, a$lab)), 5), collapse = ", "))
    o <- match(a$lab, b$lab); b$c <- b$c[o]; b$h <- b$h[o]; b$lab <- b$lab[o]
  } else {
    if (length(a$c) != length(b$c)) stop("'x' and 'y' have different lengths")
    warning("'x' or 'y' has no actor identifiers: values are matched by position")
    if (is.null(a$lab)) a$lab <- if (!is.null(b$lab)) b$lab else as.character(seq_along(a$c))
  }
  n <- length(a$c)
  if (!is.null(colour)) {
    if (!is.null(names(colour)) && !is.null(a$lab)) {
      if (!all(a$lab %in% names(colour))) stop("'colour' is named but does not cover all actors")
      colour <- colour[a$lab]
    } else if (length(colour) != n) stop("'colour' must have one value per actor")
  }
  cols <- if (is.null(colour)) rep("#0072B2", n) else {
    pal <- grDevices::colorRampPalette(c("#2166AC", "#F4A582", "#D55E00"))(100)
    pal[1 + round(99 * (colour - min(colour)) / max(diff(range(colour)), 1e-12))]
  }
  graphics::plot(NA, xlim = range(c(a$c - a$h, a$c + a$h)), ylim = range(c(b$c - b$h, b$c + b$h)),
                 xlab = xlab, ylab = ylab, main = main, las = 1, ...)
  th <- seq(0, 2 * pi, length.out = 60)
  for (i in seq_len(n)) if (a$h[i] > 0 || b$h[i] > 0)
    graphics::polygon(a$c[i] + a$h[i] * cos(th), b$c[i] + b$h[i] * sin(th), col = .alpha_col(cols[i], 0.25), border = NA)
  graphics::points(a$c, b$c, pch = 19, col = cols, cex = 0.8)
  if (labels) graphics::text(a$c, b$c, a$lab, pos = 3, cex = 0.6, col = cols)
  if (!is.null(colour)) graphics::legend("bottomright", bty = "n", cex = 0.75, title = colour_label,
                                         legend = signif(range(colour), 2), pch = 19, col = cols[c(which.min(colour), which.max(colour))])
  invisible(data.frame(actor = a$lab, x = a$c, y = b$c, x_half_width = a$h, y_half_width = b$h))
}

#' Plot the TFN centralities of Porreca, Maturo and Ventre (2025)
#' @param x An \code{"ijar_index"} object.
#' @param n Number of actors (highest scores) to draw; \code{NULL} for all.
#' @param ... Passed to \code{\link{plot.tfn}}.
#' @return \code{x}, invisibly.
#' @examples
#' plot(ijar_degree(fsn_simulated(n = 20, respondents = 10), "in"), n = 8)
#' @export
#' @section Methodology:
#' Draws the TFN centralities returned by the indices of Porreca, Maturo and Ventre (2025), as in their Figs. 5, 6, 8 and 9. See Section 3 of Porreca, Maturo and Ventre (2025).
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' A <- fsn_simulated(n = 20, respondents = 10)
#' plot(ijar_closeness(A, "out"), n = 6)
plot.ijar_index <- function(x, n = 10, ...) {
  o <- order(x$rank); if (!is.null(n)) o <- o[seq_len(min(n, length(o)))]
  plot.tfn(as_tfn(x)[o], xlab = attr(x, "index"), ...)
  invisible(x)
}

#' @export
plot.ijar_betweenness <- function(x, n = NULL, main = "", ...) {
  o <- order(-x$betweenness); if (!is.null(n)) o <- o[seq_len(min(n, length(o)))]
  graphics::plot(seq_along(o), x$betweenness[o], type = "h", lwd = 2, col = "#0072B2", xaxt = "n",
                 xlab = "Actors", ylab = "Betweenness on best paths", main = main, las = 1, ...)
  graphics::axis(1, at = seq_along(o), labels = x$actor[o], las = 2, cex.axis = max(0.4, min(0.8, 25 / length(o))))
  invisible(x)
}

# small person icon centred at (x0, y0) with height s
.draw_person <- function(x0, y0, s, col) {
  th <- seq(0, 2 * pi, length.out = 40)
  hr <- 0.17 * s
  graphics::polygon(x0 + hr * cos(th), y0 + 0.22 * s + hr * sin(th), col = col, border = "white", lwd = 0.8)
  tb <- seq(0, pi, length.out = 30)
  bx <- c(x0 + 0.32 * s * cos(tb), x0 - 0.32 * s, x0 + 0.32 * s)
  by <- c(y0 - 0.02 * s + 0.30 * s * sin(tb), y0 - 0.42 * s, y0 - 0.42 * s)
  graphics::polygon(c(x0 + 0.32 * s * cos(tb), x0 - 0.32 * s, x0 + 0.32 * s),
                    c(y0 - 0.12 * s + 0.26 * s * sin(tb), y0 - 0.40 * s, y0 - 0.40 * s),
                    col = col, border = "white", lwd = 0.8)
  invisible(NULL)
}

#' Network plot with person icons and triangular fuzzy ties
#'
#' Draws a small fuzzy social network as in Fig. 3 of Porreca, Maturo and
#' Ventre (2025): actors are person icons and every tie is an arrow carrying a
#' small panel with the membership function of its triangular fuzzy number on
#' the scale of the ties. Reciprocated ties are drawn as two curved arrows, so
#' that the asymmetry of the ratings is visible. Unknown ties
#' (\code{\link{set_unknown}}) are drawn as dashed arrows with a question mark.
#'
#' @param x An \code{\link{fsn}} object, preferably with at most 12-15 actors.
#' @param actors Optional subset of actors (labels or indices).
#' @param layout Optional two-column matrix of coordinates; default a circle.
#' @param min_core Ties with core below this value are omitted.
#' @param glyph_size Width of the tie panels, in units of the layout radius.
#' @param person_col Colours of the actors (default: respondents dark blue,
#'   others light blue).
#' @param fill Fill colour of the triangles.
#' @param values Logical: write \eqn{(l, m, r)} under each triangle.
#' @param digits Digits of the values.
#' @param curvature Bending of reciprocated ties.
#' @param glyph_pos Position of the tie panels along the arrows, from 0
#'   (sender) to 1 (receiver); values below 0.5 keep the panels of reciprocated
#'   ties apart.
#' @param scale Range of the tie values used for the panels.
#' @param main Title.
#' @param label_cex Size of the actor labels.
#' @return \code{x}, invisibly.
#' @section Methodology:
#' Each panel shows the possibility distribution of the latent strength of a
#' tie on the common scale \code{scale}: its position along the horizontal axis
#' shows the core and its width the vagueness of the rating; asymmetric
#' triangles reveal on which side of the final answer the respondent hesitated.
#' The plot is meant for small networks, to explain the data structure and
#' illustrative examples such as weak direct ties bypassed by strong indirect
#' routes (best paths, \code{\link{ijar_best_paths}}). For larger networks use
#' \code{\link{plot.fsn}}.
#' @references Porreca A., Maturo F. (2026). FSNA: an R package for
#'   centrality analysis in social networks with triangular fuzzy ties.
#'   Manuscript; see \code{citation("FSNA")}.
#'
#'   Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in
#'   social network analysis: Theory and application in a university department
#'   collaboration network. \emph{International Journal of Approximate
#'   Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @seealso \code{\link{plot.fsn}}, \code{\link{plot.tfn}}
#' @examples
#' # Five researchers: each tie is a fuzzy rating (l, m, r) on [0, 1]. The example
#' # contains asymmetric reciprocated ties (Anna-Bruno), a weak direct tie
#' # (Anna -> Dario) bypassed by a strong indirect route (Anna -> Bruno -> Carla
#' # -> Dario) and an asymmetric rating (Elena -> Dario).
#' ratings <- data.frame(
#'   from = c("Anna", "Bruno", "Bruno", "Carla", "Dario", "Dario", "Elena", "Anna", "Elena"),
#'   to   = c("Bruno", "Anna", "Carla", "Dario", "Carla", "Elena", "Anna", "Dario", "Dario"),
#'   l    = c(0.70, 0.35, 0.60, 0.70, 0.40, 0.55, 0.60, 0.05, 0.20),
#'   m    = c(0.85, 0.60, 0.75, 0.85, 0.50, 0.80, 0.75, 0.15, 0.30),
#'   r    = c(0.95, 0.70, 0.85, 0.90, 0.80, 0.90, 0.95, 0.35, 0.70))
#' x5 <- fsn_from_ratings(ratings, actors = c("Anna", "Bruno", "Carla", "Dario", "Elena"),
#'                        drop_isolated = FALSE)
#' th <- pi / 2 - 2 * pi * (0:4) / 5          # regular pentagon, Anna on top
#' plot_fsn_glyphs(x5, layout = cbind(cos(th), sin(th)), values = TRUE,
#'                 glyph_size = 0.19, curvature = 0.25, glyph_pos = 0.3, label_cex = 1.1,
#'                 main = "Collaboration among five researchers")
#' # for a publication: a square device
#' \donttest{
#' pdf(file.path(tempdir(), "five_researchers.pdf"), width = 7, height = 7)
#' plot_fsn_glyphs(x5, layout = cbind(cos(th), sin(th)), values = TRUE,
#'                 glyph_size = 0.19, curvature = 0.25, glyph_pos = 0.3, label_cex = 1.1)
#' dev.off()
#' }
#'
#' plot_fsn_glyphs(fsn_example(), main = "Example network")
#' plot_fsn_glyphs(fsn_example(), values = TRUE, min_core = 0.4)
#' # three actors with a weak direct tie and a strong indirect route
#' M <- matrix(0, 3, 3, dimnames = list(c("Perry", "Smith", "Porreca"), NULL))
#' M[1, 3] <- 0.15; M[1, 2] <- 0.85; M[2, 3] <- 0.8; M[3, 2] <- 0.7
#' x <- fsn(L = pmax(0, M - 0.1) * (M > 0), M = M, R = pmin(1, M + 0.1) * (M > 0),
#'          labels = c("Perry", "Smith", "Porreca"))
#' plot_fsn_glyphs(x, layout = cbind(c(-1, 0, 1), c(0, 0.9, 0)), values = TRUE)
#' # ten actors of a simulated network
#' plot_fsn_glyphs(simulate_fsn(10, 0.25, seed = 3))
#' # a subset of a simulated survey network, with unknown ties of non-respondents
#' A <- set_unknown(fsn_simulated(n = 20, respondents = 10), "nonrespondents")
#' plot_fsn_glyphs(A, actors = c("a01", "a02", "a03", "a16", "a17", "a18"))
#' @export
plot_fsn_glyphs <- function(x, actors = NULL, layout = NULL, min_core = 0, glyph_size = 0.12,
                            person_col = NULL, fill = "#0072B2", values = FALSE, digits = 2,
                            curvature = 0.15, scale = c(0, 1), main = "", label_cex = 0.85, glyph_pos = 0.4) {
  if (!is.null(actors)) x <- x[actors]
  n <- n_actors(x)
  if (n > 15) warning("the plot is designed for small networks (n = ", n, "); consider 'actors' or plot.fsn()")
  if (is.null(layout)) {
    th <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
    layout <- cbind(cos(th), sin(th))
  }
  layout <- as.matrix(layout)
  if (is.null(person_col)) person_col <- if (is.null(x$respondents)) "#0072B2" else ifelse(x$respondents, "#0072B2", "#9ecae1")
  person_col <- rep_len(person_col, n)
  U <- if (is.null(x$unknown)) matrix(FALSE, n, n) else x$unknown
  E <- (x$M >= min_core & (x$M > 0 | x$R > 0)) | U
  diag(E) <- FALSE
  rg <- apply(layout, 2, range); span <- max(rg[2, ] - rg[1, ], 1)
  ps <- 0.16 * span; gw <- glyph_size * span; gh <- 0.6 * gw
  op <- graphics::par(mar = c(0.4, 0.4, if (nzchar(main)) 2.2 else 0.4, 0.4)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = rg[, 1] + c(-1, 1) * 0.2 * span, ylim = rg[, 2] + c(-0.3, 0.16) * span,
                 asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)
  idx <- which(E, arr.ind = TRUE)
  bez <- function(p0, p2, c1, t) (1 - t)^2 * p0 + 2 * (1 - t) * t * c1 + t^2 * p2
  glyphs <- list()
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    p0 <- layout[i, ]; p2 <- layout[j, ]; d <- p2 - p0; L <- sqrt(sum(d^2)); u <- d / L
    nrm <- c(-u[2], u[1])
    bend <- if (E[j, i]) curvature * L else 0
    c1 <- (p0 + p2) / 2 + bend * nrm
    tt <- seq(0, 1, length.out = 60)
    pts <- t(vapply(tt, function(t) bez(p0, p2, c1, t), numeric(2)))
    keep <- sqrt(colSums((t(pts) - p0)^2)) > 0.45 * ps & sqrt(colSums((t(pts) - p2)^2)) > 0.5 * ps
    pts <- pts[keep, , drop = FALSE]
    unk <- U[i, j]
    col <- if (unk) "grey55" else "grey35"
    lwd <- if (unk) 1 else 0.8 + 2.2 * x$M[i, j]
    graphics::lines(pts, col = col, lwd = lwd, lty = if (unk) 2 else 1)
    m <- nrow(pts)
    graphics::arrows(pts[m - 1, 1], pts[m - 1, 2], pts[m, 1], pts[m, 2], length = 0.09, col = col, lwd = lwd,
                     lty = 1)
    glyphs[[k]] <- list(c = bez(p0, p2, c1, glyph_pos), i = i, j = j, unk = unk)
  }
  for (p in seq_len(n)) .draw_person(layout[p, 1], layout[p, 2], ps, person_col[p])
  for (g in glyphs) {
    cx <- g$c[1]; cy <- g$c[2]
    x0 <- cx - gw / 2; x1 <- cx + gw / 2; y0 <- cy - gh / 2; y1 <- cy + gh / 2
    yt <- y0 - if (values) 0.38 * gh else 0          # room for the values inside the panel
    graphics::rect(x0, yt, x1, y1, col = "white", border = "grey70", lwd = 0.7)
    if (g$unk) { graphics::text(cx, cy, "?", cex = 1.1, col = "grey40"); next }
    sx <- function(v) x0 + 0.08 * gw + 0.84 * gw * (v - scale[1]) / diff(scale)
    by <- y0 + 0.18 * gh; ty <- y1 - 0.12 * gh
    graphics::segments(sx(scale[1]), by, sx(scale[2]), by, col = "grey50", lwd = 0.6)
    graphics::segments(sx(scale), by, sx(scale), by - 0.06 * gh, col = "grey50", lwd = 0.6)
    l <- x$L[g$i, g$j]; mm <- x$M[g$i, g$j]; r <- x$R[g$i, g$j]
    graphics::polygon(c(sx(l), sx(mm), sx(r)), c(by, ty, by), col = .alpha_col(fill, 0.35), border = fill, lwd = 1)
    if (values) graphics::text(cx, y0 - 0.17 * gh, sprintf(paste0("(%.", digits, "f, %.", digits, "f, %.", digits, "f)"), l, mm, r),
                               cex = 0.62, col = "grey20")
  }
  # actor labels on top of everything, on a white background
  ly <- layout[, 2] - 0.55 * ps
  tw <- graphics::strwidth(x$labels, cex = label_cex, font = 2); th <- graphics::strheight("Ag", cex = label_cex, font = 2)
  graphics::rect(layout[, 1] - tw / 2 - 0.15 * th, ly - 0.7 * th, layout[, 1] + tw / 2 + 0.15 * th, ly + 0.7 * th,
                 col = grDevices::adjustcolor("white", 0.85), border = NA)
  graphics::text(layout[, 1], ly, x$labels, cex = label_cex, font = 2)
  invisible(x)
}
