# ---------------------------------------------------------------------------
# Simulation study: exact cuts, approximations, rankings and certification
# ---------------------------------------------------------------------------

#' Simulation study of exact fuzzy centralities
#'
#' For each combination of \code{n}, \code{p} and \code{kappa}, simulates
#' \code{nrep} fuzzy networks with \code{\link{simulate_fsn}} and, for each
#' index, compares the exact \eqn{\alpha}-cuts with the triangular approximation
#' and with Monte Carlo inner approximations, and evaluates rankings and
#' certified orderings against the centralities of the latent network.
#'
#' Coverage of the latent centralities by the exact cuts is guaranteed only when
#' the latent network lies in the support box, which is rare in large networks;
#' it is otherwise an empirical property of the index. Certified orderings agree
#' with the CoG of the triangle by construction, so the informative quantities
#' are the proportion of certified pairs, the error rate among certified pairs and
#' the error rate of the CoG among the pairs that are not certified.
#'
#' @param n,p,kappa Vectors of network sizes, tie probabilities and spread
#'   multipliers (full factorial design).
#' @param nrep Replicates per scenario.
#' @param indices Isotone indices (see \code{\link{centrality_index}}).
#' @param B Numbers of Monte Carlo draws.
#' @param alpha Grid of levels; the levels 0, 0.5 and 1 are always added, so
#'   that the indicators declared at these levels are computed there (the level
#'   used for the indicators at 0.5 is returned in \code{alpha_eval}).
#' @param max_steps Step limit for t-norm closeness.
#' @param seed Seed; each replicate uses a seed derived from it, so that results
#'   do not depend on the order of execution.
#' @param progress Logical, print progress.
#' @param model,sigma Generating model and reliability levels passed to
#'   \code{\link{simulate_fsn}}; the defaults (\code{"logit"},
#'   \code{c(0.05, 0.20)}) are the design of Porreca and Maturo (2026), whose
#'   simulation study is reproduced with \code{nrep = 100}.
#' @param bracket Rule for the suprema of the possibility degrees used by the
#'   fuzzy set of the most central actor, see \code{\link{possibility_matrix}}:
#'   \code{"auto"} (interpolation for strength, bisection otherwise, as in
#'   Porreca and Maturo, 2026) or \code{"upper"}, faster for large designs.
#' @param ... Further arguments to \code{\link{simulate_fsn}}.
#' @return A data frame of class \code{"fsna_simulation"} with one row per
#'   replicate and index. Columns: scenario parameters; \code{ties_in_support}
#'   and \code{network_in_box}; \code{coverage_exact} and \code{coverage_mc<B>}
#'   (share of actors whose latent centrality lies in the exact or Monte Carlo
#'   0-cut); \code{mc<B>_outside} (draws outside the exact 0-cut, which must be
#'   zero); \code{width_mc<B>} (Monte Carlo width relative to the exact width);
#'   \code{tri_error} (mean over actors of the error of the triangle at
#'   \eqn{\alpha = 0.5} relative to the support width); \code{disc_cog},
#'   \code{disc_core}, \code{disc_cog_uncertified}; \code{cert0}, \code{cert05}
#'   (share of certified pairs) with the counts \code{n_cert0}, \code{wrong_cert0},
#'   \code{n_cert05}, \code{wrong_cert05}; \code{top_cog}, \code{top_in_set},
#'   \code{top_set_size}.
#' @examples
#' sim <- simulation_study(n = 12, p = 0.3, kappa = c(1, 3), nrep = 1,
#'                         indices = c("in_strength", "in_closeness_prod"), B = 20)
#' summary(sim)
#' @export
#' @section Methodology:
#' The design separates what follows from theory from what must be assessed empirically. Exactness of the cuts and agreement of certified orderings with the CoG follow from the theory; coverage when the latent network falls outside the box, the size of the error of the triangular and Monte Carlo approximations, the proportion of certifiable orderings and the error rates among certified and uncertified pairs do not. The results depend on the generating model, in particular on the informativeness of the spreads and on \eqn{\kappa}.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' sim <- simulation_study(n = 12, p = 0.3, kappa = 2, nrep = 1,
#'                         indices = "in_closeness_min", B = 20)
#' cols <- c("index", "coverage_exact", "coverage_mc20", "tri_error", "cert0")
#' as.data.frame(unclass(sim))[, cols]
simulation_study <- function(n = c(30, 60), p = c(0.1, 0.3), kappa = 1:3, nrep = 20,
                             indices = c("in_strength", "in_closeness_min", "in_closeness_prod", "in_katz"),
                             B = c(50, 500), alpha = seq(0, 1, by = 0.05), max_steps = Inf, seed = 1,
                             progress = interactive(), model = "logit", sigma = c(0.05, 0.20),
                             bracket = "auto", ...) {
  design <- expand.grid(n = n, p = p, kappa = kappa, rep = seq_len(nrep))
  # the levels 0, 0.5 and 1 are always computed: the indicators declared at
  # these levels (support, cert05, tri_error, top set, cores) use them exactly
  alpha <- sort(unique(c(alpha, 0, 0.5, 1)))
  g5 <- match(0.5, alpha); out <- vector("list", nrow(design))
  for (s in seq_len(nrow(design))) {
    d <- design[s, ]
    set.seed(seed * 100003L + s)
    x <- simulate_fsn(d$n, d$p, d$kappa, model = model, sigma = sigma, ...)
    inside <- latent_in_support(x)
    rows <- list()
    for (ix in indices) {
      args <- list(x = x, index = ix, alpha = alpha)
      if (grepl("closeness", ix)) args$max_steps <- max_steps
      fc <- do.call(fuzzy_centrality, args)
      pr <- fc$index$params; pr$label <- NULL
      truth <- do.call(latent_centrality, c(list(x = x, index = ix), pr))
      r <- data.frame(n = d$n, p = d$p, kappa = d$kappa, rep = d$rep, index = ix, alpha_eval = alpha[g5],
                      ties_in_support = mean(inside), network_in_box = all(inside),
                      coverage_exact = mean(coverage(fc, truth)))
      for (b in B) {
        mc <- do.call(mc_centrality, c(list(x = x, index = ix, B = b), pr))
        w <- fc$upper[, 1] - fc$lower[, 1]
        r[[paste0("coverage_mc", b)]] <- mean(truth >= mc$lower - 1e-10 & truth <= mc$upper + 1e-10)
        r[[paste0("mc", b, "_outside")]] <- sum(mc$lower < fc$lower[, 1] - 1e-9 | mc$upper > fc$upper[, 1] + 1e-9)
        r[[paste0("width_mc", b)]] <- mean(((mc$upper - mc$lower) / w)[w > 1e-12])
      }
      ta <- triangular_approx(fc); r$tri_error <- mean(ta$error[, g5])
      cog <- tfn_score(fc$triangle, "cog"); core <- (fc$lower[, length(alpha)] + fc$upper[, length(alpha)]) / 2
      r$disc_cog <- discordance(cog, truth); r$disc_core <- discordance(core, truth)
      C0 <- certified_pairs(fc, 0); C5 <- certified_pairs(fc, alpha[g5])
      np <- d$n * (d$n - 1) / 2
      r$disc_cog_uncertified <- discordance(cog, truth, pairs = !(C0 | t(C0)) & row(C0) != col(C0))
      wrong <- function(C) sum(C & outer(truth, truth, "<"))
      r$cert0 <- sum(C0) / np; r$n_cert0 <- sum(C0); r$wrong_cert0 <- wrong(C0)
      r$cert05 <- sum(C5) / np; r$n_cert05 <- sum(C5); r$wrong_cert05 <- wrong(C5)
      tl <- which.max(truth)
      ts <- top_set(fc, bracket = bracket, alpha = alpha[g5])
      r$top_cog <- which.max(cog) == tl; r$top_in_set <- x$labels[tl] %in% ts; r$top_set_size <- length(ts)
      rows[[ix]] <- r
    }
    out[[s]] <- do.call(rbind, rows)
    if (progress && s %% 10 == 0) message(s, " of ", nrow(design), " networks")
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  class(res) <- c("fsna_simulation", "data.frame"); res
}

#' @export
summary.fsna_simulation <- function(object, by = c("index", "kappa"), ...) {
  d <- as.data.frame(unclass(object))
  key <- interaction(d[, by], drop = TRUE, lex.order = TRUE)
  num <- setdiff(names(d)[vapply(d, function(z) is.numeric(z) || is.logical(z), TRUE)],
                 c("n", "p", "kappa", "rep", "n_cert0", "wrong_cert0", "n_cert05", "wrong_cert05"))
  s <- lapply(split(d, key), function(z) {
    m <- as.data.frame(lapply(z[, num], function(v) mean(as.numeric(v), na.rm = TRUE)))
    cbind(z[1, by, drop = FALSE], replicates = nrow(z), m,
          wrong_among_cert0 = if (sum(z$n_cert0)) sum(z$wrong_cert0) / sum(z$n_cert0) else NA,
          wrong_among_cert05 = if (sum(z$n_cert05)) sum(z$wrong_cert05) / sum(z$n_cert05) else NA)
  })
  s <- do.call(rbind, s); rownames(s) <- NULL; s
}

#' @export
plot.fsna_simulation <- function(x, what = c("coverage_exact", "disc_cog", "cert0", "tri_error"), ...) {
  s <- summary(x)
  ix <- unique(s$index); cols <- fsna_palette(length(ix))
  op <- graphics::par(mfrow = c(1, length(what)), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
  for (w in what) {
    graphics::plot(NA, xlim = range(s$kappa), ylim = range(c(0, s[[w]]), na.rm = TRUE), xlab = expression(kappa),
                   ylab = w, main = w, las = 1, xaxt = "n")
    graphics::axis(1, at = unique(s$kappa))
    for (k in seq_along(ix)) { z <- s[s$index == ix[k], ]; graphics::lines(z$kappa, z[[w]], type = "b", col = cols[k], pch = 19) }
  }
  graphics::legend("topright", legend = ix, col = cols, lwd = 1, pch = 19, bty = "n", cex = 0.7)
  invisible(x)
}
