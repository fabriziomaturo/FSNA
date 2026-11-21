# ---------------------------------------------------------------------------
# Simulation of fuzzy social networks and of questionnaire ratings
# ---------------------------------------------------------------------------

#' Simulate a fuzzy social network with a latent crisp network
#'
#' Generates a latent crisp weighted network \eqn{W^*} and a fuzzy network whose
#' ties are vague, possibly erroneous, statements about it. Ties exist with
#' probability \code{p} (Erdos-Renyi model) or with block-dependent
#' probabilities (stochastic block model). For an existing tie, the latent
#' strength is uniform on \code{strength} (\code{model = "er"} or
#' \code{"sbm"}), or, with \code{model = "logit"}, the design of Porreca and
#' Maturo (2026): ties exist with probability \code{p} and
#' \eqn{W^*_{ij} = \mathrm{logit}^{-1}(a_i + b_j + e_{ij})}, with sender,
#' receiver and tie effects normal with standard deviations \code{effects_sd}
#' (default 0.5, 1, 0.5), so that actors differ in centrality. Each tie has a reliability
#' \eqn{\sigma_{ij}}, equal to \code{sigma[2]} with probability \code{pi_vague}
#' and to \code{sigma[1]} otherwise; the core is
#' \eqn{m_{ij} = W^*_{ij} + \sigma_{ij}\varepsilon_{ij}}, truncated to
#' [\code{core_min}, 1], with \eqn{\varepsilon_{ij}} standard normal; the
#' support is \eqn{[m_{ij} - \kappa\sigma_{ij} s^-_{ij},\, m_{ij} + \kappa\sigma_{ij} s^+_{ij}]}
#' truncated to \eqn{[0, 1]}, where \eqn{s^\pm} are asymmetry factors with mean 1.
#'
#' With symmetric spreads (\code{asymmetry = 0}) the latent strength lies in the
#' support with probability \eqn{2\Phi(\kappa) - 1} before truncation: about 68,
#' 95 and 99.7 per cent for \eqn{\kappa = 1, 2, 3}. The spreads are informative
#' about the error of the cores when \code{informative = TRUE}; otherwise
#' \eqn{\sigma_{ij}} used for the spreads is drawn independently of the one used
#' for the error.
#'
#' @param n Number of actors.
#' @param p Tie probability (Erdos-Renyi model).
#' @param kappa Spread multiplier.
#' @param model \code{"er"}, \code{"sbm"} or \code{"logit"}.
#' @param effects_sd Standard deviations of the sender, receiver and tie effects
#'   on the logit scale (\code{model = "logit"}).
#' @param K,p_in,p_out Number of blocks and within/between block tie
#'   probabilities for \code{model = "sbm"}.
#' @param strength Range of the latent strengths (for the SBM, of within-block
#'   ties; between-block ties use \code{strength_out}).
#' @param strength_out Range of the latent strengths of between-block ties.
#' @param sigma Reliability levels (precise, vague).
#' @param pi_vague Probability that a tie is vague.
#' @param asymmetry In \eqn{[0, 1)}: 0 gives symmetric supports; larger values
#'   draw the share of the support on the left from a Beta distribution with
#'   mean 1/2 and increasing dispersion.
#' @param informative Logical, see Details.
#' @param core_min Minimum core of an existing tie.
#' @param seed Optional random seed.
#' @return An \code{\link{fsn}} object with additional components
#'   \code{latent} (the crisp network \eqn{W^*}), \code{sigma} (reliabilities) and
#'   \code{blocks} (block labels, SBM only), and the simulation parameters in
#'   \code{info}.
#' @seealso \code{\link{latent_centrality}}, \code{\link{simulation_study}}
#' @examples
#' set.seed(1)
#' x <- simulate_fsn(30, p = 0.2, kappa = 2)
#' x
#' mean(latent_in_support(x))
#' y <- simulate_fsn(30, model = "sbm", K = 3, p_in = 0.4, p_out = 0.05)
#' table(y$blocks)
#' # design of the simulation study of Porreca and Maturo (2026)
#' z <- simulate_fsn(30, p = 0.1, kappa = 2, model = "logit", sigma = c(0.05, 0.20), seed = 4)
#' summary(colSums(z$latent))
#' @export
#' @section Methodology:
#' The simulator separates a latent crisp network, which plays the role of the true collaboration intensities, from the fuzzy statements about it. The cores contain measurement error with tie-specific reliability \eqn{\sigma_{ij}}, and the spreads are proportional to \eqn{\kappa\sigma_{ij}}: with informative spreads the vaguer ratings are also the less accurate ones. The parameter \eqn{\kappa} controls how often the latent strength falls outside the support: with symmetric spreads about 32, 5 and 0.3 per cent of the ties for \eqn{\kappa = 1, 2, 3}. This allows the coverage of latent centralities by the exact \eqn{\alpha}-cuts and the accuracy of rankings to be assessed (\code{\link{simulation_study}}).
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' # share of latent ties inside the supports for increasing kappa
#' sapply(1:3, function(k) mean(latent_in_support(simulate_fsn(40, 0.3, kappa = k, seed = 1))))
#' # asymmetric ratings and non-informative spreads
#' x <- simulate_fsn(30, 0.2, kappa = 2, asymmetry = 0.6, informative = FALSE, seed = 2)
#' summary(tfn_spread(tfn(fsn_ties(x)$l, fsn_ties(x)$m, fsn_ties(x)$r))$asymmetry)
#' # stochastic block model
#' y <- simulate_fsn(45, model = "sbm", K = 3, p_in = 0.5, p_out = 0.05, seed = 3)
#' plot(y, vertex_col = fsna_palette(3)[y$blocks])
simulate_fsn <- function(n = 30, p = 0.2, kappa = 2, model = c("er", "sbm", "logit"), K = 3,
                         p_in = 0.3, p_out = 0.05, strength = c(0.2, 1), strength_out = c(0.1, 0.5),
                         sigma = c(0.03, 0.10), pi_vague = 0.3, asymmetry = 0, informative = TRUE,
                         core_min = 0.05, effects_sd = c(sender = 0.5, receiver = 1, tie = 0.5),
                         seed = NULL) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)
  if (asymmetry < 0 || asymmetry >= 1) stop("'asymmetry' must lie in [0, 1)")
  nn <- n * n
  blocks <- NULL
  if (model == "er") {
    E <- matrix(stats::runif(nn) < p, n, n)
    Wst <- matrix(stats::runif(nn, strength[1], strength[2]), n, n)
  } else if (model == "logit") {
    E <- matrix(stats::runif(nn) < p, n, n)
    a <- stats::rnorm(n, 0, effects_sd[1]); b <- stats::rnorm(n, 0, effects_sd[2])
    Wst <- stats::plogis(outer(a, b, "+") + matrix(stats::rnorm(nn, 0, effects_sd[3]), n, n))
  } else {
    blocks <- sample(rep(seq_len(K), length.out = n))
    same <- outer(blocks, blocks, "==")
    E <- matrix(stats::runif(nn), n, n) < ifelse(same, p_in, p_out)
    Wst <- matrix(ifelse(same, stats::runif(nn, strength[1], strength[2]),
                         stats::runif(nn, strength_out[1], strength_out[2])), n, n)
  }
  diag(E) <- FALSE
  sig <- matrix(ifelse(stats::runif(nn) < pi_vague, sigma[2], sigma[1]), n, n)
  sig_sp <- if (informative) sig else matrix(ifelse(stats::runif(nn) < pi_vague, sigma[2], sigma[1]), n, n)
  M <- Wst + sig * matrix(stats::rnorm(nn), n, n)
  M[] <- pmin(1, pmax(core_min, M))
  if (asymmetry > 0) {
    conc <- 2 * (1 - asymmetry) / asymmetry
    sh <- matrix(stats::rbeta(nn, conc, conc), n, n)
  } else sh <- matrix(0.5, n, n)
  L <- M - kappa * sig_sp * 2 * sh; R <- M + kappa * sig_sp * 2 * (1 - sh)
  L[] <- pmax(0, L); R[] <- pmin(1, R)
  L[!E] <- 0; M[!E] <- 0; R[!E] <- 0; Wst[!E] <- 0
  lab <- sprintf("a%02d", seq_len(n))
  x <- fsn(L, M, R, labels = lab,
           info = list(description = sprintf("Simulated network (%s, n = %d, kappa = %g)", model, n, kappa),
                       model = model, n = n, p = p, kappa = kappa, K = K, p_in = p_in, p_out = p_out,
                       sigma = sigma, pi_vague = pi_vague, asymmetry = asymmetry, informative = informative))
  dimnames(Wst) <- dimnames(sig) <- list(lab, lab)
  x$latent <- Wst; x$sigma <- sig; x$blocks <- blocks
  x
}

#' Latent ties inside the supports of the fuzzy ties
#'
#' @param x A simulated network from \code{\link{simulate_fsn}}.
#' @param alpha Level of the cut.
#' @return Logical vector over the existing ties: \code{TRUE} if the latent
#'   strength lies in the \eqn{\alpha}-cut. The whole latent network lies in the
#'   box \eqn{\mathcal W_\alpha} when all values are \code{TRUE}.
#' @examples
#' x <- simulate_fsn(20, 0.3, kappa = 1, seed = 2)
#' mean(latent_in_support(x))
#' @export
#' @section Methodology:
#' Proposition 2 of Porreca and Maturo (2026) guarantees coverage of the latent centrality by the exact \eqn{\alpha}-cut only when the whole latent network lies in the box \eqn{\mathcal W_\alpha}, an event that is rare in large networks even when most ties are covered.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' x <- simulate_fsn(60, 0.3, kappa = 3, seed = 4)
#' mean(latent_in_support(x))   # share of covered ties
#' all(latent_in_support(x))    # whole network in the box
latent_in_support <- function(x, alpha = 0) {
  if (is.null(x$latent)) stop("'x' has no latent network")
  ct <- fsn_cut(x, alpha); E <- x$M > 0
  (x$latent >= ct$lower - 1e-12 & x$latent <= ct$upper + 1e-12)[E]
}

#' Simulate mouse-tracking questionnaire ratings
#'
#' Generates ratings in the format of \code{\link{fsna_questionnaire}}: a subset of
#' actors respond, and each respondent rates a subset of colleagues on a
#' continuous scale \eqn{[1, s]}, \eqn{s} = \code{scale_max}. The pointer moves from the
#' centre of the scale toward the final position with random hesitations; the
#' rating is summarised by (minimum, final position, maximum) of the path,
#' divided by \code{scale_max} (\code{\link{tfn_from_trajectory}}). Some pairs
#' are rated twice, to illustrate the \code{keep} rule of
#' \code{\link{fsn_from_ratings}}.
#'
#' @param n Number of members.
#' @param response_rate Share of members who respond.
#' @param rate_prob Probability that a respondent rates a given colleague.
#' @param n_groups Number of latent working groups: ratings are higher within
#'   groups.
#' @param hesitation Standard deviation of the hesitations of the pointer, on
#'   the rating scale.
#' @param n_steps Number of recorded positions per rating.
#' @param repeat_prob Probability that a rated pair is rated a second time.
#' @param scale_max Maximum of the scale.
#' @param keep_paths Logical, return the simulated paths in a list column.
#' @param seed Optional seed.
#' @return A data frame with columns \code{order}, \code{rater}, \code{rated},
#'   \code{l}, \code{m}, \code{r}, \code{l_raw}, \code{m_raw}, \code{r_raw} and,
#'   optionally, \code{path}; the latent groups are in the attribute
#'   \code{"groups"}.
#' @examples
#' sim <- simulate_ratings(40, seed = 3)
#' head(sim)
#' x <- fsn_from_ratings(sim, from = "rater", to = "rated")
#' x
#' @export
#' @section Methodology:
#' Generates data with the features of the departmental questionnaire: nonresponse, partial rating of colleagues, working groups with higher ratings within groups, hesitating pointer paths and occasional repeated ratings. Useful to test import and analysis pipelines before collecting real data.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
#' @examples
#' s <- simulate_ratings(30, response_rate = 0.4, keep_paths = TRUE, seed = 5)
#' plot(s$path[[1]], type = "b", ylim = c(1, 5),
#'      ylab = "position on the scale", xlab = "recorded point")
#' x <- fsn_from_ratings(s, from = "rater", to = "rated", order = "order")
#' table(respondent = x$respondents)
simulate_ratings <- function(n = 40, response_rate = 0.5, rate_prob = 0.35, n_groups = 3,
                             hesitation = 0.35, n_steps = 12, repeat_prob = 0.02, scale_max = 5,
                             keep_paths = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- sample(rep(seq_len(n_groups), length.out = n))
  resp <- sort(sample.int(n, max(1, round(response_rate * n))))
  rows <- list(); k <- 0L
  one_path <- function(final) {
    tt <- seq(0, 1, length.out = n_steps)
    target <- 3 + (final - 3) * tt^0.6
    noise <- cumsum(stats::rnorm(n_steps, 0, hesitation / sqrt(n_steps))) * (1 - tt^2)
    path <- pmin(scale_max, pmax(1, target + noise)); path[n_steps] <- final; path
  }
  for (i in resp) {
    rated <- setdiff(seq_len(n), i)
    p_rate <- ifelse(g[rated] == g[i], min(1, 2 * rate_prob), rate_prob / 2)
    rated <- rated[stats::runif(length(rated)) < p_rate]
    for (j in rated) {
      mu <- if (g[j] == g[i]) 3.8 else 2.2
      final <- min(scale_max, max(1, stats::rnorm(1, mu, 0.8)))
      reps <- if (stats::runif(1) < repeat_prob) 2 else 1
      for (rr in seq_len(reps)) {
        path <- one_path(if (rr == 1) final else min(scale_max, max(1, final + stats::rnorm(1, 0, 0.5))))
        k <- k + 1L
        rows[[k]] <- data.frame(order = k, rater = i, rated = j, l_raw = min(path),
                                m_raw = path[n_steps], r_raw = max(path))
        if (keep_paths) rows[[k]]$path <- I(list(path))
      }
    }
  }
  d <- do.call(rbind, rows)
  if (is.null(d)) {
    # no rating generated: an empty data frame with the documented columns
    d <- data.frame(order = integer(0), rater = integer(0), rated = integer(0),
                    l_raw = numeric(0), m_raw = numeric(0), r_raw = numeric(0))
    if (keep_paths) d$path <- list()
  }
  d$l <- d$l_raw / scale_max; d$m <- d$m_raw / scale_max; d$r <- d$r_raw / scale_max
  d <- d[, c("order", "rater", "rated", "l", "m", "r", "l_raw", "m_raw", "r_raw", if (keep_paths) "path")]
  attr(d, "groups") <- g
  d
}

#' A simulated network with the structure of a questionnaire survey
#'
#' Simulates a directed fuzzy network in which only a subset of the actors (the
#' respondents) send ties, as in a collaboration survey where each respondent
#' rates some colleagues and the non-respondents are only rated. The ties are
#' generated by the latent model of \code{\link{simulate_fsn}} with
#' \code{model = "logit"}. The network replaces in the examples the
#' questionnaire data of the two departments, which are not distributed with
#' the package (\code{\link{fsna_questionnaire}}).
#' @param n Number of actors.
#' @param respondents Number of respondents (the first \code{respondents}
#'   actors send ties).
#' @param p Probability that a respondent rates a given colleague.
#' @param kappa Spread multiplier of the supports.
#' @param sigma Noise levels of the ties, see \code{\link{simulate_fsn}}.
#' @param seed Random seed.
#' @param ... Further arguments to \code{\link{simulate_fsn}} (not \code{model},
#'   which is \code{"logit"}).
#' @return An \code{\link{fsn}} object with the latent network in
#'   \code{$latent} and the respondents in \code{$respondents}.
#' @examples
#' x <- fsn_simulated()
#' x
#' sum(x$respondents)
#' fc <- fuzzy_centrality(x, "in_strength")
#' sort(top_set(fc), decreasing = TRUE)[1:5]
#' @export
#' @section Methodology:
#' Only respondents send ties, so every actor is rated by the same set of potential raters and paths can pass only through respondents, as in the networks analysed in Section 4 of Porreca and Maturo (2026). The sender, receiver and tie effects make actors differ in centrality; the spreads of the supports are proportional to the noise of each tie, as in the simulation design of Section 3.1 of that paper.
#' @references Porreca A., Maturo F. (2026). Exact fuzzy centrality in directed fuzzy social networks: alpha-cuts by monotonicity and possibilistic rankings of actors. Manuscript submitted for publication; preprint forthcoming on arXiv.
fsn_simulated <- function(n = 50, respondents = 25, p = 0.25, kappa = 2, sigma = c(0.05, 0.20),
                          seed = 1, ...) {
  if (respondents < 1 || respondents > n) stop("'respondents' must be between 1 and n")
  x <- simulate_fsn(n, p = p, kappa = kappa, model = "logit", sigma = sigma, seed = seed, ...)
  nr <- seq_len(n) > respondents
  for (k in c("L", "M", "R", "latent")) x[[k]][nr, ] <- 0
  x$respondents <- stats::setNames(!nr, x$labels)
  x$info$description <- sprintf("Simulated survey network (n = %d, %d respondents, kappa = %g)", n, respondents, kappa)
  x
}
