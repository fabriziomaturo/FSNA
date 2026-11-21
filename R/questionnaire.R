# ---------------------------------------------------------------------------
# Mouse-tracking questionnaire: parsing of the raw export and conversion of
# trajectories into triangular fuzzy ratings
# ---------------------------------------------------------------------------

# angle (degrees) at p1 from the direction p1 -> p2 to the direction p1 -> p3,
# truncated at zero; screen coordinates (y grows downwards).
.angle_deg <- function(p2, p1, p3) {
  max(0, (atan2(p3[2] - p1[2], p3[1] - p1[1]) - atan2(p2[2] - p1[2], p2[1] - p1[1])) * 180 / pi)
}

#' Convert a mouse trajectory on the semicircular scale into a rating path
#'
#' The questionnaire shows five buttons on a semicircle around the starting
#' position of the pointer. Each recorded position of the pointer is converted
#' into a continuous position on the rating scale \eqn{[1, 5]} through the angle
#' that it forms, at the starting position, with the direction of the first
#' button: the angle is divided by the angle of the fifth button, multiplied by
#' 4 and shifted by 1. Positions beyond the fifth button but above the
#' horizontal line through the starting point are set to 5; positions below it
#' are set to 1.
#'
#' @param x,y Numeric vectors, screen coordinates of the recorded positions
#'   (absolute, y growing downwards).
#' @param start Numeric vector \code{c(x, y)}, the centre of the semicircle.
#' @param buttons Numeric matrix with five rows and columns \code{x}, \code{y}:
#'   the positions of the five buttons, from rating 1 to rating 5.
#' @return Numeric vector of positions on \eqn{[1, 5]}, one per recorded point.
#' @seealso \code{\link{tfn_from_trajectory}}, \code{\link{read_fuzzy_questionnaire}}
#' @examples
#' buttons <- cbind(x = c(414, 522, 709, 896, 1004), y = c(497, 326, 264, 325, 495))
#' trajectory_to_scale(x = c(700, 850, 890), y = c(400, 330, 326),
#'                     start = c(709, 557), buttons = buttons)
#' @export
#' @section Methodology:
#' The semicircular layout gives the pointer more freedom than a linear scale; the angle with the first button, rescaled to \eqn{[1, 5]}, is the position on the scale. Positions beyond the fifth button but above the horizontal line through the starting point are truncated to 5, and those below it to 1.
#' @examples
#' buttons <- cbind(x = c(414, 522, 709, 896, 1004), y = c(497, 326, 264, 325, 495))
#' th <- seq(pi, 0, length.out = 9)
#' round(trajectory_to_scale(709 + 290 * cos(th), 557 - 290 * sin(th), c(709, 557), buttons), 2)
trajectory_to_scale <- function(x, y, start, buttons) {
  stopifnot(length(x) == length(y), length(start) == 2, nrow(buttons) == 5)
  b <- as.matrix(buttons)
  a5 <- .angle_deg(b[1, ], start, b[5, ])
  divisor <- .angle_deg(b[1, ], start, c(start[1], start[2] + 20))
  ang <- vapply(seq_along(x), function(k) .angle_deg(b[1, ], start, c(x[k], y[k])), numeric(1))
  ang[ang > a5 & ang <= divisor] <- a5
  ang[ang > divisor] <- 0
  ang / a5 * 4 + 1
}

#' Triangular fuzzy rating from a trajectory on the rating scale
#'
#' Following the Dynamic Fuzzy Rating Tracker idea of Calcagni and Lombardi
#' (2014), a rating is summarised by the triangular fuzzy number
#' (minimum, final position, maximum) of the path of the pointer on the rating
#' scale, divided by the maximum of the scale.
#'
#' @param path Numeric vector, positions of the pointer on the rating scale in
#'   temporal order (the last one is the final position).
#' @param scale_max Maximum of the rating scale.
#' @param scale_min Value mapped to 0. The TFN is
#'   \eqn{(x - s_{min})/(s_{max} - s_{min})}: the default \code{scale_min = 0}
#'   divides by the maximum, as in Porreca, Maturo and Ventre (2025);
#'   \code{scale_min = 1} maps the scale \eqn{[1, 5]} onto \eqn{[0, 1]}, so
#'   that the lowest answer is the absence of a tie, as in Porreca and Maturo
#'   (2026). Use \code{scale_max = 1, scale_min = 0} to keep the original scale.
#' @return A \code{\link{tfn}} of length one.
#' @examples
#' tfn_from_trajectory(c(3.1, 3.8, 4.4, 4.1))
#' @export
#' @section Methodology:
#' The minimum and maximum of the path give the support and the final position gives the core: the support collects the values the respondent considered while hesitating. This construction captures hesitation, not necessarily uncertainty about the tie: a respondent who moves the pointer quickly and directly produces a narrow TFN irrespective of how sure he or she is.
#' @examples
#' tfn_from_trajectory(c(3, 3.5, 4.2, 3.9))
#' tfn_from_trajectory(c(3, 3.9), scale_max = 1)     # keep the 1-5 scale
#' tfn_from_trajectory(c(3, 3.9), scale_min = 1)     # (x - 1) / 4
tfn_from_trajectory <- function(path, scale_max = 5, scale_min = 0) {
  path <- as.numeric(path)
  if (!length(path) || anyNA(path)) stop("'path' must be a non-empty numeric vector without NA")
  if (scale_max <= scale_min) stop("'scale_max' must exceed 'scale_min'")
  s <- scale_max - scale_min
  tfn(min(path) - scale_min, path[length(path)] - scale_min, max(path) - scale_min) / s
}

#' Read the raw export of the fuzzy questionnaire
#'
#' Parses the semicolon-separated export of the web questionnaire used in
#' Porreca, Maturo and Ventre (2025): one row per submitted questionnaire, with
#' the department, the date and a JSON field with, for each rated colleague,
#' the positions of the buttons, the starting point of the pointer, the recorded
#' positions and the clicked button. Each trajectory is converted with
#' \code{\link{trajectory_to_scale}} and summarised by
#' \code{\link{tfn_from_trajectory}}.
#'
#' The JSON field must contain, as its first element, the identifier of the
#' respondent (with a null value) and then one element per rated colleague.
#' Elements without positions, empty questionnaires and respondents with
#' identifier above \code{max_id} (test accounts) are skipped. When a respondent rates the same colleague more than once, all
#' ratings are returned with their submission order; \code{\link{fsn_from_ratings}}
#' keeps the last one by default.
#'
#' @param file Path to the CSV export.
#' @param sep Field separator.
#' @param max_id Respondents with a larger identifier are discarded.
#' @param scale_max Maximum of the rating scale (5 in the original study).
#' @param scale_min Value mapped to 0 (see \code{\link{tfn_from_trajectory}}):
#'   0 divides by \code{scale_max}, as in the original study; 1 gives
#'   \eqn{(x - 1)/4}.
#' @return A data frame with one row per rating and columns \code{response_id},
#'   \code{department}, \code{date}, \code{order} (submission order),
#'   \code{rater}, \code{rated}, \code{l}, \code{m}, \code{r} (TFN on
#'   \eqn{[0, 1]}), \code{l_raw}, \code{m_raw}, \code{r_raw} (on
#'   \eqn{[1, s]}, with \eqn{s} = \code{scale_max}), \code{button} (clicked button) and
#'   \code{n_points} (number of recorded positions).
#' @note Requires the \pkg{jsonlite} package.
#' @seealso \code{\link{fsn_from_ratings}}
#' @examples
#' # requires the questionnaire export, which is not distributed
#' if (fsna_data_available()) {
#'   raw <- read_fuzzy_questionnaire(fsna_data_path())
#'   head(raw)
#' }
#' @export
#' @section Methodology:
#' For each rated colleague, only the movements after the last start of the pointer are used, each position is converted with \code{\link{trajectory_to_scale}}, and the rating is (minimum, final position, maximum) divided by 5.
#' @examples
#' # the questionnaire data of the two departments are not distributed with
#' # the package (see ?fsna_questionnaire)
#' if (fsna_data_available()) {
#'   raw <- read_fuzzy_questionnaire(fsna_data_path())
#'   print(table(raw$department))
#'   print(summary(raw$n_points))
#'   # the final position agrees with the clicked button
#'   mean(round(raw$m_raw) == raw$button)
#' }
read_fuzzy_questionnaire <- function(file, sep = ";", max_id = 1000, scale_max = 5, scale_min = 0) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("package 'jsonlite' is required to parse the questionnaire")
  raw <- utils::read.csv(file, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  names(raw) <- c("response_id", "department", "date", "answer")[seq_len(ncol(raw))]
  out <- vector("list", nrow(raw)); ord <- 0L
  for (row in seq_len(nrow(raw))) {
    ans <- jsonlite::fromJSON(as.character(raw$answer[row]), simplifyVector = FALSE)
    if (!length(ans)) next
    person <- suppressWarnings(as.numeric(names(ans)[1]))
    if (is.na(person) || person > max_id || length(ans) < 2) next
    recs <- list()
    for (id in 2:length(ans)) {
      el <- ans[[id]]
      if (!length(el) || is.null(el[[1]]$positions$start$left)) next
      pos <- el[[1]]$positions
      start <- c(pos$start$left, pos$start$top)
      buttons <- t(vapply(1:5, function(b) {
        p <- pos[[paste0("button", b)]]; c(p$left, p$top)
      }, numeric(2)))
      nm <- vapply(el, function(e) names(e)[1], character(1))
      last <- length(el) - 1L
      st <- max(which(nm[seq_len(last)] == "start"))
      click <- c(el[[st]]$start$x, el[[st]]$start$y)
      if (st >= last) next
      pts <- t(vapply(el[(st + 1):last], function(e) c(e$info$x, e$info$y), numeric(2)))
      pts <- sweep(pts, 2, click, "+")
      path <- trajectory_to_scale(pts[, 1], pts[, 2], start, buttons)
      ord <- ord + 1L
      btn <- suppressWarnings(as.numeric(el[[length(el)]]$resp))
      recs[[length(recs) + 1L]] <- data.frame(
        response_id = raw$response_id[row], department = raw$department[row], date = raw$date[row],
        order = ord, rater = person, rated = as.numeric(names(ans)[id]),
        l_raw = min(path), m_raw = path[length(path)], r_raw = max(path),
        button = if (length(btn)) btn else NA_real_, n_points = length(path))
    }
    out[[row]] <- do.call(rbind, recs)
  }
  res <- do.call(rbind, out)
  sc <- function(v) (v - scale_min) / (scale_max - scale_min)
  res$l <- sc(res$l_raw); res$m <- sc(res$m_raw); res$r <- sc(res$r_raw)
  rownames(res) <- NULL
  res[, c("response_id", "department", "date", "order", "rater", "rated",
          "l", "m", "r", "l_raw", "m_raw", "r_raw", "button", "n_points")]
}
