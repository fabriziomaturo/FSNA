# ---------------------------------------------------------------------------
# Elicitation of fuzzy ratings without mouse tracking; reading rating files
# ---------------------------------------------------------------------------

#' Default linguistic scale for collaboration ratings
#'
#' Five linguistic terms with triangular meanings on \eqn{[0, 1]} forming a
#' Ruspini partition (memberships sum to one).
#' @return A data frame with columns \code{label}, \code{l}, \code{m}, \code{r}.
#' @examples
#' linguistic_scale()
#' @export
#' @section Methodology:
#' A Ruspini partition of \eqn{[0, 1]}: each value belongs to at most two adjacent terms with memberships summing to one. Replace it with a term set suited to the question.
#' @examples
#' ls <- linguistic_scale()
#' plot(tfn(ls$l, ls$m, ls$r, names = ls$label), main = "Linguistic scale")
linguistic_scale <- function() {
  data.frame(label = c("none", "occasional", "moderate", "strong", "very strong"),
             l = c(0, 0, 0.25, 0.5, 0.75), m = c(0, 0.25, 0.5, 0.75, 1), r = c(0.25, 0.5, 0.75, 1, 1),
             stringsAsFactors = FALSE)
}

#' Convert questionnaire responses into triangular fuzzy ratings
#'
#' Converts responses collected without mouse tracking into TFNs on
#' \eqn{[0, 1]}. Four elicitation formats are supported:
#' \describe{
#'   \item{\code{"three_point"}}{the respondent gives the lowest plausible, the
#'     most plausible and the highest plausible value (columns \code{low},
#'     \code{best}, \code{high}): the TFN is \eqn{(low, best, high)}. This is the
#'     most direct elicitation of a possibility distribution and the one
#'     recommended for new data.}
#'   \item{\code{"interval"}}{a range slider for the plausible values and the
#'     most plausible value inside it: same columns as \code{"three_point"}.}
#'   \item{\code{"confidence"}}{a single value \code{best} and a confidence
#'     \code{confidence} in \eqn{[0, 1]} (or on a Likert scale, rescaled): the
#'     support is \eqn{best \pm (1 - confidence)\, s_{\max}}, truncated to the
#'     scale, with \eqn{s_{\max}} = \code{max_spread}.}
#'   \item{\code{"linguistic"}}{a linguistic label (column \code{label}) mapped to
#'     a TFN through \code{dictionary} (default \code{\link{linguistic_scale}}); a
#'     second, optional label column \code{label2} expresses hesitation between two
#'     adjacent terms, in which case the TFN spans both.}
#' }
#' Values on a scale \eqn{[a, b]} are mapped to \eqn{[0, 1]} by
#' \eqn{(x - a)/(b - a)} (\code{rescale = "minmax"}) or by \eqn{x / b}
#' (\code{"divide"}, as in the original study, which maps the lowest rating to
#' \eqn{a/b} rather than to zero).
#'
#' @param data A data frame with one row per response.
#' @param method \code{"three_point"}, \code{"interval"}, \code{"confidence"} or
#'   \code{"linguistic"}.
#' @param from,to Columns with rater and rated actor.
#' @param low,best,high,confidence,label,label2 Column names.
#' @param scale Range \eqn{[a, b]} of the numeric responses.
#' @param rescale \code{"minmax"} or \code{"divide"}.
#' @param confidence_scale Range of the confidence column.
#' @param max_spread Half-width of the support at zero confidence, on
#'   \eqn{[0, 1]}.
#' @param dictionary Data frame with columns \code{label}, \code{l}, \code{m},
#'   \code{r}.
#' @param fix_order Logical: if \code{low > best} or \code{best > high} the values
#'   are sorted (with a warning) instead of raising an error.
#' @return A data frame with columns \code{from}, \code{to}, \code{l}, \code{m},
#'   \code{r}, ready for \code{\link{fsn_from_ratings}}, plus the original
#'   columns.
#' @examples
#' d <- data.frame(rater = c("a", "a", "b"), rated = c("b", "c", "a"),
#'                 low = c(3, 1, 5), best = c(5, 2, 7), high = c(6, 4, 7))
#' fuzzify_responses(d, "three_point", scale = c(1, 7))
#' d2 <- data.frame(rater = "a", rated = c("b", "c"), best = c(4, 2), confidence = c(5, 2))
#' fuzzify_responses(d2, "confidence", scale = c(1, 5), confidence_scale = c(1, 5))
#' d3 <- data.frame(rater = "a", rated = c("b", "c"), label = c("strong", "moderate"),
#'                  label2 = c(NA, "strong"))
#' fuzzify_responses(d3, "linguistic")
#' @export
#' @section Methodology:
#' Mouse tracking captures hesitation indirectly. The formats implemented here elicit vagueness directly. The three-point format asks for the lowest plausible, most plausible and highest plausible value, a direct statement of a triangular possibility distribution; it is the format of three-point estimation in project management (Malcolm et al., 1959) and of expert elicitation (Garthwaite, Kadane and O'Hagan, 2005). The confidence format derives the spread from a self-assessed confidence and imposes symmetry. The linguistic format maps verbal labels to fuzzy numbers and allows hesitation between adjacent labels. The choice of format and scale defines what the spreads mean, and results obtained with different formats are not directly comparable.
#' @examples
#' # interval responses on a 1-7 scale, divided by the maximum as in the original study
#' d <- data.frame(rater = "a", rated = c("b", "c"),
#'                 low = c(2, 5), best = c(3, 6), high = c(5, 7))
#' fuzzify_responses(d, "interval", scale = c(1, 7), rescale = "divide")
#' # unordered answers can be sorted with a warning
#' d$low[1] <- 4
#' suppressWarnings(fuzzify_responses(d, "three_point", scale = c(1, 7), fix_order = TRUE))
fuzzify_responses <- function(data, method = c("three_point", "interval", "confidence", "linguistic"),
                              from = "rater", to = "rated", low = "low", best = "best", high = "high",
                              confidence = "confidence", label = "label", label2 = "label2",
                              scale = c(0, 1), rescale = c("minmax", "divide"), confidence_scale = c(0, 1),
                              max_spread = 0.25, dictionary = linguistic_scale(), fix_order = FALSE) {
  method <- match.arg(method); rescale <- match.arg(rescale); d <- as.data.frame(data)
  sc <- function(v) { v <- as.numeric(v); if (rescale == "minmax") (v - scale[1]) / (scale[2] - scale[1]) else v / scale[2] }
  if (method %in% c("three_point", "interval")) {
    l <- sc(d[[low]]); m <- sc(d[[best]]); r <- sc(d[[high]])
    bad <- which(l > m | m > r)
    if (length(bad)) {
      if (!fix_order) stop("responses with low > best or best > high in rows ", paste(utils::head(bad, 10), collapse = ", "),
                           "; use fix_order = TRUE to sort them")
      warning(length(bad), " responses were not ordered and have been sorted")
      s <- t(apply(cbind(l, m, r)[bad, , drop = FALSE], 1, sort)); l[bad] <- s[, 1]; m[bad] <- s[, 2]; r[bad] <- s[, 3]
    }
  } else if (method == "confidence") {
    m <- sc(d[[best]])
    cf <- (as.numeric(d[[confidence]]) - confidence_scale[1]) / diff(confidence_scale)
    if (any(cf < 0 | cf > 1, na.rm = TRUE)) stop("confidence outside 'confidence_scale'")
    lo_b <- if (rescale == "minmax") 0 else scale[1] / scale[2]
    l <- pmax(lo_b, m - (1 - cf) * max_spread); r <- pmin(1, m + (1 - cf) * max_spread)
  } else {
    lab1 <- trimws(as.character(d[[label]])); empty <- is.na(lab1) | lab1 == ""
    i1 <- match(tolower(lab1), tolower(dictionary$label))
    if (anyNA(i1[!empty])) stop("unknown labels: ", paste(unique(d[[label]][is.na(i1) & !empty]), collapse = ", "))
    l <- dictionary$l[i1]; m <- dictionary$m[i1]; r <- dictionary$r[i1]
    if (label2 %in% names(d)) {
      i2 <- match(tolower(trimws(d[[label2]])), tolower(dictionary$label))
      h <- !is.na(i2) & !is.na(i1)
      if (any(h)) {
        l[h] <- pmin(dictionary$l[i1[h]], dictionary$l[i2[h]])
        r[h] <- pmax(dictionary$r[i1[h]], dictionary$r[i2[h]])
        m[h] <- (dictionary$m[i1[h]] + dictionary$m[i2[h]]) / 2
      }
    }
  }
  out <- data.frame(from = d[[from]], to = d[[to]], l = l, m = m, r = r, stringsAsFactors = FALSE)
  extra <- setdiff(names(d), c(from, to, "l", "m", "r"))
  cbind(out, d[, extra, drop = FALSE])
}

#' Read a file of ratings
#'
#' Reads ratings exported by a data collection tool: CSV (comma or semicolon
#' separated, detected automatically), TSV or Excel (\code{.xlsx}, \code{.xls},
#' requires \pkg{readxl}). A URL is accepted for CSV files, e.g. the CSV export
#' of a published Google Sheet. Column names are converted to lower case and
#' common synonyms are recognised (\code{rater}/\code{from}/\code{respondent},
#' \code{rated}/\code{to}/\code{colleague}, \code{min}/\code{low},
#' \code{mode}/\code{best}/\code{most_likely}, \code{max}/\code{high}).
#' @param file Path or URL.
#' @param sheet Sheet of an Excel file.
#' @return A data frame with standardised column names.
#' @examples
#' f <- system.file("extdata", "three_point_example.csv", package = "FSNA")
#' d <- read_ratings(f)
#' head(d)
#' x <- fsn_from_ratings(fuzzify_responses(d, "three_point", scale = c(0, 10)),
#'                       drop_isolated = FALSE)
#' x
#' @export
#' @section Methodology:
#' Reads exports of spreadsheets and forms and standardises column names, so that the same pipeline (\code{\link{fuzzify_responses}}, \code{\link{fsn_from_ratings}}, \code{\link{fsna_analyze}}) applies to data collected with \code{\link{make_collector}}, Google Forms or other tools.
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' raw <- data.frame(Respondent = "a", Colleague = c("b", "c"),
#'                   Min = c(2, 4), Mode = c(3, 6), Max = c(5, 7))
#' write.csv(raw, tmp, row.names = FALSE)
#' read_ratings(tmp)
read_ratings <- function(file, sheet = 1) {
  ext <- tolower(tools::file_ext(sub("\\?.*$", "", file)))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("package 'readxl' is required to read Excel files")
    d <- as.data.frame(readxl::read_excel(file, sheet = sheet))
  } else {
    first <- readLines(file, n = 1, warn = FALSE)
    sep <- if (ext == "tsv" || grepl("\t", first)) "\t" else if (lengths(regmatches(first, gregexpr(";", first))) >
                                                             lengths(regmatches(first, gregexpr(",", first)))) ";" else ","
    d <- utils::read.csv(file, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  }
  nm <- tolower(trimws(names(d))); nm <- gsub("[^a-z0-9]+", "_", nm)
  syn <- c(from = "rater", respondent = "rater", rater_id = "rater", to = "rated", colleague = "rated",
           rated_id = "rated", min = "low", minimum = "low", lowest = "low", mode = "best",
           most_likely = "best", most_plausible = "best", value = "best", max = "high", maximum = "high",
           highest = "high")
  hit <- nm %in% names(syn); nm[hit] <- syn[nm[hit]]
  names(d) <- nm; d
}
