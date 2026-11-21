# ---------------------------------------------------------------------------
# Complete analysis of a fuzzy social network with figures, tables and report
# ---------------------------------------------------------------------------

.html_table <- function(d, digits = 3, caption = NULL) {
  d <- as.data.frame(d)
  for (k in seq_along(d)) if (is.numeric(d[[k]])) d[[k]] <- formatC(d[[k]], digits = digits, format = "fg", flag = "#")
  esc <- function(s) gsub("<", "&lt;", gsub("&", "&amp;", s))
  paste0("<table>", if (!is.null(caption)) paste0("<caption>", esc(caption), "</caption>"),
         "<tr>", paste0("<th>", esc(names(d)), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")), collapse = ""),
         "</table>")
}

#' Complete analysis of a fuzzy social network
#'
#' Runs, in one call, the analyses implemented in the package and saves figures
#' (PNG), tables (CSV) and an HTML report:
#' \enumerate{
#'   \item description of the network and of the ties (size, density,
#'     respondents, reciprocity, widths and asymmetry of the TFNs), network plot;
#'   \item the indices of Porreca, Maturo and Ventre (2025): fuzzy in- and
#'     out-degree, in- and out-closeness and betweenness on best paths, with the
#'     fuzzy scatter plots of their Figs. 10 and 11;
#'   \item for each isotone index, the exact fuzzy centrality: membership
#'     functions of the most central actors with the triangular approximation,
#'     interval plot, rank intervals, fuzzy set of the most central actor,
#'     proportion of certified orderings, necessity heat map, error of the
#'     triangular approximation and comparison of nine ranking methods.
#' }
#' @param x An \code{\link{fsn}} object, a data frame of ratings with columns
#'   \code{from}/\code{to} (or \code{rater}/\code{rated}) and \code{l},
#'   \code{m}, \code{r}, or the path of a file readable by
#'   \code{\link{read_ratings}} containing such columns.
#' @param output_dir Directory for the outputs (created if needed); \code{NULL}
#'   for no files.
#' @param indices Isotone indices for the exact analysis (see
#'   \code{\link{centrality_index}}).
#' @param ijar Logical, compute the indices of Porreca, Maturo and Ventre (2025).
#' @param ijar_variant \code{"definition"} or \code{"alternative"}, see
#'   \code{\link{ijar_closeness}}.
#' @param max_steps Maximum number of steps of best paths (IJAR indices and
#'   t-norm closeness).
#' @param alpha Grid of levels.
#' @param cut Level used for cuts, rank intervals and top sets in tables.
#' @param rank_methods Ranking methods compared.
#' @param top Number of actors drawn in membership plots.
#' @param title Title of the report.
#' @param open Logical, open the report in the browser.
#' @param pdf Logical, also write all the figures in a single file
#'   \code{figures.pdf} in \code{output_dir}.
#' @param person_col Optional colours of the person icons in the plot of
#'   \code{\link{plot_fsn_glyphs}} (drawn for networks with at most 15 actors).
#' @return Invisibly, a list with all the results (network, IJAR indices, exact
#'   fuzzy centralities, summaries, comparisons, certification) and the paths of
#'   the files written.
#' @examples
#' res <- fsna_analyze(fsn_example(), output_dir = file.path(tempdir(), "fsna_example"))
#' names(res)
#' \donttest{
#' res <- fsna_analyze(fsn_simulated(n = 20, respondents = 10),
#'                     output_dir = file.path(tempdir(), "fsna_sim"))
#' res$files$report
#' }
#' @export
#' @section Methodology:
#' The function chains the steps of the analysis described in the article: description of the ties, the indices of Porreca, Maturo and Ventre (2025), exact fuzzy centralities, comparisons and rankings. The outputs are descriptive; they do not provide inferential statements about a population.
#' @references Porreca A., Maturo F., Ventre V. (2025). Fuzzy centrality measures in social network analysis: Theory and application in a university department collaboration network. \emph{International Journal of Approximate Reasoning} 176, 109319. \doi{10.1016/j.ijar.2024.109319}
#' @examples
#' res <- fsna_analyze(fsn_example(), output_dir = NULL, indices = "in_closeness_prod")
#' res$certification_summary
#' names(res$exact)
fsna_analyze <- function(x, output_dir = NULL, indices = c("in_strength", "out_strength", "in_closeness_prod",
                                                           "out_closeness_prod", "in_katz"),
                         ijar = TRUE, ijar_variant = c("definition", "alternative"), max_steps = 2,
                         alpha = seq(0, 1, by = 0.05), cut = 0.5,
                         rank_methods = c("cog", "core", "yager", "graded_mean", "liou_wang", "possibility",
                                          "necessity", "certified", "top"),
                         top = 6, title = "Fuzzy social network analysis", open = FALSE, pdf = FALSE,
                         person_col = NULL) {
  ijar_variant <- match.arg(ijar_variant)
  if (is.character(x) && length(x) == 1) x <- read_ratings(x)
  if (is.data.frame(x)) {
    f <- if ("from" %in% names(x)) "from" else "rater"; t <- if ("to" %in% names(x)) "to" else "rated"
    x <- fsn_from_ratings(x, from = f, to = t)
  }
  if (!is.fsn(x)) stop("'x' must be an fsn object, a data frame of ratings or a file")
  files <- list(); fig <- character(0); tabs <- list()
  save <- !is.null(output_dir)
  if (save) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_dev <- NULL
  if (save && pdf) {
    grDevices::pdf(file.path(output_dir, "figures.pdf"), width = 8, height = 6.5, onefile = TRUE)
    pdf_dev <- grDevices::dev.cur()
    on.exit(if (!is.null(pdf_dev) && pdf_dev %in% grDevices::dev.list()) grDevices::dev.off(pdf_dev), add = TRUE)
  }
  png_fig <- function(name, expr, w = 8, h = 5.5) {
    if (!save) return(invisible(NULL))
    e <- substitute(expr); env <- parent.frame()
    f <- file.path(output_dir, paste0(name, ".png"))
    # a figure that cannot be drawn (e.g. no ties) is replaced by a note; the
    # device is closed also on error
    draw <- function() {
      ok <- tryCatch({ eval(e, env); TRUE }, error = function(err) {
        graphics::par(mfrow = c(1, 1)); graphics::plot.new()
        graphics::text(0.5, 0.5, paste("Figure not available:", conditionMessage(err)), cex = 0.8)
        warning(sprintf("figure '%s' not drawn: %s", name, conditionMessage(err)), call. = FALSE)
        FALSE
      })
      ok
    }
    grDevices::png(f, width = w, height = h, units = "in", res = 150)
    dev_png <- grDevices::dev.cur()
    local({
      on.exit(if (dev_png %in% grDevices::dev.list()) grDevices::dev.off(dev_png))
      graphics::par(mar = c(4.2, 4.2, 2, 1)); draw()
    })
    if (!is.null(pdf_dev)) {
      grDevices::dev.set(pdf_dev)
      graphics::par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))   # one figure per page
      draw()
      graphics::par(mfrow = c(1, 1))
    }
    fig[name] <<- basename(f)
  }
  csv <- function(d, name) {
    tabs[[name]] <<- d
    if (save) utils::write.csv(d, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  res <- list(network = x, summary = summary(x))
  # 1. network
  png_fig("Fig01_network", plot(x, vertex_size = colSums(x$M), main = "Network (arrow width = core)"), 7, 7)
  glyph_min_core <- if (n_actors(x) > 6) 0.4 else 0
  if (n_actors(x) <= 15)
    png_fig("Fig01b_glyphs", plot_fsn_glyphs(x, values = n_actors(x) <= 6, min_core = glyph_min_core,
                                             glyph_size = if (n_actors(x) <= 6) 0.19 else 0.13, curvature = 0.25,
                                             glyph_pos = 0.35, person_col = person_col,
                                             main = if (glyph_min_core > 0) sprintf("Fuzzy ties with core >= %.1f", glyph_min_core)
                                                    else "Fuzzy ties"), 7, 7)
  ties_df <- fsn_ties(x)
  if (nrow(ties_df)) {
  sp <- tfn_spread(tfn(ties_df$l, ties_df$m, ties_df$r))
  png_fig("Fig02_tie_widths", {
    graphics::par(mfrow = c(1, 2))
    graphics::hist(sp$width, breaks = 20, col = "#56B4E9", border = "white", main = "Support width r - l", xlab = "")
    graphics::hist(sp$asymmetry, breaks = 20, col = "#E69F00", border = "white", main = "Asymmetry (r-m)-(m-l)", xlab = "")
  }, 9, 4)
  }
  # 2. IJAR indices
  if (ijar) {
    paths <- ijar_best_paths(x, max_steps, "cog", .ijar_variant(ijar_variant)$path_tfn, .ijar_variant(ijar_variant)$tie_break)
    ij <- list(in_degree = ijar_degree(x, "in"), out_degree = ijar_degree(x, "out"),
               in_closeness = ijar_closeness(x, "in", max_steps, ijar_variant, paths = paths),
               out_closeness = ijar_closeness(x, "out", max_steps, ijar_variant, paths = paths),
               betweenness = ijar_betweenness(x, max_steps, if (max_steps == 2) ijar_variant else "definition"))
    res$ijar <- ij
    tab <- data.frame(actor = x$labels)
    for (nm in names(ij)[1:4]) {
      d <- ij[[nm]]; tab[[paste0(nm, "_l")]] <- d$l; tab[[paste0(nm, "_m")]] <- d$m
      tab[[paste0(nm, "_r")]] <- d$r; tab[[paste0(nm, "_cog")]] <- d$score
    }
    tab$betweenness <- ij$betweenness$betweenness
    csv(tab[order(-tab$in_degree_cog), ], "Tab01_IJAR_indices")
    png_fig("Fig03_IJAR_in_degree", plot(ij$in_degree, n = 10, main = "Fuzzy in-degree (FTOWA, top 10 by CoG)"))
    png_fig("Fig04_IJAR_out_degree", plot(ij$out_degree, n = 10, main = "Fuzzy out-degree (FTOWA, top 10 by CoG)"))
    png_fig("Fig05_IJAR_betweenness", plot(ij$betweenness, n = min(30, n_actors(x)), main = "Betweenness on best paths"))
    png_fig("Fig06_IJAR_scatter_degree", plot_fuzzy_scatter(ij$in_degree, ij$out_degree, colour = ij$betweenness$betweenness,
            xlab = "Fuzzy in-degree", ylab = "Fuzzy out-degree", colour_label = "Betweenness"), 8, 6.5)
    png_fig("Fig07_IJAR_scatter_closeness", plot_fuzzy_scatter(ij$in_closeness, ij$out_closeness, colour = ij$betweenness$betweenness,
            xlab = "Fuzzy in-closeness", ylab = "Fuzzy out-closeness", colour_label = "Betweenness"), 8, 6.5)
  }
  # 3. exact fuzzy centralities
  res$exact <- list(); res$comparison <- list(); res$certification <- list(); k <- 7
  cert_rows <- list()
  for (ix in indices) {
    args <- list(x = x, index = ix, alpha = alpha)
    if (grepl("closeness", ix)) args$max_steps <- max_steps
    fc <- do.call(fuzzy_centrality, args)
    res$exact[[ix]] <- fc
    sm <- summary(fc, alpha = cut); csv(as.data.frame(unclass(sm)), paste0("Tab_exact_", ix))
    cmp <- compare_rankings(fc, rank_methods); res$comparison[[ix]] <- cmp
    csv(data.frame(actor = rownames(cmp$ranks), cmp$ranks, check.names = FALSE), paste0("Tab_ranks_", ix))
    ce <- certification(fc, c(0, cut, 1)); res$certification[[ix]] <- ce
    ta <- triangular_approx(fc)
    cert_rows[[ix]] <- data.frame(index = fc$index$label, certified_0 = ce$proportion_certified[1],
                                  certified_cut = ce$proportion_certified[2], N1_not_separated = ce$N1_not_separated,
                                  N1_not_separated_grid = ce$N1_not_separated_grid,
                                  top_set_cut_size = length(top_set(fc, alpha = cut)),
                                  max_triangle_error = max(ta$max_error),
                                  kendall_cog_necessity = if (all(c("cog", "necessity") %in% rank_methods))
                                    cmp$correlation["cog", "necessity"] else NA_real_)
    k <- k + 1; png_fig(sprintf("Fig%02d_exact_membership_%s", k, ix), plot(fc, actors = top, main = fc$index$label))
    k <- k + 1; png_fig(sprintf("Fig%02d_exact_intervals_%s", k, ix), plot(fc, "intervals", alpha = cut, main = fc$index$label),
                        8, max(5, 0.16 * n_actors(x)))
    k <- k + 1; png_fig(sprintf("Fig%02d_necessity_%s", k, ix), plot(possibility_matrix(fc), main = paste("Necessity,", fc$index$label)), 8, 8)
  }
  cert <- do.call(rbind, cert_rows); rownames(cert) <- NULL; csv(cert, "Tab02_certification_summary")
  res$certification_summary <- cert
  # 4. report
  if (save) {
    html <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><title>", title, "</title><style>",
              "body{font-family:Georgia,serif;max-width:1000px;margin:2em auto;color:#222;line-height:1.45}",
              "h1,h2,h3{font-family:Helvetica,Arial,sans-serif}table{border-collapse:collapse;font-size:12px;margin:1em 0}",
              "td,th{border-bottom:1px solid #ddd;padding:3px 8px;text-align:right}caption{text-align:left;font-weight:bold}",
              "img{max-width:100%;border:1px solid #eee;margin:0.5em 0}</style></head><body>",
              "<h1>", title, "</h1>",
              sprintf("<p>%d actors, %d non-null ties, density %.3f%s. Generated by FSNA %s on %s.</p>", n_actors(x),
                      res$summary$n_ties, res$summary$density,
                      if (!is.na(res$summary$respondents)) sprintf(", %d respondents", res$summary$respondents) else "",
                      utils::packageVersion("FSNA"), format(Sys.time(), "%Y-%m-%d %H:%M")),
              "<h2>1. Network and ties</h2>", sprintf("<img src='%s'>", fig[intersect(c("Fig01_network", "Fig01b_glyphs", "Fig02_tie_widths"), names(fig))]))
    if (ijar) html <- c(html, "<h2>2. Indices of Porreca, Maturo and Ventre (2025)</h2>",
                        sprintf("<p>Variant: %s; best paths with %s; ranking by centre of gravity.</p>", ijar_variant,
                                if (is.finite(max_steps)) sprintf("at most %d steps", as.integer(max_steps)) else "any number of steps"),
                        .html_table(utils::head(tabs$Tab01_IJAR_indices[, c("actor", "in_degree_l", "in_degree_m", "in_degree_r",
                          "out_degree_m", "in_closeness_m", "out_closeness_m", "betweenness")], 15), caption = "First 15 actors by fuzzy in-degree (CoG)"),
                        sprintf("<img src='%s'>", fig[grepl("IJAR", names(fig))]))
    html <- c(html, "<h2>3. Exact fuzzy centralities</h2>",
              "<p>Alpha-cuts computed exactly from two crisp evaluations per level (extension principle). Certified: share of pairs of actors whose supports (alpha = 0) or alpha-cuts are strictly separated.</p>",
              .html_table(cert, caption = "Certification and approximation by index"))
    for (ix in indices) {
      html <- c(html, sprintf("<h3>%s</h3>", res$exact[[ix]]$index$label),
                .html_table(utils::head(tabs[[paste0("Tab_exact_", ix)]], 15), caption = sprintf("First 15 actors; cuts at alpha = %s", cut)),
                sprintf("<img src='%s'>", fig[grepl(paste0("_", ix, "$"), names(fig))]))
    }
    html <- c(html, "</body></html>")
    rep <- file.path(output_dir, "report.html"); writeLines(html, rep)
    if (!is.null(pdf_dev)) { grDevices::dev.off(pdf_dev); pdf_dev <- NULL }
    files <- list(report = rep, figures = file.path(output_dir, fig),
                  pdf = if (pdf) file.path(output_dir, "figures.pdf") else NULL,
                  tables = file.path(output_dir, paste0(names(tabs), ".csv")))
    if (open) utils::browseURL(rep)
  }
  res$files <- files
  invisible(res)
}


#' Demonstration: simulate a small network and run all the analyses
#'
#' Simulates a small fuzzy social network with named actors in two working
#' groups (stochastic block model of \code{\link{simulate_fsn}}, with vague and
#' asymmetric ratings), runs \code{\link{fsna_analyze}} on it and opens
#' the HTML report in the browser. It is the quickest way to see what the
#' package produces.
#' @param n Number of actors (default 8; the plot with person icons is
#'   included up to 15 actors).
#' @param p Tie probability: within groups 2p, between groups p/3.
#' @param kappa Spread multiplier of the simulated ratings.
#' @param seed Random seed.
#' @param output_dir Directory of the outputs (default a temporary directory).
#' @param open Logical, open the report in the browser (default in interactive
#'   sessions).
#' @param pdf Logical, also write all the figures in \code{figures.pdf}.
#' @param colours Colours of the person icons, one per actor (default: a
#'   different colour for each actor).
#' @param ... Further arguments to \code{\link{fsna_analyze}}.
#' @return Invisibly, the result of \code{\link{fsna_analyze}}, with the
#'   simulated network in \code{network}.
#' @section Methodology:
#' The network is simulated with a latent crisp network and vague, possibly
#' erroneous, fuzzy ratings of it (\code{\link{simulate_fsn}}), and is analysed
#' with the indices of Porreca, Maturo and Ventre (2025) and with the exact
#' fuzzy centralities.
#' @references Porreca A., Maturo F. (2026). FSNA: an R package for centrality
#'   analysis in social networks with triangular fuzzy ties. Manuscript; see
#'   \code{citation("FSNA")}.
#' @examples
#' res <- fsna_demo(n = 7, open = FALSE)
#' res$files$report
#' \donttest{
#' res <- fsna_demo(n = 6, open = FALSE, pdf = FALSE,
#'                  output_dir = file.path(tempdir(), "fsna_demo"))
#' res$files$report
#' }
#' if (interactive()) fsna_demo()                # opens the report in the browser
#' @export
fsna_demo <- function(n = 8, p = 0.3, kappa = 2, seed = 1, output_dir = file.path(tempdir(), "fsna_demo"),
                      open = interactive(), pdf = TRUE, colours = NULL, ...) {
  nm <- c("Anna", "Bruno", "Carla", "Dario", "Elena", "Fabio", "Giulia", "Luca", "Marta", "Nicola",
          "Olga", "Paolo", "Rita", "Sergio", "Teresa")
  if (n < 3) stop("n must be at least 3")
  # two working groups: dense, strong ties within groups, few weak ties between;
  # vague and asymmetric ratings, so that the triangles are visible
  x <- simulate_fsn(n, model = "sbm", K = 2, p_in = min(1, 2 * p), p_out = p / 3, kappa = kappa,
                    sigma = c(0.05, 0.12), asymmetry = 0.5, seed = seed)
  o <- order(x$blocks)                       # actors of the same group are adjacent
  x <- x[o]
  lab <- if (n <= length(nm)) nm[seq_len(n)] else sprintf("A%02d", seq_len(n))
  x <- fsn(x$L, x$M, x$R, labels = lab,
           info = list(description = sprintf("Simulated example (n = %d, p = %g, kappa = %g, seed = %d)", n, p, kappa, seed)))
  if (is.null(colours)) {
    base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#999999")
    colours <- if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base[1:6])(n)
  }
  res <- fsna_analyze(x, output_dir = output_dir, open = open, pdf = pdf, person_col = colours,
                      title = sprintf("FSNA demonstration: simulated network of %d actors", n), ...)
  if (!is.null(res$files$report)) message("Report: ", res$files$report)
  invisible(res)
}
