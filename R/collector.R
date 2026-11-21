# ---------------------------------------------------------------------------
# Web form for the collection of fuzzy ratings
# ---------------------------------------------------------------------------

#' Create a web form for collecting fuzzy collaboration ratings
#'
#' Writes a self-contained HTML page in which each respondent selects his or her
#' name and, for each colleague with whom he or she collaborates, gives the
#' lowest plausible, the most plausible and the highest plausible value of the
#' collaboration on a numeric scale (three-point elicitation of a triangular
#' possibility distribution). The triangle is drawn while the sliders move.
#' Answers are sent to a Google Sheet through a Google Apps Script web app
#' (\code{endpoint}; the script is in
#' \code{system.file("collector", "google_apps_script.gs", package = "FSNA")}),
#' or downloaded by the respondent as a CSV file when no endpoint is given.
#'
#' The page can be opened locally, hosted on any web server, or embedded in a
#' Google Site (Insert > Embed > Embed code, pasting the content of the file).
#' The collected sheet, downloaded as Excel or CSV, is read with
#' \code{\link{read_ratings}}, converted with \code{\link{fuzzify_responses}}
#' (\code{method = "three_point"}) and analysed with \code{\link{fsna_analyze}}.
#'
#' @param members Character vector with the names (or codes) of the members.
#' @param file Output HTML file.
#' @param title,question Texts shown on the page.
#' @param scale Numeric range of the sliders.
#' @param step Slider step.
#' @param low_label,high_label Meaning of the ends of the scale.
#' @param endpoint URL of the Google Apps Script web app (empty: CSV download only).
#' @param group Identifier of the group or department, stored with the answers.
#' @param language Language of the interface: \code{"en"} or \code{"it"}.
#' @param open Logical, open the page in the browser.
#' @return The path of the file, invisibly.
#' @examples
#' f <- make_collector(c("Anna", "Bruno", "Carla", "Dario"),
#'                     file = file.path(tempdir(), "form.html"))
#' f
#' @export
#' @section Methodology:
#' The form implements three-point elicitation with sliders constrained to \eqn{low \le best \le high} and a live drawing of the triangle, so that respondents see the possibility distribution they state. Respondents tick only the colleagues they collaborate with; unticked colleagues produce no record, so the analyst must decide whether they are null or unknown ties (\code{\link{set_unknown}}).
#' @examples
#' f <- make_collector(paste("Member", 1:8), file = file.path(tempdir(), "form8.html"),
#'                     title = "Research collaboration", scale = c(1, 5), step = 0.25,
#'                     low_label = "never", high_label = "every week")
#' file.exists(f)
make_collector <- function(members, file = "fsna_form.html", title = "Collaboration questionnaire",
                           question = "How much do you collaborate with each colleague?", scale = c(0, 10),
                           step = 0.5, low_label = "no collaboration", high_label = "very strong collaboration",
                           endpoint = "", group = "G1", language = c("en", "it"), open = FALSE) {
  language <- match.arg(language)
  tpl <- system.file("collector", "collector_template.html", package = "FSNA")
  if (!nzchar(tpl)) stop("template not found")
  h <- paste(readLines(tpl, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  esc <- function(s) paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", s)), '"')
  ui <- if (language == "it") c(who = "Chi sei?", sel = "-- seleziona il tuo nome --",
    priv = "Le risposte sono collegate al tuo nome solo per costruire la rete; i risultati sono riportati senza nomi.",
    howto = "Per ogni collega con cui collabori, spunta <b>Valuta</b> e muovi i tre cursori: il valore <b>pi\u00f9 plausibile</b> della collaborazione e i valori <b>minimo</b> e <b>massimo</b> che consideri ancora plausibili. Pi\u00f9 l'intervallo \u00e8 ampio, meno sei sicuro. Il triangolo a destra mostra la tua risposta. Lascia senza spunta i colleghi con cui non collabori.",
    rate = "Valuta", low = "Minimo plausibile", best = "Pi\u00f9 plausibile", high = "Massimo plausibile", send = "Invia",
    csv = "Scarica le mie risposte (CSV)", first = "Seleziona prima il tuo nome.",
    none = "Non hai valutato alcun collega. Inviare comunque?",
    noend = "Nessun server di raccolta configurato: usa 'Scarica le mie risposte' e invia il file all'organizzatore.",
    thanks = "Grazie: ricezione confermata dal server. Valutazioni: ",
    unconf = "Ricezione non confermata: le risposte potrebbero essere state registrate oppure no. Puoi premere di nuovo Invia (il server ignora i duplicati) e scaricare le risposte. Valutazioni: ",
    fail = "Il server ha segnalato un errore: le risposte potrebbero non essere state registrate. Premi di nuovo Invia (il server ignora i duplicati) oppure scarica le risposte e invia il file.") else NULL
  uij <- if (is.null(ui)) "" else paste0(',"ui":{', paste(sprintf('"%s":%s', names(ui), vapply(ui, esc, "")), collapse = ","), "}")
  cfg <- sprintf('{"title":%s,"question":%s,"members":[%s],"scale":[%s,%s],"step":%s,"low_label":%s,"high_label":%s,"endpoint":%s,"group":%s%s}',
                 esc(title), esc(question), paste(esc(as.character(members)), collapse = ","), scale[1], scale[2], step,
                 esc(low_label), esc(high_label), esc(endpoint), esc(group), uij)
  h <- sub("/\\*CONFIG\\*/.*?/\\*END\\*/", paste0("/*CONFIG*/", gsub("\\\\", "\\\\\\\\", cfg), "/*END*/"), h, perl = TRUE)
  writeLines(h, file, useBytes = TRUE)
  if (open) utils::browseURL(file)
  invisible(normalizePath(file))
}
