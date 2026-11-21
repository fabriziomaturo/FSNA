/**
 * Google Apps Script backend for the FSNA collection form.
 *
 * 1. Create a Google Sheet; in Extensions > Apps Script paste this code.
 * 2. Deploy > New deployment > type "Web app"; Execute as: Me;
 *    Who has access: Anyone. Copy the web app URL.
 * 3. In R: make_collector(members, endpoint = "<web app URL>", file = "form.html")
 * 4. Publish the HTML (e.g. Google Sites > Insert > Embed > Embed code) and send
 *    the link to the respondents.
 * 5. Download the sheet as .xlsx or .csv (File > Download) and run
 *    d <- read_ratings("file.xlsx"); fsna_analyze(fuzzify_responses(d, "three_point", scale = c(0, 10)))
 */
function doPost(e) {
  var lock = LockService.getScriptLock();
  lock.waitLock(20000);
  try {
    var ss = SpreadsheetApp.getActiveSpreadsheet();
    var sheet = ss.getSheetByName("ratings") || ss.insertSheet("ratings");
    if (sheet.getLastRow() === 0) {
      sheet.appendRow(["group", "timestamp", "submission_id", "rater", "rated", "low", "best", "high"]);
    }
    var data = JSON.parse(e.postData.contents);
    var id = String(data.submission_id || "");
    // idempotency: the submission identifier is stored in every row written,
    // in a single write. A retry is recognised from the rows actually stored
    // (column C), so that it is not written again even if an earlier request
    // failed after writing. The check and the write happen under the same lock.
    if (id && sheet.getLastRow() > 1 &&
        sheet.getRange(2, 3, sheet.getLastRow() - 1, 1).createTextFinder(id).matchEntireCell(true).findNext()) {
      return ContentService.createTextOutput("duplicate");
    }
    var out = [];
    if (!data.rows || data.rows.length === 0) {
      // record participation of a respondent who rated nobody
      out.push([data.group || "", new Date().toISOString(), id, data.rater, "", "", "", ""]);
    } else {
      data.rows.forEach(function (r) {
        out.push([r.group, r.timestamp || new Date().toISOString(), id, r.rater, r.rated, r.low, r.best, r.high]);
      });
    }
    sheet.getRange(sheet.getLastRow() + 1, 1, out.length, out[0].length).setValues(out);
    return ContentService.createTextOutput("ok");
  } catch (err) {
    // the outcome may be uncertain: the client tells the respondent so
    return ContentService.createTextOutput("error");
  } finally {
    lock.releaseLock();
  }
}
