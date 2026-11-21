// Mock test of the Apps Script backend (no Google services).
// Usage: node tools/test_backend.cjs inst/collector/google_apps_script.gs
const fs = require('fs'), vm = require('vm');
function makeSheet(failAfterWrite) {
  const sh = { data: [], failNext: failAfterWrite,
    getLastRow() { return this.data.length; },
    appendRow(row) { this.data.push(row); },
    getRange(row, col, n, m) { const self = this; return {
      setValues(values) { for (let i = 0; i < n; i++) self.data[row - 1 + i] = values[i];
        if (self.failNext) { self.failNext = false; throw new Error('failure after the rows were stored'); } },
      createTextFinder(id) { return { matchEntireCell() { return this; },
        findNext() { for (let i = 0; i < n; i++) if (String(self.data[row - 1 + i][col - 1]) === id) return {}; return null; } }; } }; } };
  return sh;
}
function backend(failAfterWrite) {
  const ratings = makeSheet(failAfterWrite);
  const ctx = { LockService: { getScriptLock: () => ({ waitLock() {}, releaseLock() {} }) },
    SpreadsheetApp: { getActiveSpreadsheet: () => ({ getSheetByName: n => (n === 'ratings' ? ratings : null), insertSheet: () => ratings }) },
    ContentService: { createTextOutput: x => x }, Date, JSON, String };
  vm.createContext(ctx); vm.runInContext(fs.readFileSync(process.argv[2], 'utf8'), ctx);
  return { ctx, ratings };
}
const ev = id => ({ postData: { contents: JSON.stringify({ submission_id: id, group: 'G1', rater: 'A',
  rows: [{ group: 'G1', timestamp: 't0', rater: 'A', rated: 'B', low: 1, best: 2, high: 3 },
         { group: 'G1', timestamp: 't0', rater: 'A', rated: 'C', low: 4, best: 5, high: 6 }] }) } });
let fail = 0; const check = (c, m) => { if (!c) { fail++; console.log('FAIL:', m); } else console.log('PASS:', m); };
let b = backend(false);
check(b.ctx.doPost(ev('id1')) === 'ok' && b.ratings.data.length === 3, 'first submission recorded (header + 2 rows)');
check(b.ctx.doPost(ev('id1')) === 'duplicate' && b.ratings.data.length === 3, 'same id: duplicate, nothing written');
check(b.ctx.doPost(ev('id2')) === 'ok' && b.ratings.data.length === 5, 'new id: recorded');
b = backend(true);
const r1 = b.ctx.doPost(ev('id9')), n1 = b.ratings.data.length;
const r2 = b.ctx.doPost(ev('id9')), n2 = b.ratings.data.length;
check(r1 === 'error' && n1 === 3, 'error after the write: reported as error, rows stored once');
check(r2 === 'duplicate' && n2 === 3, 'retry after a partial failure is recognised from the stored rows');
process.exit(fail ? 1 : 0);
