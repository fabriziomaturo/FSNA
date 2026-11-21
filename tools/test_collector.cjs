// Mock tests of the submission protocol of the FSNA collector (no network).
// Usage: node tools/test_collector.cjs inst/collector/collector_template.html
// The real rows() function of the template is used, with a clock that advances
// between clicks, so that timestamps differ across retries.
const fs = require('fs');
const s = fs.readFileSync(process.argv[2], 'utf8');
const rowsSrc = s.slice(s.indexOf('function rows() {'), s.indexOf('function toCSV('));
const sendSrc = s.slice(s.indexOf('let pending = null;'), s.indexOf('</script>', s.indexOf('let pending = null;')));
function run(replies) {
  const calls = [], msgs = []; let handler = null; let clock = Date.UTC(2026, 10, 20, 10, 0, 0);
  class FakeDate extends Date { constructor(...a) { super(...(a.length ? a : [clock])); clock += 1000; } static now() { return clock; } }
  const document = { getElementById: () => ({ addEventListener: (_, h) => { handler = h; } }) };
  const fetch = async (url, opts) => { calls.push(JSON.parse(opts.body)); const r = replies[calls.length - 1];
    if (r === 'throw') throw new Error('unreadable'); return { ok: r !== 'http500', text: async () => r }; };
  const state = { B: { on: true, low: 1, best: 2, high: 3 }, C: { on: true, low: 4, best: 5, high: 6 } };
  const fn = new Function('me', 'state', 'UI', 'confirm', 'CONFIG', 'fetch', 'say', 'document', 'window', 'Date',
    rowsSrc + '\n' + sendSrc + '\nreturn null;');
  fn({ value: 'A' }, state, { unconf: 'UNCONF ', thanks: 'THANKS ', fail: 'FAIL ', first: '', none: '', noend: '' }, () => true,
     { endpoint: 'https://mock.invalid', group: 'G1' }, fetch, (t, ok) => msgs.push(t.split(' ')[0]), document, { crypto: null }, FakeDate);
  return { calls, msgs, state, click: () => handler() };
}
(async () => {
  let fail = 0; const check = (c, m) => { if (!c) { fail++; console.log('FAIL:', m); } else console.log('PASS:', m); };
  let t = run(['throw', 'duplicate']);
  await t.click(); check(t.calls.length === 1 && t.msgs[0] === 'UNCONF', 'unreadable reply: one POST, receipt uncertain');
  await t.click();
  check(t.calls.length === 2 && t.calls[0].submission_id === t.calls[1].submission_id,
    'retry after uncertain receipt reuses the submission_id although rows() gives new timestamps');
  check(JSON.stringify(t.calls[0].rows) === JSON.stringify(t.calls[1].rows), 'retry resends the original rows (same timestamps)');
  check(t.msgs[1] === 'THANKS', 'duplicate answered as recorded');
  t = run(['error', 'ok']); await t.click();
  check(t.calls.length === 1 && t.msgs[0] === 'FAIL', 'backend error: one POST, uncertain-failure message');
  await t.click(); check(t.calls[0].submission_id === t.calls[1].submission_id, 'retry after backend error reuses the id');
  t = run(['throw', 'ok']); await t.click(); t.state.B.best = 2.5; await t.click();
  check(t.calls[0].submission_id !== t.calls[1].submission_id, 'changed ratings give a new submission');
  t = run(['ok', 'ok']); await t.click(); await t.click();
  check(t.calls[0].submission_id !== t.calls[1].submission_id, 'after a confirmed submission a new click is a new submission');
  process.exit(fail ? 1 : 0);
})();
