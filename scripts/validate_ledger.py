#!/usr/bin/env python3
"""Ledger integrity gate (added round-11, J3 repair — M3).

Asserts, on theorem_ledger_v2.json:
  L1  all depends_on endpoints exist (round-11 found nine dangling);
  L2  every entry carries status_class ∈ the declared status_classes
      (free-form status text is historical record, never rewritten);
  L3  entry IDs unique; count printed (NOT pinned — pinned counts are the
      job of round-pinned auditor snapshots, see the expected-red list);
  L4  the round-11 honesty edits hold: RECURSION-TERMINATION is 'broken'
      until M1 lands; the rev-3.5 LEAFKILL-CLASSIFICATION is 'superseded';
      TYPED-THEOREM-GRAPH's verification names 71 math edges.
"""
import json, sys

led = json.load(open('theorem_ledger_v2.json'))
entries = led['entries']
ids = [e['id'] for e in entries]
by = {e['id']: e for e in entries}
fails = []

def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + ' ' + name)
    if not ok:
        fails.append(name)

check('L3 ids unique', len(ids) == len(set(ids)))
print(f'      ({len(ids)} entries)')
dangling = [(e['id'], d) for e in entries for d in e.get('depends_on', []) if d not in by]
check('L1 zero dangling depends_on', not dangling)
if dangling:
    print('      dangling:', dangling)
classes = set(led.get('status_classes', []))
bad = [e['id'] for e in entries if e.get('status_class') not in classes]
check('L2 every entry has status_class in declared vocabulary', not bad)
if bad:
    print('      missing/invalid:', bad[:10])
check('L4a RECURSION-TERMINATION is proved (lexicographic v2, round-13 VERIFIED)',
      by['RECURSION-TERMINATION']['status_class'] == 'proved'
      and 'ROUND-13 VERIFIED' in by['RECURSION-TERMINATION']['status'])
check('L4b rev-3.5 LEAFKILL-CLASSIFICATION is superseded',
      by['LEAFKILL-CLASSIFICATION']['status_class'] == 'superseded')
check('L4c TYPED-THEOREM-GRAPH verification names 71 math edges',
      '71 math edges' in by['TYPED-THEOREM-GRAPH']['verification'])

print(('\nLEDGER GATE: ALL GREEN' if not fails else f'\nLEDGER GATE: {len(fails)} FAILURES'))
sys.exit(1 if fails else 0)
