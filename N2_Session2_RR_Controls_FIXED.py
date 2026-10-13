#!/usr/bin/env python3
r"""N2 Session 2: frozen annotation normalization and synthetic RR controls.
No primary real-data RR, selection, MAE or bias is computed here.
Python 3.9+; standard library only.
Usage: py -3 .\N2_Session2_RR_Controls.py --input-dir .\N2_inputs_v1.0.0
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

HEADERS = ('breaths ann1 [signal sample no]', 'breaths ann2 [signal sample no]')
FIRST, LAST, STEP, WIDTH = 0, 149, 3, 32

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(65536), b''):
            h.update(block)
    return h.hexdigest()

def in_window(sample, start_s, end_s):
    """Sample ID is 1-based and positive; [start,end) exact INTEGER test."""
    return 125 * start_s + 1 <= sample < 125 * end_s + 1

def rate(samples, start_s, end_s):
    """Samples must be normalized increasing integers. Window-exclusive inputs ignored."""
    s = [x for x in samples if in_window(x, start_s, end_s)]
    if len(s) < 2:
        return None
    if any(b <= a for a, b in zip(s, s[1:])):
        raise ValueError('non-increasing sample IDs')
    return 60.0 * 125.0 * (len(s) - 1) / (s[-1] - s[0])

def reference(rr1, rr2):
    """Return agreement flag and selected RR; NEVER threshold monitor-reference."""
    if rr1 is None or rr2 is None:
        return None, None
    if not (math.isfinite(rr1) and math.isfinite(rr2)):
        raise ValueError('nonfinite annotation RR')
    passed = abs(rr1 - rr2) <= 2.0
    return passed, (rr1 + rr2) / 2.0 if passed else None

def normalize(cells, record_id, annotator):
    raw = []
    last_nonpadding = -1
    for i, value in enumerate(cells, start=2):  # CSV line number, header line 1
        v = value.strip()
        if v and v.lower() != 'nan':
            last_nonpadding = i
    for line, value in enumerate(cells, start=2):
        v = value.strip()
        if not v or v.lower() == 'nan':
            if line <= last_nonpadding:
                raise ValueError(f'{record_id} ann{annotator}: internal NaN/padding, CSV line {line}')
            continue
        try:
            number = float(v)
        except ValueError as exc:
            raise ValueError(f'{record_id} ann{annotator}: invalid sample at line {line}') from exc
        if not math.isfinite(number) or not number.is_integer() or not 1 <= number <= 60001:
            raise ValueError(f'{record_id} ann{annotator}: out-of-domain sample at line {line}: {v}')
        raw.append((int(number), line, v))
    ordered = sorted(raw, key=lambda z: z[0])
    sorted_samples = [x[0] for x in ordered]
    unique = sorted(set(sorted_samples))
    duplicate = len(raw) - len(unique)
    inversions = sum(raw[i][0] > raw[i+1][0] for i in range(len(raw)-1))
    changed = duplicate > 0 or inversions > 0
    if changed:
        affected = []
        # Track any window that contains a sample from a duplicate group or a storage-order inversion pair.
        duplicate_samples = {x for x, ct in Counter(x[0] for x in raw).items() if ct > 1}
        implicated = set(duplicate_samples)
        for a, b in zip(raw, raw[1:]):
            if a[0] > b[0]:
                implicated.add(a[0]); implicated.add(b[0])
        for k in range(FIRST, LAST+1):
            start = k*STEP
            if any(in_window(x, start, start+WIDTH) for x in implicated):
                affected.append(k)
        entry = {
            'record_id': record_id, 'annotator': annotator,
            'raw_event_count': len(raw), 'normalized_event_count': len(unique),
            'exact_duplicate_instances_removed': duplicate,
            'raw_adjacent_order_inversions': inversions,
            'raw_rows_samples': [{'csv_line':ln, 'raw':val, 'sample':sm} for sm,ln,val in raw],
            'operations': (['sort_by_sample'] if inversions else []) + (['collapse_identical_sample'] if duplicate else []),
            'implicated_samples': sorted(implicated),
            'potentially_affected_primary_window_k': affected,
            'potentially_affected_primary_window_count': len(affected),
            'affected_windows_interpretation': 'potential membership, not actual usable RR changes'
        }
        return unique, entry
    return unique, None

def controls():
    count = 0
    def ck(test, desc):
        nonlocal count
        count += 1
        if not test:
            raise AssertionError(desc)
    ck(rate([], 0, 32) is None, 'n0')
    ck(rate([1], 0, 32) is None, 'n1')
    ck(math.isclose(rate([1, 126], 0, 32), 60.0), 'n2: one second')
    ck(math.isclose(rate([1, 126, 376], 0, 32), 40.0), 'mean intervals: 1,2 seconds')
    ck(in_window(1, 0, 32) and in_window(4000, 0, 32), 'last sample 31.992')
    ck(not in_window(4001, 0, 32), 'right boundary excluded')
    ck(not in_window(1, 3, 35) and in_window(376, 3, 35), 'left boundary inclusive')
    ck(rate([1, 126, 4001], 0, 32) == 60.0, 'outside events ignored')
    ck(rate([1, 4001], 0, 32) is None, 'no boundary-crossing interval')
    ck(math.isclose((376-1)/(125*2), 1.5), 'interval telescope')
    ck(reference(12.0,14.0) == (True,13.0), 'agreement exact equality 2')
    ck(reference(12.0,14.0000001) == (False,None), 'unrounded agreement')
    ck(reference(None,14.0) == (None,None), 'not evaluable')
    ck(reference(12.0,11.0) == (True,11.5), 'mean of RR values')
    ck(len(range(0, 480-32+1,3)) == 150, 'primary window count')
    ck(len(range(0, 480-64+1,6)) == 70, 'sensitivity window count')
    u, op = normalize(['126', '1', '126', 'nan'], 'test', 1)
    ck(u == [1,126] and op['exact_duplicate_instances_removed'] == 1 and op['raw_adjacent_order_inversions'] == 1, 'normalization')
    return count

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', default='N2_inputs_v1.0.0')
    p.add_argument('--verification', default=None, help='Session1 verification JSON (defaults inside input-dir)')
    p.add_argument('--output-dir', default='N2_Session2_output')
    a = p.parse_args()
    base = Path(a.input_dir)
    verification = Path(a.verification) if a.verification else base/'N2_Session1_Verification.json'
    if not verification.is_file():
        raise FileNotFoundError(f'Missing verification JSON: {verification}')
    manifest = json.loads(verification.read_text(encoding='utf-8'))
    if manifest.get('status') != 'PASS' or manifest.get('passed_files') != 106 or manifest.get('record_count') != 53:
        raise ValueError('Session 1 verification failed/incomplete')
    expected = {row['file']:row for row in manifest['files']}
    names = {f'bidmc_{rid:02d}_{kind}.csv' for rid in range(1,54) for kind in ('Breaths','Numerics')}
    if set(expected) != names or len(expected) != 106:
        raise ValueError('Expected CSV filename set mismatch')
    for name in sorted(names):
        file = base/name
        if not file.is_file() or sha256(file) != expected[name]['sha256'] or file.stat().st_size != expected[name]['bytes']:
            raise ValueError(f'Input missing or hash/size changed: {file}')
    ntests = controls()
    log = []
    normalized_total = 0
    for rid in range(1,54):
        with (base/f'bidmc_{rid:02d}_Breaths.csv').open(newline='',encoding='utf-8-sig') as fh:
            reader = csv.DictReader(fh)
            # Official BIDMC CSV has a space after the delimiter in its header.
            # Match the pre-protocol audit: strip header whitespace only.
            if reader.fieldnames is None:
                raise ValueError(f'Missing Breaths headers: {rid}')
            reader.fieldnames = [name.strip() for name in reader.fieldnames]
            if reader.fieldnames != list(HEADERS):
                raise ValueError(f'Unexpected Breaths headers: {rid}: {reader.fieldnames!r}')
            rows = list(reader)
        for j, header in enumerate(HEADERS, 1):
            unique, change = normalize([row[header] for row in rows], f'{rid:02d}', j)
            normalized_total += len(unique)
            if change: log.append(change)
    actual = {(row['record_id'],row['annotator']): (row['exact_duplicate_instances_removed'],row['raw_adjacent_order_inversions']) for row in log}
    expected_changes = {('06',1):(1,0),('20',1):(1,0),('09',2):(0,1),('19',2):(0,1)}
    if actual != expected_changes:
        raise ValueError(f'Unexpected normalization anomalies: {actual} (expected {expected_changes})')
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    report = {'stage':'Session 2 ONLY: annotation normalization and synthetic controls; NO real-data RR computed',
        'protocol':'N2 Minimum Research Protocol v1.0 FROZEN 2026-10-12',
        'python':sys.version, 'script_sha256':sha256(Path(__file__)),
        'input_verification_sha256':sha256(verification), 'raw_csv_hash_verified':106,
        'records':53,'annotator_streams':106,'normalized_total_event_count':normalized_total,
        'synthetic_controls_passed':ntests,'normalization_changed_streams':len(log),
        'duplicate_instances_removed':sum(e['exact_duplicate_instances_removed'] for e in log),
        'raw_adjacent_order_inversions':sum(e['raw_adjacent_order_inversions'] for e in log),
        'primary_rr_computed':False, 'status':'PASS'}
    (out/'N2_Session2_Controls.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    (out/'N2_Session2_Normalization_Ledger.json').write_text(json.dumps(log,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('SESSION 2: PASS synthetic_controls=',ntests,'records=53 input_hashes=106 normalization_streams=',len(log))
    print('Output:',out.resolve())

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('SESSION 2: FAIL',type(exc).__name__,str(exc),file=sys.stderr)
        sys.exit(1)
