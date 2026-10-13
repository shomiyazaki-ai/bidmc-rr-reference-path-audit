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


def fmt(v):
    return '' if v is None else repr(v)

def annotation_window(samples, start, end):
    selected = [x for x in samples if in_window(x,start,end)]
    n = len(selected)
    if not n:
        return {'n':0,'first_sample':None,'last_sample':None,'first_s':None,'last_s':None,
                'support_s':None,'support_fraction':None,'lead_gap_s':None,'trail_gap_s':None,'rr':None}
    first=(selected[0]-1)/125.0; last=(selected[-1]-1)/125.0
    span=(selected[-1]-selected[0])/125.0 if n>=2 else 0.0
    rr = rate(selected,start,end)
    return {'n':n,'first_sample':selected[0],'last_sample':selected[-1],
        'first_s':first,'last_s':last,'support_s':span,'support_fraction':span/(end-start),
        'lead_gap_s':first-start,'trail_gap_s':end-last,'rr':rr}

def monitor_window(tmap,start,width):
    expected=list(range(start,start+width))
    missing=[t for t in expected if t not in tmap]
    bad=[t for t in expected if t in tmap and not math.isfinite(tmap[t])]
    finite=[tmap[t] for t in expected if t in tmap and math.isfinite(tmap[t])]
    complete=not missing and not bad
    zero_count=sum(v==0 for v in finite)
    zero_state=('zero_all' if zero_count==width else 'zero_none' if zero_count==0 else 'zero_partial') if complete else 'not_complete'
    return {'valid':complete,'value':sum(finite)/width if complete else None,
        'finite_count':len(finite),'missing_slots':missing,'nonfinite_slots':bad,
        'zero_count':zero_count,'zero_state':zero_state}

def load_numerics(path,rid):
    with path.open(newline='', encoding='utf-8-sig') as fh:
        rdr=csv.DictReader(fh)
        if rdr.fieldnames is None: raise ValueError(f'{rid}: no Numerics header')
        rdr.fieldnames=[x.strip() for x in rdr.fieldnames]
        if rdr.fieldnames != ['Time [s]','HR','PULSE','RESP','SpO2']:
            raise ValueError(f'{rid}: unexpected Numerics headers {rdr.fieldnames}')
        rows=list(rdr)
    tm={}
    for row in rows:
        t=float(row['Time [s]'])
        if not math.isfinite(t) or not t.is_integer() or int(t) in tm:
            raise ValueError(f'{rid}: non-integer, duplicate or nonfinite Numerics time')
        v=row['RESP'].strip()
        tm[int(t)]=float(v) if v else float('nan')
    actual=sorted(tm)
    expected=list(range(-1,480)) if rid==23 else list(range(481))
    if actual != expected: raise ValueError(f'{rid}: stored Time grid differs from frozen input audit')
    return tm

def load_annotations(path,rid):
    with path.open(newline='',encoding='utf-8-sig') as fh:
        rdr=csv.DictReader(fh)
        if rdr.fieldnames is None: raise ValueError(f'{rid}: no Breaths header')
        rdr.fieldnames=[x.strip() for x in rdr.fieldnames]
        if rdr.fieldnames != list(HEADERS): raise ValueError(f'{rid}: Breaths headers wrong')
        rows=list(rdr)
    ret=[]; changes=[]
    for j,key in enumerate(HEADERS,1):
        unique,change=normalize([row[key] for row in rows],f'{rid:02d}',j)
        ret.append(unique)
        if change: changes.append(change)
    return ret,changes

def test_s3():
    checks=0
    def ck(x,s):
        nonlocal checks
        checks+=1
        if not x: raise AssertionError(s)
    a=annotation_window([1,126,376,4001],0,32)
    ck(a['n']==3 and math.isclose(a['rr'],40), 'window rr')
    ck(a['support_s']==3 and a['lead_gap_s']==0 and a['trail_gap_s']==29,'window support')
    ck(annotation_window([4001],0,32)['rr'] is None,'right boundary excludes')
    m=monitor_window({t:10 for t in range(32)},0,32)
    ck(m['valid'] and m['value']==10 and m['zero_state']=='zero_none','complete monitor')
    m=monitor_window({t:0 for t in range(32)},0,32)
    ck(m['valid'] and m['value']==0 and m['zero_state']=='zero_all','zeros valid')
    m=monitor_window({t:(0 if t==1 else 10) for t in range(32)},0,32)
    ck(m['valid'] and m['zero_state']=='zero_partial','partial zeros')
    m=monitor_window({t:10 for t in range(31)},0,32)
    ck(not m['valid'] and m['missing_slots']==[31] and m['value'] is None,'missing slot')
    m=monitor_window({t:(float('nan') if t==12 else 10) for t in range(32)},0,32)
    ck(not m['valid'] and m['nonfinite_slots']==[12] and m['value'] is None,'nan invalid')
    ck(reference(12,14)==(True,13),'agreement boundary')
    ck(reference(12,14.000001)==(False,None),'no pre-round')
    return checks

def main3():
    ap=argparse.ArgumentParser(description='N2 frozen 32s/3s per-window primary computation only; no pooled metrics.')
    ap.add_argument('--input-dir',default='N2_inputs_v1.0.0')
    ap.add_argument('--verification',default=None)
    ap.add_argument('--session2-ledger',default='N2_Session2_output/N2_Session2_Normalization_Ledger.json')
    ap.add_argument('--output-dir',default='N2_Session3_output')
    args=ap.parse_args()
    base=Path(args.input_dir)
    ver=Path(args.verification) if args.verification else base/'N2_Session1_Verification.json'
    info=json.loads(ver.read_text(encoding='utf-8'))
    if info.get('status')!='PASS' or info.get('passed_files')!=106 or info.get('record_count')!=53:
        raise ValueError('Session 1 verification not PASS')
    expected={r['file']:r for r in info['files']}
    names={f'bidmc_{rid:02d}_{kind}.csv' for rid in range(1,54) for kind in ('Breaths','Numerics')}
    if set(expected)!=names:raise ValueError('106-file name set incorrect')
    for name in sorted(names):
        p=base/name
        if not p.is_file() or sha256(p)!=expected[name]['sha256'] or p.stat().st_size!=expected[name]['bytes']:
            raise ValueError(f'Input hash/size mismatch: {name}')
    assert controls()==17
    assert test_s3()==10
    previous=json.loads(Path(args.session2_ledger).read_text(encoding='utf-8'))
    outputs=[]; changes=[]; counters=Counter(); count_by_record={}
    for rid in range(1,54):
        tmap=load_numerics(base/f'bidmc_{rid:02d}_Numerics.csv',rid)
        anns,change=load_annotations(base/f'bidmc_{rid:02d}_Breaths.csv',rid)
        changes.extend(change)
        byrec=Counter()
        affected={j:set() for j in (1,2)}
        for e in change:affected[e['annotator']]=set(e['potentially_affected_primary_window_k'])
        for k in range(150):
            start=k*3;end=start+32
            monitor=monitor_window(tmap,start,32)
            a1=annotation_window(anns[0],start,end);a2=annotation_window(anns[1],start,end)
            agree,ref=reference(a1['rr'],a2['rr'])
            both=a1['rr'] is not None and a2['rr'] is not None
            eligible=monitor['valid'] and agree is True
            diff=monitor['value']-ref if eligible else None
            reasons=[]
            if not monitor['valid']:
                if monitor['missing_slots']:reasons.append('monitor_missing_slot')
                if monitor['nonfinite_slots']:reasons.append('monitor_nonfinite_RESP')
            if a1['rr'] is None:reasons.append('ann1_insufficient_events')
            if a2['rr'] is None:reasons.append('ann2_insufficient_events')
            if agree is False:reasons.append('interannotator_disagreement_gt_2')
            row={'record_id':f'{rid:02d}','window_k':k,'start_s':start,'end_s':end,
                 'numerics_sha256':expected[f'bidmc_{rid:02d}_Numerics.csv']['sha256'],
                 'breaths_sha256':expected[f'bidmc_{rid:02d}_Breaths.csv']['sha256'],
                 'monitor_valid':monitor['valid'],'monitor_finite_count':monitor['finite_count'],
                 'monitor_missing_slot_count':len(monitor['missing_slots']),
                 'monitor_missing_slots_s':';'.join(map(str,monitor['missing_slots'])),
                 'monitor_nonfinite_slots_s':';'.join(map(str,monitor['nonfinite_slots'])),
                 'monitor_zero_count_finite':monitor['zero_count'], 'monitor_zero_state':monitor['zero_state'],
                 'monitor_RR_mean':fmt(monitor['value'])}
            for j,a in [(1,a1),(2,a2)]:
                for key in ['n','first_sample','last_sample','first_s','last_s','support_s','support_fraction','lead_gap_s','trail_gap_s','rr']:
                    row[f'ann{j}_{key}']=fmt(a[key])
                row[f'ann{j}_valid']=a['rr'] is not None
                row[f'ann{j}_normalization_affected_candidate']=k in affected[j]
            row.update(both_annotation_computable=both,agreement_evaluable=both,
                 interannotator_abs_diff=fmt(abs(a1['rr']-a2['rr']) if both else None),
                 agreement_pass=fmt(agree), RR_ref=fmt(ref),
                 primary_comparison_eligible=eligible, signed_difference=fmt(diff),
                 absolute_difference=fmt(abs(diff) if diff is not None else None),
                 exclusion_reasons=';'.join(reasons))
            outputs.append(row)
            for key,val in [('base_windows',True),('monitor_complete',monitor['valid']),
              ('ann1_computable',a1['rr'] is not None),('ann2_computable',a2['rr'] is not None),
              ('both_computable',both),('agreement_pass',agree is True),('agreement_fail',agree is False),
              ('agreement_not_evaluable',agree is None),('pre_screen_monitor_intersection',monitor['valid'] and both),
              ('eligible',eligible),('zero_all_eligible',eligible and monitor['zero_state']=='zero_all')]:
                if val:counters[key]+=1;byrec[key]+=1
        assert byrec['base_windows']==150
        count_by_record[f'{rid:02d}']=dict(byrec)
    new={(e['record_id'],e['annotator']):(e['exact_duplicate_instances_removed'],e['raw_adjacent_order_inversions'],e['implicated_samples']) for e in changes}
    old={(e['record_id'],e['annotator']):(e['exact_duplicate_instances_removed'],e['raw_adjacent_order_inversions'],e['implicated_samples']) for e in previous}
    if new!=old or len(new)!=4:raise ValueError('Session 2 ledger disagreement')
    if len(outputs)!=7950:raise AssertionError('wrong primary row count')
    if counters['agreement_pass']+counters['agreement_fail']+counters['agreement_not_evaluable']!=7950:raise AssertionError('agreement partition')
    if counters['eligible']>counters['pre_screen_monitor_intersection']:raise AssertionError('eligibility set')
    if any((r['primary_comparison_eligible'] is False and r['signed_difference']!='') or
           (r['primary_comparison_eligible'] is True and (r['signed_difference']=='' or r['RR_ref']=='')) for r in outputs):
        raise AssertionError('difference eligibility')
    norm_ids=sorted({(e['record_id'],k) for e in changes for k in e['potentially_affected_primary_window_k']})
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    csvpath=out/'N2_Session3_Window_Results.csv'
    with csvpath.open('w',newline='',encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(outputs[0]));writer.writeheader();writer.writerows(outputs)
    report={'stage':'Session 3 Primary window computation; NO pooled MAE/bias summary',
      'protocol':'N2 Minimum Research Protocol v1.0 FROZEN 2026-10-12',
      'python':sys.version,'script_sha256':sha256(Path(__file__)),
      'session1_verification_sha256':sha256(ver),'session2_ledger_sha256':sha256(Path(args.session2_ledger)),
      'input_hash_verified':106,'records':53,'primary_windows':len(outputs),
      'synthetic_controls_session2':17,'synthetic_controls_session3':10,
      'counts':dict(counters),'count_by_record':count_by_record,
      'changed_normalization_streams':len(changes),
      'normalization_candidate_windows_unique_record_k':len(norm_ids),
      'normalization_candidate_windows_sum_streams':sum(len(e['potentially_affected_primary_window_k']) for e in changes),
      'normalization_40_vs_47_note':'47 candidate record/windows under Session 2 implicated-sample rule; pre-protocol audit reported 40. Reasons for this discrepancy are not established; no biological-error or primary-impact conclusion is drawn.',
      'window_csv_sha256':sha256(csvpath),'status':'PASS'}
    (out/'N2_Session3_Controls.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('SESSION 3: PASS windows=',len(outputs),'records=53 input_hashes=106 eligible=',counters['eligible'])
    print('Output:',out.resolve())

if __name__=='__main__':
    try:main3()
    except Exception as e:
        print('SESSION 3: FAIL',type(e).__name__,str(e),file=sys.stderr)
        sys.exit(1)
