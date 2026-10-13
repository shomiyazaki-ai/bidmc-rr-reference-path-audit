#!/usr/bin/env python3
"""N2 FROZEN v1.0: ONE 64 s / 6 s sensitivity grid, no parameter sweep.
Run beside N2_Session3_Primary_Windows.py. Python 3.9+, stdlib only.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
import N2_Session3_Primary_Windows as core

WIDTH, STEP, PER_RECORD = 64, 6, 70

def check(condition, message):
    if not condition: raise AssertionError(message)

def test_sensitivity():
    n = 0
    def t(ok, label):
        nonlocal n
        n += 1
        check(ok, label)
    t(len(range(0, 480-WIDTH+1, STEP)) == PER_RECORD, '70 windows')
    t((PER_RECORD-1)*STEP == 414 and 414+WIDTH == 478, 'last [414,478)')
    t(core.in_window(1, 0,64), 'first sample inclusive')
    t(core.in_window(8000, 0,64) and not core.in_window(8001,0,64), 'right exclusive')
    t(core.rate([1,126,251],0,64)==60, '60 breaths/min')
    t(core.rate([1,8001],0,64) is None, 'no external interval')
    m=core.monitor_window({i:0 for i in range(64)},0,64)
    t(m['valid'] and m['zero_state']=='zero_all' and m['value']==0, 'all zero valid')
    m=core.monitor_window({i:(0 if i==3 else 12) for i in range(64)},0,64)
    t(m['valid'] and m['zero_state']=='zero_partial', 'partial zero')
    m=core.monitor_window({i:12 for i in range(63)},0,64)
    t(not m['valid'] and m['value'] is None and m['missing_slots']==[63], '63/64 invalid')
    m=core.monitor_window({i:(float('nan') if i==3 else 12) for i in range(64)},0,64)
    t(not m['valid'] and m['value'] is None and m['nonfinite_slots']==[3], 'NaN invalid')
    t(core.reference(12,14)==(True,13) and core.reference(12,14.000001)==(False,None), 'selection <=2 before rounding')
    return n

def mean_or_none(values):
    return math.fsum(values)/len(values) if values else None

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir',default='N2_inputs_v1.0.0')
    ap.add_argument('--session2-ledger',default='N2_Session2_output/N2_Session2_Normalization_Ledger.json')
    ap.add_argument('--session3-controls',default='N2_Session3_output/N2_Session3_Controls.json')
    ap.add_argument('--session4-summary',default='N2_Session4_output/N2_Primary_Summary.json')
    ap.add_argument('--output-dir',default='N2_Session5_output')
    args=ap.parse_args()
    base=Path(args.input_dir)
    ver=base/'N2_Session1_Verification.json'
    verification=json.loads(ver.read_text(encoding='utf-8'))
    check(verification['status']=='PASS' and verification['passed_files']==106 and verification['record_count']==53,'S1 verification')
    expected={d['file']:d for d in verification['files']}
    names={f'bidmc_{rid:02d}_{kind}.csv' for rid in range(1,54) for kind in ('Breaths','Numerics')}
    check(set(expected)==names,'exact 106 file list')
    for name in sorted(names):
        f=base/name
        check(f.is_file() and f.stat().st_size==expected[name]['bytes'] and core.sha256(f)==expected[name]['sha256'],f'input integrity {name}')
    ctl3=json.loads(Path(args.session3_controls).read_text(encoding='utf-8'))
    check(ctl3['status']=='PASS' and ctl3['window_csv_sha256'] and ctl3['input_hash_verified']==106, 'S3 control')
    summary4_path=Path(args.session4_summary)
    if summary4_path.is_file():
        primary_summary=json.loads(summary4_path.read_text(encoding='utf-8'))
        check(primary_summary.get('source_window_csv_sha256')==ctl3['window_csv_sha256'],'S4/S3 lineage')
    check(core.controls()==17 and core.test_s3()==10 and test_sensitivity()==11, 'synthetic tests')
    previous=json.loads(Path(args.session2_ledger).read_text(encoding='utf-8'))
    expected_changes={(e['record_id'],e['annotator']):(e['exact_duplicate_instances_removed'],e['raw_adjacent_order_inversions'],e['implicated_samples']) for e in previous}
    rows=[]; by_record={}; counters=Counter(); changed=[]
    for rid in range(1,54):
        ridstr=f'{rid:02d}'
        tmap=core.load_numerics(base/f'bidmc_{ridstr}_Numerics.csv',rid)
        anns,changes=core.load_annotations(base/f'bidmc_{ridstr}_Breaths.csv',rid)
        changed.extend(changes)
        record_count=Counter()
        # Membership in 64s sensitivity windows based on implicated unchanged raw samples.
        impacted={e['annotator']:{k for k in range(PER_RECORD) if any(core.in_window(s,k*STEP,k*STEP+WIDTH) for s in e['implicated_samples'])} for e in changes}
        for k in range(PER_RECORD):
            start=k*STEP; end=start+WIDTH
            m=core.monitor_window(tmap,start,WIDTH)
            aa=[core.annotation_window(a,start,end) for a in anns]
            a1,a2=aa
            agree,ref=core.reference(a1['rr'],a2['rr'])
            both=a1['rr'] is not None and a2['rr'] is not None
            eligible=m['valid'] and agree is True
            delta=m['value']-ref if eligible else None
            reasons=[]
            if m['missing_slots']: reasons.append('monitor_missing_slot')
            if m['nonfinite_slots']: reasons.append('monitor_nonfinite_RESP')
            if a1['rr'] is None:reasons.append('ann1_insufficient_events')
            if a2['rr'] is None:reasons.append('ann2_insufficient_events')
            if agree is False:reasons.append('interannotator_disagreement_gt_2')
            row={'record_id':ridstr,'window_k':k,'start_s':start,'end_s':end,
                'numerics_sha256':expected[f'bidmc_{ridstr}_Numerics.csv']['sha256'],
                'breaths_sha256':expected[f'bidmc_{ridstr}_Breaths.csv']['sha256'],
                'monitor_valid':m['valid'],'monitor_finite_count':m['finite_count'],
                'monitor_missing_slot_count':len(m['missing_slots']),
                'monitor_missing_slots_s':';'.join(map(str,m['missing_slots'])),
                'monitor_nonfinite_slots_s':';'.join(map(str,m['nonfinite_slots'])),
                'monitor_zero_count_finite':m['zero_count'],'monitor_zero_state':m['zero_state'],
                'monitor_RR_mean':core.fmt(m['value'])}
            for j,a in enumerate(aa,1):
                for key in ('n','first_sample','last_sample','first_s','last_s','support_s','support_fraction','lead_gap_s','trail_gap_s','rr'):
                    row[f'ann{j}_{key}']=core.fmt(a[key])
                row[f'ann{j}_valid']=a['rr'] is not None
                row[f'ann{j}_normalization_affected_candidate']=k in impacted.get(j,set())
            row.update(both_annotation_computable=both, agreement_evaluable=both,
                interannotator_abs_diff=core.fmt(abs(a1['rr']-a2['rr']) if both else None),
                agreement_pass=core.fmt(agree),RR_ref=core.fmt(ref),
                primary_comparison_eligible=eligible,signed_difference=core.fmt(delta),
                absolute_difference=core.fmt(abs(delta) if delta is not None else None),
                exclusion_reasons=';'.join(reasons))
            rows.append(row)
            conditions={'base_windows':True,'monitor_complete':m['valid'],'ann1_computable':a1['rr'] is not None,
                'ann2_computable':a2['rr'] is not None,'both_computable':both,
                'agreement_pass':agree is True,'agreement_fail':agree is False,
                'agreement_not_evaluable':agree is None,
                'pre_screen_monitor_intersection':m['valid'] and both,
                'eligible':eligible,'zero_all_eligible':eligible and m['zero_state']=='zero_all'}
            for name,condition in conditions.items():
                if condition: counters[name]+=1;record_count[name]+=1
        check(record_count['base_windows']==PER_RECORD,f'record {rid} window count')
        by_record[ridstr]=dict(record_count)
    new_changes={(e['record_id'],e['annotator']):(e['exact_duplicate_instances_removed'],e['raw_adjacent_order_inversions'],e['implicated_samples']) for e in changed}
    check(new_changes==expected_changes and len(changed)==4,'S2 normalization conservation')
    check(len(rows)==3710 and counters['base_windows']==3710,'3710 windows')
    check(counters['agreement_pass']+counters['agreement_fail']+counters['agreement_not_evaluable']==3710,'agreement partition')
    check(counters['both_computable']==counters['agreement_pass']+counters['agreement_fail'],'both valid partition')
    check(counters['eligible']<=counters['pre_screen_monitor_intersection'],'eligible subset')
    eligible=[r for r in rows if r['primary_comparison_eligible']]
    check(len(eligible)==counters['eligible'] and all(r['signed_difference']!='' and r['RR_ref']!='' for r in eligible),'delta existence')
    check(all(r['signed_difference']=='' and r['absolute_difference']=='' for r in rows if not r['primary_comparison_eligible']),'no ineligible differences')
    errors=[float(r['absolute_difference']) for r in eligible]; signed=[float(r['signed_difference']) for r in eligible]
    zero={state:[float(r['absolute_difference']) for r in eligible if r['monitor_zero_state']==state] for state in ('zero_none','zero_partial','zero_all')}
    record_rows=[]
    for rid in range(1,54):
        rr=[r for r in eligible if r['record_id']==f'{rid:02d}']
        record_rows.append({'record_id':f'{rid:02d}','base_windows':PER_RECORD,
            **{name:by_record[f'{rid:02d}'].get(name,0) for name in ('monitor_complete','ann1_computable','ann2_computable','both_computable','agreement_pass','agreement_fail','agreement_not_evaluable','pre_screen_monitor_intersection','eligible','zero_all_eligible')},
            'MAE':core.fmt(mean_or_none([float(r['absolute_difference']) for r in rr])),
            'signed_bias':core.fmt(mean_or_none([float(r['signed_difference']) for r in rr]))})
    valid_record_maes=[float(r['MAE']) for r in record_rows if r['MAE']!='']
    summary={'protocol':'N2 Minimum Research Protocol v1.0 FROZEN 2026-10-12',
        'stage':'Session 5 ONLY, prespecified 64s / 6s sensitivity; NOT a second primary',
        'grid':{'width_s':WIDTH,'step_s':STEP,'per_record':PER_RECORD,'last_window':[414,478]},
        'counts':dict(counters),'records_with_comparisons':len(valid_record_maes),
        'MAE_window_pooled':mean_or_none(errors),'bias_window_pooled':mean_or_none(signed),
        'MAE_record_equal_weight':mean_or_none(valid_record_maes),
        'zero_states':{name:{'n':len(vals),'MAE':mean_or_none(vals)} for name,vals in zero.items()},
        'synthetic_controls':{'session2':17,'session3':10,'sensitivity':11},
        'input_hash_verified':106,'changed_normalization_streams':4,
        'python_version':sys.version,'script_sha256':core.sha256(Path(__file__)),
        'session1_verification_sha256':core.sha256(ver),'session2_ledger_sha256':core.sha256(Path(args.session2_ledger)),
        'session3_controls_sha256':core.sha256(Path(args.session3_controls)),
        'session4_summary_sha256':core.sha256(summary4_path) if summary4_path.is_file() else None,
        'status':'PASS'}
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    windowpath=out/'N2_Session5_Window_Results.csv'
    recordpath=out/'N2_Session5_Record_Summary.csv'
    for path,records in ((windowpath,rows),(recordpath,record_rows)):
        with path.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    summary['window_csv_sha256']=core.sha256(windowpath)
    summary['record_csv_sha256']=core.sha256(recordpath)
    (out/'N2_Session5_Sensitivity_Summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print(f"SESSION 5: PASS windows= {len(rows)} records=53 input_hashes=106 eligible= {counters['eligible']}")
    print('MAE=',summary['MAE_window_pooled'],'bias=',summary['bias_window_pooled'])
    print('Output:',out.resolve())

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print('SESSION 5: FAIL',type(exc).__name__,str(exc),file=sys.stderr)
        sys.exit(1)
