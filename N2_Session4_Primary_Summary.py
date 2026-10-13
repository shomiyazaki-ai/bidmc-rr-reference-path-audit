#!/usr/bin/env python3
"""N2 Session 4: deterministic descriptive summary of frozen Session 3 window CSV.
No modification to primary methods, no sensitivity or inferential tests.
"""
import csv, json, hashlib, math, statistics, collections, pathlib, sys
ROOT=pathlib.Path(__file__).resolve().parent
SRC=ROOT/'N2_Session3_Window_Results.csv'
CTL=ROOT/'N2_Session3_Controls.json'
OUT=ROOT/'N2_Session4_output'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def mean(xs): return sum(xs)/len(xs) if xs else None
def fmt(x): return 'NA' if x is None else f'{x:.9f}'
def write_csv(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader();w.writerows(rows)
def flag(row,k):
    s=row[k]; assert s in ('True','False') or (k=='agreement_pass' and s=='' and row['agreement_evaluable']=='False'),(k,s);return s=='True'
def num(row,k):return float(row[k]) if row[k] else None

def summarize(rows):
    n=len(rows); eligible=[r for r in rows if flag(r,'primary_comparison_eligible')]
    for r in rows:
        if flag(r,'primary_comparison_eligible'):
            assert flag(r,'monitor_valid') and flag(r,'both_annotation_computable') and flag(r,'agreement_pass')
            assert all(r[k] for k in ['monitor_RR_mean','RR_ref','signed_difference','absolute_difference'])
            assert math.isclose(num(r,'signed_difference'),num(r,'monitor_RR_mean')-num(r,'RR_ref'),abs_tol=1e-9)
            assert math.isclose(num(r,'absolute_difference'),abs(num(r,'signed_difference')),abs_tol=1e-9)
        else:
            assert not r['signed_difference'] and not r['absolute_difference']
    groups={k: [r for r in rows if r['monitor_zero_state']==k and flag(r,'primary_comparison_eligible')] for k in ('zero_none','zero_partial','zero_all')}
    d=[num(r,'signed_difference') for r in eligible]; a=[abs(x) for x in d]
    counts={
        'base_windows':n,
        'monitor_complete':sum(flag(r,'monitor_valid') for r in rows),
        'ann1_computable':sum(flag(r,'ann1_valid') for r in rows),
        'ann2_computable':sum(flag(r,'ann2_valid') for r in rows),
        'both_computable':sum(flag(r,'both_annotation_computable') for r in rows),
        'agreement_pass':sum(flag(r,'agreement_pass') for r in rows),
        'agreement_fail':sum(flag(r,'agreement_evaluable') and not flag(r,'agreement_pass') for r in rows),
        'agreement_not_evaluable':sum(not flag(r,'agreement_evaluable') for r in rows),
        'pre_screen_monitor_intersection':sum(flag(r,'monitor_valid') and flag(r,'both_annotation_computable') for r in rows),
        'eligible':len(eligible),
        'monitor_incomplete':sum(not flag(r,'monitor_valid') for r in rows),
        'both_not_computable':sum(not flag(r,'both_annotation_computable') for r in rows),
        'pre_intersection_screen_fail':sum(flag(r,'monitor_valid') and flag(r,'both_annotation_computable') and not flag(r,'agreement_pass') for r in rows),
        'zero_none_eligible':len(groups['zero_none']),
        'zero_partial_eligible':len(groups['zero_partial']),
        'zero_all_eligible':len(groups['zero_all'])
    }
    assert counts['agreement_pass']+counts['agreement_fail']+counts['agreement_not_evaluable']==n
    assert counts['pre_screen_monitor_intersection']-counts['pre_intersection_screen_fail']==counts['eligible']
    assert counts['eligible']==sum(len(x) for x in groups.values())
    reasons=collections.Counter()
    for r in rows:
        if r['exclusion_reasons']:
            for reason in r['exclusion_reasons'].split(';'):
                reasons[reason.strip()]+=1
    return {'counts':counts,'mae':mean(a),'bias':mean(d),'zero':{k:{'eligible':len(v),'mae':mean([num(r,'absolute_difference') for r in v])} for k,v in groups.items()},'reasons':dict(reasons)}

def main():
    with SRC.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    control=json.loads(CTL.read_text(encoding='utf-8'))
    assert sha(SRC)==control['window_csv_sha256'],'Session 3 source hash mismatch'
    assert len(rows)==7950
    by=collections.defaultdict(list)
    for r in rows: by[r['record_id']].append(r)
    assert sorted(by)==[f'{i:02d}' for i in range(1,54)]
    for rid,group in by.items():
        assert len(group)==150 and [int(x['window_k']) for x in group]==list(range(150)),rid
    summary=summarize(rows)
    for k,v in control['counts'].items():
        if k in summary['counts']: assert v==summary['counts'][k],(k,v,summary['counts'][k])
    recs=[]
    for rid,group in sorted(by.items()):
        x=summarize(group); c=x['counts'];recs.append({'record_id':rid,**c,'MAE_bpm':fmt(x['mae']),'signed_bias_bpm':fmt(x['bias']),**{f'{s}_MAE_bpm':fmt(x['zero'][s]['mae']) for s in x['zero']}})
    available=[x['mae'] for x in (summarize(g) for g in by.values()) if x['mae'] is not None]
    assert len(available)==sum(1 for r in recs if int(r['eligible'])>0)
    assert sum(int(r['eligible']) for r in recs)==summary['counts']['eligible']
    eq=mean(available)
    overview={'protocol':'N2 Minimum Research Protocol v1.0 FROZEN 2026-10-12','stage':'Session 4 descriptive primary summary (no sensitivity, MSC or independent reproducibility claimed)', 'source_window_csv_sha256':sha(SRC),'source_control_sha256':sha(CTL),'python_version':sys.version,'counts':summary['counts'],'window_pooled_MAE_bpm':summary['mae'],'window_pooled_signed_bias_bpm':summary['bias'],'equal_record_MAE_bpm':eq,'equal_record_R':len(available),'records_total':53,'zero_state':summary['zero'],'exclusion_reason_occurrences_nonexclusive':summary['reasons'],'limitations':['Overlapping windows are not independent observations','The monitor and manual annotation pathways refer to the same underlying impedance respiration signal','MAE here is numerical discordance, not prediction error or clinical accuracy','Two bpm screening applies to annotators, not monitor-reference differences','Finite monitor zero retained; no physiological apnoea inference','Reasons overlap and do not partition ineligible windows']}
    OUT.mkdir(exist_ok=True)
    write_csv(OUT/'N2_Record_Summary.csv',recs)
    (OUT/'N2_Primary_Summary.json').write_text(json.dumps(overview,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    text=['N2 PRIMARY SUMMARY — Session 4','Protocol: v1.0 FROZEN (2026-10-12)',f'Input Window CSV SHA256: {sha(SRC)}','','PRIMARY — all eligible overlapping windows',f'Eligible: {summary["counts"]["eligible"]} / {len(rows)}',f'Window-pooled MAE (bpm): {fmt(summary["mae"])}',f'Window-pooled signed bias, monitor - RR_ref (bpm): {fmt(summary["bias"])}','','SUPPLEMENTARY',f'Equal-record MAE (bpm): {fmt(eq)} (R={len(available)} / 53 records with >=1 eligible window)']
    for k,v in summary['zero'].items():text.append(f'{k}: N={v["eligible"]}, MAE_bpm={fmt(v["mae"])}')
    text+=['','AVAILABILITY (window counts; overlapping conditions do not partition):']+[f'{k}: {v}' for k,v in summary['counts'].items()]+['','EXCLUSION REASONS (nonexclusive event counts):']+[f'{k}: {v}' for k,v in sorted(summary['reasons'].items())]+['','Limitations:']+['- '+s for s in overview['limitations']]
    (OUT/'N2_Primary_Summary.txt').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(f'SESSION 4: PASS records={len(recs)} eligible={summary["counts"]["eligible"]} pooled_MAE={summary["mae"]:.9f} bias={summary["bias"]:.9f} equal_record_MAE={eq:.9f}')
    print('Output:',OUT)
if __name__=='__main__': main()
