#!/usr/bin/env python3
"""N2 Session 6 deterministic replay: existing FROZEN pipelines, no new scientific analysis.
Place alongside Session3 and Session5 scripts, run from project root with original folders intact.
Creates N2_Session6_output with verification JSON and replay outputs only.
"""
import csv,hashlib,json,math,subprocess,sys
from pathlib import Path
from collections import Counter,defaultdict

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(131072),b''):h.update(b)
 return h.hexdigest()
def req(cond,msg):
 if not cond:raise AssertionError(msg)
def verify_csv(p, width, step, expected, expected_eligible):
 with p.open(newline='',encoding='utf-8') as f: rows=list(csv.DictReader(f))
 req(len(rows)==expected,f'{p} rowcount')
 ids=Counter(r['record_id'] for r in rows)
 req(len(ids)==53 and all(ids[f'{i:02d}']==expected//53 for i in range(1,54)),'record counts')
 elig=[]
 for r in rows:
  k=int(r['window_k']); start=int(r['start_s']); end=int(r['end_s'])
  req(start==k*step and end==start+width,'window grid')
  valid=r['monitor_valid']=='True'; agree=r['agreement_pass']=='True'; selected=r['primary_comparison_eligible']=='True'
  req(selected==(valid and agree),'eligibility')
  if selected:
   delta=float(r['monitor_RR_mean'])-float(r['RR_ref'])
   req(math.isclose(delta,float(r['signed_difference']),rel_tol=1e-12,abs_tol=1e-10),'delta')
   req(math.isclose(abs(delta),float(r['absolute_difference']),rel_tol=1e-12,abs_tol=1e-10),'abs delta')
   elig.append((r,delta))
 req(len(elig)==expected_eligible,'eligible count')
 return {'eligible':len(elig),'mae':math.fsum(abs(d) for _,d in elig)/len(elig),'bias':math.fsum(d for _,d in elig)/len(elig)}
def main():
 root=Path.cwd()
 s3=root/'N2_Session3_output';s5=root/'N2_Session5_output'
 ctr3=json.loads((s3/'N2_Session3_Controls.json').read_text(encoding='utf-8'))
 ctr5=json.loads((s5/'N2_Session5_Sensitivity_Summary.json').read_text(encoding='utf-8'))
 summary4=json.loads((root/'N2_Session4_output/N2_Primary_Summary.json').read_text(encoding='utf-8'))
 orig3=s3/'N2_Session3_Window_Results.csv';orig5=s5/'N2_Session5_Window_Results.csv';orig5r=s5/'N2_Session5_Record_Summary.csv'
 for p in [orig3,orig5,orig5r]:req(p.exists(),f'missing {p}')
 for p,h in [(orig3,ctr3['window_csv_sha256']),(orig5,ctr5['window_csv_sha256']),(orig5r,ctr5['record_csv_sha256'])]:req(sha(p)==h,f'original hash differs {p}')
 a=verify_csv(orig3,32,3,7950,7271);b=verify_csv(orig5,64,6,3710,3503)
 for res,exp,key in [(a,summary4,'window_pooled_MAE_bpm'),(b,ctr5,'MAE_window_pooled')]:req(math.isclose(res['mae'],exp[key],abs_tol=1e-10),'MAE summary')
 req(math.isclose(a['bias'],summary4['window_pooled_signed_bias_bpm'],abs_tol=1e-10),'primary bias')
 req(math.isclose(b['bias'],ctr5['bias_window_pooled'],abs_tol=1e-10),'sensitivity bias')
 out=root/'N2_Session6_output';out.mkdir(exist_ok=True)
 replay3=out/'replay_primary';replay5=out/'replay_sensitivity'
 for p in (replay3,replay5):
  req(not p.exists(),f'{p} already exists; rename/remove before rerun (do not overwrite prior replay)')
 script3=root/'N2_Session3_Primary_Windows.py'; script5=root/'N2_Session5_Sensitivity.py'
 commands=[([sys.executable,str(script3),'--output-dir',str(replay3)],'primary'),([sys.executable,str(script5),'--session3-controls',str(replay3/'N2_Session3_Controls.json'),'--output-dir',str(replay5)],'sensitivity')]
 logs={}
 for cmd,name in commands:
  p=subprocess.run(cmd,cwd=root,capture_output=True,text=True)
  logs[name]={'returncode':p.returncode,'stdout':p.stdout.strip(),'stderr':p.stderr.strip()}
  req(p.returncode==0,f'{name} replay failed: {p.stderr}')
 paths=[('primary_window',orig3,replay3/'N2_Session3_Window_Results.csv'),('sensitivity_window',orig5,replay5/'N2_Session5_Window_Results.csv'),('sensitivity_record',orig5r,replay5/'N2_Session5_Record_Summary.csv')]
 comparisons={}
 for label,original,replay in paths:
  comparisons[label]={'original_sha256':sha(original),'replay_sha256':sha(replay),'identical':sha(original)==sha(replay)}
  req(comparisons[label]['identical'],f'{label} deterministic mismatch')
 report={'protocol':'N2 Minimum Research Protocol v1.0 FROZEN','stage':'Session 6 deterministic replay','status':'PASS','python':sys.version,'script3_sha256':sha(script3),'script5_sha256':sha(script5),'input_hashes_verified_by_replay':106,'original_summary_recalculated':{'primary':a,'sensitivity':b},'replay_comparisons':comparisons,'logs':logs,'note':'Replay shares same implementation as original; not independent raw-to-RR reimplementation.'}
 (out/'N2_Session6_Verification.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
 print('SESSION 6: PASS replay primary=7950 sensitivity=3710 exact_csv_hashes=3')
 print('Output:',out.resolve())
if __name__=='__main__':
 try:main()
 except Exception as e:
  print('SESSION 6: FAIL',type(e).__name__,e,file=sys.stderr);sys.exit(1)
