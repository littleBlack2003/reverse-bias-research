from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def load(name):return json.loads((ROOT/'data'/name).read_text())
def runfile(t):return f'T{t}_N321_I1_blocking_bn0_bp0'+('_lightonly' if t==100 else'')
def main():
 temps=[100,150,200,250,300];metrics=[]
 fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
 for T,col in zip(temps,['C0','C1','C2','C3','C4']):
  tag=runfile(T);a=load(tag+'.json');rr=sorted(a['curves']['light'],key=lambda r:r['V']);V=[r['V'] for r in rr];J=[r['J_Acm2']*1000 for r in rr]
  axs[0,0].plot(V,J,label=f'{T} K',color=col)
  axs[0,1].semilogy(V,np.maximum(abs(np.array(J)),1e-25),label=f'{T} K light',color=col)
  darkfile=ROOT/'data'/('T100_N321_dark_partial_snapshot.json' if T==100 else f'T{T}_N321_I1_blocking_bn0_bp0.json')
  if darkfile.exists():
   aa=json.loads(darkfile.read_text());dr=sorted(aa['curves'].get('dark',[]),key=lambda r:r['V']);dr=[r for r in dr if r['V']>0 and r['gate_passed']]
   axs[0,1].semilogy([r['V'] for r in dr],[abs(r['J_Acm2']*1000) for r in dr],'--',color=col,alpha=.75)
  metrics.append(load(tag+'_metrics.json'))
 axs[0,0].set(xlabel='Applied voltage (V)',ylabel='Signed current (mA/cm²)',ylim=(-25,15),xlim=(0,1.15),title='Conditional 1-sun-equivalent JV');axs[0,0].axhline(0,c='k',lw=.5);axs[0,0].legend(fontsize=8)
 axs[0,1].set(xlabel='Applied voltage (V)',ylabel='Current magnitude (mA/cm²)',ylim=(1e-12,1e3),xlim=(0,1.15),title='Light solid; dark dashed; 100 K extrapolated');axs[0,1].legend(fontsize=8)
 axs[1,0].semilogy(temps,[abs(m['Jsc_Acm2'])*1000 for m in metrics],'o-',label='Conditional model Jsc')
 axs[1,0].scatter([100,300],[17.9,29.2],marker='s',label='S15 |J(0.04 V)| endpoints')
 axs[1,0].set(xlabel='Temperature (K)',ylabel='Photocurrent (mA/cm²)',title='Low-T transport extrapolation is incompatible');axs[1,0].legend(fontsize=8)
 axs[1,1].plot(temps,[m['Voc_V'] for m in metrics],'o-',label='Conditional full-device Voc')
 audit=load('calibration_audit.json');axs[1,1].plot([r['T_K'] for r in audit['rows']],[r['Voc_fixed142_V'] for r in audit['rows']],'--',label='Bulk QFLS with same Eg1.42 eV')
 axs[1,1].scatter([r['T_K'] for r in audit['rows']],[r['Voc_readout_V'] for r in audit['rows']],c='k',label='Fig5a measured read-offs')
 axs[1,1].set(xlabel='Temperature (K)',ylabel='Voltage (V)',title='Bulk QFLS is not automatically terminal Voc');axs[1,1].legend(fontsize=8)
 fig.savefig(ROOT/'figures/conditional_jv.png',dpi=200);plt.close(fig)
 simple=[{k:v for k,v in m.items() if k!='refined_rows'} for m in metrics]
 (ROOT/'data/metrics_summary.json').write_text(json.dumps(simple,indent=2))
 with (ROOT/'data/metrics_summary.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['T_K','Jsc_mAcm2','Voc_V','FF_percent','Vmp_V','Jmp_mAcm2','interpretation']);w.writeheader()
  for m in metrics:w.writerow(dict(T_K=m['T'],Jsc_mAcm2=abs(m['Jsc_Acm2'])*1000,Voc_V=m['Voc_V'],FF_percent=100*m['FF'],Vmp_V=m['Vmp_V'],Jmp_mAcm2=m['Jmp_Acm2']*1000,interpretation='unvalidated contact/constant-G scenario; mobility extrapolated below223K'))
 print('summary done')
if __name__=='__main__':main()
