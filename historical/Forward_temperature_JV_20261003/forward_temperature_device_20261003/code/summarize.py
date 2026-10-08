from pathlib import Path
import numpy as np,json,csv
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];data={}
for b in [.4,.5]:
 for T in [300,340,380,420]:
  data[b,T]=json.loads((R/'data'/f'b{b:g}_T{T}_N161.json').read_text())['rows']
rows=sum(data.values(),[])
with (R/'forward_JV_all_points.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(11,4.9),sharey=True)
for a,b in zip(ax,[.4,.5]):
 for T,col in zip([300,340,380,420],['#235789','#3b9b70','#d99823','#b63d38']):
  rr=data[b,T];a.semilogy([x['V'] for x in rr[1:]],[x['J_Acm2']*1000 for x in rr[1:]],'o-',markersize=3,color=col,label=f'{T} K')
 a.set(xlabel='Forward voltage (V)',title=f'Contact barrier b = {b:.1f} eV',xlim=(0,1.22),ylim=(1e-4,3e6));a.grid(alpha=.2);a.legend()
ax[0].set_ylabel('Dark current density (mA cm$^{-2}$)')
fig.suptitle('Conditional forward J–V: same parameters at every temperature',fontsize=13)
fig.text(.5,.02,'Isothermal model; no series resistance or self-heating. High-current branch is a diagnostic, not an operating prediction.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.055,1,.95));fig.savefig(R/'forward_temperature_JV.png',dpi=220);fig.savefig(R/'forward_temperature_JV.pdf');plt.close(fig)
summary={'points':len(rows),'temperatures_K':[300,340,380,420],'barriers_eV':[.4,.5],'nodes_curves':161,'maxima':{},'representative_currents':[],'grid_checks':[]}
for k in ['relative_error','charge_relative_error','energy_relative_error','max_LDB_cycle_error','max_density_fraction','max_reduced_field','fraction_faces_low_density']:
 summary['maxima'][k]=max(x[k] for x in rows)
for b in [.4,.5]:
 for V in [.1,.3,.5,1.,1.2]:
  vals=[next(x['J_Acm2'] for x in data[b,T] if x['V']==V) for T in [300,340,380,420]]
  ea=-np.polyfit(1/np.array([300,340,380,420]),np.log(vals),1)[0]*8.617333262145e-5
  summary['representative_currents'].append(dict(barrier=b,V=V,J_Acm2=vals,ratio420_300=vals[-1]/vals[0],apparent_Ea_eV=ea))
for f in sorted((R/'data').glob('refine*.json')):
 p=json.loads(f.read_text());b=p['barrier_eV'];T=p['T'];v=p['V'];base=next(x['J_Acm2'] for x in data[b,T] if x['V']==v)
 summary['grid_checks'].append(dict(b=b,T=T,V=v,N=p['nodes'],J_Acm2=p['J_Acm2'],relative_to_161=p['J_Acm2']/base-1))
(R/'SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
