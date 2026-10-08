"""Typography-only replay of frozen plotting code and numeric sources. No device fitting or scan rerun."""
from pathlib import Path
from types import SimpleNamespace
import numpy as np,csv,json,hashlib,textwrap,os
P=Path(__file__).resolve().parents[1];os.environ['MPLCONFIGDIR']=str(P/'qa/mpl')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.font_manager import fontManager
for f in (P/'assets').glob('*.ttf'):
 if 'variable' not in f.name:fontManager.addfont(str(f))
plt.rcParams.update({'font.family':['Liberation Serif','Noto Serif SC','Noto Serif'],'mathtext.fontset':'stix','font.size':10,'figure.dpi':160,'svg.fonttype':'path'})
D=P/'data/research';I=P/'typography_inputs';F=P/'figures';seen=[];outputs=[]
map_names={'branch_and_mesh':'branch_mesh','spatial_balance':'spatial_balance','newton_failure':'newton_failure','recovered_branch_and_gate':'qf_recovery','cancellation_confidence':'lowT_cancellation','cycle_vs_field':'graph_field','failed_naive_reverse':'graph_wrong','coulomb_extraction_feedback':'graph_feedback','temperature_and_pair_accumulation':'reservoir_temperature','same_budget_spectral_support_control':'reservoir_support','quantum_and_balance':'quantum_kernel','dos_and_extraction':'dos_extract'}
orig_save=Figure.savefig

def save(self,path,*args,**kw):
 key=Path(path).stem
 if key not in map_names:return
 name=map_names[key];dst=F/(name+'.png')
 if name=='graph_wrong':
  from matplotlib.ticker import FixedLocator
  ax=self.axes[0];ax.yaxis.set_major_locator(FixedLocator([v for v in ax.get_yticks() if v==0 or abs(v)>=1e-18]))
 # Preserve curve values, labels, axes scales and limits independently of font/render changes.
 axes=[]
 for ax in self.axes:
  curves=[]
  for l in ax.lines:
   xy=np.asarray(l.get_xydata(),dtype=float);curves.append({'label':l.get_label(),'xy_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'n':len(xy)})
  axes.append({'xlabel':ax.get_xlabel(),'ylabel':ax.get_ylabel(),'title':ax.get_title(),'xscale':ax.get_xscale(),'yscale':ax.get_yscale(),'xlim':list(ax.get_xlim()),'ylim':list(ax.get_ylim()),'curves':curves})
 kw['dpi']=250
 orig_save(self,dst,*args,**kw);orig_save(self,dst.with_suffix('.svg'),*args,**{k:v for k,v in kw.items() if k!='dpi'})
 outputs.append({'name':name,'axes':axes,'png_sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'svg_sha256':hashlib.sha256(dst.with_suffix('.svg').read_bytes()).hexdigest()})
Figure.savefig=save

def note(f):seen.append({'path':str(f.relative_to(P)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()});return f

def read(name,base):
 out=list(csv.DictReader(note(base/name).open()))
 for r in out:
  for k,v in r.items():
   if k in ['state','source','destination','process']:continue
   try:r[k]=float(v)
   except (ValueError,TypeError):pass
 return out

def run(name,start=None,end=None,extra=None):
 f=note(I/name);s=f.read_text()
 if start:s=s[s.index(start):]
 if end:s=s[:s.index(end)]
 ns=globals().copy()
 if extra:ns.update(extra)
 exec(compile(textwrap.dedent(s),str(f),'exec'),ns)

def savefig(name):plt.tight_layout();plt.savefig(F/name,dpi=250);plt.close()
# Direct plotting modules, replacing only their import of solver machinery and output destinations.
root=D/'revalidated_rate_audit_20261002'
s=(I/'rate_plot.py').read_text().replace('from rate_audit import *','')
s=s.replace("ROOT/'data/baseline_161_-20.npz'", "I/'baseline_161_-20.npz'")
ns=globals().copy();ns.update(ROOT=root,device=lambda n:SimpleNamespace(x=np.linspace(0,1,n)));exec(compile(s,'rate_plot.py','exec'),ns)
run('solver_plot.py',start='import matplotlib',extra={'ROOT':D/'asymmetric_solver_diagnostics_20261002'})
run('confidence_plot.py',start='import matplotlib',end='fig,axes=plt.subplots(1,2,figsize=(10.5',extra={'ROOT':D/'qf_temperature_confidence_20261002'})
# The selected graph figures use original plotting blocks with already frozen rows.
b=D/'quantum_cycle_audit_20261002/results';rows=read('full_reversible_cycle_scan.csv',b);bad=read('counterexamples.csv',b)
run('cycle_plot_source.py',start=' fig,ax=plt.subplots(1,2,figsize=(11,4.4),sharey=True)',end=' fig,ax=plt.subplots(1,2,figsize=(11,4.3))')
run('cycle_plot_source.py',start=' fig,ax=plt.subplots(figsize=(7,4.3))',end=' # Hard checks')
rows=read('extraction_coulomb_feedback_scan.csv',b)
run('feedback_plot_source.py',start='fig,axes=plt.subplots(1,2,figsize=(10.8,4.3))',end='print(json.dumps(summary',extra={'FIG':F})
b=D/'quantum_reservoir_audit_20261002/results';rows=read('matched_boundary_full_scan.csv',b);states=read('representative_state_probabilities.csv',b);supportrows=read('one_sided_vs_broad_contact_support.csv',b);reps=read('representative_main_points.csv',b);Ts=sorted({r['T_K'] for r in rows})
run('reservoir_plot_source.py',start=' fig,ax=plt.subplots(1,2,figsize=(11,4.4))',end=' fig,axes=plt.subplots(1,2,figsize=(11,4.4))',extra={'save':savefig})
run('reservoir_plot_source.py',start=' fig,axes=plt.subplots(1,2,figsize=(11,4.4))\n for T,ax in zip([80,300],axes):\n  for support',end=' assert tests',extra={'save':savefig})
# Elementary analytic kernels replay the same pure formula functions at the original 501 points.
from scipy.special import logsumexp,gammaln
KB=8.617333262145e-5;HBAR=6.582119569e-16
params=json.loads(note(I/'kernel_parameters.json').read_text());source=note(I/'kernel_source.py').read_text();functions=source[source.index('def classical_lograte'):source.index('def cycle(')]
ns=globals().copy();ns['P']=params;exec(compile(functions,'kernel_functions','exec'),ns)
colors={80:'#244c8c',150:'#5aa6a8',300:'#d97a22',330:'#a74c65'}
block=source[source.index('    fig,axs=plt.subplots(1,2,figsize=(10,4))'):source.index('    fig,axs=plt.subplots(1,2,figsize=(10,4))',source.index('    fig,axs=plt.subplots(1,2,figsize=(10,4))')+1)]
ns.update(FIG=F,colors=colors);exec(compile(textwrap.dedent(block),'kernel_plot','exec'),ns)
b=D/'temperature_mechanism_audit_20261002/results';rows=read('transport_capacity.csv',b)
a=source.index('    fig,axs=plt.subplots(1,2,figsize=(10,4))',source.index('    fig,axs=plt.subplots(1,2,figsize=(10,4))')+1);block=source[a:source.index("    tests['all_assertions_passed']",a)]
ns=globals().copy();ns.update(R=b,FIG=F,frac=1e-4);exec(compile(textwrap.dedent(block),'dos_plot','exec'),ns)
assert len(outputs)==12,len(outputs)
(P/'qa/frozen_plot_font_replay.json').write_text(json.dumps({'scope':'font-only replay of original plotting blocks; pure analytic kernel evaluated at original same501 points; all other scientific values loaded from frozen results','inputs':seen,'outputs':outputs},ensure_ascii=False,indent=2))
print('12 frozen plots restyled with original data and plotting blocks')
