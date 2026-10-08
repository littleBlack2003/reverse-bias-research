from pathlib import Path
import csv,json,sys,importlib.util,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import loggamma,gammaln
import kernel_benchmark as b
root=Path(__file__).resolve().parent;out=root/'results'
checks={}
# Independent spectral integration versus closed Gamma characteristic function.
errs=[]
for T in [80,300]:
 for Ec in [.001,.005,.02]:
  beta=1/(b.KB*T);a=1/(beta*Ec);alpha=.05/Ec
  for u in [5,35,100]:
   def realint(v):
    if v==0:return 0.
    return alpha*np.exp(-v)/v*(2*np.sin(Ec*v*u/2)**2)/np.tanh(beta*Ec*v/2)
   re=quad(realint,0,100,epsabs=1e-10,epsrel=1e-10)[0]
   im=quad(lambda v:alpha*np.exp(-v)*np.sin(Ec*v*u)/v if v else alpha*Ec*u,0,100,epsabs=1e-10,epsrel=1e-10)[0]
   lc=alpha*(loggamma(a+1j*u/beta)+loggamma(1+a-1j*u/beta)-np.log(beta*Ec)-2*gammaln(1+a))
   errs.append(abs(lc+re+1j*im))
checks['spectral_quadrature_vs_gamma_max_log_C_error']=max(errs)
# Low temperature gamma density limit independently known for Ohmic exponential cutoff.
errs=[]
for E in [.01,.05,.2]:
 Ec=.005;alpha=.05/Ec
 lp=(alpha-1)*np.log(E)-E/Ec-gammaln(alpha)-alpha*np.log(Ec)
 errs.append(abs(b.ohmic_logp(E,.0125,Ec)-lp))
checks['T_00125K_vs_exact_zero_T_Gamma_max_log_P_difference']=max(errs)
# Exact upstream code crosscheck.
spec=importlib.util.spec_from_file_location('upstream',root/'inputs/old_cycle/model.py');upstream=importlib.util.module_from_spec(spec);spec.loader.exec_module(upstream)
checks['old_kernel_vs_current_baseline_max_log_error']=max(abs(float(upstream.logkernel(dg,T))-b.quantum_logk(dg,T)) for T in [80,150,300] for dg in [-.9,-.2,0,.2,.9])
checks['two_mode_L18_L24_max_log_error']=max(abs(b.quantum_logk(dg,T,.005,L=18,modes=((.1,.75),(.2,.375)))-b.quantum_logk(dg,T,.005,L=24,modes=((.1,.75),(.2,.375)))) for dg in [-.9,.65] for T in [80,300])
# Existing graph times: all exits, not only net flux. Read previous verified edge ledger.
p=root/'inputs/representative_edges.csv'
rows=list(csv.DictReader(p.open()));groups={}
for r in rows:
 key=(r['T_K'],r['F_MVcm'],r['variant'],r['s_nm'])
 groups.setdefault(key,[]).append(r)
ans=[]
for (T,F,variant,s),rr in groups.items():
 exits={}
 for r in rr:
  exits[r['source']]=exits.get(r['source'],0)+float(r['forward_rate_s'])
  exits[r['destination']]=exits.get(r['destination'],0)+float(r['reverse_rate_s'])
 k=max(exits.values());tau=.01/k
 ans.append(dict(T_K=T,F_MVcm=F,variant=variant,s_nm=s,max_total_exit_rate_s=k,min_state_dwell_ns=1e9/k,tau_relax_for_one_percent_reset_ns=tau*1e9,largest_edge_gross_s=max(float(r['gross_flux_s']) for r in rr),scope='existing MLJ graph only; relaxation criterion is a conditional 1% diagnostic'))
b.savecsv('prior_graph_thermalization_scales.csv',ans)
# Focus actual hop energies without interpreting them as full graph-current ratios.
ans=[]
rr=[r for r in rows if r['T_K']=='80' and r['F_MVcm']=='1.5' and r['variant']=='explicit' and r['s_nm']=='6.0']
seen=set()
for r in sorted(rr,key=lambda r:float(r['gross_flux_s']),reverse=True):
 if r['process'].startswith('bath'):continue
 dg=float(r['energy_change_eV']);key=(r['process'],round(dg,12))
 if key in seen:continue
 seen.add(key)
 for Ec in [.001,.005,.02]:
  lk=b.quantum_logk(dg,80,Ec);lold=b.quantum_logk(dg,80)
  ans.append(dict(process=r['process'],source=r['source'],destination=r['destination'],delta_G_eV=dg,Ec_meV=Ec*1000,quantum_over_existing=np.exp(lk-lold),scope='single edge kernel only; not current ratio'))
 if len(seen)>=8:break
b.savecsv('representative_edge_kernel_sensitivity.csv',ans)
assert checks['spectral_quadrature_vs_gamma_max_log_C_error']<1e-8
assert checks['T_00125K_vs_exact_zero_T_Gamma_max_log_P_difference']<1e-5
assert checks['old_kernel_vs_current_baseline_max_log_error']<1e-10
assert checks['two_mode_L18_L24_max_log_error']<1e-8
checks['all_assertions_passed']=True
(out/'focused_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
