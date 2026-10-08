"""Transport-only equation benchmark, no reaction, no material fit.
Reproduces published parameterization curves, NOT the original master-equation data.
"""
import json,numpy as np
from pathlib import Path
from gaussian_transport import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
rows=[];sig=.1;a=1e-7
for T in [300,340,380,420]:
 m=EGDM(sig,a,1.)
 for n in [0,1e16]:
  v=m.evaluate(T,n,1e5);rows.append(dict(T=T,n=n,F=1e5,mu=float(v['mu']),s=v['s'],reduced_field=float(v['reduced_field'])))
fig,axs=plt.subplots(1,3,figsize=(13,3.6))
for s in [2,3,4,5,6]:
 T=sig/(KB*s);m=EGDM(sig,a,1.8e-9)
 c=np.logspace(-6,-2,160);mu=m.evaluate(T,c/a**3,0)['mu'];axs[0].loglog(c,mu,label=str(s))
 u=np.linspace(0,3,150);mu=m.evaluate(T,1e-5/a**3,u*sig/a)['mu'];axs[1].semilogy(u,mu,label=str(s))
axs[0].set(xlabel='Occupation n a³',ylabel='Mobility / intrinsic μ₀',title='Pasveer Eq. (3), F=0')
axs[1].set(xlabel='Reduced field qFa/σ',ylabel='Mobility / intrinsic μ₀',title='Pasveer Eqs. (3–5), n a³=10⁻⁵')
T=np.linspace(300,420,121);m=EGDM(sig,a,1)
ratio=np.array([m.evaluate(t,1e16,1e5)['mu']/m.evaluate(300,1e16,1e5)['mu'] for t in T]);axs[2].plot(T,ratio)
axs[2].set(xlabel='Temperature (K)',ylabel='μ(T) / μ(300 K)',title='Illustrative fixed n,F; σ=0.1 eV')
for ax in axs:ax.grid(alpha=.2)
axs[0].legend(title='σ / kT');axs[1].legend(title='σ / kT')
fig.tight_layout();fig.savefig(ROOT/'data/transport_only_benchmark.png',dpi=180)
T=np.array([300,340,380,420]);mu=np.array([m.evaluate(t,1e16,1e5)['mu'] for t in T]);Ea=-KB*np.polyfit(1/T,np.log(mu),1)[0]
(ROOT/'data/transport_only_benchmark.json').write_text(json.dumps(dict(rows=rows,illustrative_mu420_over_300=float(mu[-1]/mu[0]),illustrative_apparent_Ea_eV=float(Ea),not_material_fit=True,reproduces_equations_not_master_equation_data=True),indent=2))
print('mu420/mu300',mu[-1]/mu[0],'apparent Ea',Ea)
