"""Independent read-only audit of coupled source; no production source edits."""
from pathlib import Path
import sys,json,hashlib,numpy as np
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from coupled import Coupled
sys.path.insert(0,str(ROOT.parent/'nonlocal_adapter_300K_20261003'))
from adapter import Mesh,Reaction,SpectralRates,Adapter
from rates import KB
out={'source_sha256':hashlib.sha256((ROOT/'code/coupled.py').read_bytes()).hexdigest(),'snapshots':[]}
for N in (81,161,321):
 for distance in (0,1):
  d=Coupled(N,distance); dx=d.h*d.p.d;x=d.x*d.p.d
  faces=np.r_[x[0]-dx/2,(x[:-1]+x[1:])/2,x[-1]+dx/2]
  kernel=SpectralRates(np.array([0.]),np.ones(1),np.array([d.p.Eg]),np.ones(1),d.p.Eg/2,
   dict(t0=d.t0,xi=d.xi,lam=d.lam),dict(t0=d.t0,xi=d.xi,lam=d.lam),temperature=d.vt/KB)
  # Matching numerical kBT isolates assembly from inherited rounded Boltzmann constant.
  aa=Adapter(Mesh(faces,x),[Reaction(h,t,l,a,kernel,.5) for h,t,l,a in zip(d.xH,d.xt,d.xL,d.a)])
  for V in (0,-1,-5,-15):
   path=ROOT/'data'/f'd{distance}_N{N}_V{V}.npz'
   if not path.exists():continue
   z=np.load(path)['z'];o=d.evaluate(z,V)
   state=np.r_[z[:,0]*d.vt,o['muH'],o['muL'],np.full(len(d.a),.5)]
   r=aa.evaluate(state,steady=True)
   gross=(o['aH']*o['bL']+o['aL']*o['bH'])/(o['aH']+o['bH']+o['aL']+o['bL'])
   scale=max(np.max(d.L.T@(d.a*gross)/dx),np.max(d.H.T@(d.a*gross)/dx),1e-300)
   rate_error=max(np.max(abs(r['Sn']+o['Rn'])),np.max(abs(r['Sp']+o['Rp'])))/scale
   charge_error=np.max(abs(r['rho']/d.q-o['rhot']))/max(np.max(abs(o['rhot'])),d.p.Nt*1e-10)
   ferror=np.max(abs(r['f']-o['actual_f']))
   jinternal=np.r_[0.,-d.q*dx*np.cumsum((o['Rn']-o['Rp'])[1:-1])]
   current_error=np.max(abs(r['Jtransfer'][1:-1]-jinternal))/(d.q*dx*scale)
   res=d.residual(z,V).reshape(N,3)
   rp=np.diff(z[:,0],n=2)/d.h**2+(o['ps']-o['ns']+o['rhot']/d.n0)[1:-1]
   poisson_error=np.max(abs(res[1:-1,0]-rp))
   l,_=d.ledger(z,V)
   row=dict(N=N,distance_nm=distance,V=V,adapter_source_scaled_error=float(rate_error),adapter_charge_relative_error=float(charge_error),adapter_f_absolute_error=float(ferror),adapter_current_scaled_error=float(current_error),poisson_assembly_absolute_error=float(poisson_error),energy_relative_error=l['energy_relative_error'],charge_relative_error=l['charge_relative_error'],current_relative_error=l['relative_error'],trap_areal_cm2=float(sum(d.a)),boundary_source_max=float(max(abs(o['Rn'][[0,-1]]).max(),abs(o['Rp'][[0,-1]]).max())))
   if distance==0:
    fn=expit(o['y'][:,1]+np.log(d.n0/d.p.Nc));fp=expit(o['y'][:,2]+np.log(d.n0/d.p.Nc))
    kc=d.kernel(d.p.Eg/2);kr=d.kernel(-d.p.Eg/2)
    denom=kc*(2-fn-fp)+kr*(fn+fp)
    ff=(kc*(1-fp)+kr*fn)/denom
    jj=(kc*kc*(1-fn)*(1-fp)-kr*kr*fn*fp)/denom
    density=np.asarray(d.T.T@d.a/dx)
    expected=density*jj
    row['analytic_local_source_scaled_error']=float(max(np.max(abs(expected+o['Rn'])),np.max(abs(expected+o['Rp'])))/scale)
    row['analytic_local_charge_relative_error']=float(np.max(abs(density*(.5-ff)-o['rhot']))/d.p.Nt)
   out['snapshots'].append(row)
# Independent one-column complex steps detect incorrect sparsity/coloring, and
# a real five-point directional derivative checks complex-step calculus itself.
d=Coupled(81,1);V=-5.;z=np.load(ROOT/'data/d1_N81_V-5.npz')['z'];J=d.jacobian(z,V,0).toarray();direct=np.zeros_like(J)
for k in range(z.size):
 zz=z.astype(complex).ravel();zz[k]+=1e-26j;direct[:,k]=d.residual(zz.reshape(z.shape),V).imag/1e-26
scale=np.maximum(np.max(abs(direct),axis=1),1e-300)
out['colored_vs_column_complex_step_scaled_error']=float(np.max(np.abs(J-direct)/scale[:,None]))
rng=np.random.default_rng(803);v=rng.normal(size=z.shape);v/=np.max(abs(v));h=1e-4
rs=[d.residual(z+m*h*v,V) for m in (-2,-1,1,2)]
fd=(rs[0]-8*rs[1]+8*rs[2]-rs[3])/(12*h)
out['five_point_vs_jacobian_direction_scaled_error']=float(np.max(abs(fd-J@v.ravel())/scale))
out['jacobian_nnz']=int(np.count_nonzero(J))
# Directional currents have the proper thermodynamic sign.
o=d.evaluate(z,V);out['min_edge_transport_dissipation_Wcm2']=float(min(np.min(o['Jn']*np.diff(o['muL'])),np.min(o['Jp']*np.diff(o['muH']))))
out['min_reaction_dissipation_pertrap_eVs']=float(np.min(o['j']*(o['mh']-o['ml'])))
out['min_bim_dissipation_eVcm3s']=float(np.min(o['bim']*(o['muL']-o['muH'])))
(ROOT/'independent_audit/results.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='snapshots'},indent=2))
print('Snapshot count:',len(out['snapshots']))
for key in out['snapshots'][0]:
 if key not in ('N','distance_nm','V'):print(key,max(row.get(key,0) for row in out['snapshots']))
assert len(out['snapshots'])==24
assert out['colored_vs_column_complex_step_scaled_error']<1e-12
assert out['five_point_vs_jacobian_direction_scaled_error']<2e-8
assert max(x['adapter_source_scaled_error'] for x in out['snapshots'])<2e-10
assert max(x['adapter_charge_relative_error'] for x in out['snapshots'])<2e-11
assert max(x.get('analytic_local_source_scaled_error',0) for x in out['snapshots'])<2e-10
assert max(x['boundary_source_max'] for x in out['snapshots'])==0
