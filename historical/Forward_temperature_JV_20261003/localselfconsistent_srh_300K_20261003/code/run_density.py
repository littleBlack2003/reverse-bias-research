"""One-factor Nt study; stable QF/cycle arithmetic, unchanged physical equations."""
from pathlib import Path
import sys,json,time,traceback,numpy as np
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference/code'))
from stable_cycle_solver import StableCycleQF,PFParams,solve_stable
BASE=json.loads((ROOT/'baseline_ff78.json').read_text())['full_local_parameters']
def save(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False,default=float),encoding='utf-8')
class DensityDevice(StableCycleQF):
 def __init__(self,*a,light=0,**kw):
  self.light=light;super().__init__(*a,**kw)
  if self.p.nonlocal_fraction==0:
   rr=[];cc=[]
   for i in range(self.nodes):
    for a in range(3):
     for k in range(max(0,i-1),min(self.nodes,i+2)):
      for b in range(3):rr.append(3*i+a);cc.append(3*k+b)
   self.rr=np.array(rr);self.cc=np.array(cc);self.ncolors=9;self.colors=self.cc%9
 def residual(self,z,V,light=0):return super().residual(z,V,self.light)
 def physical(self,z,V):
  o=self.evaluate(z,V);dx=self.h*self.p.d;Jn=o['Jn'];Jp=o['Jp'];photo=self.light*self.p.Jgen*(self.nodes-2)/(self.nodes-1)
  pair=self.q*self.p.Nt*dx*np.sum(o['pair'][1:-1]);bim=self.q*dx*np.sum(o['bim'][1:-1]);expected=pair+photo-bim
  nout=Jn[0]-Jn[-1];pout=Jp[-1]-Jp[0];jnl=np.r_[0.,-self.q*dx*np.cumsum((o['Rn']-o['Rp'])[1:-1])];jt=Jn+Jp+jnl
  error=max(abs(nout-expected),abs(pout-expected),np.ptp(jt),abs(jnl[-1]));scale=max(abs(pair),abs(bim),photo,*abs(Jn[[0,-1]]),*abs(Jp[[0,-1]]),1e-300)
  budget=(1e-17 if V==0 and self.light==0 else 1e-280)+1e-8*scale
  return dict(J_Acm2=float((Jn[0]+Jp[0]+Jn[-1]+Jp[-1])/2),pair_source_Acm2=float(pair),photo_Acm2=float(photo),bim_Acm2=float(bim),expected_particle_output_Acm2=float(expected),electron_out_Acm2=float(nout),hole_out_Acm2=float(pout),electron_cycle_error_Acm2=float(nout-expected),hole_cycle_error_Acm2=float(pout-expected),Jn_left_Acm2=float(Jn[0]),Jn_right_Acm2=float(Jn[-1]),Jp_left_Acm2=float(Jp[0]),Jp_right_Acm2=float(Jp[-1]),nonlocal_mean_Acm2=float(np.mean(jnl)),nonlocal_max_Acm2=float(max(abs(jnl))),total_current_spread_Acm2=float(np.ptp(jt)),nonlocal_boundary_error_Acm2=float(abs(jnl[-1])),budget_Acm2=float(budget),error_Acm2=float(error),gate_passed=bool(error<=budget),relative_error=float(error/scale))
def details(d,z,V,info):
 o=d.evaluate(z,V);n=o['ns']*d.n0;p=o['ps']*d.n0;f=o['f'];phi=z[:,0]*d.vt;field=-np.gradient(phi,d.h*d.p.d);Ef=-np.diff(phi)/(d.h*d.p.d);sl=slice(1,-1);bulk=(d.x>=.1)&(d.x<=.9);dx=d.h*d.p.d
 rho=p-n+d.p.Nt*(.5-f);field_jump=d.eps*(Ef[-1]-Ef[0]);charge=d.q*dx*sum(rho[sl]);gross=d.q*dx*sum((p+n+d.p.Nt*abs(.5-f))[sl])
 out=dict(V=float(V),U=float(-V),light=d.light,Nt_cm3=d.p.Nt,nodes=d.nodes,**{k:v for k,v in info.items() if k!='history'},field_mean_bulk_Vcm=float(np.mean(abs(field[bulk]))),field_center_Vcm=float(abs(field[d.nodes//2])),field_min_bulk_Vcm=float(min(abs(field[bulk]))),field_max_bulk_Vcm=float(max(abs(field[bulk]))),field_max_Vcm=float(max(abs(field))),f_mean=float(np.mean(f[sl])),f_min=float(min(f[sl])),f_max=float(max(f[sl])),f_center=float(f[d.nodes//2]),n_mean_cm3=float(np.mean(n[sl])),p_mean_cm3=float(np.mean(p[sl])),trap_charge_Ccm2=float(d.q*dx*d.p.Nt*sum((.5-f)[sl])),trap_abs_charge_Ccm2=float(d.q*dx*d.p.Nt*sum(abs(.5-f)[sl])),carrier_charge_Ccm2=float(d.q*dx*sum((p-n)[sl])),carrier_abs_charge_Ccm2=float(d.q*dx*sum(abs(p-n)[sl])),poisson_charge_Ccm2=float(charge),field_jump_charge_Ccm2=float(field_jump),charge_relative_error=float(abs(charge-field_jump)/max(gross,1e-300)))
 return out,dict(phi=phi,n=n,p=p,f=f,field=field,pair_rate_s=o['pair'],Jn=o['Jn'],Jp=o['Jp'],Rn=o['Rn'],Rp=o['Rp'],bim=o['bim'])
def advance(d,z,oldV,V,trace,depth=0):
 guess=z.copy();shift=-(V-oldV)/d.vt;guess[:,0]+=shift*d.x;guess[:,1]+=shift*(1-d.x);guess[:,2]+=shift*d.x
 zz,info=solve_stable(d,V,guess);trace.append(dict(V=V,oldV=oldV,depth=depth,**{k:v for k,v in info.items() if k!='history'}))
 if info.get('status'):
  if depth>=6:raise RuntimeError(f'Unresolved at {V}: {info.get("status")}')
  mid=(oldV+V)/2;half,_=advance(d,z,oldV,mid,trace,depth+1);return advance(d,half,mid,V,trace,depth+1)
 return zz,info
def run(mode,Nt,light,nodes=81,step=.25,maxU=30,forward=True):
 tag=f'{mode}_Nt{Nt:.0e}_L{light:g}_N{nodes}_h{step:g}';p=dict(BASE,Nt=Nt,nonlocal_fraction=.5 if mode=='PF' else 0.);d=DensityDevice(PFParams(**p),nodes,trap_e_depth=.65);start=time.monotonic();rows=[];states=[];arrays=[];trace=[];fw=[];failure=None
 try:
  z=np.c_[d.x*d.vbi/d.vt,np.zeros(nodes),np.zeros(nodes)];z,info=solve_stable(d,0,z)
  if info.get('status'):raise RuntimeError('dark equilibrium '+info['status'])
  equilibrium=dict(**{k:v for k,v in info.items() if k!='history'},max_qf=float(np.max(abs(z[:,1:]))))
  for L in np.linspace(0,light,21)[1:] if light else []:
   d.light=float(L);z,info=solve_stable(d,0,z)
   if info.get('status'):raise RuntimeError('Light ramp '+info['status'])
  d.light=light;z0=z.copy();row,a=details(d,z,0,info);rows.append(row);states.append(z.copy());arrays.append(a)
  prev=0.
  for U in np.arange(step,maxU+step/2,step):
   V=-float(U);z,info=advance(d,z,prev,V,trace);prev=V;row,a=details(d,z,V,info);rows.append(row);states.append(z.copy());arrays.append(a)
   if abs(U-round(U/5)*5)<1e-8:print(tag,'reverse',V,row['J_Acm2'],row['relative_error'],time.monotonic()-start,flush=True)
  if forward and light:
   z=z0.copy();prev=0.;row=rows[0].copy();fw.append(row)
   for V in np.arange(.01,1.071,.01):
    z,info=advance(d,z,prev,float(V),trace);prev=float(V);row,a=details(d,z,V,info);fw.append(row)
   print(tag,'forward complete',time.monotonic()-start,flush=True)
 except Exception as exc:failure=dict(message=str(exc),traceback=traceback.format_exc());print(tag,'FAIL',failure['message'],flush=True)
 np.savez_compressed(ROOT/'data'/f'{tag}.npz',V=np.array([r['V'] for r in rows]),z=np.array(states),**{k:np.array([a[k] for a in arrays]) for k in arrays[0]})
 save(ROOT/'data'/f'{tag}.json',dict(tag=tag,parameters=asdict(d.p),light=light,nodes=nodes,voltage_step=step,maximum_reverse_voltage=maxU,equilibrium=equilibrium,rows=rows,forward_rows=fw,trace=trace,failure=failure,seconds=time.monotonic()-start))
 return tag
if __name__=='__main__':run(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]) if len(sys.argv)>4 else 81,float(sys.argv[5]) if len(sys.argv)>5 else .25,float(sys.argv[6]) if len(sys.argv)>6 else 30,len(sys.argv)<8 or sys.argv[7]!='noforward')
