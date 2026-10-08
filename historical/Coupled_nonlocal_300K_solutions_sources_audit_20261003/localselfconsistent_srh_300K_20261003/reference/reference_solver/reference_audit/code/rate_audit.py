"""Independent event/cycle audit and controlled solvers; source files untouched.
Units cm,s,V,A,eV,K; positive R is recombination, negative J extraction.
"""
from pathlib import Path
import sys,json,hashlib,time,traceback,csv
import numpy as np
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'source_snapshot'))
from solvers.pf_solver import PFDevice,PFParams
from solvers.base_dd_solver import Device,Params
BASE=json.loads((ROOT/'source_snapshot/model_comparison/baseline_ff78.json').read_text())['full_local_parameters']

def device(nodes=81,**kw):
 p=dict(BASE);p.update(kw);return PFDevice(PFParams(**p),nodes)

class ControlledPF(PFDevice):
 """Diagnostic variant: turn PF lowering off per branch, preserving e/c.
 trap_e_depth controls n depth; p depth=Eg-n depth. No actual material fit.
 """
 def __init__(self,p,nodes=81,branch_n=1.,branch_p=1.,trap_e_depth=None):
  self.branch_n,self.branch_p=branch_n,branch_p
  self.trap_e_depth=p.Eg/2 if trap_e_depth is None else trap_e_depth
  super().__init__(p,nodes)
 def trap_rates(self,y):
  p=self.p;phi=y[:,0]*self.vt;n=np.exp(y[:,1])*self.n0;holes=np.exp(y[:,2])*self.n0
  e=-np.gradient(phi,self.h)/p.d;F=np.sqrt(e*e+1.)
  delta=p.pf_strength*np.sqrt(self.q*F/(np.pi*self.eps))
  dn=self.trap_e_depth;dp=p.Eg-dn;c0=p.capture*(1-p.nonlocal_fraction)
  an=c0*n;ap=c0*holes;en=c0*p.Nc*np.exp(-dn/self.vt);ep=c0*p.Nc*np.exp(-dp/self.vt)
  fill=an+ep;empty=ap+en;links=[];pref=p.capture*p.Nc*p.nonlocal_fraction/2;w=p.barrier_smoothing_eV
  for t,end in self.pairs:
   du=phi[end]-phi[t]
   for carrier,energy,conc,depth,mult in [('n',dn-du,n[end],dn,self.branch_n),('p',dp+du,holes[end],dp,self.branch_p)]:
    base=self.smoothmax(0.,depth-mult*delta,w)
    barrier=self.smoothmax(base[t],energy,w)
    emission=pref*np.exp(-barrier/self.vt);capture=(pref/p.Nc)*np.exp(-(barrier-energy)/self.vt)
    inc=capture*conc
    if carrier=='n':fill[t]+=inc;empty[t]+=emission
    else:fill[t]+=emission;empty[t]+=inc
    links.append((carrier,t,end,inc,emission,energy,barrier))
  f=fill/(fill+empty);fb=empty/(fill+empty)
  rn=p.Nt*(an*fb-en*f);rp=p.Nt*(ap*f-ep*fb);rn[[0,-1]]=0;rp[[0,-1]]=0
  for carrier,t,end,inc,em,energy,barrier in links:
   net=p.Nt*(inc*fb[t]-em*f[t]) if carrier=='n' else p.Nt*(inc*f[t]-em*fb[t])
   np.add.at(rn if carrier=='n' else rp,end,net)
  return rn,rp,f,e,delta,links,fb

def channels(dev,y):
 """Build per-trap channel arrays from emitted event metadata, not Rn/Rp."""
 n=np.exp(y[:,1])*dev.n0;p=np.exp(y[:,2])*dev.n0
 rates=dev.trap_rates(y);f,fb=rates[2],rates[-1]
 dn=getattr(dev,'trap_e_depth',dev.p.Eg/2);dp=dev.p.Eg-dn
 c0=dev.p.capture*(1-dev.p.nonlocal_fraction)
 ch={'n':[], 'p':[]}
 for i in range(dev.nodes):
  ch['n'].append([(i,c0*n[i],c0*dev.p.Nc*np.exp(-dn/dev.vt))])
  ch['p'].append([(i,c0*p[i],c0*dev.p.Nc*np.exp(-dp/dev.vt))])
 for carrier,t,end,inc,em,energy,barrier in rates[5]:
  for a,b,c,d in zip(t,end,inc,em):ch[carrier][a].append((int(b),float(c),float(d)))
 return ch,rates

def audit(dev,y,V,light):
 ch,rates=channels(dev,y);rn,rp,f,field,delta,links,fb=rates
 n,p,jn,jp,bim=dev.transport(y);jn=jn*dev.j0;jp=jp*dev.j0
 phi=y[:,0]*dev.vt;dx=dev.h*dev.p.d;fac=dev.p.Nt*dx
 independent_n=np.zeros(dev.nodes);independent_p=np.zeros(dev.nodes);Jnl=np.zeros(dev.nodes-1)
 gen=np.zeros(dev.nodes);rec=np.zeros(dev.nodes);tn=np.zeros(dev.nodes);tp=np.zeros(dev.nodes)
 work_pair=0.;work_transport=0.;energy_pair=0.;balance=[]
 # Reconstruction via complete cycles; no observed rn/rp inputs used here.
 for t in range(1,dev.nodes-1):
  cn=ch['n'][t];cp=ch['p'][t]
  An=sum(z[1] for z in cn);En=sum(z[2] for z in cn);Ap=sum(z[1] for z in cp);Ep=sum(z[2] for z in cp)
  S=An+En+Ap+Ep
  gen[t]=En*Ep/S;rec[t]=An*Ap/S
  balance.append((An+Ep)*fb[t]-(Ap+En)*f[t])
  for e,a,em in cn:
   for h,b,ep in cp:
    g=em*ep/S;r=a*b/S;net=g-r
    independent_n[e]-=dev.p.Nt*net;independent_p[h]-=dev.p.Nt*net
    work_pair+=dev.q*fac*net*(phi[e]-phi[h])
    energy_pair+=dev.q*fac*net*(dev.p.Eg-phi[e]+phi[h])
   for s,a2,e2 in cn:
    flow=a*e2/S # capture at e then emit at s
    tn[t]+=flow
    independent_n[e]+=dev.p.Nt*flow;independent_n[s]-=dev.p.Nt*flow
    work_transport+=dev.q*fac*flow*(phi[s]-phi[e])
  for h,b,ep in cp:
   for s,b2,e2 in cp:
    flow=b*e2/S
    tp[t]+=flow
    independent_p[h]+=dev.p.Nt*flow;independent_p[s]-=dev.p.Nt*flow
    work_transport+=dev.q*fac*flow*(phi[h]-phi[s])
 # Event charge displacement integrated over faces, independent of cumulative-source formula.
 work_events=0.
 for carrier,t,end,inc,em,energy,barrier in links:
  for a,b,c,d in zip(t,end,inc,em):
   net=c*fb[a]-d*f[a] if carrier=='n' else c*f[a]-d*fb[a]
   sign=1 if carrier=='n' else -1
   Jnl[min(a,b):max(a,b)]+=sign*dev.q*fac*net*np.sign(b-a)
   work_events-=sign*dev.q*fac*net*(phi[b]-phi[a])
 o=dev.observe(y);G=light*dev.p.Jgen/(dev.q*dev.p.d)
 pairnet=fac*np.sum(gen-rec);bim_areal=dx*np.sum(bim[1:-1]);photo_areal=dx*(dev.nodes-2)*G
 n_out=(jn[0]-jn[-1])/dev.q;p_out=(jp[-1]-jp[0])/dev.q
 scale=max(1.,np.max(np.abs(rn)),np.max(np.abs(rp)))
 gross_scale=max(1.,dev.p.Nt*max(sum(z[1]+z[2] for z in ch['n'][t])+sum(z[1]+z[2] for z in ch['p'][t]) for t in range(1,dev.nodes-1)))
 out=dict(V=V,light=light,nodes=dev.nodes,escape_actual_nm=dev.actual_escape_nm,J_mAcm2=o['J_Acm2']*1000,
 terminal_difference_Acm2=o['terminal_difference_Acm2'],current_spread_Acm2=o['spread_Acm2'],nonlocal_boundary_Acm2=o['nonlocal_boundary_error_Acm2'],
 pair_generation_mAcm2=dev.q*fac*np.sum(gen)*1000,pair_recombination_mAcm2=dev.q*fac*np.sum(rec)*1000,pair_net_mAcm2=dev.q*pairnet*1000,
 electron_transport_gross_mAcm2=dev.q*fac*np.sum(tn)*1000,hole_transport_gross_mAcm2=dev.q*fac*np.sum(tp)*1000,
 photo_mAcm2=dev.q*photo_areal*1000,bim_net_mAcm2=dev.q*bim_areal*1000,
 electron_out_mAcm2=dev.q*n_out*1000,hole_out_mAcm2=dev.q*p_out*1000,
 Jn_left_mAcm2=jn[0]*1000,Jn_right_mAcm2=jn[-1]*1000,Jp_left_mAcm2=jp[0]*1000,Jp_right_mAcm2=jp[-1]*1000,
 electron_continuity_error_Acm2=dev.q*(n_out-(pairnet+photo_areal-bim_areal)),hole_continuity_error_Acm2=dev.q*(p_out-(pairnet+photo_areal-bim_areal)),
 cycle_source_absolute_error_cm3s=max(float(np.max(abs(rn-independent_n))),float(np.max(abs(rp-independent_p)))),
 cycle_source_error_over_gross=max(float(np.max(abs(rn-independent_n))),float(np.max(abs(rp-independent_p))))/gross_scale,
 cycle_source_relative_error=max(np.max(abs(rn-independent_n)),np.max(abs(rp-independent_p)))/scale,
 displacement_current_max_error_Acm2=float(np.max(abs(Jnl-o['Jnonlocal']))),
 fieldwork_event_Wcm2=work_events,fieldwork_face_Wcm2=float(-np.dot(Jnl,np.diff(phi))),fieldwork_pair_Wcm2=work_pair,fieldwork_transport_Wcm2=work_transport,
 endpoint_pair_energy_Wcm2=energy_pair,bandgap_pair_energy_Wcm2=dev.q*dev.p.Eg*pairnet,
 total_fieldwork_Wcm2=float(-np.dot(jn+jp+Jnl,np.diff(phi))),dd_fieldwork_Wcm2=float(-np.dot(jn+jp,np.diff(phi))),
 port_input_power_Wcm2=V*o['J_Acm2'],trap_balance_max_s=float(max(abs(np.array(balance)))),
 min_f=float(min(f[1:-1])),min_fbar=float(min(fb[1:-1])),max_field_Vcm=float(max(abs(field))))
 arrays=dict(phi=phi,n=n*dev.n0,p=p*dev.n0,f=f,fb=fb,En=np.array([sum(z[2] for z in x) for x in ch['n']]),Ep=np.array([sum(z[2] for z in x) for x in ch['p']]),pair_gen=gen,pair_rec=rec,Rn=rn,Rp=rp,Jnl_independent=Jnl,Jn=jn,Jp=jp,y=y)
 return out,arrays

def savejson(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')

def run_case(name,dev,targets=(-5,-10,-15,-20),dv=.5,light=0):
 rows=[];trace=[];failures=[];start=time.monotonic();yp=None;prev=0.
 try:
  y,diag=dev.solve(0,0);o,a=audit(dev,y,0,0);o.update(diag);rows.append(o);np.savez_compressed(ROOT/'data'/f'{name}_0.npz',**a)
  for target in targets:
   steps=int(np.ceil(abs(target-prev)/dv));vstart=prev
   for j,V in enumerate(np.linspace(prev,target,steps+1)[1:]):
    before=time.monotonic();y,diag=dev.solve(float(V),light,y);trace.append(dict(V=float(V),light=light,seconds=time.monotonic()-before,**diag))
   prev=target;o,a=audit(dev,y,target,light);o.update(diag);rows.append(o);np.savez_compressed(ROOT/'data'/f'{name}_{target}.npz',**a)
   print(name,target,o['J_mAcm2'],o['pair_net_mAcm2'],o['current_spread_Acm2'],flush=True)
   savejson(ROOT/'data'/f'{name}.json',dict(parameters=asdict(dev.p),variant={k:getattr(dev,k,None) for k in ['branch_n','branch_p','trap_e_depth']},seconds=time.monotonic()-start,rows=rows,trace=trace,failures=failures))
 except Exception as exc:
  failures.append(dict(type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),last_successful_V=trace[-1]['V'] if trace else None));print('FAIL',name,failures[-1],flush=True)
 savejson(ROOT/'data'/f'{name}.json',dict(parameters=asdict(dev.p),variant={k:getattr(dev,k,None) for k in ['branch_n','branch_p','trap_e_depth']},seconds=time.monotonic()-start,rows=rows,trace=trace,failures=failures))
 return rows

if __name__=='__main__':
 import argparse
 pa=argparse.ArgumentParser();pa.add_argument('case');args=pa.parse_args();name=args.case
 if name.startswith('baseline'):
  nodes=int(name.split('_')[1]);dev=device(nodes);dv=.25 if name.endswith('fine') else .5
 elif name in ['n_only','p_only','neither']:
  dev=ControlledPF(device().p,81,branch_n=float(name=='n_only'),branch_p=float(name=='p_only'));dv=.5
 else:raise ValueError(name)
 run_case(name,dev,dv=dv)
