from run_density import *
import csv
from scipy.signal import savgol_filter,find_peaks
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq,minimize_scalar

def csvout(name,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (ROOT/'data'/name).open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def read(mode,Nt,L,N=81,h=.25):return json.loads((ROOT/'data'/f'{mode}_Nt{Nt:.0e}_L{L:g}_N{N}_h{h:g}.json').read_text())
def crossing(x,y,value):
 if y[0]>=value:return dict(U=None,status='already_at_or_above_at_zero')
 for i in range(1,len(x)):
  if y[i]>=value and y[i-1]<value:return dict(U=float(x[i-1]+(x[i]-x[i-1])*(value-y[i-1])/(y[i]-y[i-1])),status='linear_first_upcrossing')
 return dict(U=None,status='not_crossed_in_range')
def shape(x,y,w):
 h=x[1]-x[0];d1=savgol_filter(y,w,min(w-1,3),deriv=1,delta=h);d2=savgol_filter(y,w,min(w-1,3),deriv=2,delta=h);peaks,meta=find_peaks(d2,prominence=0);c=[]
 for ii,k in enumerate(peaks):
  if not 2<=x[k]<=x[-1]-2 or d2[k]<=0:continue
  a=np.argmin(abs(x-(x[k]-2)));b=np.argmin(abs(x-(x[k]+2)));ratio=float(d1[b]/d1[a]) if d1[a]>0 else None
  c.append(dict(U=float(x[k]),curvature=float(d2[k]),normalized_curvature=float(d2[k]/max(abs(y))),prominence_fraction=float((d2[k]-max(d2[a],d2[b]))/d2[k]),scipy_local_prominence_fraction=float(meta['prominences'][ii]/d2[k]),slope_before=float(d1[a]),slope_after=float(d1[b]),slope_ratio=ratio,passes=bool(ratio is not None and ratio>=4 and (d2[k]-max(d2[a],d2[b]))>=.25*d2[k])))
 c.sort(key=lambda z:z['curvature'],reverse=True);chosen=c[0] if c else None
 return dict(window_samples=w,sample_footprint_V=float(w*h),actual_support_V=float((w-1)*h),candidates=c,largest_positive_curvature=chosen,strict_knee_U=(chosen['U'] if chosen and chosen['passes'] else None)),d1,d2

def main():
 out=[];curves=[];shape_rows=[];thresholds=[];forward=[];frozen=[];check=[]
 for L in [0,1]:
  ref=read('PF',1e15,L);rr=ref['rows'];u=np.array([r['U'] for r in rr]);refpair=np.array([r['pair_source_Acm2'] for r in rr]);refi=-np.array([r['J_Acm2'] for r in rr]);refdd=read('DD',1e15,L);refdel=refi+np.array([r['J_Acm2'] for r in refdd['rows']]);arr=np.load(ROOT/'data'/f'DD_Nt1e+15_L{L}_N81_h0.25.npz');fz=[]
  for Nt in [1e14,1e15,1e16]:
   pp=dict(BASE,Nt=Nt);d=DensityDevice(PFParams(**pp),81,light=L);this=[]
   for V,z in zip(arr['V'],arr['z']):
    e=d.evaluate(z,float(V));pair=d.q*Nt*d.h*d.p.d*sum(e['pair'][1:-1]);this.append(pair);frozen.append(dict(U=float(-V),V=float(V),Nt_cm3=Nt,light=L,frozen_pair_Acm2=float(pair),frozen_pair_per_Nt=float(pair/Nt),frozen_fmean=float(np.mean(e['f'][1:-1]))))
   fz.append(np.array(this)/Nt)
  check.append(dict(light=L,frozen_max_abs_per_center_spread=float(max(np.max(abs(x-fz[1])) for x in fz)),frozen_max_relative_spread=float(max(np.max(abs(x-fz[1])/np.maximum(abs(fz[1]),1e-200)) for x in fz))))
  for Nt in [1e14,1e15,1e16]:
   pf=read('PF',Nt,L);dd=read('DD',Nt,L);pr=pf['rows'];dr=dd['rows'];I=-np.array([r['J_Acm2'] for r in pr]);I0=-np.array([r['J_Acm2'] for r in dr]);delta=I-I0;pair=np.array([r['pair_source_Acm2'] for r in pr]);shapes={}
   for cname,y in [('terminal',I),('excess',delta)]:
    shapes[cname]=[]
    for w in [3,5,7]:
     ss,s1,s2=shape(u,y,w);shapes[cname].append(ss);c=ss['largest_positive_curvature'];shape_rows.append(dict(Nt_cm3=Nt,light=L,curve=cname,window_samples=w,support_V=ss['actual_support_V'],curvature_landmark_U=c['U'] if c else None,strict_knee_U=ss['strict_knee_U'],slope_ratio=c['slope_ratio'] if c else None,prominence_fraction=c['prominence_fraction'] if c else None))
   ss,s1,s2=shape(u,delta,5)
   for k,r in enumerate(pr):
    row=dict(model='PF',**r,I_Acm2=float(I[k]),DD_I_Acm2=float(I0[k]),excess_I_Acm2=float(delta[k]),I_per_Nt=float(I[k]/Nt),excess_per_Nt=float(delta[k]/Nt),pair_per_Nt=float(pair[k]/Nt),normalized_pair_vs_Nt15=float((pair[k]/Nt)/(refpair[k]/1e15)) if abs(refpair[k])>1e-30 else None,normalized_excess_vs_Nt15=float((delta[k]/Nt)/(refdel[k]/1e15)) if abs(refdel[k])>1e-20 else None,excess_slope_Acm2V=float(s1[k]),excess_curvature_Acm2V2=float(s2[k]),log_excess_slope_perV=float(s1[k]/delta[k]) if delta[k]>1e-15 else None,log_excess_curvature_perV2=float(s2[k]/delta[k]-(s1[k]/delta[k])**2) if delta[k]>1e-15 else None);curves.append(row)
   for val in [.025,.05,.1]:thresholds.append(dict(Nt_cm3=Nt,light=L,current_mAcm2=val*1000,**crossing(u,I,val)))
   idx=(u>=10)&(u<=20);nr=(pair[idx]/Nt)/(refpair[idx]/1e15);nx=(delta[idx]/Nt)/(refdel[idx]/1e15)
   out.append(dict(Nt_cm3=Nt,light=L,shapes=shapes,pair_perNt_ratio_10to20_min=float(min(nr)),pair_perNt_ratio_10to20_max=float(max(nr)),excess_perNt_ratio_10to20_min=float(min(nx)),excess_perNt_ratio_10to20_max=float(max(nx))))
   for model,run in [('PF',pf),('DD',dd)]:
    if L:
     vr=np.array([r['V'] for r in run['forward_rows']]);jr=np.array([r['J_Acm2'] for r in run['forward_rows']]);f=PchipInterpolator(vr,jr);voc=brentq(f,0,vr[-1]);opt=minimize_scalar(lambda v:v*f(v),bounds=(0,voc),method='bounded');jsc=-jr[0];power=-opt.fun;forward.append(dict(model=model,Nt_cm3=Nt,Jsc_mAcm2=jsc*1000,Voc_V=voc,FF=power/(voc*jsc),Vmpp_V=float(opt.x),Jmpp_mAcm2=float(-f(opt.x)*1000),Pmpp_mWcm2=power*1000))
 csvout('selfconsistent_curves.csv',curves);csvout('frozen_source_reference.csv',frozen);csvout('shape_landmarks.csv',shape_rows);csvout('current_crossings.csv',thresholds);csvout('forward_metrics.csv',forward);save(ROOT/'data/analysis.json',dict(runs=out,frozen_checks=check,thresholds=thresholds,forward=forward));print(json.dumps(dict(thresholds=thresholds,forward=forward,normalization=[{k:v for k,v in r.items() if k!='shapes'} for r in out]),indent=2))
if __name__=='__main__':main()
