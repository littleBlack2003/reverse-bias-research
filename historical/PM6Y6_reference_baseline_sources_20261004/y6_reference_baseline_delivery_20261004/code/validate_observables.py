from calibration import *
from device import mobility,Q,EPS0
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 audit=calibration_audit();eg=audit['Eg_single_global_fit_eV'];s16=[]
 for r in rows('figS16a_multilight_voc_holdout.csv'):
  if r['split']!='holdout':
   if 'holdout' not in r['split']:continue
  t=float(r['T_K']);I=float(r['illumination_sun_equiv']);meas=float(r['Voc_V'])
  pred=float(bulk_voc(t,I,eg));shift=pred-float(bulk_voc(t,1,eg))
  onesun=[float(x['Voc_V']) for x in rows('figS16a_multilight_voc_holdout.csv') if float(x['T_K'])==t and float(x['illumination_sun_equiv'])==1][0]
  s16.append(dict(T_K=t,I_sun=I,Voc_measured_V=meas,Voc_predicted_V=pred,residual_mV=1e3*(pred-meas),deltaVoc_predicted_V=shift,deltaVoc_measured_V=meas-onesun,delta_residual_mV=1e3*(shift-meas+onesun)))
 bace=[];br=rows('fig3a_bace_k2_validation.csv');b300=float(br[0]['k2_cm3_s'])
 for r in br:
  t=float(r['T_K']);meas=float(r['k2_cm3_s']);pred=float(beta(t))
  bace.append(dict(T_K=t,k2_bace_cm3s=meas,beta_PIA_conditional_cm3s=pred,bace_ratio_to300=meas/b300,pia_ratio_to300=pred/BETA300,normalized_ratio_mismatch=(pred/BETA300)/(meas/b300)))
 coef=[]
 for t in TP:
  mn=mobility(8.4e-4,.06,t);mp=mobility(1.3e-4,.074,t);bl=Q*(mn+mp)/(EPS0*3.5)
  coef.append(dict(T_K=float(t),mu_n=mn,mu_p=mp,beta=float(beta(t)),beta_Langevin_scenario=bl,beta_over_Langevin=float(beta(t))/bl))
 out=dict(calibration=audit,s16_holdout=s16,s16_rmse_mV=float(np.sqrt(np.mean([r['residual_mV']**2 for r in s16]))),s16_intensity_shift_rmse_mV=float(np.sqrt(np.mean([r['delta_residual_mV']**2 for r in s16]))),bace_independent_condition=bace,coefficients=coef,interpretation='S16 absolute predictions use one globally fitted Eg. Intensity shifts cancel Eg and same-temperature baseline offset. BACE absolute data are not asserted same device/optical condition. PIA n(T) was consumed as input; it cannot count as successful holdout.')
 (ROOT/'data/observable_validation.json').write_text(json.dumps(out,indent=2))
 fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
 ax=axs[0,0]
 ax.plot(TP,[r['Voc_readout_V'] for r in audit['rows']],'ko',label='Fig5a digitized (calibration)')
 ax.plot(TP,[r['Voc_global_fit_V'] for r in audit['rows']],label=f'Exact FD, one fitted Eg={eg:.4f} eV')
 ax.plot(TP,[r['Voc_fixed142_V'] for r in audit['rows']],'--',label='Published-bracket Eg=1.4200 eV')
 ax.set(xlabel='Temperature (K)',ylabel='Open-circuit voltage (V)',title='Bulk EOS consistency; PIA n is an input');ax.legend(fontsize=8)
 ax=axs[0,1]
 for I,col in zip([.5,.25,.08],['C0','C1','C2']):
  r=[r for r in s16 if r['I_sun']==I]
  ax.errorbar([q['T_K'] for q in r],[1e3*q['deltaVoc_measured_V'] for q in r],yerr=8.5,fmt='o',color=col,label=f'S16 {I:g} sun')
  ax.plot([q['T_K'] for q in r],[1e3*q['deltaVoc_predicted_V'] for q in r],color=col)
 ax.set(xlabel='Temperature (K)',ylabel='Voc(I)-Voc(1) (mV)',title='Independent intensity-shift holdout');ax.legend(fontsize=8)
 ax=axs[1,0]
 ax.semilogy(TP,[beta(t)/BETA300 for t in TP],label='PIA-derived, constant G')
 ax.errorbar([r['T_K'] for r in bace],[r['bace_ratio_to300'] for r in bace],fmt='s',label='BACE normalized (different condition)')
 ax.set(xlabel='Temperature (K)',ylabel='beta(T)/beta(300)',title='Independent recombination mismatch');ax.legend(fontsize=8)
 ax=axs[1,1]
 ax.semilogy(TP,[r['mu_n'] for r in coef],label='electron effective SCLC GDM')
 ax.semilogy(TP,[r['mu_p'] for r in coef],label='hole effective SCLC GDM')
 ax.axvspan(100,223,color='orange',alpha=.18,label='Unverified extrapolation')
 ax.set(xlabel='Temperature (K)',ylabel='Mobility (cm²/Vs)',title='Corrected-sign, separately normalized transport');ax.legend(fontsize=8)
 fig.savefig(ROOT/'figures/observable_validation.png',dpi=200);plt.close(fig)
 print(json.dumps({k:v for k,v in out.items() if k in ['s16_rmse_mV','s16_intensity_shift_rmse_mV']},indent=2))
if __name__=='__main__':main()
