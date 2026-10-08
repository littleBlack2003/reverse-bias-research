"""One explicit conditional inference: PIA n(T) fixes beta(T)/G(T).
No temperature-dependent contact or gap fits; no interpolation beyond data bounds.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.interpolate import PchipInterpolator
from fd_eos import GaussianDOS,KB
ROOT=Path(__file__).resolve().parents[1]
N0=2.4e20;BETA300=8e-12
DATA=ROOT/'reference'
def rows(name):return list(csv.DictReader((DATA/name).open()))
PIA=rows('fig5b_pia_noc_calibration.csv')
TP=np.array([float(r['T_K']) for r in PIA]);NP=np.array([float(r['n_oc_cm3']) for r in PIA])
N300=NP[-1];G0=BETA300*N300**2
_inter=PchipInterpolator(TP,np.log(NP),extrapolate=False)
def noc(T):
    result=np.exp(_inter(T))
    if not np.all(np.isfinite(result)):raise ValueError('PIA input requires100≤T≤300K')
    return result
def beta(T):return BETA300*(N300/noc(T))**2
def bulk_voc(T,I=1.,Eg=1.42):
    c=noc(T)*np.sqrt(I)/N0
    return Eg+KB*T*(GaussianDOS(.06,T).eta_from_c(c)+GaussianDOS(.074,T).eta_from_c(c))
def calibration_audit():
    voc=rows('fig5a_voc_calibration.csv');data=[]
    for r in voc:
        T=float(r['T_K'])
        if T<100:continue
        measured=float(r['Voc_V']);offset=float(bulk_voc(T,Eg=0))
        data.append(dict(T_K=T,n_pia_cm3=float(noc(T)),beta_inferred_cm3s=float(beta(T)),Voc_readout_V=measured,offset_V=offset,Eg_individual_diagnostic_eV=measured-offset))
    eg=np.mean([r['Eg_individual_diagnostic_eV'] for r in data])
    for r in data:
        r.update(Voc_fixed142_V=1.42+r['offset_V'],Voc_global_fit_V=eg+r['offset_V'],residual_global_fit_mV=1000*(eg+r['offset_V']-r['Voc_readout_V']))
    result=dict(material='PM6:Y6',N0_cm3=N0,beta300_cm3s=BETA300,n300_cm3=N300,G0_conditional_cm3s=G0,Jgeneration_110nm_mAcm2=1.602176634e-19*G0*110e-7*1000,Eg_single_global_fit_eV=float(eg),global_fit_rmse_mV=float(np.sqrt(np.mean([r['residual_global_fit_mV']**2 for r in data]))),rows=data,limitations=['PIA is an INPUT not validation','G(T)=G300 is an explicit unverified assumption','beta300 is printed approximate independent RT anchor; common sample/optical condition unproven','Eg142 chosen from published fitted bracket, so Fig5a is reproduction/consistency rather than independent prediction','Only S16 lower intensities and separately conditioned BACE data are holdouts','No temperature-by-temperature fitted Eg used'])
    (ROOT/'data/calibration_audit.json').write_text(json.dumps(result,indent=2))
    return result
if __name__=='__main__':print(json.dumps(calibration_audit(),indent=2))
