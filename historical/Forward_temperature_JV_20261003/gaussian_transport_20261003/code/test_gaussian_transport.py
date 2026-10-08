import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from gaussian_transport import *
checks=[]
def check(name,value,limit):
    assert value<limit,(name,value,limit)
    checks.append(dict(name=name,error=float(value),limit=limit))
# EOS independent adaptive integration, s up to literature region endpoint.
for s in [0,2,4,6]:
    eta=np.array([-60.,-20.,-5.,0.,5.,20.]);c,dc,g=gaussian_eos(eta,s)
    cr=np.array([quad(lambda x:np.exp(-x*x/2)/np.sqrt(2*np.pi)*expit(y-s*x),-14,14,epsabs=1e-28,epsrel=2e-11)[0] for y in eta])
    dr=np.array([quad(lambda x:np.exp(-x*x/2)/np.sqrt(2*np.pi)*expit(y-s*x)*expit(-y+s*x),-14,14,epsabs=1e-28,epsrel=2e-11)[0] for y in eta])
    check('EOS quadrature s='+str(s),max(abs(c-cr)/np.maximum(cr,1e-30)),2e-5)
    check('EOS derivative s='+str(s),max(abs(dc-dr)/np.maximum(dr,1e-30)),2e-5)
    assert np.all(dc>0) and np.all(g>=1-1e-14)
    check('particle hole symmetry '+str(s),max(abs(gaussian_eos(-eta,s)[0]+c-1)),1e-14)
# Gaussian Boltzmann asymptote, zero disorder Fermi limit, derivative finite differences
for s in [0,2,4,6]:
    c,dc,g=gaussian_eos(-120,s)
    check('Boltzmann asymptote '+str(s),abs(c/np.exp(-120+s*s/2)-1),1e-10)
    check('Einstein dilute '+str(s),abs(g-1),1e-10)
c,dc,g=gaussian_eos(np.linspace(-10,10,21),0)
check('single level EOS',max(abs(c-expit(np.linspace(-10,10,21)))),1e-14)
# mobility exact low density T law and field symmetry, monotone thermal trend in benchmark
m=EGDM(.1,1e-7,1.)
T=np.array([300,340,380,420]);mu=np.array([m.evaluate(t,0,0)['mu'] for t in T])
check('Pasveer T law',abs(np.polyfit(1/T**2,np.log(mu),1)[0]/(-.42*(.1/KB)**2)-1),1e-12)
assert np.all(np.diff(mu)>0)
for T in [300,420]:
    a=m.evaluate(T,1e16,5e5);b=m.evaluate(T,1e16,-5e5)
    check('field reversal '+str(T),float(abs(a['mu']/b['mu']-1)),1e-14)
assert not m.evaluate(300,1e16,4e6)['field_plot_coverage']
try:m.evaluate(300,1e16,4e6,strict=True)
except ValueError:pass
else:raise AssertionError('strict domain guard failed')
# exact common-QF zero current despite bent bands; current sign and dilute SG recovery
psi=np.linspace(-3,3,31);dx=1e-7
for car in ['n','p']:
    eta=(-20+psi) if car=='n' else (-20-psi)
    j=inverse_activity_flux(eta,psi,300,1e21,.1,dx,1e-4,car,affinity=np.zeros(30))
    check('zero equilibrium '+car,max(abs(j)),1e-30)
    j=inverse_activity_flux(eta,psi,300,1e21,.1,dx,1e-4,car,affinity=np.ones(30)*.01)
    assert np.all(j>0) if car=='n' else np.all(j<0)
# Continuum homogeneous constitutive consistency: derivative flux approaches n*grad EF,
# giving D/mu=vt*c/dc rather than classical vt.
for s in [0,4,6]:
    h=1e-5;eta=np.array([-2-h/2,-2+h/2]);psi=np.zeros(2);c,dc,g=gaussian_eos(-2,s)
    j=inverse_activity_flux(eta,psi,300,1e21,s*KB*300,1e-7,1e-4,'n')
    expect=Q*1e-4*KB*300*1e21*c*h/1e-7
    check('continuum Einstein '+str(s),abs(j[0]/expect-1),1e-9)
import json,pathlib
pathlib.Path(__file__).parents[1].joinpath('data/tests.json').write_text(json.dumps(checks,indent=2))
print(len(checks),'numerical checks and additional domain/sign assertions passed')
