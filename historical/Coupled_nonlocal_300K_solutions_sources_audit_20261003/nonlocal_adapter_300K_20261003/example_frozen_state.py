"""Numerical wiring example only: no self-consistent device or material fit."""
import numpy as np
from adapter import Mesh,Reaction,Adapter,FixedRates
from rates import from_positive,Q
m=Mesh(np.geomspace(1,11,42)*1e-6-1e-6)
k=FixedRates(from_positive(2.,7.,3.,11.))
a=Adapter(m,[Reaction(2.31e-6,4.62e-6,7.74e-6,3e7,k)])
s=np.r_[np.zeros(3*m.size),.37]
out=a.evaluate(s)
terms=a.device_terms(s)
scale=max(abs(Q*out['Sn']).max(),abs(Q*out['Sp']).max(),abs(out['drho']).max())
print('transient normalized local charge defect:',max(abs(a.continuity_defect(out)))/scale)
print('transient sparse reaction Jacobian:',a.jacobian(s).shape)
steady,j=a.condensed_steady(s[:3*m.size])
print('steady reduced vector/Jacobian:',steady['vector'].shape,j.shape)
print('endpoint distance cm:',a.reactions[0].xL-a.reactions[0].xH)
print('No Poisson/transport solve or J-V prediction performed.')
