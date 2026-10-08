"""Source-motivated Murgatroyd/Gill local-field extension of fixed zero-field mobility.
mu(E)=mu0 exp(gamma sqrt(|E|)); gamma units (cm/V)^1/2.
Same positive factor multiplies drift AND FD generalized diffusion on each edge.
A declared field cap must bound extrapolation beyond the source extraction range.
No density dependence is added. Carrier zero-field mobilities remain source inputs.
"""
import numpy as np
from device import Device

class FieldDevice(Device):
    def __init__(self,p,nodes=321,mesh_strength=5.,gamma_n=0.,gamma_p=0.,field_cap=1e6,field_regularization=1.):
        if not np.isfinite([gamma_n,gamma_p,field_cap,field_regularization]).all() or min(gamma_n,gamma_p)<0 or min(field_cap,field_regularization)<=0:raise ValueError('Nonnegative finite gamma and positive finite field limits required')
        self.gamma_n=float(gamma_n);self.gamma_p=float(gamma_p);self.field_cap=float(field_cap);self.field_regularization=float(field_regularization)
        super().__init__(p,nodes,mesh_strength)
    def evaluate(self,z,V):
        o=super().evaluate(z,V)
        E=np.abs(np.diff(z[:,0])*self.vt/self.dx);Ec=np.minimum(E,self.field_cap)
        # Subtracted smoothing preserves mu(E=0)=mu0 exactly, while removing cusp.
        sq=(Ec**2+self.field_regularization**2)**.25-np.sqrt(self.field_regularization)
        fn=np.exp(self.gamma_n*sq);fp=np.exp(self.gamma_p*sq)
        o['Jn']*=fn;o['Jp']*=fp
        o.update(field_abs_Vcm=E,mobility_factor_n=fn,mobility_factor_p=fp)
        return o
    def ledger(self,z,V,light):
        r=super().ledger(z,V,light);o=self.evaluate(z,V)
        r.update(gamma_n_sqrtcmV=self.gamma_n,gamma_p_sqrtcmV=self.gamma_p,field_cap_Vcm=self.field_cap,field_regularization_Vcm=self.field_regularization,
                 field_min_Vcm=float(o['field_abs_Vcm'].min()),field_max_Vcm=float(o['field_abs_Vcm'].max()),
                 mobility_factor_n_min=float(o['mobility_factor_n'].min()),mobility_factor_n_max=float(o['mobility_factor_n'].max()),
                 mobility_factor_p_min=float(o['mobility_factor_p'].min()),mobility_factor_p_max=float(o['mobility_factor_p'].max()),
                 field_capped_length_fraction=float(np.sum(self.dx*(o['field_abs_Vcm']>self.field_cap))/self.p.d_cm),
                 field_law_status='source-motivated local-field extension; gamma provenance and field-cap sensitivity must be separately checked')
        return r
