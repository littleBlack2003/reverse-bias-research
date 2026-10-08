"""closure_v2: density-dependent mobility + optional Langevin-tied beta, as a SUBCLASS.
No original file is modified. Everything is a conditional diagnostic, not a validated model.

mu_i(c,F,T) = amp_i * mu_i,dilute(T) * g1_i(c_edge,T)**lam_d * exp(gamma_i*sq(F))
  g1: Pasveer PRL 94, 206601 (2005) density factor, sigma_hat = sigma/kT.
      Fit range of the source is 2..6; we EXTRAPOLATE (flagged, see `sigma_hat_max`).
      c > c_cut is saturated (same convention as the candidate code).
  Same positive factor multiplies drift AND FD generalized diffusion (as in FieldDevice).
beta_loc = beta0 * (beta_L/beta0)**lam_b,  beta_L = zeta*q*(mu_n_loc+mu_p_loc)/eps
  (lam_b=0: independent beta(T) of the original model; lam_b=1: Langevin-tied)
"""
import os
import numpy as np
from device import Device, KB, Q
from field_device import FieldDevice
from extend_device import parameters, ANCHOR


def g1(c, sigma, T, c_cut=0.1):
    s = sigma / (KB * T)
    delta = 2 * (np.log(s * s - s) - np.log(np.log(4))) / (s * s)
    return np.exp(0.5 * (s * s - s) * (2 * np.minimum(np.clip(c, 0, 1), c_cut)) ** delta)


class ClosureDevice(FieldDevice):
    def __init__(self, p, N=321, amp_n=1., amp_p=1., lam_d=0., lam_b=0., zeta=None, c_cut=None, dilute='bassler'):
        c_cut = float(os.environ.get('CLOSURE_CCUT', .1)) if c_cut is None else c_cut
        _, gn, gp = parameters(p.T)
        self.amp_n, self.amp_p, self.lam_d, self.lam_b, self.zeta, self.c_cut = amp_n, amp_p, lam_d, lam_b, zeta, c_cut
        super().__init__(p, N, gamma_n=gn, gamma_p=gp, field_cap=ANCHOR['field_cap'])
        self.mn *= amp_n
        self.mp *= amp_p

    def _factors(self, z, o):
        E = np.abs(np.diff(z[:, 0]) * self.vt / self.dx)
        Ec = np.minimum(E, self.field_cap); r = self.field_regularization
        sq = (Ec ** 2 + r ** 2) ** .25 - np.sqrt(r)
        ffn, ffp = np.exp(self.gamma_n * sq), np.exp(self.gamma_p * sq)
        cn = np.exp(self.en.logc((o['eta_n'][:-1] + o['eta_n'][1:]) / 2))
        cp = np.exp(self.ep.logc((o['eta_p'][:-1] + o['eta_p'][1:]) / 2))
        gdn = g1(cn, self.p.sigma_n_eV, self.p.T, self.c_cut) ** self.lam_d
        gdp = g1(cp, self.p.sigma_p_eV, self.p.T, self.c_cut) ** self.lam_d
        return E, ffn * gdn, ffp * gdp, cn, cp

    def mu_node(self, z, V):
        o = Device.evaluate(self, z, V)
        E, fn, fp, cn, cp = self._factors(z, o)
        node = lambda a: np.r_[a[0], .5 * (a[:-1] + a[1:]), a[-1]]
        return self.mn * node(fn), self.mp * node(fp)

    def evaluate(self, z, V):
        o = Device.evaluate(self, z, V)
        E, fn, fp, cn, cp = self._factors(z, o)
        o['Jn'] *= fn; o['Jp'] *= fp
        if self.lam_b > 0:
            node = lambda a: np.r_[a[0], .5 * (a[:-1] + a[1:]), a[-1]]
            mu = self.mn * node(fn) + self.mp * node(fp)
            bl = self.zeta * Q * mu / self.eps
            o['R'] = o['R'] * (bl / self.p.beta_cm3s) ** self.lam_b
        o.update(field_abs_Vcm=E, mobility_factor_n=fn, mobility_factor_p=fp, c_edge_n=cn, c_edge_p=cp)
        return o
