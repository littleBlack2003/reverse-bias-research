#!/usr/bin/env python3
"""NEW portable review wrapper. The adopted code/ files are unmodified.
No network access. No full temperature scan. Saved-state checks are not a fresh solve.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
import argparse
import hashlib
import json
from dataclasses import asdict
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from device import solve, advance
from extend_device import make_device
from warm_device import observable, RSH


def check_core_identity():
    source = json.loads((ROOT / 'SOURCE_PROVENANCE.json').read_text())
    for name, item in source['core_and_inputs'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != item['sha256']:
            raise AssertionError('Adopted source/input differs: ' + name)
    if RSH is not None:
        raise AssertionError('External shunt must be disabled')
    return len(source['core_and_inputs'])


def check_manifest():
    count = 0
    for line in (ROOT / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1)
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise AssertionError('Package hash mismatch: ' + name)
        count += 1
    return count


def load_case(T):
    if T not in (100, 300):
        raise ValueError('The compact review package contains only 100 K and 300 K')
    d = make_device(T, 321)
    expected = json.loads((ROOT / 'reference' / 'selected_state_reference.json').read_text())[str(T)]
    if asdict(d.p) != expected['parameters']:
        raise AssertionError('Frozen parameters differ')
    for key in ('gamma_n', 'gamma_p', 'field_cap'):
        if getattr(d, key) != expected[key]:
            raise AssertionError('Frozen field parameter differs: ' + key)
    with np.load(ROOT / 'reference' / f'critical_states_T{T}_N321.npz', allow_pickle=False) as f:
        states = {str(name): (float(label[0]), float(label[1]), z.copy())
                  for name, label, z in zip(f['names'], f['labels'], f['z'])}
        np.testing.assert_allclose(f['x_cm'], d.x*d.p.d_cm, rtol=0, atol=0)
    return d, states, expected


def audit_state(d, z, V, L, zeq):
    row = observable(d, z, V, L, zeq)
    jac = d.jacobian(z, V, L)
    scale = np.maximum(np.asarray(abs(jac).max(axis=1).toarray()).ravel(), 1e-250)
    norm = float(np.max(np.abs(d.residual(z, V, L))/scale))
    if not row['gate_passed'] or norm >= 3e-10:
        raise AssertionError(f'Saved-state accuracy gate failed: T={d.p.T}, V={V}, norm={norm}')
    if row['J_shunt_Acm2'] != 0 or row['external_shunt_heat_Wcm2'] != 0:
        raise AssertionError('Nonzero external shunt')
    if row['J_terminal_Acm2'] != row['J_intrinsic_Acm2']:
        raise AssertionError('Terminal/intrinsic current mismatch')
    return dict(T_K=d.p.T, V_V=V, light=L, J_Acm2=row['J_terminal_Acm2'],
                scaled_residual=norm, gate_passed=True,
                current_closure_uncertainty_Acm2=row['current_closure_uncertainty_Acm2'],
                charge_relative=row['charge_relative'], J_shunt_Acm2=0.0)


def smoke():
    rows = []
    for T in (100, 300):
        d, states, expected = load_case(T)
        zeq = states['equilibrium'][2]
        for name, (L, V, z) in states.items():
            row = audit_state(d, z, V, L, zeq)
            old = expected['rows'][expected['state_names'].index(name)]
            tolerance = max(1e-20, abs(old['J_terminal_Acm2'])*1e-10,
                            old['current_closure_uncertainty_Acm2']*2)
            if abs(row['J_Acm2'] - old['J_terminal_Acm2']) > tolerance:
                raise AssertionError('Reference current differs: ' + name)
            rows.append(dict(state=name, **row))
        if not (rows[-3]['J_Acm2'] < 0 < rows[-2]['J_Acm2']):
            raise AssertionError('Voc bracket lost its signed crossing')
    return dict(kind='saved-state residual re-evaluation, not a fresh PDE solve',
                unchanged_core_and_input_files=check_core_identity(),
                checked_package_files=check_manifest(), passed=True, states=rows)


def solve_one(T, name):
    d, states, _ = load_case(T)
    if name not in ('sc', 'mpp'):
        raise ValueError('Use sc or mpp for the bounded fresh-solve smoke test')
    L, V, z = states[name]
    guess = z.copy()
    perturb = 1e-5*np.sin(np.pi*d.x)
    guess[:, 0] += perturb
    guess[:, 1] += perturb
    guess[:, 2] -= perturb
    zz, info = solve(d, V, guess, L, maxiter=80)
    if not info['success']:
        raise RuntimeError(info.get('reason', 'Fresh solve failed'))
    row = audit_state(d, zz, V, L, states['equilibrium'][2])
    reference_current = float(d.ledger(z, V, L)['J_Acm2'])
    if abs(row['J_Acm2'] - reference_current) > max(1e-16, abs(reference_current)*1e-5):
        raise AssertionError('Fresh solution current disagrees with saved solution')
    return dict(kind='fresh nonlinear solve from a perturbed archived seed; not seed-free',
                state=name, Newton_iterations=info['iterations'],
                relative_current_change=(row['J_Acm2']/reference_current-1), passed=True, **row)


def metrics(T):
    """Optional narrow root/MPP recalculation, not a complete JV or temperature scan."""
    d, states, _ = load_case(T)
    zeq = states['equilibrium'][2]
    cache = {}
    for name, (L, V, z) in states.items():
        audit_state(d, z, V, L, zeq)
        if L == 1:
            cache[V] = (float(d.ledger(z, V, L)['J_Acm2']), z)
    solved = 0
    def evaluate(V):
        nonlocal solved
        V = float(V)
        if V not in cache:
            V0 = min(cache, key=lambda v: abs(v-V))
            zz, info = advance(d, cache[V0][1], V0, V, 1., 1.)
            row = audit_state(d, zz, V, 1., zeq)
            cache[V] = (row['J_Acm2'], zz)
            solved += 1
        return cache[V][0]
    lo, hi = states['bracket_low'][1], states['bracket_high'][1]
    voc = float(brentq(evaluate, lo, hi, xtol=2e-10))
    opt = minimize_scalar(lambda V: V*evaluate(V), bounds=(0., voc), method='bounded',
                          options={'xatol': 2e-8})
    if not opt.success:
        raise RuntimeError('MPP optimization did not converge')
    vmp = float(opt.x); jsc = evaluate(0.); jmp = evaluate(vmp)
    fresh = dict(T_K=T, N=321, Jsc_mAcm2=-1000*jsc, Voc_V=voc,
                 FF=vmp*jmp/(voc*jsc), Pmax_mWcm2=-1000*vmp*jmp,
                 Vmp_V=vmp, Jmp_mAcm2=-1000*jmp)
    expected = json.loads((ROOT/'reference'/'expected_metrics.json').read_text())[str(T)]
    return dict(kind='narrow root and MPP recalculation using archived seeds',
                newly_solved_bias_points=solved, metrics=fresh,
                difference_from_saved={k: fresh[k]-expected[k] for k in fresh if k not in ('T_K','N')},
                note='Solver tolerance is not mesh uncertainty or physical validation')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['smoke', 'solve', 'metrics', 'hashes'])
    ap.add_argument('--temperature', type=int, choices=[100, 300], default=300)
    ap.add_argument('--state', choices=['sc', 'mpp'], default='sc')
    a = ap.parse_args()
    check_core_identity()
    if a.command == 'smoke': result = smoke()
    elif a.command == 'solve': result = solve_one(a.temperature, a.state)
    elif a.command == 'metrics': result = metrics(a.temperature)
    else: result = dict(files=check_manifest(), passed=True)
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))


if __name__ == '__main__':
    main()
