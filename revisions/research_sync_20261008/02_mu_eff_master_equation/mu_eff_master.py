#!/usr/bin/env python3
"""
mu_eff_master.py
高斯能量无序格点上的单载流子输运：用主方程（线性方程组）直接求
  (1) 平均首达时间 MFPT（平衡起点，用于与 Bässler 公式比较并标定 ν0）；
  (2) 带一阶损失 k_r 的收集概率 η（热载流子均匀产生，用于与 Jsc 比较）。
不使用 KMC，所以深陷阱下没有采样噪声；k_r>0 时矩阵严格对角占优，数值稳定。

模型（与 Bässler 的蒙特卡洛模型同类）：
  - 简单立方格点，晶格常数 a，x 方向 Lx = d/a 个格点，横向 Ly×Lz 周期；
  - 格点能量 E_i ~ N(0, σ²)，独立（无空间关联）；
  - 仅最近邻 Miller–Abrahams 跃迁：w_ij = ν0·exp(-max(ΔU,0)/kT)，
    ΔU = (E_j-E_i) - q·F·a·Δx，电场沿 +x 推动载流子；
  - x=0 面反射，x=d 面吸收（电极收集，接触能量取 0）；
  - 单载流子（低密度），不含库仑相互作用、复合以外的相互作用、极化子、变程跃迁。

用法：
  python3 mu_eff_master.py validate --out opv_mu      # 验证 Bässler 系数并标定 ν0
  python3 mu_eff_master.py sweep    --out opv_mu      # 75–300 K 收集概率
"""
import argparse
import csv
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

K_B = 8.617333262e-5   # eV/K
SEED0 = 12345
OFFS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

# 图 S12 的 PM6:Y6 参数（σ 为图中标注；μ0(300 K) 由图中直线读出，约 ±30%）
CARRIERS = {
    "hole":     dict(sigma=0.074, mu300=1.4e-4),
    "electron": dict(sigma=0.060, mu300=8.8e-4),
}


def make_energies(sigma, seed, shape):
    return np.random.default_rng(seed).normal(0.0, sigma, size=shape)


def assemble(E, F, T, nu0, a_cm):
    """返回 Woff（w_ij 稀疏矩阵）、每个格点的总出射率 Wsum、到吸收面的速率 wabs。"""
    N = E.size
    kT = K_B * T
    qFa = F * a_cm                      # 每步沿 +x 的场功 (eV)
    idx = np.arange(N).reshape(E.shape)
    rows, cols, vals = [], [], []
    Wsum = np.zeros(E.shape)
    wabs = np.zeros(E.shape)
    for dx, dy, dz in OFFS:
        Ej = np.roll(E, (-dx, -dy, -dz), axis=(0, 1, 2))
        jj = np.roll(idx, (-dx, -dy, -dz), axis=(0, 1, 2))
        dU = (Ej - E) - qFa * dx
        w = nu0 * np.exp(-np.maximum(dU, 0.0) / kT)
        mask = np.ones(E.shape, bool)
        if dx == 1:
            mask[-1] = False            # x=d 面的前向跃迁由吸收项处理
        if dx == -1:
            mask[0] = False             # x=0 面反射
        w = np.where(mask, w, 0.0)
        Wsum += w
        rows.append(idx[mask]); cols.append(jj[mask]); vals.append(w[mask])
    dUa = (0.0 - E[-1]) - qFa
    wa = nu0 * np.exp(-np.maximum(dUa, 0.0) / kT)
    wabs[-1] = wa
    Wsum[-1] += wa
    Woff = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(N, N))
    return Woff, Wsum.ravel(), wabs.ravel()


def collect_prob(Woff, Wsum, wabs, kr):
    """p_i：从格点 i 出发，在以速率 kr 损失之前到达 x=d 的概率。"""
    A = (sp.diags(Wsum + kr) - Woff).tocsc()
    return spsolve(A, wabs)


def mfpt(Woff, Wsum):
    """平均首达时间 t_i（无损失）。"""
    A = (sp.diags(Wsum) - Woff).tocsc()
    return spsolve(A, np.ones_like(Wsum))


def mu_from_t(t, F, T, Lx, a_cm):
    """由平均首达时间反推迁移率（均匀介质的漂移-扩散解，x=0 反射、x=L 吸收）：
       t = [L - (kT/qF)(1-exp(-qFL/kT))]/(μ F)。"""
    L = Lx * a_cm
    Vt = K_B * T
    return (L - Vt / F * (1.0 - np.exp(-F * L / Vt))) / (F * t)


# ---------------- 任务函数（供多进程调用） ----------------
def task_validate(a):
    name, sigma, T, r, F, nu0, Lx, Ly, Lz, a_nm = a
    a_cm = a_nm * 1e-7
    E = make_energies(sigma, SEED0 + r, (Lx, Ly, Lz))
    Woff, Wsum, wabs = assemble(E, F, T, nu0, a_cm)
    t = mfpt(Woff, Wsum).reshape(E.shape)
    w0 = np.exp(-E[0] / (K_B * T))      # 平衡起点：x=0 面上按玻尔兹曼权重
    return name, T, r, float((w0 * t[0]).sum() / w0.sum())


def task_sweep(a):
    name, sigma, T, r, F, nu0, krs, Lx, Ly, Lz, a_nm = a
    a_cm = a_nm * 1e-7
    E = make_energies(sigma, SEED0 + r, (Lx, Ly, Lz))
    Woff, Wsum, wabs = assemble(E, F, T, nu0, a_cm)
    out = []
    for kr in krs:
        p = collect_prob(Woff, Wsum, wabs, kr)
        out.append(float(p.mean()))     # 热载流子，沿膜厚均匀产生
    return name, T, F, r, out


def run_pool(fn, tasks, workers):
    if workers <= 1:
        return [fn(t) for t in tasks]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn, tasks, chunksize=1))


def cmd_validate(args):
    os.makedirs(args.out, exist_ok=True)
    Ts = [330, 300, 275, 250, 225]
    nu_base = 1e10
    a_cm = args.a * 1e-7
    res, calib = {}, {}
    rows = []
    for name, p in CARRIERS.items():
        tasks = [(name, p["sigma"], T, r, args.Fcal, nu_base, args.Lx, args.Ly, args.Ly, args.a)
                 for T in Ts for r in range(args.R)]
        t0 = time.time()
        out = run_pool(task_validate, tasks, args.workers)
        print(f"[validate] {name}: {len(tasks)} solves in {time.time()-t0:.1f}s")
        tm = {T: np.mean([o[3] for o in out if o[1] == T]) for T in Ts}
        mu = {T: mu_from_t(tm[T], args.Fcal, T, args.Lx, a_cm) for T in Ts}
        x = np.array([1.0 / T ** 2 for T in Ts])
        y = np.log([mu[T] for T in Ts])
        slope = np.polyfit(x, y, 1)[0]
        sig_fit = 1.5 * K_B * np.sqrt(-slope)             # slope = -(2σ/3k)^2
        # 用 300 K 标定 ν0（μ 与 ν0 成正比）
        nu0 = nu_base * p["mu300"] / mu[300]
        calib[name] = dict(sigma_eV=p["sigma"], sigma_fit_eV=float(sig_fit), nu0=float(nu0),
                           mu300_target=p["mu300"], a_nm=args.a, Lx=args.Lx)
        print(f"  σ_input={p['sigma']*1e3:.0f} meV, σ_fit(ln μ vs T^-2, 225–330 K)={sig_fit*1e3:.1f} meV, "
              f"ν0(calibrated)={nu0:.3e} s^-1")
        for T in Ts:
            rows.append([name, T, mu[T] * nu0 / nu_base])
    with open(os.path.join(args.out, "calibration.json"), "w") as f:
        json.dump(calib, f, indent=2)
    with open(os.path.join(args.out, "validate_mu.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["carrier", "T_K", "mu_cm2_Vs_calibrated"]); w.writerows(rows)


def cmd_sweep(args):
    calib = json.load(open(os.path.join(args.out, "calibration.json")))
    Ts = [300, 275, 250, 225, 200, 175, 150, 125, 100, 75]
    krs = [float(k) for k in args.kr.split(",")]
    Fs = [float(f) for f in args.F.split(",")]
    tasks = []
    for name, p in CARRIERS.items():
        for F in Fs:
            for T in Ts:
                for r in range(args.R):
                    tasks.append((name, p["sigma"], T, r, F, calib[name]["nu0"], krs,
                                  args.Lx, args.Ly, args.Ly, args.a))
    t0 = time.time()
    out = run_pool(task_sweep, tasks, args.workers)
    print(f"[sweep] {len(tasks)} tasks in {time.time()-t0:.1f}s")
    rows = []
    for name in CARRIERS:
        for F in Fs:
            for T in Ts:
                sel = [o[4] for o in out if o[0] == name and o[1] == T and o[2] == F]
                arr = np.array(sel)
                for i, kr in enumerate(krs):
                    rows.append([name, T, F, kr, arr[:, i].mean(), arr[:, i].std()])
    with open(os.path.join(args.out, "eta_sweep.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["carrier", "T_K", "F_V_cm", "kr_s", "eta_mean", "eta_std_over_realizations"])
        w.writerows(rows)
    print("T(K) " + " ".join(f"{n[:4]}@kr={k:.0e}" for n in CARRIERS for k in krs))
    for T in Ts:
        vals = []
        for name in CARRIERS:
            for kr in krs:
                m = [r for r in rows if r[0] == name and r[1] == T and r[2] == Fs[0] and r[3] == kr][0]
                vals.append(f"{m[4]:.3e}")
        print(f"{T:4d} " + " ".join(vals))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["validate", "sweep"])
    ap.add_argument("--out", default="opv_mu")
    ap.add_argument("--a", type=float, default=1.0, help="晶格常数 nm")
    ap.add_argument("--Lx", type=int, default=100, help="x 方向格点数 (d = Lx*a)")
    ap.add_argument("--Ly", type=int, default=10, help="横向周期格点数")
    ap.add_argument("--R", type=int, default=4, help="无序实现个数")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--Fcal", type=float, default=1e4, help="标定用低场 V/cm")
    ap.add_argument("--F", default="1e5", help="扫描用场 V/cm，逗号分隔")
    ap.add_argument("--kr", default="1e4,1e5", help="一阶损失速率 s^-1，逗号分隔")
    args = ap.parse_args()
    {"validate": cmd_validate, "sweep": cmd_sweep}[args.cmd](args)
