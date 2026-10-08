#!/usr/bin/env python3
"""
opv_reverse_bias.py
一维稳态漂移-扩散模型：有机太阳能电池/光电二极管的反偏行为。

机制（均可单独开关，见 Params）：
  1. 光生激子 -> CT 态 -> 自由载流子，解离概率 P(F) 用 Onsager-Braun 场依赖（或常数）
  2. 双分子复合（Langevin × 压低因子 zeta）
  3. 体内 SRH 陷阱辅助复合/热生成（反偏暗电流的来源之一）
  4. 电极注入：热发射 + 镜像力势垒降低（Schottky），接触密度随接触处局部电场变化
  5. 局域陷阱层：陷阱占据由稳态速率方程给出，陷阱电荷进入 Poisson 方程
     （可复现陷阱电荷 -> 接触场增强 -> 注入增益 的"光电倍增"正反馈）

约定：x=0 为阳极（高功函数），x=d 为阴极。外加电压 V 加在阳极，V<0 为反偏。
电流为传统电流，+x 方向为正；光电流在反偏下为负。
单位：cm, V, s, A, cm^-3。
数值方法：Scharfetter-Gummel 离散 + Gummel 迭代，非均匀网格（两端加密）。
"""
import argparse
import csv
import os
import time
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import ive

_trapz = getattr(np, 'trapezoid', None) or np.trapz
q = 1.602176634e-19
kB = 8.617333262e-5      # eV/K
eps0 = 8.8541878128e-14  # F/cm


@dataclass
class Trap:
    kind: str = "e"      # 'e' 电子陷阱（空时中性，占据时带-1）；'h' 空穴陷阱（占据时带+1）
    x1: float = 0.0      # 区域起点 cm
    x2: float = 0.0      # 区域终点 cm
    Nt: float = 1e17     # 陷阱密度 cm^-3
    cn: float = 1e-9     # 电子俘获系数 cm^3/s
    cp: float = 1e-9     # 空穴俘获系数 cm^3/s
    Et: float = 0.5      # 能级距 LUMO 的深度 eV


@dataclass
class Params:
    # 几何与网格
    d: float = 100e-7
    N: int = 241
    stretch: float = 2.5
    # 材料
    eps_r: float = 3.5
    T: float = 300.0
    Nc: float = 1e21
    Nv: float = 1e21
    chi: float = 3.9          # 电子亲和能 eV
    Eg: float = 1.4           # 有效带隙 eV
    WF_a: float = 5.1         # 阳极功函数 eV
    WF_c: float = 4.3         # 阴极功函数 eV
    mu_n: float = 1e-3        # cm^2/Vs
    mu_p: float = 1e-3
    pf_gamma: float = 0.0     # Poole-Frenkel: mu = mu0 exp(gamma sqrt(F)), F in V/cm
    # 光学
    Jmax: float = 20e-3       # A/cm^2，qGd：全部激子产生对应的电流上限
    alpha: float = 0.0        # 吸收系数 cm^-1，0 = 均匀产生
    # CT 态解离
    ct_field_dep: bool = True
    ct_a: float = 1.2e-7      # 电子-空穴间距 cm
    ct_kF: float = 1e6        # CT 态衰减速率 s^-1
    ct_P0: float = 0.9        # 不用场依赖时的常数解离概率
    # 复合与体陷阱
    zeta: float = 0.1         # Langevin 压低因子
    tau_srh: float = 1e-6     # SRH 寿命 s；0 表示关闭
    # 注入
    schottky: bool = True
    # 局域陷阱
    traps: list = field(default_factory=list)


def bern(x):
    x = np.clip(x, -600.0, 600.0)
    small = np.abs(x) < 1e-6
    xs = np.where(small, 1.0, x)
    return np.where(small, 1.0 - x / 2.0, xs / np.expm1(xs))


def ob_f(b):
    """Onsager-Braun 场因子 f(b) = I1(2*sqrt(2b)) / sqrt(2b)"""
    b = np.maximum(b, 1e-12)
    z = np.minimum(2.0 * np.sqrt(2.0 * b), 300.0)
    return ive(1, z) * np.exp(z) / (z / 2.0)


class Device:
    def __init__(self, P: Params):
        self.P = P
        N = P.N
        u = np.linspace(0.0, 1.0, N)
        b = P.stretch
        s = 0.5 * (1.0 + np.tanh(b * (2 * u - 1)) / np.tanh(b))
        self.x = P.d * s
        self.h = np.diff(self.x)
        hc = np.empty(N)
        hc[1:-1] = 0.5 * (self.h[1:] + self.h[:-1])
        hc[0] = self.h[0] / 2
        hc[-1] = self.h[-1] / 2
        self.hc = hc

        self.eps = P.eps_r * eps0
        self.Vt = kB * P.T
        self.IP = P.chi + P.Eg
        self.Vbi = P.WF_a - P.WF_c
        self.ni2 = P.Nc * P.Nv * np.exp(-P.Eg / self.Vt)
        self.ni = np.sqrt(self.ni2)
        self.gamL = q * (P.mu_n + P.mu_p) / self.eps   # Langevin 系数 cm^3/s

        self.phi = dict(
            n_a=max(P.WF_a - P.chi, 0.0), p_a=max(self.IP - P.WF_a, 0.0),
            n_c=max(P.WF_c - P.chi, 0.0), p_c=max(self.IP - P.WF_c, 0.0))

        # 激子产生剖面 (cm^-3 s^-1 @ 1 sun)
        w = np.exp(-P.alpha * self.x) if P.alpha > 0 else np.ones(N)
        self.Gexc = P.Jmax / q * w / _trapz(w, self.x)

        # CT 参数
        a = P.ct_a
        self.Eb = q / (4 * np.pi * self.eps * a)               # V (= eV)
        self.rc = q / (4 * np.pi * self.eps * self.Vt)         # Onsager 半径 cm
        self.kd0 = 3 * self.gamL / (4 * np.pi * a ** 3) * np.exp(-self.Eb / self.Vt)

        # 陷阱区域掩码
        self.masks = [((self.x >= t.x1) & (self.x <= t.x2)).astype(float) for t in P.traps]

    # ---------------- 物理子模块 ----------------
    def Pdiss(self, F):
        P = self.P
        if not P.ct_field_dep:
            return np.full_like(F, P.ct_P0)
        bb = self.rc * np.abs(F) / (2 * self.Vt)  # b = rc*F/(2Vt)
        kd = self.kd0 * ob_f(bb)
        return kd / (kd + P.ct_kF)

    def mobility(self, Fcell):
        P = self.P
        f = np.exp(np.minimum(P.pf_gamma * np.sqrt(np.abs(Fcell)), 30.0))
        return P.mu_n * f, P.mu_p * f

    def recomb(self, n, p):
        """净复合率 R(n,p)（cm^-3 s^-1）与固定(陷阱)电荷数密度 rho_num"""
        P = self.P
        np_ = n * p - self.ni2
        R = P.zeta * self.gamL * np_
        if P.tau_srh and P.tau_srh > 0:
            R = R + np_ / (P.tau_srh * (n + self.ni) + P.tau_srh * (p + self.ni))
        rho = np.zeros_like(n)
        for t, m in zip(P.traps, self.masks):
            en = t.cn * P.Nc * np.exp(-t.Et / self.Vt)
            ep = t.cp * P.Nv * np.exp(-(P.Eg - t.Et) / self.Vt)
            if t.kind == "e":
                f = (t.cn * n + ep) / (t.cn * n + en + t.cp * p + ep)
                Rt = t.cn * n * (1 - f) - en * f
                rho = rho - t.Nt * f * m
            else:
                g = (t.cp * p + en) / (t.cp * p + ep + t.cn * n + en)
                Rt = t.cp * p * (1 - g) - ep * g
                rho = rho + t.Nt * g * m
            R = R + t.Nt * Rt * m
        return R, rho

    def contacts(self, psi):
        """由接触处局部电场（镜像力势垒降低）给出四个接触密度"""
        P = self.P
        F0 = (psi[1] - psi[0]) / self.h[0]
        FL = (psi[-1] - psi[-2]) / self.h[-1]

        def lower(Finj):
            if not P.schottky:
                return 0.0
            return np.sqrt(q * max(Finj, 0.0) / (4 * np.pi * self.eps))

        ph, Vt = self.phi, self.Vt
        n0 = P.Nc * np.exp(-max(ph["n_a"] - lower(F0), 0.0) / Vt)
        p0 = P.Nv * np.exp(-max(ph["p_a"] - lower(-F0), 0.0) / Vt)
        nL = P.Nc * np.exp(-max(ph["n_c"] - lower(-FL), 0.0) / Vt)
        pL = P.Nv * np.exp(-max(ph["p_c"] - lower(FL), 0.0) / Vt)
        return n0, p0, nL, pL

    # ---------------- 求解 ----------------
    def solve(self, V, light=1.0, state=None, maxit=6000, tol=1e-8,
              w_c=0.5, verbose=False):
        P = self.P
        N, x, h, hc, Vt, eps = P.N, self.x, self.h, self.hc, self.Vt, self.eps

        if state is None:
            psi = V + (self.Vbi - V) * x / P.d
            n0, p0, nL, pL = self.contacts(psi)
            n = np.exp(np.interp(x, [0, P.d], [np.log(n0), np.log(nL)]))
            p = np.exp(np.interp(x, [0, P.d], [np.log(p0), np.log(pL)]))
        else:
            psi = state["psi"] + (V - state["V"]) * (1 - x / P.d)
            n, p = state["n"].copy(), state["p"].copy()
        psi[0], psi[-1] = V, self.Vbi

        Jprev, converged = np.inf, False
        for it in range(maxit):
            # 1) 接触密度（对数空间阻尼）
            tgt = self.contacts(psi)
            n[0] = np.exp((1 - w_c) * np.log(n[0]) + w_c * np.log(tgt[0]))
            p[0] = np.exp((1 - w_c) * np.log(p[0]) + w_c * np.log(tgt[1]))
            n[-1] = np.exp((1 - w_c) * np.log(n[-1]) + w_c * np.log(tgt[2]))
            p[-1] = np.exp((1 - w_c) * np.log(p[-1]) + w_c * np.log(tgt[3]))

            # 2) Poisson 修正（载流子 Boltzmann 线性化，陷阱电荷显式）
            R0, rho = self.recomb(n, p)
            dps = np.diff(psi) / h
            res = eps * (dps[1:] - dps[:-1]) / hc[1:-1] + q * (p[1:-1] - n[1:-1] + rho[1:-1])
            lo = eps / (h[:-1] * hc[1:-1])
            up = eps / (h[1:] * hc[1:-1])
            dg = -(lo + up) - q * (n[1:-1] + p[1:-1]) / Vt
            if P.traps:
                # 陷阱电荷对电势的响应（n~exp(+psi/Vt), p~exp(-psi/Vt)）也隐式放入 Poisson 线性化，抑制振荡
                et = 1e-3
                rho_t = self.recomb(n * np.exp(et), p * np.exp(-et))[1]
                dg = dg + q * ((rho_t - rho) / (et * Vt))[1:-1]
            ab = np.zeros((3, N - 2))
            ab[0, 1:] = up[:-1]
            ab[1, :] = dg
            ab[2, :-1] = lo[1:]
            dpsi = solve_banded((1, 1), ab, -res)
            dpsi = np.clip(dpsi, -2.0, 2.0)
            psi[1:-1] += dpsi
            dmax = np.max(np.abs(dpsi))

            # 3) 连续性方程
            dps = np.diff(psi)
            Dl = dps / Vt
            Fc = -dps / h
            Bp, Bm = bern(Dl), bern(-Dl)
            mun, mup = self.mobility(Fc)
            Dn, Dp = mun * Vt, mup * Vt
            Fnode = np.empty(N)
            Fnode[1:-1] = 0.5 * (np.abs(Fc[1:]) + np.abs(Fc[:-1]))
            Fnode[0], Fnode[-1] = abs(Fc[0]), abs(Fc[-1])
            G = light * self.Gexc * self.Pdiss(Fnode)

            hci = hc[1:-1]
            # --- 电子 ---
            R0, _ = self.recomb(n, p)
            e = 1e-4
            an = np.maximum((self.recomb(n * (1 + e), p)[0] - R0) / (n * e), 0.0)[1:-1]
            cu, cl = Dn[1:] / h[1:], Dn[:-1] / h[:-1]
            Cn = cu * Bp[1:] / hci
            Bn = (-cu * Bm[1:] - cl * Bp[:-1]) / hci - an
            An = cl * Bm[:-1] / hci
            rhs = R0[1:-1] - an * n[1:-1] - G[1:-1]
            rhs[0] -= An[0] * n[0]
            rhs[-1] -= Cn[-1] * n[-1]
            ab = np.zeros((3, N - 2))
            ab[0, 1:] = Cn[:-1]
            ab[1, :] = Bn
            ab[2, :-1] = An[1:]
            n[1:-1] = np.maximum(solve_banded((1, 1), ab, rhs), 1e-30)

            # --- 空穴 ---
            R0, _ = self.recomb(n, p)
            ap = np.maximum((self.recomb(n, p * (1 + e))[0] - R0) / (p * e), 0.0)[1:-1]
            cu, cl = Dp[1:] / h[1:], Dp[:-1] / h[:-1]
            Cp = -cu * Bm[1:] / hci
            Bpp = (cu * Bp[1:] + cl * Bm[:-1]) / hci + ap
            Ap = -cl * Bp[:-1] / hci
            rhs = G[1:-1] - R0[1:-1] + ap * p[1:-1]
            rhs[0] -= Ap[0] * p[0]
            rhs[-1] -= Cp[-1] * p[-1]
            ab = np.zeros((3, N - 2))
            ab[0, 1:] = Cp[:-1]
            ab[1, :] = Bpp
            ab[2, :-1] = Ap[1:]
            p[1:-1] = np.maximum(solve_banded((1, 1), ab, rhs), 1e-30)

            # 4) 电流与收敛判据
            Jn = q * Dn / h * (n[1:] * Bp - n[:-1] * Bm)
            Jp = q * Dp / h * (p[:-1] * Bp - p[1:] * Bm)
            Jc = Jn + Jp
            J = np.average(Jc, weights=h)
            # 绝对容差 1e-15 A/cm^2：低于此量级的电流（如无生成时的热注入）只是数值噪声
            if it >= 8 and dmax < tol and abs(J - Jprev) <= 1e-6 * abs(J) + 1e-15:
                converged = True
                break
            Jprev = J
            if verbose and it % 200 == 0:
                print(f"   it={it} dpsi={dmax:.2e} J={J:.4e}")

        R, rho = self.recomb(n, p)
        Fc = -np.diff(psi) / h
        return dict(
            V=V, J=J, Jn=np.average(Jn, weights=h), Jp=np.average(Jp, weights=h),
            spread=(Jc.max() - Jc.min()) / (abs(J) + 1e-30), it=it, converged=converged,
            x=x.copy(), psi=psi.copy(), n=n.copy(), p=p.copy(), rho=rho, G=G, Fc=Fc,
            state=dict(V=V, psi=psi.copy(), n=n.copy(), p=p.copy()))


def sweep(P, voltages, light, verbose=False):
    dev = Device(P)
    state, out = None, []
    for V in voltages:
        r = dev.solve(V, light=light, state=state, verbose=verbose)
        state = r["state"]
        out.append(r)
        if not r["converged"]:
            print(f"   [warn] V={V:+.2f} 未收敛 (it={r['it']}, spread={r['spread']:.1e})")
    return dev, out


# ---------------- 情景与作图 ----------------
def make_scenarios(base: Params):
    d = base.d
    trap = Trap(kind="e", x1=d - 10e-7, x2=d, Nt=1e18, cn=1e-9, cp=1e-9, Et=0.5)
    s0 = replace(base, ct_field_dep=False, schottky=False, tau_srh=0.0, traps=[])
    s1 = replace(s0, ct_field_dep=True)
    s2 = replace(s1, tau_srh=base.tau_srh)
    s3 = replace(s2, schottky=True)
    s4 = replace(s3, traps=[trap])
    return {
        "S0 const P, Langevin only": s0,
        "S1 + field-dep CT dissociation": s1,
        "S2 + bulk SRH generation": s2,
        "S3 + Schottky injection": s3,
        "S4 + cathode e-trap layer": s4,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vmax", type=float, default=12.0, help="最大反偏 |V|")
    ap.add_argument("--dv", type=float, default=0.5)
    ap.add_argument("--vstart", type=float, default=0.0)
    ap.add_argument("--out", default="opv_out")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--profile_v", type=float, default=-6.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    base = Params(N=121 if args.quick else 241)
    volts = np.arange(args.vstart, -args.vmax - 1e-9, -args.dv)
    scen = make_scenarios(base)

    rows, results = [], {}
    for name, P in scen.items():
        t0 = time.time()
        devd, dark = sweep(P, volts, 0.0, args.verbose)
        devl, lite = sweep(P, volts, 1.0, args.verbose)
        results[name] = (devl, dark, lite)
        for V, rd, rl in zip(volts, dark, lite):
            rows.append([name, V, rd["J"], rl["J"], (rl["J"] - rd["J"]) / (-P.Jmax)])
        print(f"{name}: {time.time() - t0:.1f}s")

    with open(os.path.join(args.out, "jv_scenarios.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "V", "J_dark_A_cm2", "J_light_A_cm2", "Jph_over_qGd"])
        w.writerows(rows)

    # 汇总表
    print("\nV(V)   " + "  ".join(f"{k.split()[0]:>10}" for k in scen))
    for i, V in enumerate(volts):
        if abs(V * 2 - round(V * 2)) < 1e-9 and (abs(V) % 2 < 1e-9):
            vals = []
            for name in scen:
                _, dk, lt = results[name]
                vals.append((lt[i]["J"] - dk[i]["J"]) / (-base.Jmax))
            print(f"{V:+5.1f}  " + "  ".join(f"{v:10.3f}" for v in vals))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(scen)
    cols = plt.cm.viridis(np.linspace(0.05, 0.9, len(names)))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for c, name in zip(cols, names):
        _, dk, lt = results[name]
        Jd = np.array([r["J"] for r in dk])
        Jl = np.array([r["J"] for r in lt])
        ax[0].semilogy(-volts, np.abs(Jd) + 1e-30, color=c, label=name)
        ax[1].plot(-volts, (Jl - Jd) / (-base.Jmax), color=c, label=name)
        ax[2].semilogy(-volts, np.abs(Jl) + 1e-30, color=c, label=name)
    ax[1].axhline(1.0, color="k", ls="--", lw=0.8)
    ax[1].text(0.2, 1.02, "qGd (optical limit)", fontsize=8)
    ax[0].set_title("(a) dark |J|")
    ax[1].set_title("(b) photocurrent / qGd  (>1 => injection gain)")
    ax[2].set_title("(c) light |J|")
    for a in ax:
        a.set_xlabel("reverse bias  -V (V)")
        a.grid(alpha=0.3)
    ax[0].set_ylabel("|J| (A/cm$^2$)")
    ax[1].set_ylim(0, None)
    ax[1].legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "fig1_jv_scenarios.png"), dpi=150)

    # 剖面图：S3 与 S4，暗/亮
    iV = int(np.argmin(np.abs(volts - args.profile_v)))
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    for row, name in enumerate([names[3], names[4]]):
        dev, dk, lt = results[name]
        P = dev.P
        for tag, r, ls in (("dark", dk[iV], "--"), ("light", lt[iV], "-")):
            x = r["x"] * 1e7
            Ec = -P.chi - r["psi"]
            Ev = -(P.chi + P.Eg) - r["psi"]
            axs[row, 0].plot(x, Ec, ls=ls, color="C0", label=f"Ec {tag}")
            axs[row, 0].plot(x, Ev, ls=ls, color="C3", label=f"Ev {tag}")
            axs[row, 1].semilogy(x, r["n"], ls=ls, color="C0", label=f"n {tag}")
            axs[row, 1].semilogy(x, r["p"], ls=ls, color="C3", label=f"p {tag}")
            xc = 0.5 * (x[1:] + x[:-1])
            axs[row, 2].plot(xc, r["Fc"] / 1e6, ls=ls, color="k", label=f"F {tag}")
        axs[row, 0].set_title(f"{name}\nbands @ V={volts[iV]:+.1f} V")
        axs[row, 1].set_title("carrier densities")
        axs[row, 1].set_ylim(1e-2, 1e22)
        axs[row, 2].set_title("field (MV/cm)")
        for a in axs[row]:
            a.set_xlabel("x (nm)")
            a.grid(alpha=0.3)
            a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "fig2_profiles.png"), dpi=150)
    print("done ->", args.out)


if __name__ == "__main__":
    main()
