#!/usr/bin/env python3
"""plot_mu_eff.py : 读取 opv_mu/eta_sweep.csv、validate_mu.csv，作两联图。"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K_B = 8.617333262e-5
out = sys.argv[1] if len(sys.argv) > 1 else "opv_mu"
F0 = 1e5                 # V/cm
L = 100e-7               # cm
SIG = {"hole": 0.074, "electron": 0.060}
MU300 = {"hole": 1.4e-4, "electron": 8.8e-4}

# 2020 批次 Jsc（mA/cm2）
T_d = np.array([300, 280, 250, 230, 200, 175, 150, 125, 100, 75.])
J_d = np.array([22.80, 22.12, 21.46, 20.82, 19.60, 16.85, 12.83, 8.66, 5.50, 3.94])

rows = list(csv.DictReader(open(os.path.join(out, "eta_sweep.csv"))))
def eta(name, kr, F=F0):
    sel = sorted([(int(r["T_K"]), float(r["eta_mean"])) for r in rows
                  if r["carrier"] == name and float(r["kr_s"]) == kr and float(r["F_V_cm"]) == F])
    return np.array([s[0] for s in sel], float), np.array([s[1] for s in sel])

def mu_ext(name, T):                      # Bässler 零场外推（图中直线）
    f = lambda T: (2 * SIG[name] / (3 * K_B * T)) ** 2
    return MU300[name] * np.exp(-(f(T) - f(300.0)))

krs = [1e3, 1e4, 1e5]
ink, muted = "#0b0b0b", "#898781"
blue_ramp = {1e3: "#86b6ef", 1e4: "#3987e5", 1e5: "#184f95"}
c_h, c_e = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#c3c2b7",
                     "axes.linewidth": 0.8, "axes.facecolor": "#fcfcfb", "figure.facecolor": "#fcfcfb"})

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 5.2))

# ---- (a) 收集概率 vs T ----
for kr in krs:
    T, eh = eta("hole", kr); _, ee = eta("electron", kr)
    ax.semilogy(T, (eh + ee) / 2, "-", color=blue_ramp[kr], lw=2,
                label=f"master eq., 1/k$_r$ = {1e6/kr:.0f} µs")
Tt = np.linspace(75, 300, 200)
ext = 0.5 * (np.exp(-1e3 * L / (mu_ext("hole", Tt) * F0)) + np.exp(-1e3 * L / (mu_ext("electron", Tt) * F0)))
ax.semilogy(Tt, np.maximum(ext, 1e-4), "--", color=muted, lw=2,
            label="extrapolated Bässler $\\mu_0$ (1/k$_r$ = 1 ms)")
ax.semilogy(T_d, J_d / J_d[0], "o", color=ink, mfc="#fcfcfb", mew=1.5, ms=8, label="measured $J_{sc}$ / $J_{sc}$(300 K), 2020 batch")
ax.set_xlim(70, 305); ax.set_ylim(2e-3, 1.4)
ax.set_xlabel("Temperature (K)"); ax.set_ylabel("collected fraction  (a.u. of $J_{sc}$ proxy)")
ax.set_title("(a) Hot-carrier collection, d = 100 nm, F = 10$^5$ V/cm", loc="left", fontsize=11, color=ink)
ax.grid(alpha=0.25, lw=0.5); ax.legend(frameon=False, fontsize=8.5, loc="lower right")

# ---- (b) 有效迁移率 vs (1000/T)^2 ----
xg = np.linspace(8, 180, 200)
for name, c in (("hole", c_h), ("electron", c_e)):
    Tg = 1000.0 / np.sqrt(xg)
    bx.semilogy(xg, mu_ext(name, Tg), "--", color=c, lw=1.5, alpha=0.9)
    T, e4 = eta(name, 1e4)
    xs = (1000.0 / T) ** 2
    mu_eff = {kr: -kr * L / (F0 * np.log(eta(name, kr)[1])) for kr in krs}
    bx.fill_between(xs, mu_eff[1e5], mu_eff[1e3], color=c, alpha=0.15, lw=0)
    bx.semilogy(xs, mu_eff[1e4], "-o", color=c, lw=2, ms=6, label=f"{name}: from collection probability (band: 1/k$_r$ = 10 µs–1 ms)")
# 平衡验证点
val = list(csv.DictReader(open(os.path.join(out, "validate_mu.csv"))))
for name, c in (("hole", c_h), ("electron", c_e)):
    pts = [(float(r["T_K"]), float(r["mu_cm2_Vs_calibrated"])) for r in val if r["carrier"] == name]
    bx.semilogy([(1000 / t) ** 2 for t, _ in pts], [m for _, m in pts], "s", color=c, mfc="#fcfcfb", mew=1.5, ms=7)
bx.plot([], [], "--", color=muted, label="Bässler line, extrapolated (Fig. S12 $\\sigma$)")
bx.plot([], [], "s", color=muted, mfc="#fcfcfb", mew=1.5, label="solver, equilibrium start, 225–330 K")
bx.set_xlim(8, 182); bx.set_ylim(1e-18, 1e-2)
bx.set_xlabel("(1000/T)$^2$ (K$^{-2}$)   [75 K ≈ 178]"); bx.set_ylabel("effective mobility (cm$^2$/Vs)")
bx.set_title("(b) Effective mobility vs extrapolation", loc="left", fontsize=11, color=ink)
bx.grid(alpha=0.25, lw=0.5); bx.legend(frameon=False, fontsize=8, loc="lower left")

fig.suptitle("PM6:Y6 Gaussian disorder (hole 74 meV, electron 60 meV): master equation vs extrapolation", fontsize=11.5, color=ink, y=0.995)
fig.tight_layout()
fig.savefig(os.path.join(out, "fig_mu_eff.png"), dpi=150)
print("saved", os.path.join(out, "fig_mu_eff.png"))
