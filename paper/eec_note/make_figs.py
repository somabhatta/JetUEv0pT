#!/usr/bin/env python3
"""Schematic figures for the note on the mathematical relation between
energy-energy correlators and v0(pT)^jet-UE. No data and no model input.
Run: python3 make_figs.py   (writes figs/*.pdf)
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(OUT, exist_ok=True)

BLUE, BLUE2 = "#0F4D92", "#3775BA"
RED, RED1, RED2 = "#B64342", "#F6CFCB", "#E9A6A1"
GREEN, GREEN1 = "#8BCF8B", "#DDF3DE"
GREY, GREY_D = "#CFCECE", "#4D4D4D"
TEAL, VIOLET, GOLD = "#42949E", "#9A4D8E", "#E0B000"

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "sans-serif"],
    "font.size": 12, "axes.linewidth": 1.2,
    "axes.spines.right": False, "axes.spines.top": False,
    "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 1
def fig_geometry():
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.5), gridspec_kw={"width_ratios": [1.0, 2.9]})
    # (a) in-jet EEC: pairs inside one jet, pair angle R_L
    ax = axs[0]
    ax.set_aspect("equal")
    ax.add_patch(Circle((0, 0), 0.4, fc="#E8EFF8", ec=BLUE, lw=2))
    rng = np.random.default_rng(3)
    r = 0.38 * np.sqrt(rng.random(14)); a = 2 * np.pi * rng.random(14)
    px, py = r * np.cos(a), r * np.sin(a)
    sizes = 20 + 160 * rng.random(14) ** 3
    ax.scatter(px, py, s=sizes, color=BLUE, zorder=3)
    i, j = np.argmax(sizes), np.argsort(sizes)[-3]
    ax.annotate("", xy=(px[j], py[j]), xytext=(px[i], py[i]),
                arrowprops=dict(arrowstyle="<->", color=RED, lw=2))
    ax.text(0.03, 0.42, r"$R_L=\theta_{12}$", color=RED, fontsize=11)
    ax.text(0, -0.52, "jet cone, radius $R$", ha="center", color=BLUE)
    ax.set_xlim(-0.6, 0.6); ax.set_ylim(-0.62, 0.58)
    ax.set_xlabel(r"$\Delta\varphi$ to jet axis"); ax.set_ylabel(r"$\Delta\eta$")
    ax.set_title("(a) in-jet EEC:\nweight $E_1E_2$, angle $R_L$", fontsize=10.5, loc="left")

    # (b) jet-UE: x1 in the cone, x2 in the reference region (eta-phi plane)
    ax = axs[1]
    ax.add_patch(Rectangle((-np.pi / 2, -0.9), 2 * np.pi, 1.8, fc=GREEN1, ec="none", zorder=0))
    R = 0.4
    jet = (0.0, 0.1)
    for dphi, lab in [(np.pi / 2, r"perp. cone ($+\pi/2$)"), (-np.pi / 2, None)]:
        ax.add_patch(Circle((jet[0] + dphi, jet[1]), R, fc="white", ec=GREY_D, ls="--", lw=1.4, zorder=1))
    ax.add_patch(Circle((jet[0] + 3 * np.pi / 2, jet[1]), R, fc="white", ec=GREY_D, ls="--", lw=1.4, zorder=1))
    ax.text(np.pi / 2, jet[1] + 0.48, "perp. cone (excluded)", ha="center", fontsize=9, color=GREY_D)
    ax.text(-np.pi / 2 + 0.05, jet[1] + 0.48, "perp. cone", ha="left", fontsize=9, color=GREY_D)
    ax.add_patch(Circle(jet, R, fc="#E8EFF8", ec=BLUE, lw=2, zorder=2))
    ax.text(jet[0], jet[1] + 0.48, "jet cone $R$: $f(p_T)$", ha="center", fontsize=9, color=BLUE)
    ax.text(2.35, 0.74, r"reference region: $0.2<p_T<3$ GeV/$c$, outside jet and perp. cones $\to[p_T]_{\rm ref}$",
            fontsize=8.5, color="#2E7D32", ha="center")
    x1 = (jet[0] + 0.15, jet[1] + 0.12); x2 = (0.95, -0.62)
    ax.scatter(*x1, s=70, color=BLUE, zorder=4); ax.scatter(*x2, s=40, color=TEAL, zorder=4)
    ax.annotate("", xy=x2, xytext=x1, arrowprops=dict(arrowstyle="->", color=RED, lw=2), zorder=4)
    ax.text(0.62, -0.12, r"$\theta(x_2)>R$", color=RED, fontsize=12)
    ax.text(x1[0] - 0.1, x1[1] + 0.13, "$x_1$", color=BLUE, fontsize=10)
    ax.text(x2[0] + 0.07, x2[1] - 0.08, "$x_2$", color=TEAL, fontsize=12)
    ax.set_xlim(-np.pi / 2, 3 * np.pi / 2); ax.set_ylim(-0.9, 0.9); ax.set_aspect("equal")
    ax.set_xticks([-np.pi / 2, 0, np.pi / 2, np.pi, 3 * np.pi / 2])
    ax.set_xticklabels([r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$"])
    ax.set_xlabel(r"$\Delta\varphi$ to jet axis"); ax.set_ylabel(r"$\eta$")
    ax.set_title(r"(b) jet–UE pairs: $x_1$ in the cone, $x_2$ in the reference region" "\n"
                 r"($v_0^{\rm jet\text{-}UE}$ and jet–UE EEC; angle measured from the jet axis)", fontsize=10.5, loc="left")
    fig.tight_layout(pad=0.8)
    save(fig, "fig1_geometry")


# ---------------------------------------------------------------- Fig. 2
def fig_weights():
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.7))
    # (a) jet side, EEC: energy weight
    ax = axs[0]
    pt = np.array([0.35, 0.75, 1.5, 2.5, 4.5, 9, 21, 65])
    lab = ["0.2-0.5", "0.5-1", "1-2", "2-3", "3-6", "6-12", "12-30", "30-100"]
    x = np.arange(len(pt))
    w = pt / 50.0
    ax.bar(x, w, color=BLUE2, edgecolor="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, fontsize=8.5)
    ax.set_xlabel(r"in-cone $p_{T,1}$ (GeV/$c$)"); ax.set_ylabel(r"weight $p_{T,1}/p_T^{\rm jet}$")
    ax.set_title("(a) EEC, jet side: energy weight\n(dominated by hard constituents)", fontsize=10.5, loc="left")

    # (b) jet side, v0: shape projector; two in-cone spectra and their difference
    ax = axs[1]
    f0 = np.array([0.20, 0.22, 0.20, 0.11, 0.12, 0.08, 0.05, 0.02]); f0 /= f0.sum()
    f1 = f0 * np.array([1.18, 1.15, 1.08, 1.0, 0.92, 0.8, 0.7, 0.6]); f1 /= f1.sum()
    d = f1 - f0
    ax.bar(x, d, color=[GREEN if v > 0 else RED2 for v in d], edgecolor="black", lw=0.8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, fontsize=8.5)
    ax.set_xlabel(r"in-cone $p_{T,1}$ (GeV/$c$)"); ax.set_ylabel(r"$\delta f(p_T)$ (schematic)")
    ax.text(0.97, 0.93, r"$\sum_{p_T}\delta f=0$", transform=ax.transAxes, ha="right", fontsize=11)
    ax.set_title("(b) $v_0$, jet side: projector $\\chi_{p_T}-\\langle f(p_T)\\rangle$\n(only the shape of $f(p_T)$ enters)", fontsize=10.5, loc="left")

    # (c) reference side: weights vs p_T2 with schematic spectra
    ax = axs[2]
    p = np.linspace(0.2, 3.0, 300)
    mu = 0.70
    ax.plot(p, p, color=BLUE, lw=2.2, label=r"EEC: $w_E\propto p_{T,2}$")
    ax.plot(p, p - mu, color=RED, lw=2.2, label=r"$v_0$: $w_R=p_{T,2}-\mu$")
    ax.axhline(0, color="black", lw=0.8); ax.axvline(mu, color=GREY_D, ls=":", lw=1.2)
    ax.text(mu + 0.05, -0.55, r"$\mu=\langle[p_T]_{\rm ref}\rangle$", fontsize=9.5, color=GREY_D)
    ax.fill_between(p, p - mu, 0, where=p < mu, color=RED1, alpha=0.8)
    ax.fill_between(p, p - mu, 0, where=p > mu, color=GREEN1, alpha=0.9)
    ax.text(0.22, -0.8, r"$w_R<0$", color=RED, fontsize=10)
    ax.text(2.3, 1.1, r"$w_R>0$", color="#2E7D32", fontsize=10)
    ax.set_ylim(-1.05, 3.0); ax.set_xlim(0.2, 3.0)
    ax.set_xlabel(r"reference $p_{T,2}$ (GeV/$c$)"); ax.set_ylabel("weight")
    ax.legend(loc="upper left", frameon=False, fontsize=9.5)
    ax.set_title("(c) reference side: $w_E\\geq0$ vs. signed $w_R$", fontsize=10.5, loc="left")
    fig.tight_layout(pad=0.8)
    save(fig, "fig2_weights")


# ---------------------------------------------------------------- Fig. 3
def fig_angle():
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.6))
    th = np.linspace(0.0, 1.6, 400)
    G = 1.6 * np.exp(-th / 0.12) + 0.16 * np.exp(-((th - 0.55) / 0.4) ** 2) + 0.03
    ax = axs[0]
    ax.plot(th, G, color=BLUE, lw=2)
    for t0 in [0.1, 0.3, 0.6, 1.0]:
        g0 = np.interp(t0, th, G)
        ax.plot([t0, t0], [0, g0], color=BLUE2, lw=1, ls=":"); ax.scatter([t0], [g0], color=BLUE2, zorder=3, s=25)
    ax.text(0.65, 0.85, "EEC: one value per angle\n" r"$\Sigma^c_{\rm JU}(\theta)\propto\int dp\,p\,\Gamma_1(p,\theta)$", color=BLUE, fontsize=10)
    ax.axvline(0.4, color="black", lw=1, ls=":"); ax.text(0.41, 1.12, "$R$", fontsize=11)
    ax.set_ylim(0, 1.2); ax.set_xlim(0, 1.6)
    ax.set_xlabel(r"angle to jet axis $\theta$"); ax.set_ylabel(r"$\Gamma_k(p,\theta)$ (schematic)")
    ax.set_title(r"(a) EEC: $\Gamma$ resolved in angle", fontsize=10.5, loc="left")

    ax = axs[1]
    ax.plot(th, G, color=RED, lw=2)
    R1, R2 = 0.2, 0.4
    ax.fill_between(th, 0, G, where=th > R2, color=GREEN1, label=r"$v_0$ at radius $R_2$: $\int_{\theta>R_2}d\theta$")
    ax.fill_between(th, 0, G, where=(th > R1) & (th <= R2), color=GOLD, alpha=0.5,
                    label=r"extra region at $R_1$: $R_1<\theta<R_2$")
    for r_, l_ in [(R1, "$R_1$"), (R2, "$R_2$")]:
        ax.axvline(r_, color="black", lw=1, ls=":"); ax.text(r_ + 0.01, 1.12, l_, fontsize=11)
    ax.set_ylim(0, 1.2); ax.set_xlim(0, 1.6)
    ax.set_xlabel(r"angle to jet axis $\theta$"); ax.set_ylabel(r"$\Gamma_k(p,\theta)$ (schematic)")
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    ax.set_title(r"(b) $v_0$: $\Gamma$ integrated over the reference region $\theta>R$", fontsize=10.5, loc="left")
    fig.tight_layout(pad=0.8)
    save(fig, "fig3_angle")


# ---------------------------------------------------------------- Fig. 4
def fig_flow():
    fig, ax = plt.subplots(figsize=(10.0, 3.7))
    ax.set_xlim(0, 12.5); ax.set_ylim(0, 4.6); ax.axis("off")

    def box(x, y, w, h, txt, fc, ec, fs=11):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=1.6))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs)

    def arrow(a, b, col=GREY_D):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=16, color=col, lw=1.6))

    box(0.15, 1.6, 2.6, 1.4, "connected density\n" r"$C(x_1,x_2)=\rho_2-\rho_1\rho_1$", "#F3F3F3", GREY_D)
    box(3.6, 1.35, 3.5, 1.9, "two angular kernels ($k=0,1$)\n"
        r"$\Gamma_k(p,\theta)=\int_{\rm cone}\!dx_1\,\delta(p-p_{T,1})$" "\n"
        r"$\times\int_{\rm ref}\!dx_2\,p_{T,2}^{\,k}\,\delta(\theta-\theta_2)\,C$", "#E8EFF8", BLUE)
    box(8.6, 2.75, 3.75, 1.55, "jet–UE EEC (connected)\n"
        r"$\Sigma^c_{\rm JU}(\theta)=\frac{1}{(p_T^{\rm jet})^2}\int\!dp\,p\,\Gamma_1(p,\theta)$", GREEN1, GREEN)
    box(8.6, 0.3, 3.75, 1.75, r"numerator of $v_0^{\rm jet\text{-}UE}(p_T)$" "\n"
        r"$\frac{1}{\bar N\bar N_{\rm ref}}\int_{\theta>R}\!d\theta\,[\int_{\rm bin}-\langle f\rangle\int_{\rm all}]\,dp$" "\n"
        r"$\times\,[\Gamma_1-\mu\,\Gamma_0](p,\theta)$", RED1, RED)
    arrow((2.75, 2.3), (3.6, 2.3))
    arrow((7.1, 2.6), (8.6, 3.5)); arrow((7.1, 2.0), (8.6, 1.2))
    ax.text(7.0, 3.55, r"weight $p$, keep $\theta$", fontsize=9.5, color="#2E7D32")
    ax.text(7.0, 0.75, "signed, projected,\n" r"integrate $\theta$", fontsize=9.5, color=RED)
    fig.tight_layout(pad=0.3)
    save(fig, "fig4_flow")


if __name__ == "__main__":
    fig_geometry(); fig_weights(); fig_angle(); fig_flow()
    print("figures in", OUT)
