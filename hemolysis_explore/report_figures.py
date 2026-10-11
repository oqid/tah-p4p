"""Curated, report-sized figures (16 cm wide) telling the 50 cc hemolysis story.

    python report_figures.py ..\\z_inputs                (50 cc, 120 bpm)
    python report_figures.py ..\\z_inputs --bpm 120 --outdir report_figs
    python report_figures.py ..\\z_inputs --only F2 F3        # rebuild just these
    python report_figures.py ..\\z_inputs --formats png svg    # other output types

Styling lives in p4p_report.mplstyle (fonts, sizes, grid) and report_style.py
(colours, markers, figure width). Each figure is one function below (f1_… f6_)
so titles, labels, limits and annotations are edited there.

Storyline (see FIGURE_GUIDE.md for captions and the sentence each one supports):
  F1  the model through systole: diaphragm dome + inlet speed, valve disc angle (to scale)
  F2  mechanism: HI3 falls with flow but tracks valve opening; stress term dominates
  F3  spatial: path stress (top) and where the dose is deposited (bottom)
  F4  timing: dose is deposited in the last tens of ms -> quasi-steady defence
  F5  path completion: paths stopping at the valve ring
  F6  35 cc hypothesis from similarity scaling
  F7  worst-case stress–duration envelopes vs approximate NIH lysis curves
PNG (300 dpi) + PDF are written for each.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

import explore_core as ec
import explore_hemolysis as eh
import operating_point as op
import empirical_envelope as ee
import report_style as rs
import scaling

HERE = Path(__file__).resolve().parent
rs.apply()                              # look & feel: p4p_report.mplstyle + report_style.py
W = rs.WIDTH
PC, PM = rs.PHASE_COLOURS, rs.PHASE_MARKERS
BLUE, PINK = rs.DEVICE_COLOURS["50 cc"], rs.DEVICE_COLOURS["35 cc"]
INK, MUTED = rs.INK, rs.MUTED
FORMATS = list(rs.FORMATS)


def save(fig, out, name):
    rs.save(fig, out, name, FORMATS)


def label_points(ax, x, y, labels, dx=4, dy=3, offsets=None):
    """Label markers; a white box keeps any line underneath from covering the text."""
    for k, (xi, yi, l) in enumerate(zip(x, y, labels)):
        ox, oy = offsets[k] if offsets is not None else (dx, dy)
        ax.annotate(l, (xi, yi), xytext=(ox, oy), textcoords="offset points", fontsize=6.5, color=INK,
                    ha="left" if ox >= 0 else "right", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9))


def smooth_path(ax, x, y, logy=True, n=40, **kw):
    """One smooth curve through the points in order (centripetal Catmull-Rom,
    built in axis-normalised coordinates so the shape suits the plot scales).
    A small arrowhead at the end shows the direction of travel."""
    x = np.asarray(x, float)
    yt = np.log10(y) if logy else np.asarray(y, float)
    sx, sy = np.ptp(x) or 1.0, np.ptp(yt) or 1.0
    P = np.c_[x / sx, yt / sy]
    P = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])   # phantom end points
    pts = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1:i + 3]
        t0 = 0.0
        t1 = t0 + np.linalg.norm(p1 - p0) ** 0.5
        t2 = t1 + np.linalg.norm(p2 - p1) ** 0.5
        t3 = t2 + np.linalg.norm(p3 - p2) ** 0.5
        t = np.linspace(t1, t2, n)[:, None]
        A1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1
        A2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2
        A3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3
        B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
        B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
        pts.append((t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2)
    C = np.vstack(pts)
    cx, cy = C[:, 0] * sx, C[:, 1] * sy
    cy = 10 ** cy if logy else cy
    style = dict(color=MUTED, lw=1.0, zorder=1)
    style.update(kw)
    ax.plot(cx, cy, **style)
    k = int(len(cx) - 0.5 * n)   # arrowhead half-way along the last segment
    ax.annotate("", (cx[k + 1], cy[k + 1]), (cx[k - 2], cy[k - 2]),
                arrowprops=dict(arrowstyle="-|>", color=style["color"], lw=0, mutation_scale=11), zorder=1)


def logfmt(ax, axis="y"):
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_minor_formatter(mt.NullFormatter())


# ---------------------------------------------------------------------------

def f1_model_states(case, opt, out):
    """Schematic of the model changing through systole, drawn to scale from the
    CAD parameters: diaphragm dome (inlet) and outlet-valve disc angle."""
    g = op.geometry(1.0)
    L, S = g["L"] * 1e3, g["S"] * 1e3            # mm
    R = L / 2
    D_valve = 18.0                               # orifice diameter [mm]
    n = len(opt)
    fig, axs = plt.subplots(2, n, figsize=(W, 2.55), constrained_layout=True,
                            gridspec_kw=dict(height_ratios=[1.0, 1.0]))
    r = np.linspace(-R, R, 200)
    vmax = opt.diaphragm_speed_m_s.max()
    for j, row in opt.iterrows():
        f = row.phase_pct / 100
        col = PC[j]
        # --- diaphragm (velocity inlet) ---
        ax = axs[0, j]
        for k, rr in opt.iterrows():
            ax.plot(r, S * rr.phase_pct / 100 * np.cos(np.pi * r / L), color="#D0D0D0", lw=0.6, zorder=1)
        ax.plot([-R, R], [0, 0], color=MUTED, lw=0.8, ls=":", zorder=1)
        y = S * f * np.cos(np.pi * r / L)
        ax.fill_between(r, 0, y, color=col, alpha=0.15, lw=0, zorder=2)
        ax.plot(r, y, color=col, lw=1.8, zorder=3)
        # inlet velocity arrows, normal to the surface, length ∝ speed
        alen = 9.0 * row.diaphragm_speed_m_s / vmax
        for x0 in np.linspace(-0.75 * R, 0.75 * R, 5):
            y0 = S * f * np.cos(np.pi * x0 / L)
            dy = -S * f * np.pi / L * np.sin(np.pi * x0 / L)
            nx, ny = -dy, 1.0
            nn = np.hypot(nx, ny)
            ax.annotate("", (x0 + alen * nx / nn, y0 + alen * ny / nn), (x0, y0),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=0.9, mutation_scale=6), zorder=4)
        ax.set_xlim(-R - 1, R + 1); ax.set_ylim(-9, S + 9)
        ax.set_title(f"{row.phase_pct:g}% systole", color=INK)
        ax.text(0, -1.5, f"h = {S * f:.1f} mm\nv = {row.diaphragm_speed_m_s * 1e3:.0f} mm/s",
                ha="center", va="top", fontsize=6.3, color=INK)
        # --- outlet valve (section through the orifice, flow left -> right) ---
        ax = axs[1, j]
        hw = D_valve / 2
        for sgn in (1, -1):
            ax.add_patch(plt.Rectangle((-2.5, sgn * hw if sgn > 0 else -hw - 2.5), 5, 2.5,
                                       color="#9E9E9E", lw=0))
        th = np.radians(row.valve_angle_deg)
        dl = 0.47 * D_valve
        ax.plot([-dl * np.sin(th), dl * np.sin(th)], [dl * np.cos(th), -dl * np.cos(th)],
                color=col, lw=3, solid_capstyle="round", zorder=3)
        ax.plot(0, 0, "o", color="white", mec=col, ms=3, zorder=4)
        tip = dl * np.cos(th)
        yy = 0.5 * (tip + hw)                      # centre of the gap between disc tip and housing
        for sgn in (1, -1):
            ax.annotate("", (12.5, sgn * yy), (-12.5, sgn * yy), arrowprops=dict(arrowstyle="-|>", color=MUTED,
                        lw=0.7, mutation_scale=6))
        ax.set_xlim(-14, 14); ax.set_ylim(-hw - 4, hw + 4)
        ax.text(0, -hw - 3.3, f"θ = {row.valve_angle_deg:g}°", ha="center", va="top", fontsize=6.8, color=INK)
    for a in axs.ravel():
        a.set_aspect("equal"); a.axis("off")
    axs[0, 0].text(-0.08, 0.5, "Diaphragm\n(velocity inlet)", transform=axs[0, 0].transAxes, rotation=90,
                   ha="right", va="center", fontsize=7, color=INK)
    axs[1, 0].text(-0.08, 0.5, "Outlet valve\n(18 mm orifice)", transform=axs[1, 0].transAxes, rotation=90,
                   ha="right", va="center", fontsize=7, color=INK)
    save(fig, out, "F1_model_through_systole")


def f2_mechanism(case, opt, out):
    t = case.table
    ph = t.phase_pct.to_numpy()
    fig, axs = plt.subplots(2, 2, figsize=(W, 4.9), constrained_layout=True)
    ax = axs[0, 0]
    y = t.GW_HI3_mean.to_numpy()
    err = np.vstack([y - t.GW_HI3_mean_ci_lo, t.GW_HI3_mean_ci_hi - y])
    ax.errorbar(ph, y, yerr=err, marker="o", color=BLUE, capsize=2.5, label="all exported paths")
    ax.plot(ph, t.GW_HI3_mean_completed_paths, marker="o", mfc="white", color=BLUE, lw=1,
            label="outlet-reaching paths only")
    ax.set_yscale("log"); logfmt(ax)
    ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel("Mean per-path HI3 [%] (GW)")
    ax.set_title("a) HI3 is U-shaped through systole")
    ax.legend(loc="upper center")
    ax.set_ylim(1.5e-5, 6e-4)
    ax.set_yticks([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], ["2e-5", "5e-5", "1e-4", "2e-4", "5e-4"])

    for ax, x, xl, ttl in ((axs[0, 1], opt.Q_ml_s.to_numpy(), "Instantaneous inflow Q [mL/s]",
                            "b) HI3 falls as flow rises…"),
                           (axs[1, 0], opt.valve_angle_deg.to_numpy(), "Outlet-valve opening angle θ [°]",
                            "c) …and tracks how far the valve is open")):
        smooth_path(ax, x, y)
        for i in range(len(ph)):
            ax.plot(x[i], y[i], marker=PM[i], color=PC[i], ms=6.5, ls="none")
        offs = [(6, -6), (6, -7), (6, 6), (6, 7), (-7, 4) if x.max() > 100 else (2, 10)]
        label_points(ax, x, y, [f"{p:g}%" for p in ph], offsets=offs)
        ax.text(0.96, 0.95, "line: phase order\n16.67 → 83%", transform=ax.transAxes,
                ha="right", va="top", fontsize=6.5, color=MUTED)
        ax.set_yscale("log"); logfmt(ax)
        ax.set_xlabel(xl); ax.set_ylabel("Mean HI3 [%] (GW)"); ax.set_title(ttl)
        ax.set_ylim(3e-5, 4e-4)
        yt = [3e-5, 1e-4, 3e-4]
        ax.set_yticks(yt, ["3e-5", "1e-4", "3e-4"])
        pad = 0.08 * (x.max() - x.min())
        ax.set_xlim(x.min() - 1.6 * pad, x.max() + 2.2 * pad)
    ax = axs[1, 1]
    ref_i = int(np.argmin(np.abs(ph - 50)))
    dec = [eh.log_decomposition(case.snaps[ref_i].paths, s.paths) for s in case.snaps]
    x = np.arange(len(ph))
    ax.bar(x - 0.19, [d["time"] for d in dec], 0.36, color=rs.TIME_TERM, label="transit-time term β·Δ⟨ln T⟩")
    ax.bar(x + 0.19, [d["stress"] for d in dec], 0.36, color=rs.STRESS_TERM, label="stress term α·Δ⟨ln τ_eff⟩")
    ax.plot(x, [d["total_geo"] for d in dec], "_", color="k", ms=14, mew=1.6, label="total")
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.set_xticks(x, [f"{p:g}" for p in ph]); ax.set_xlabel("Systolic phase [%]")
    ax.set_ylabel(f"log₁₀ change vs {ph[ref_i]:g}%")
    sec = ax.secondary_yaxis("right", functions=(lambda l: 10 ** l, np.log10))
    sec.set_ylabel("× factor"); sec.set_yticks([1, 2, 4, 8, 16])
    sec.yaxis.set_major_formatter(mt.FormatStrFormatter("%g"))
    ax.set_title("d) Stress, not residence time, drives it")
    ax.legend(loc="upper center", fontsize=6.3)
    save(fig, out, "F2_mechanism")


def _path_stress(ax, s, n=200, seed=0):
    rng = np.random.default_rng(seed)
    pick = rng.choice(s.n_paths, size=min(n, s.n_paths), replace=False)
    mask = np.isin(s.pid, pick)
    xyz = s.xyz[mask] * 1e3
    tau = np.maximum(s.tau[mask], 1e-3)
    o = np.argsort(tau)
    return ax.scatter(xyz[o, 0], xyz[o, 2], c=tau[o], s=0.6, cmap="magma", norm=LogNorm(1e-2, 1e2),
                      rasterized=True, linewidths=0)


def f3_spatial(case, opt, out, phases=(16.67, 50.0, 83.0)):
    snaps = [s for s in case.all_snaps if s.phase in phases]
    allxyz = np.vstack([s.xyz for s in case.all_snaps]) * 1e3
    lo, hi = allxyz.min(0) - 1, allxyz.max(0) + 1
    edges = [np.arange(lo[k], hi[k] + 1.0, 1.0) for k in range(3)]
    fig, axs = plt.subplots(2, len(snaps), figsize=(W, 3.75), constrained_layout=True, sharex=True, sharey=True)
    cmap = plt.get_cmap("magma_r").copy(); cmap.set_under("#EDEDED")
    for j, s in enumerate(snaps):
        r = opt[np.isclose(opt.phase_pct, s.phase)].iloc[0]
        ax = axs[0, j]
        m1 = _path_stress(ax, s)
        ax.set_title(f"{s.phase:g}% systole · valve {r.valve_angle_deg:g}°")
        ax = axs[1, j]
        xyz = s.xyz * 1e3
        mid = xyz.copy(); mid[1:] = 0.5 * (xyz[1:] + xyz[:-1])
        inner = s.inner("GW")
        H, _, _ = np.histogram2d(mid[:, 0], mid[:, 2], bins=[edges[0], edges[2]], weights=inner)
        H /= inner.sum()
        occ, _, _ = np.histogram2d(xyz[:, 0], xyz[:, 2], bins=[edges[0], edges[2]])
        H = np.where(occ > 0, np.maximum(H, 1e-9), np.nan)
        m2 = ax.pcolormesh(edges[0], edges[2], H.T, cmap=cmap, norm=LogNorm(1e-5, 0.1), rasterized=True)
        H3, _ = np.histogramdd(mid, bins=edges, weights=inner)
        occ3, _ = np.histogramdd(xyz, bins=edges)
        vals = np.sort(H3[occ3 > 0])[::-1] / inner.sum()
        k50 = np.searchsorted(np.cumsum(vals), 0.5) + 1
        ax.text(0.03, 0.03, f"50% of dose in\n{100 * k50 / len(vals):.2f}% of volume", transform=ax.transAxes,
                fontsize=6.3, color=INK, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        ax.set_xlabel("X [mm]")
    for a in axs.ravel():
        a.set_aspect("equal"); a.grid(False)
    axs[0, 0].set_ylabel("Z [mm]"); axs[1, 0].set_ylabel("Z [mm]")
    fig.colorbar(m1, ax=axs[0, :].tolist(), shrink=0.85, label="Local scalar stress [Pa]")
    fig.colorbar(m2, ax=axs[1, :].tolist(), shrink=0.85, extend="min", label="Share of linear HI3 dose\nper 1 mm column")
    save(fig, out, "F3_spatial_stress_and_dose")


def f4_timing(case, out):
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.5), constrained_layout=True)
    w = np.logspace(0, np.log10(3000), 120)
    ws = np.r_[1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000]
    for i, s in enumerate(case.all_snaps):
        s = s.subset(s.paths.completed.to_numpy())   # outlet-reaching paths only
        inner, tte = s.inner("GW"), s.time_to_exit() * 1e3
        o = np.argsort(tte); cum = np.cumsum(inner[o]) / inner.sum()
        rem = case.Ts * 1e3 * (1 - s.phase / 100)
        axs[0].plot(w, np.interp(w, tte[o], cum), color=PC[i], label=f"{s.phase:g}%")
        axs[0].plot(rem, np.interp(rem, tte[o], cum), marker=PM[i], color=PC[i], ms=5)
        base = s.paths.GW_HI3.mean()
        axs[1].plot(ws, [eh.truncated_hi3(s, x / 1e3).mean() / base for x in ws], color=PC[i])
        axs[1].plot(rem, eh.truncated_hi3(s, rem / 1e3).mean() / base, marker=PM[i], color=PC[i], ms=5)
    axs[0].set_title("a) Dose deposited before exit")
    axs[0].set_xlabel("Time before path exit [ms]"); axs[0].set_ylabel("Cumulative share of dose")
    axs[0].legend(title="phase", loc="upper left", fontsize=6.3, title_fontsize=6.5)
    axs[1].set_title("b) HI3 retained in a final window")
    axs[1].set_xlabel("Exposure window before exit [ms]"); axs[1].set_ylabel("HI3(window) / HI3(full)")
    for a in axs:
        a.set_xscale("log"); a.set_ylim(0, 1.03)
    save(fig, out, "F4_damage_timing")


def f5_completion(case, out):
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.6), constrained_layout=True,
                            gridspec_kw=dict(width_ratios=[1.15, 1]))
    fr, sh, ph = [], [], []
    for i, s in enumerate(case.all_snaps):
        p = s.paths; inc = ~p.completed.to_numpy()
        E = p[["end_x", "end_z"]].to_numpy() * 1e3
        axs[0].scatter(E[~inc, 0], E[~inc, 1], s=1, color="#C8C8C8", rasterized=True,
                       label="ends on outlet" if i == 0 else None)
        axs[0].scatter(E[inc, 0], E[inc, 1], s=5 + 60 * p.GW_HI3.to_numpy()[inc] / p.GW_HI3.max(),
                       marker=PM[i], facecolor="none", edgecolor=PC[i], lw=0.6, label=f"{s.phase:g}% stops early")
        fr.append(inc.mean() * 100); sh.append(p.GW_HI3[inc].sum() / p.GW_HI3.sum() * 100); ph.append(s.phase)
    axs[0].set_aspect("equal"); axs[0].set_xlabel("X [mm]"); axs[0].set_ylabel("Z [mm]")
    axs[0].set_title("a) Path end points (size ∝ HI3)")
    axs[0].legend(fontsize=5.8, loc="upper left", markerscale=1.0, handletextpad=0.2, labelspacing=0.3)
    axs[0].set_xlim(-38, 32); axs[0].set_ylim(-8, 30)
    x = np.arange(len(ph))
    axs[1].bar(x - 0.2, fr, 0.38, color="#BDBDBD", label="% of paths")
    axs[1].bar(x + 0.2, sh, 0.38, color=BLUE, label="% of summed HI3")
    for xi, (a, b) in enumerate(zip(fr, sh)):
        axs[1].text(xi - 0.2, a + 1, f"{a:.1f}", ha="center", fontsize=6)
        axs[1].text(xi + 0.2, b + 1, f"{b:.0f}", ha="center", fontsize=6)
    axs[1].set_xticks(x, [f"{p:g}" for p in ph]); axs[1].set_xlabel("Systolic phase [%]")
    axs[1].set_ylabel("%"); axs[1].set_title("b) Early-stopping paths vs their HI3")
    axs[1].legend(loc="upper left", ncol=2); axs[1].set_ylim(0, 62)
    save(fig, out, "F5_path_completion")


def f6_scaling(case, out):
    t = case.table
    ph = t.phase_pct.to_numpy()
    sc = scaling.scenarios(case.bpm)
    cols = dict(zip(sc, ("#5E8C31", "#B8860B")))
    regs = list(scaling.REGIMES)
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.6), constrained_layout=True)
    ax = axs[0]
    ax.plot(ph, t.GW_HI3_mean, marker="o", color=BLUE, label="50 cc CFD")
    for scen, r in sc.items():
        g = [r[f"HI3_ratio|GW|{k}"] for k in regs]
        ax.fill_between(ph, t.GW_HI3_mean * min(g), t.GW_HI3_mean * max(g), color=cols[scen], alpha=0.3, lw=0,
                        label=f"35 cc, {scen}")
    ax.set_yscale("log"); logfmt(ax)
    ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel("Mean HI3 [%] (GW)")
    ax.set_title("a) Predicted 35 cc band")
    ax.legend(loc="upper center", fontsize=6.3)
    ax = axs[1]
    short = ["visc.", "BL", "inert."]
    xt, xl = [], []
    for j, (scen, r) in enumerate(sc.items()):
        for gi, g in enumerate(regs):
            xp = j * 4 + gi; xt.append(xp); xl.append(short[gi])
            for k, mk in (("GW", "o"), ("HO", "s"), ("TZ", "^")):
                ax.plot(xp, r[f"HI3_ratio|{k}|{g}"], marker=mk, ms=5.5, ls="none", color=cols[scen],
                        mfc=cols[scen] if k == "GW" else "white", label=k if (j == 0 and gi == 0) else None)
        ax.text(j * 4 + 1, 0.43, scen, ha="center", fontsize=6.5, color=cols[scen], fontweight="bold")
    ax.axhline(1, color=MUTED, lw=0.8)
    ax.axhline(1.56, color="#7A5195", ls=":", lw=1)
    ax.text(6.6, 1.47, "Marom 2014 SA ratio", fontsize=5.8, color="#7A5195", ha="right", va="top")
    ax.set_yscale("log"); ax.set_ylim(0.4, 3); logfmt(ax)
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3], ["0.5", "0.75", "1", "1.5", "2", "3"])
    ax.set_xticks(xt, xl)
    ax.set_ylabel("Predicted HI3 ratio 35/50")
    ax.set_title("b) Ratio by stress scaling and model")
    leg = ax.legend(loc="upper left", fontsize=6.3)
    for h in leg.legend_handles:
        h.set_color(MUTED)
    save(fig, out, "F6_35cc_hypothesis")


def f7_empirical(case, out):
    """Worst-case stress–duration envelope per phase vs approximate NIH lysis traces
    (same definition as haem5_empiricalthresholds / haem5_snapshots)."""
    curves = ee.lysis_curves()
    LS = ["-", (0, (5, 1.5)), (0, (4, 1.5, 1, 1.5)), (0, (1, 1.2)), (0, (7, 2))]
    envs, rows = [], []
    for i, s in enumerate(case.all_snaps):
        env = ee.worst_envelope(s)
        sf = ee.safety_factors(env, curves)
        envs.append(env)
        rows.append({"phase_pct": s.phase, "min_safety_factor_RBC": sf.get("RBC"),
                     "min_safety_factor_platelet": sf.get("Platelet"),
                     "peak_stress_Pa": env.stress_Pa.max(), "longest_exposure_ms": env.duration_s.max() * 1e3})
        env.assign(phase_pct=s.phase).to_csv(out / f"F7_envelope_{s.phase:g}.csv", index=False)
    pd.DataFrame(rows).to_csv(out / "F7_lysis_safety_factors.csv", index=False)
    allenv = pd.concat(envs)
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.85), constrained_layout=True)
    for j, ax in enumerate(axs):
        for name, ls, lab in (("RBC", "--", "Red-cell lysis"), ("Platelet", ":", "Platelet lysis")):
            c = curves[name]
            ax.plot(c.time_s * 1e3, c.stress_Pa, color=INK, ls=ls, lw=1.2, label=lab if j == 1 else None)
        for i, (s, env) in enumerate(zip(case.all_snaps, envs)):
            e = env.sort_values(["duration_s", "stress_Pa"], ascending=[True, False])
            r = rows[i]
            ax.plot(e.duration_s * 1e3, e.stress_Pa, color=PC[i], ls=LS[i], lw=1.4,
                    marker=PM[i], markevery=0.3, ms=3.5,
                    label=(f"{s.phase:g}%  (SF {r['min_safety_factor_RBC']:.1f} RBC, "
                           f"{r['min_safety_factor_platelet']:.1f} plt)") if j == 1 else None)
        ax.set_xscale("log"); ax.set_yscale("log"); logfmt(ax); logfmt(ax, "x")
        ax.set_xlabel("Longest contiguous exposure [ms]")
    axs[0].set_ylabel("Shear stress [Pa]")
    axs[0].set_xlim(1e-3, 1e5); axs[0].set_ylim(min(1e-2, allenv.stress_Pa.min() / 2), 2e5)
    axs[0].set_title("a) Full lysis-reference range")
    # direct labels on the reference curves (panel a)
    axs[0].text(2.5e-3, 3.0e2, "red-cell lysis", fontsize=6.5, color=INK)        # below the dashed trace
    axs[0].text(1.5e3, 4.0e0, "platelet lysis", fontsize=6.5, color=INK)        # beside the dotted trace
    axs[1].set_xlim(allenv.duration_s.min() * 1e3 / 1.5, allenv.duration_s.max() * 1e3 * 1.5)
    axs[1].set_ylim(allenv.stress_Pa.min() / 1.5, 5e2)
    axs[1].set_title("b) Snapshot detail")
    axs[1].legend(loc="lower left", fontsize=5.6, handlelength=2.6, labelspacing=0.25)
    save(fig, out, "F7_empirical_stress_duration")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("--label", default="50 cc")
    ap.add_argument("--bpm", type=float, default=120.0)
    ap.add_argument("--outdir", type=Path, default=HERE / "report_figs")
    ap.add_argument("--only", nargs="+", metavar="F", help="rebuild only these, e.g. --only F2 F3")
    ap.add_argument("--formats", nargs="+", default=list(rs.FORMATS), help="e.g. --formats png svg")
    a = ap.parse_args(argv)
    FORMATS[:] = a.formats
    a.outdir.mkdir(parents=True, exist_ok=True)
    case = eh.Case(a.label, a.folder, a.bpm, 0.5, eh.CASE_COLOURS[0])
    opt = op.table([s.phase for s in case.all_snaps], bpm=a.bpm)
    opt["cfd_speed_p99_m_s"] = [np.percentile(s.vel, 99) for s in case.all_snaps]
    opt["GW_HI3_mean"] = case.table.GW_HI3_mean.to_numpy()
    opt["GW_HI3_mean_completed"] = case.table.GW_HI3_mean_completed_paths.to_numpy()
    opt.to_csv(a.outdir / "operating_points.csv", index=False)
    figures = {
        "F1": lambda: f1_model_states(case, opt, a.outdir),
        "F2": lambda: f2_mechanism(case, opt, a.outdir),
        "F3": lambda: f3_spatial(case, opt, a.outdir),
        "F4": lambda: f4_timing(case, a.outdir),
        "F5": lambda: f5_completion(case, a.outdir),
        "F6": lambda: f6_scaling(case, a.outdir),
        "F7": lambda: f7_empirical(case, a.outdir),
    }
    for key in (a.only or figures):
        print(f"  {key}")
        figures[key.upper()]()
    print(f"done -> {a.outdir}")


if __name__ == "__main__":
    main()
