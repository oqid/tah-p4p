"""Exploratory hemolysis figures for the 50 cc (and later 35 cc) SynCardia snapshots.

Builds on haem5 (identical per-path numerics, verified in explore_core) but
keeps every SEGMENT, which allows analyses haem5 / haem5_paper_comparison
cannot do:

  01_phase_mechanism      HI3 vs phase + bootstrap CI, HI3 vs inlet speed,
                          log-decomposition of the phase trend into
                          transit-time and effective-stress contributions
  02_stress_time_map      per-path effective stress tau_eff vs transit T with
                          GW iso-HI3 lines (HI3 = A T^b tau_eff^a exactly) and
                          predicted 35 cc shifts
  03_damage_timing        when along each path the dose is deposited
                          (time-before-exit), and what that means for the
                          quasi-steady (frozen-flow) assumption
  04_damage_map           where the dose is deposited: 1 mm projected
                          dose-density maps + dose concentration curve,
                          hotspot voxel table for CFD-Post cross-checking
  05_numerical_audit      path completion, velocity floor, stress sampling
                          resolution in the high-dose segments
  06_scaling_hypothesis   similarity-scaling predictions for 35 cc at the
                          same BPM and at matched cardiac output
  07_size_comparison      (only when >= 2 cases) per-phase ratios with
                          bootstrap CI, ratio decomposition, ECDFs

Usage (from this folder, Windows PowerShell):
    python explore_hemolysis.py --case "50 cc" ..\\z_inputs 120
    python explore_hemolysis.py --case "50 cc" ..\\z_inputs 120 --case "35 cc" ..\\z_inputs_35cc 120
    python explore_hemolysis.py --case "50 cc" ..\\z_inputs 120 --case "35 cc (matched CO)" ..\\z_inputs_35cc_171 171

Parsed snapshots are cached in ./cache (keyed on file hash), so re-runs take
seconds. Requires numpy, pandas, matplotlib only.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

import explore_core as ec
import scaling

HERE = Path(__file__).resolve().parent

# --- Styling --------------------------------------------------------------
# Phases are ordinal: lightness-monotone ramp (bi/trans-flag inspired purple
# -> pink -> peach) plus distinct markers and direct labels, so colour is
# never the only cue.
PHASE_COLOURS = ["#2D1E5C", "#5B3A8E", "#9B4F96", "#D0668F", "#EE9A6B"]
PHASE_MARKERS = ["o", "s", "D", "^", "v"]
# Device identity: validated pair (CVD dE 8.5, normal dE 23.6)
CASE_COLOURS = ["#2E6DB4", "#C2477F", "#5E8C31", "#7A5195"]
INK, MUTED, GRID = "#1F1F1F", "#6B6B6B", "#E6E6E6"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 200, "font.size": 9.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
    "legend.frameon": False, "lines.linewidth": 2.0,
})

GW = ec.CONSTANTS["GW"]


def phase_style(i):
    return dict(color=PHASE_COLOURS[i % 5], marker=PHASE_MARKERS[i % 5])


def note(fig, text, bottom=0.07):
    """Footnote under the axes; reserves space in the constrained layout."""
    eng = fig.get_layout_engine()
    if eng is not None:
        eng.set(rect=(0, bottom, 1, 1 - bottom))
    fig.text(0.01, 0.005, text, fontsize=7.5, color=MUTED, ha="left", va="bottom", wrap=True)


def boot_mean_ci(x, n=2000, seed=0, q=(2.5, 97.5)):
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    m = x[idx].mean(axis=1)
    return np.percentile(m, q)


# ---------------------------------------------------------------------------
# Case loading & per-phase table
# ---------------------------------------------------------------------------

class Case:
    def __init__(self, label, folder, bpm, systolic_fraction, colour, completed_only=False):
        self.label, self.folder, self.bpm = label, Path(folder), float(bpm)
        self.colour = colour
        self.Ts = 60.0 / self.bpm * systolic_fraction       # systole [s]
        files = sorted(self.folder.glob("*.csv"), key=lambda p: (ec.phase_from_name(p) or 0, p.name))
        if not files:
            raise SystemExit(f"No CSV exports in {self.folder}")
        self.snaps = []
        for f in files:
            print(f"  loading {f.name} ...", flush=True)
            snap = ec.load_snapshot(f, HERE / "cache")
            self.all_snaps = getattr(self, "all_snaps", []) + [snap]
            self.snaps.append(snap.subset(snap.paths.completed.to_numpy()) if completed_only else snap)
        if any(s.phase is None for s in self.snaps):
            raise SystemExit("Could not read phase from file names (expect e.g. _1667%.csv)")
        self.table = self._table()

    def _table(self):
        rows = []
        for s, s_all in zip(self.snaps, self.all_snaps):
            p = s.paths
            pa = s_all.paths
            r = {"case": self.label, "bpm": self.bpm, "phase_pct": s.phase,
                 "systole_ms": self.Ts * 1e3,
                 "remaining_systole_ms": self.Ts * 1e3 * (1 - s.phase / 100),
                 "n_paths": len(p), "inlet_speed_mean_m_s": p.v_inlet.mean()}
            for k in ec.CONSTANTS:
                h = p[f"{k}_HI3"]
                r[f"{k}_HI3_mean"] = h.mean()
                lo, hi = boot_mean_ci(h.to_numpy())
                r[f"{k}_HI3_mean_ci_lo"], r[f"{k}_HI3_mean_ci_hi"] = lo, hi
                r[f"{k}_HI3_geomean"] = np.exp(np.log(h).mean())
                r[f"{k}_HI3_p50"], r[f"{k}_HI3_p95"], r[f"{k}_HI3_p99"] = np.percentile(h, [50, 95, 99])
                r[f"{k}_NIH_mg_100L_from_mean"] = h.mean() * (1 - ec.HCT) * ec.HB * 1000
            r["frac_paths_incomplete"] = float((~pa.completed).mean())
            r["GW_HI3_share_from_incomplete"] = float(pa.GW_HI3[~pa.completed].sum() / pa.GW_HI3.sum())
            r["GW_HI3_mean_all_paths"] = pa.GW_HI3.mean()
            r["GW_HI3_mean_completed_paths"] = pa.GW_HI3[pa.completed].mean()
            r["GW_HI3_mean_vweighted"] = np.average(p.GW_HI3, weights=p.v_inlet)
            r["transit_ms_p50"], r["transit_ms_p95"] = np.percentile(p.transit_ms, [50, 95])
            r["transit_ms_geomean"] = np.exp(np.log(p.transit_ms).mean())
            r["frac_transit_gt_remaining_systole"] = float((p.transit_ms > r["remaining_systole_ms"]).mean())
            r["frac_transit_gt_systole"] = float((p.transit_ms > self.Ts * 1e3).mean())
            r["GW_tau_eff_p50_Pa"] = p.GW_tau_eff.median()
            r["GW_tau_eff_geomean_Pa"] = np.exp(np.log(p.GW_tau_eff).mean())
            r["tau_peak_p50_Pa"], r["tau_peak_p99_Pa"] = np.percentile(p.tau_peak, [50, 99])
            r["SA_p50"], r["SA_p95"] = np.percentile(p.SA, [50, 95])
            inner, tte = s.inner("GW"), s.time_to_exit()
            for w in (10, 25, 50):
                r[f"dose_frac_last_{w}ms"] = inner[tte <= w / 1e3].sum() / inner.sum()
            # HI3 retained if only exposure inside the remaining-systole window counts
            r["GW_HI3_mean_within_remaining_systole"] = truncated_hi3(s, r["remaining_systole_ms"] / 1e3).mean()
            r["velocity_floor_time_frac"] = s.dt[s.floored].sum() / s.dt.sum()
            rows.append(r)
        t = pd.DataFrame(rows)
        # ejection-weighted mean over sampled phases (equal phase spacing assumed;
        # Q proxy = mean inlet (diaphragm) seed speed, constant inlet area)
        w = t.inlet_speed_mean_m_s / t.inlet_speed_mean_m_s.sum()
        for k in ec.CONSTANTS:
            t[f"{k}_HI3_ejection_weighted_over_phases"] = float((w * t[f"{k}_HI3_mean"]).sum())
        return t


def truncated_hi3(s, window_s, const="GW"):
    """Per-path HI3 counting only dose deposited within `window_s` before exit."""
    c = ec.CONSTANTS[const]
    inner = s.inner(const)
    keep = s.time_to_exit() <= window_s
    lin = np.bincount(s.pid, weights=np.where(keep, inner, 0.0), minlength=s.n_paths)
    return lin ** c["beta"]


def log_decomposition(p_a, p_b, const="GW"):
    """Exact split of ln(geomean HI3_b / geomean HI3_a) into time and stress terms.

    ln HI3 = ln A + beta ln T + alpha ln tau_eff for every path, so the change in
    the MEAN of ln HI3 (= ln of the geometric mean) splits exactly into
    beta * d<ln T> + alpha * d<ln tau_eff>. Returned in log10 units.
    """
    c = ec.CONSTANTS[const]
    dT = np.log(p_b.T_s).mean() - np.log(p_a.T_s).mean()
    dtau = np.log(p_b[f"{const}_tau_eff"]).mean() - np.log(p_a[f"{const}_tau_eff"]).mean()
    l10 = np.log(10)
    return {"time": c["beta"] * dT / l10, "stress": c["alpha"] * dtau / l10,
            "total_geo": (c["beta"] * dT + c["alpha"] * dtau) / l10,
            "total_arith": np.log10(p_b[f"{const}_HI3"].mean() / p_a[f"{const}_HI3"].mean())}


# ---------------------------------------------------------------------------
# 01 Phase mechanism
# ---------------------------------------------------------------------------

def fig_phase_mechanism(case: Case, out: Path):
    t = case.table
    ph = t.phase_pct.to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True)

    ax = axs[0]
    for k, ls, mk in (("GW", "-", "o"), ("HO", "--", "s"), ("TZ", ":", "^")):
        y = t[f"{k}_HI3_mean"].to_numpy()
        err = np.vstack([y - t[f"{k}_HI3_mean_ci_lo"], t[f"{k}_HI3_mean_ci_hi"] - y])
        ax.errorbar(ph, y, yerr=err, ls=ls, marker=mk, color=case.colour, ms=6,
                    capsize=3, lw=1.6, label=f"{k} mean (95% bootstrap CI)")
    ax.plot(ph, t.GW_HI3_mean_completed_paths, ls="-", marker="o", mfc="white", color=case.colour,
            lw=1.2, label="GW mean, outlet-reaching paths only")
    ax.plot(ph, t.GW_HI3_geomean, ls="-", marker="o", mfc="white", color=MUTED,
            lw=1.2, label="GW geometric mean (typical path)")
    ax.set_yscale("log")
    ax.set_xlabel("Systolic phase [%]")
    ax.set_ylabel("Per-path HI3 [%]")
    ax.set_title("a) Predicted damage per pass through systole")
    ew = t.GW_HI3_ejection_weighted_over_phases.iloc[0]
    ax.axhline(ew, color=case.colour, lw=0.8, ls=(0, (1, 2)),
               label=f"GW ejection-weighted over phases = {ew:.2e}%")
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_, hi_ * 12)
    ax.legend(fontsize=7.5, loc="upper center", ncol=2)

    ax = axs[1]
    v = t.inlet_speed_mean_m_s.to_numpy()
    y = t.GW_HI3_mean.to_numpy()
    for i in range(len(t)):
        ax.plot(v[i], y[i], ms=9, ls="none", **phase_style(i))
        ax.annotate(f"{ph[i]:g}%", (v[i], y[i]), xytext=(6, 2), textcoords="offset points",
                    fontsize=8, color=INK)
    sl, ic = np.polyfit(np.log(v), np.log(y), 1)
    vv = np.linspace(v.min() * 0.95, v.max() * 1.05, 50)
    ax.plot(vv, np.exp(ic) * vv ** sl, color=MUTED, lw=1, ls="--")
    ax.text(0.97, 0.95, f"log–log slope ≈ {sl:.1f}\n(descriptive, 5 points;\nvalve state also changes)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=MUTED)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Mean diaphragm-inlet speed [m/s]  (∝ instantaneous flow)")
    ax.set_ylabel("Mean GW HI3 [%]")
    ax.set_title("b) Lower instantaneous flow → higher HI3")

    ax = axs[2]
    ref_i = int(np.argmin(np.abs(ph - 50)))
    ref = case.snaps[ref_i].paths
    dec = [log_decomposition(ref, s.paths) for s in case.snaps]
    x = np.arange(len(ph))
    tm = np.array([d["time"] for d in dec]); st = np.array([d["stress"] for d in dec])
    ax.bar(x - 0.18, tm, 0.34, color="#8FB8DE", label="β·Δ⟨ln T⟩  transit-time term")
    ax.bar(x + 0.18, st, 0.34, color="#D98CB3", label="α·Δ⟨ln τ_eff⟩  stress term")
    ax.plot(x, [d["total_geo"] for d in dec], "k_", ms=22, mew=2, label="total (geometric mean)")
    ax.plot(x, [d["total_arith"] for d in dec], "D", color=case.colour, ms=6, ls="none",
            label="arithmetic-mean ratio")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xticks(x, [f"{p:g}%" for p in ph])
    ax.set_xlabel("Systolic phase")
    ax.set_ylabel(f"log₁₀ change vs {ph[ref_i]:g}% phase")
    sec = ax.secondary_yaxis("right", functions=(lambda l: 10 ** l, np.log10))
    sec.set_ylabel("× factor")
    ax.set_title("c) Why: stress dominates the phase trend (GW)")
    ax.legend(fontsize=7.5, loc="upper center")
    note(fig, f"{case.label}, {case.bpm:g} bpm, systole {case.Ts*1e3:.0f} ms. Each phase is a separate steady (frozen) CFD field. "
              "Equal-path statistics. Panel c uses the exact identity HI3 = A·T^β·τ_eff^α per path.")
    fig.savefig(out / "01_phase_mechanism.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 02 Stress-time map
# ---------------------------------------------------------------------------

def _iso_lines(ax, Tlim, const="GW", levels=(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)):
    c = ec.CONSTANTS[const]
    T = np.logspace(np.log10(Tlim[0]), np.log10(Tlim[1]), 50)
    for L in levels:
        tau = (L / (c["A"] * T ** c["beta"])) ** (1 / c["alpha"])
        ax.plot(T * 1e3, tau, color=MUTED, lw=0.7, ls=(0, (4, 3)), zorder=1)
        ax.annotate(f"{L:.0e}%", (T[-1] * 1e3, tau[-1]), fontsize=6.5, color=MUTED,
                    xytext=(2, 0), textcoords="offset points", va="center")


def fig_stress_time_map(cases, out: Path):
    base = cases[0]
    n = len(base.snaps)
    fig = plt.figure(figsize=(15, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, n)
    Tlim, tlim = (40e-3, 3.0), (0.03, 30)
    xb = np.logspace(np.log10(Tlim[0] * 1e3), np.log10(Tlim[1] * 1e3), 70)
    yb = np.logspace(np.log10(tlim[0]), np.log10(tlim[1]), 70)
    for i, s in enumerate(base.snaps):
        ax = fig.add_subplot(gs[0, i])
        p = s.paths
        H, _, _ = np.histogram2d(p.transit_ms, p.GW_tau_eff, bins=[xb, yb], weights=p.GW_HI3)
        H = np.where(H > 0, H / p.GW_HI3.sum(), np.nan)
        m = ax.pcolormesh(xb, yb, H.T, cmap="magma_r", norm=LogNorm(1e-4, 0.2), rasterized=True)
        _iso_lines(ax, Tlim)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(xb[0], xb[-1]); ax.set_ylim(yb[0], yb[-1])
        ax.set_title(f"{s.phase:g}% systole")
        ax.set_xlabel("Transit T [ms]")
        if i == 0:
            ax.set_ylabel("Effective stress τ_eff [Pa]")
        ax.axvline(base.Ts * 1e3 * (1 - s.phase / 100), color="#2E6DB4", lw=1)
    cb = fig.colorbar(m, ax=fig.axes[:n], shrink=0.9, pad=0.01)
    cb.set_label("Share of phase's total HI3\n(sum over paths in bin)")

    ax = fig.add_subplot(gs[1, :])
    for ci, case in enumerate(cases):
        for i, s in enumerate(case.snaps):
            p = s.paths
            Tq = np.percentile(p.transit_ms, [25, 50, 75])
            sq = np.percentile(p.GW_tau_eff, [25, 50, 75])
            st = phase_style(i)
            ax.plot([Tq[0], Tq[2]], [sq[1]] * 2, color=st["color"], lw=1)
            ax.plot([Tq[1]] * 2, [sq[0], sq[2]], color=st["color"], lw=1)
            ax.plot(Tq[1], sq[1], ms=8, ls="none", mfc=st["color"] if ci == 0 else "white",
                    mec=st["color"], marker=st["marker"],
                    label=f"{case.label} {s.phase:g}%")
            # high-damage tail: HI3-weighted median
            w = p.GW_HI3 / p.GW_HI3.sum()
            o = np.argsort(p.transit_ms.to_numpy())
            Tw = p.transit_ms.to_numpy()[o][np.searchsorted(np.cumsum(w.to_numpy()[o]), 0.5)]
            o = np.argsort(p.GW_tau_eff.to_numpy())
            tw = p.GW_tau_eff.to_numpy()[o][np.searchsorted(np.cumsum(w.to_numpy()[o]), 0.5)]
            ax.plot(Tw, tw, marker=st["marker"], ms=8, ls="none", mfc="none",
                    mec=st["color"], mew=1.6)
            ax.plot([Tq[1], Tw], [sq[1], tw], color=st["color"], lw=0.6, ls=":")
    # predicted 35 cc arrows from the 50 % phase median (GW, viscous & inertial bounds)
    ref = base.snaps[int(np.argmin([abs(s.phase - 50) for s in base.snaps]))].paths
    T0, t0 = ref.transit_ms.median(), ref.GW_tau_eff.median()
    for (scen, r), col in zip(scaling.scenarios(base.bpm).items(), ("#5E8C31", "#B8860B")):
        for reg, ls in (("viscous (τ∝μU/L)", "-"), ("inertial (τ∝ρU²)", "--")):
            T1, t1 = T0 * r["T_ratio"], t0 * r[f"tau_ratio|{reg}"]
            ax.annotate("", (T1, t1), (T0, t0), arrowprops=dict(arrowstyle="->", color=col, lw=1.4, ls=ls))
        T1, t1 = T0 * r["T_ratio"], t0 * r["tau_ratio|inertial (τ∝ρU²)"]
        ax.annotate(f"35 cc: {scen}", (T1, t1), xytext=(-8, 10 if t1 > t0 else -14),
                    textcoords="offset points", fontsize=7.5, color=col, ha="right")
    _iso_lines(ax, (60e-3, 600e-3), levels=(1e-6, 1e-5, 1e-4))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Transit T [ms]"); ax.set_ylabel("Effective stress τ_eff [Pa]")
    ax.set_title("Phase medians (filled, with IQR) vs damage-weighted medians (hollow)")
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    note(fig, "τ_eff = [Σ τ^(α/β) dt / T]^(β/α) (GW), so per-path HI3 = A·T^β·τ_eff^α exactly; dashed grey = iso-HI3 lines. "
              "Blue line = time left in systole after the snapshot. Arrows: similarity-scaling predictions for 35 cc "
              "(solid viscous, dashed inertial bound) from the 50% median.")
    fig.savefig(out / "02_stress_time_map.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 03 Damage timing / quasi-steady validity
# ---------------------------------------------------------------------------

def fig_damage_timing(case: Case, out: Path):
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)
    w = np.logspace(0, np.log10(3000), 120)  # ms
    for i, s in enumerate(case.snaps):
        st = phase_style(i)
        inner, tte = s.inner("GW"), s.time_to_exit() * 1e3
        o = np.argsort(tte)
        cum = np.cumsum(inner[o]) / inner.sum()
        frac = np.interp(w, tte[o], cum)
        rem = case.Ts * 1e3 * (1 - s.phase / 100)
        axs[0].plot(w, frac, color=st["color"], label=f"{s.phase:g}%")
        axs[0].plot(rem, np.interp(rem, tte[o], cum), marker=st["marker"], color=st["color"], ms=7)
        base_mean = s.paths.GW_HI3.mean()
        ws = np.r_[1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000]
        ret = [truncated_hi3(s, x / 1e3).mean() / base_mean for x in ws]
        axs[1].plot(ws, ret, color=st["color"])
        axs[1].plot(rem, truncated_hi3(s, rem / 1e3).mean() / base_mean, marker=st["marker"],
                    color=st["color"], ms=7)
        T = np.sort(s.paths.transit_ms.to_numpy())
        axs[2].plot(T, np.arange(1, len(T) + 1) / len(T), color=st["color"])
        axs[2].plot(rem, np.searchsorted(T, rem) / len(T), marker=st["marker"], color=st["color"], ms=7)
    axs[0].set_title("a) Where in time the dose is deposited")
    axs[0].set_xlabel("Time before path exit [ms]")
    axs[0].set_ylabel("Cumulative share of total linear GW dose")
    axs[0].legend(title="phase", fontsize=8)
    axs[1].set_title("b) HI3 kept if only the final window counts")
    axs[1].set_xlabel("Exposure window before exit [ms]")
    axs[1].set_ylabel("Mean HI3 (window) / mean HI3 (full path)")
    axs[2].set_title("c) Transit time vs time left in systole")
    axs[2].set_xlabel("Transit time T [ms]")
    axs[2].set_ylabel("Cumulative fraction of paths")
    axs[2].axvline(case.Ts * 1e3, color=MUTED, ls="--", lw=1)
    axs[2].text(case.Ts * 1e3, 0.03, " full systole", color=MUTED, fontsize=8)
    for ax in axs:
        ax.set_xscale("log")
    note(fig, "Markers = time remaining in systole after each snapshot, Ts·(1−phase). Paths are seeded on the diaphragm inlet "
              "of a frozen field, so long transits are not physical within one systole; panels a–b show how much of the "
              "predicted damage actually depends on that early, slow part of the path.")
    fig.savefig(out / "03_damage_timing.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 04 Damage map
# ---------------------------------------------------------------------------

def fig_damage_map(case: Case, out: Path, voxel_mm=1.0, plane=("x", "z")):
    ax_idx = {"x": 0, "y": 1, "z": 2}
    i0, i1 = ax_idx[plane[0]], ax_idx[plane[1]]
    allxyz = np.vstack([s.xyz for s in case.snaps]) * 1e3
    lo, hi = allxyz.min(0) - 1, allxyz.max(0) + 1
    edges = [np.arange(lo[k], hi[k] + voxel_mm, voxel_mm) for k in range(3)]
    n = len(case.snaps)
    fig, axs = plt.subplots(2, (n + 2) // 2, figsize=(15, 9), constrained_layout=True)
    axs = axs.ravel()
    hotspots, conc = [], []
    for i, s in enumerate(case.snaps):
        ax = axs[i]
        xyz = s.xyz * 1e3
        inner = s.inner("GW")
        # midpoint of each segment carries its dose
        mid = xyz.copy()
        mid[1:] = 0.5 * (xyz[1:] + xyz[:-1])
        H2, _, _ = np.histogram2d(mid[:, i0], mid[:, i1], bins=[edges[i0], edges[i1]], weights=inner)
        H2 = H2 / inner.sum()
        occ, _, _ = np.histogram2d(xyz[:, i0], xyz[:, i1], bins=[edges[i0], edges[i1]])
        ax.pcolormesh(edges[i0], edges[i1], np.where(occ.T > 0, 1, np.nan), cmap="Greys",
                      vmin=0, vmax=8, rasterized=True)
        cmap = plt.get_cmap("magma_r").copy()
        cmap.set_under("#EDEDED")
        m = ax.pcolormesh(edges[i0], edges[i1], np.where(H2.T > 0, H2.T, np.nan), cmap=cmap,
                          norm=LogNorm(1e-5, 0.1), rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"{s.phase:g}% systole")
        ax.set_xlabel(f"{plane[0].upper()} [mm]"); ax.set_ylabel(f"{plane[1].upper()} [mm]")
        # 3D voxel concentration
        H3, _ = np.histogramdd(mid, bins=edges, weights=inner)
        occ3, _ = np.histogramdd(xyz, bins=edges)
        vals = np.sort(H3[occ3 > 0])[::-1] / inner.sum()
        cum = np.cumsum(vals)
        conc.append((s.phase, np.arange(1, len(vals) + 1) / len(vals), cum))
        k50 = np.searchsorted(cum, 0.5) + 1
        ax.text(0.02, 0.02, f"50% of dose in {k50} voxels\n({100*k50/len(vals):.2f}% of visited volume)",
                transform=ax.transAxes, fontsize=7.5, color=INK,
                bbox=dict(fc="white", ec="none", alpha=0.8))
        flat = np.argsort(H3.ravel())[::-1][:20]
        for rank, f in enumerate(flat):
            ix = np.unravel_index(f, H3.shape)
            c = [0.5 * (edges[k][ix[k]] + edges[k][ix[k] + 1]) for k in range(3)]
            hotspots.append({"case": case.label, "phase_pct": s.phase, "rank": rank + 1,
                             "x_mm": c[0], "y_mm": c[1], "z_mm": c[2],
                             "dose_share": H3.ravel()[f] / inner.sum()})
        top = hotspots[-20]
        ax.plot(top[f"{plane[0]}_mm"], top[f"{plane[1]}_mm"], "o", mfc="none", mec="#2E6DB4", ms=10, mew=1.5)
    ax = axs[n]
    for i, (ph, xf, cum) in enumerate(conc):
        st = phase_style(i)
        ax.plot(xf * 100, cum, color=st["color"], label=f"{ph:g}%")
    ax.set_xscale("log")
    ax.set_xlabel(f"Visited {voxel_mm:g} mm voxels, ranked by dose [% of visited volume]")
    ax.set_ylabel("Cumulative share of GW dose")
    ax.set_title("Dose concentration")
    ax.legend(title="phase", fontsize=8)
    for a in axs[n + 1:]:
        a.set_visible(False)
    fig.colorbar(m, ax=axs[:n].tolist(), shrink=0.8, extend="min", label=f"Share of total linear GW dose per {voxel_mm:g} mm column")
    note(fig, f"{case.label}: linear HI3 dose A^(1/β)·τ^(α/β)·dt deposited at segment midpoints, projected on the "
              f"{plane[0].upper()}–{plane[1].upper()} plane (same view as 07_spatial_path_stress). Light grey = visited, dose share <1e-5. "
              "Blue circle = highest-dose 3D voxel (see hotspots.csv to probe in CFD-Post). Linear dose shares, not final HI3 shares.")
    fig.savefig(out / f"04_damage_map_{plane[0]}{plane[1]}.png")
    plt.close(fig)
    return pd.DataFrame(hotspots)


# ---------------------------------------------------------------------------
# 05 Numerical audit
# ---------------------------------------------------------------------------

def fig_numerical_audit(case: Case, out: Path):
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.2), constrained_layout=True)
    axs = axs.ravel()
    rows = []
    snaps = case.all_snaps          # audit always looks at every exported path
    ph = [s.phase for s in snaps]
    inc_frac, inc_share = [], []
    for i, s in enumerate(snaps):
        st = phase_style(i)
        p = s.paths
        inc = ~p.completed.to_numpy()
        inc_frac.append(inc.mean() * 100)
        inc_share.append(p.GW_HI3[inc].sum() / p.GW_HI3.sum() * 100)
        E = p[["end_x", "end_z"]].to_numpy() * 1e3
        if i == 0:
            axs[0].scatter(E[~inc, 0], E[~inc, 1], s=2, color="#BDBDBD", label="ends on outlet plane", rasterized=True)
        else:
            axs[0].scatter(E[~inc, 0], E[~inc, 1], s=2, color="#BDBDBD", rasterized=True)
        hi3 = p.GW_HI3.to_numpy()[inc]
        axs[0].scatter(E[inc, 0], E[inc, 1], s=8 + 120 * hi3 / p.GW_HI3.max(), marker=st["marker"],
                       facecolor="none", edgecolor=st["color"], lw=0.8, label=f"{s.phase:g}% stops inside domain")
        # stress-sampling resolution: change of ln tau across each segment
        first = np.zeros(len(s.tau), bool); first[s.offsets[:-1]] = True
        dl = np.zeros(len(s.tau))
        dl[1:] = np.abs(np.log(np.maximum(s.tau[1:], 1e-9) / np.maximum(s.tau[:-1], 1e-9)))
        dl[first] = 0
        inner = s.inner("GW")
        o = np.argsort(dl)
        cum = np.cumsum(inner[o][::-1])[::-1] / inner.sum()
        axs[2].plot(np.exp(dl[o]), cum, color=st["color"], label=f"{s.phase:g}%")
        seg = np.zeros(len(s.tau)); seg[1:] = np.linalg.norm(np.diff(s.xyz, axis=0), axis=1) * 1e3
        ok = ~first
        hi = ok & (inner > np.percentile(inner[ok], 99))
        bins = np.logspace(-3, 1, 50)
        axs[3].hist(seg[ok], bins=bins, histtype="step", color=st["color"], lw=0.8, ls=":", density=True)
        axs[3].hist(seg[hi], bins=bins, histtype="step", color=st["color"], lw=1.6, density=True)
        last = s.offsets[1:] - 1
        rows.append({"case": case.label, "phase_pct": s.phase, "n_paths": len(p),
                     "n_incomplete": int(inc.sum()), "frac_incomplete": float(inc.mean()),
                     "GW_HI3_share_incomplete": float(p.GW_HI3[inc].sum() / p.GW_HI3.sum()),
                     "GW_HI3_mean_all": p.GW_HI3.mean(), "GW_HI3_mean_completed": p.GW_HI3[~inc].mean(),
                     "end_speed_median_incomplete": float(np.median(s.vel[last][inc])) if inc.any() else np.nan,
                     "end_speed_median_completed": float(np.median(s.vel[last][~inc])),
                     "end_stress_median_incomplete_Pa": float(np.median(s.tau[last][inc])) if inc.any() else np.nan,
                     "end_stress_median_completed_Pa": float(np.median(s.tau[last][~inc])),
                     "velocity_floor_time_frac": s.dt[s.floored].sum() / s.dt.sum(),
                     "dose_share_segments_tau_jump_gt_2x": float(inner[dl > np.log(2)].sum() / inner.sum()),
                     "dose_share_segments_tau_jump_gt_5x": float(inner[dl > np.log(5)].sum() / inner.sum()),
                     "segment_len_median_mm": float(np.median(seg[ok])),
                     "segment_len_median_top1pct_dose_mm": float(np.median(seg[hi]))})
    axs[0].set_aspect("equal")
    axs[0].set_title("a) Path end points (X–Z): marker size ∝ path HI3")
    axs[0].set_xlabel("X [mm]"); axs[0].set_ylabel("Z [mm]")
    axs[0].legend(fontsize=7, loc="lower left", markerscale=1.5)
    x = np.arange(len(ph))
    axs[1].bar(x - 0.2, inc_frac, 0.38, color="#BDBDBD", label="% of paths")
    axs[1].bar(x + 0.2, inc_share, 0.38, color=case.colour, label="% of summed GW HI3")
    for xi, (f_, h_) in enumerate(zip(inc_frac, inc_share)):
        axs[1].text(xi - 0.2, f_ + 0.8, f"{f_:.1f}", ha="center", fontsize=7.5, color=INK)
        axs[1].text(xi + 0.2, h_ + 0.8, f"{h_:.0f}", ha="center", fontsize=7.5, color=INK)
    axs[1].set_xticks(x, [f"{v:g}%" for v in ph]); axs[1].set_xlabel("Systolic phase")
    axs[1].set_ylabel("%"); axs[1].set_title("b) Paths that never reach the outlet carry a large share of HI3")
    axs[1].legend(fontsize=8)
    axs[2].set_title("c) Stress resolution along paths")
    axs[2].set_xscale("log")
    axs[2].set_xlabel("Stress ratio across one segment  max(τ_i/τ_i−1, τ_i−1/τ_i)")
    axs[2].set_ylabel("Share of GW dose in segments ≥ this ratio")
    axs[2].axvline(2, color=MUTED, ls="--", lw=1)
    axs[2].legend(title="phase", fontsize=8)
    axs[3].set_title("d) Segment length: all (dotted) vs top-1% dose segments (solid)")
    axs[3].set_xscale("log")
    axs[3].set_xlabel("Segment length [mm]"); axs[3].set_ylabel("Density")
    note(fig, "a–b) outlet = plane fitted to path end points; paths ending >1 mm off it stopped inside the domain "
              "(typically slow, high-stress near-wall/valve-gap points). Rerun with --completed-only to exclude them. "
              "c) dose in segments where τ jumps >2× is sensitive to end-point stress sampling (export more segments). "
              "d) spatial resolution where damage is accumulated.", bottom=0.05)
    fig.savefig(out / "05_numerical_audit.png")
    plt.close(fig)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 06 Scaling hypothesis
# ---------------------------------------------------------------------------

def fig_scaling(case: Case, out: Path, cases=()):
    sc = scaling.scenarios(case.bpm)
    t = case.table
    ph = t.phase_pct.to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True,
                            gridspec_kw=dict(width_ratios=[1.2, 1.1, 1.15]))
    regs = list(scaling.REGIMES)
    cols = {list(sc)[0]: "#5E8C31", list(sc)[1]: "#B8860B"}
    ax = axs[0]
    for k, mk in (("GW", "o"), ("HO", "s")):
        ax.plot(ph, t[f"{k}_HI3_mean"], marker=mk, color=case.colour, label=f"{case.label} {k} (CFD)")
        for scen, r in sc.items():
            lo = min(r[f"HI3_ratio|{k}|{g}"] for g in regs)
            hi = max(r[f"HI3_ratio|{k}|{g}"] for g in regs)
            ax.fill_between(ph, t[f"{k}_HI3_mean"] * lo, t[f"{k}_HI3_mean"] * hi, color=cols[scen],
                            alpha=0.22 if k == "GW" else 0.12, lw=0,
                            label=f"35 cc predicted: {scen} ({k})")
    for c2 in cases[1:]:
        ax.plot(c2.table.phase_pct, c2.table.GW_HI3_mean, marker="o", color=c2.colour, label=f"{c2.label} GW (CFD)")
    ax.set_yscale("log"); ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel("Mean per-path HI3 [%]")
    ax.set_title("a) What similarity scaling predicts for 35 cc")
    ax.legend(fontsize=7, loc="upper center", ncol=2)

    ax = axs[1]
    short = ["viscous", "boundary\nlayer", "inertial"]
    xt, xl = [], []
    for j, (scen, r) in enumerate(sc.items()):
        for gi, g in enumerate(regs):
            xpos = j * (len(regs) + 1) + gi
            xt.append(xpos); xl.append(short[gi])
            for k, mk in (("GW", "o"), ("HO", "s"), ("TZ", "^")):
                ax.plot(xpos, r[f"HI3_ratio|{k}|{g}"], marker=mk, ms=8, ls="none", color=cols[scen],
                        mfc=cols[scen] if k == "GW" else "white",
                        label=f"{k}" if (j == 0 and gi == 0) else None)
        mid = j * (len(regs) + 1) + (len(regs) - 1) / 2
        ax.text(mid, 0.42, scen, ha="center", fontsize=8.5, color=cols[scen], fontweight="bold")
    ax.axhline(1, color=MUTED, lw=1)
    ax.axhline(1.56, color="#7A5195", ls=":", lw=1.2)
    ax.text(xt[-1] + 0.4, 1.5, "Marom 2014 median-SA ratio\n35/50 = 1.56 (164 vs 140 bpm)", fontsize=6.5,
            color="#7A5195", ha="right", va="top")
    ax.set_yscale("log"); ax.set_ylim(0.38, 3.2)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3], ["0.5", "0.75", "1", "1.5", "2", "3"])
    ax.set_xticks(xt, xl, fontsize=8)
    ax.set_ylabel("Predicted HI3 ratio 35 cc / 50 cc")
    ax.set_title("b) Ratio by stress-scaling regime and model")
    leg = ax.legend(fontsize=8, loc="upper left", title="constants")
    for h in leg.legend_handles:
        h.set_color(MUTED)

    ax = axs[2]
    ax.axis("off")
    keys = [("length_scale", "Length scale s"), ("U_ratio", "Velocity U"), ("T_ratio", "Transit time T"),
            ("Q_ratio", "Instant. flow Q"), ("CO_ratio", "Cardiac output"), ("Re_ratio", "Reynolds no."),
            ("Wo_ratio", "Womersley no.")]
    cell = [[lab] + [f"{r[k]:.3f}" for r in sc.values()] for k, lab in keys]
    tb = ax.table(cellText=cell, colLabels=["35/50 ratio"] + list(sc), loc="center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(0.95, 1.5)
    tb.auto_set_column_width([0, 1, 2])
    ax.set_title("c) Dimensionless ratios (geometric similarity)")
    note(fig, "Assumes all lengths (incl. diaphragm stroke and valve opening) scale by (35/50)^(1/3), the same systolic "
              "fraction and the same valve state at matched phase. Bands span viscous → inertial stress scaling. "
              "A hypothesis to test against the 35 cc CFD, not a result.")
    fig.savefig(out / "06_scaling_hypothesis.png")
    plt.close(fig)
    rows = []
    for scen, r in sc.items():
        rows.append({"scenario": scen, **r})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 07 Size comparison (needs >= 2 cases)
# ---------------------------------------------------------------------------

def fig_comparison(cases, out: Path):
    ref = cases[0]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    rows = []
    for c in cases[1:]:
        common = sorted(set(ref.table.phase_pct).intersection(c.table.phase_pct))
        if not common:
            print(f"  {c.label}: no phases in common with {ref.label}; skipping comparison")
            continue
        rat, lo, hi, dec = [], [], [], []
        rng = np.random.default_rng(1)
        for ph in common:
            a = next(s for s in ref.snaps if s.phase == ph).paths
            b = next(s for s in c.snaps if s.phase == ph).paths
            r = b.GW_HI3.mean() / a.GW_HI3.mean()
            bs = []
            for _ in range(1000):
                bs.append(rng.choice(b.GW_HI3.to_numpy(), len(b)).mean() /
                          rng.choice(a.GW_HI3.to_numpy(), len(a)).mean())
            rat.append(r); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
            d = log_decomposition(a, b)
            dec.append(d)
            rows.append({"case": c.label, "ref": ref.label, "phase_pct": ph, "GW_HI3_mean_ratio": r,
                         "ci_lo": lo[-1], "ci_hi": hi[-1], "log10_time_term": d["time"],
                         "log10_stress_term": d["stress"], "log10_total_geomean": d["total_geo"]})
        common = np.array(common)
        axs[0].errorbar(common, rat, yerr=[np.subtract(rat, lo), np.subtract(hi, rat)], marker="o",
                        color=c.colour, capsize=3, label=f"{c.label} / {ref.label}")
        sc = scaling.ratios(35, 50, c.bpm, ref.bpm)
        g = [sc[f"HI3_ratio|GW|{k}"] for k in scaling.REGIMES]
        axs[0].axhspan(min(g), max(g), color=c.colour, alpha=0.12, lw=0)
        x = np.arange(len(common))
        off = (cases.index(c) - 1) * 0.12
        axs[1].bar(x - 0.18 + off, [d["time"] for d in dec], 0.16, color="#8FB8DE", edgecolor=c.colour)
        axs[1].bar(x + 0.02 + off, [d["stress"] for d in dec], 0.16, color="#D98CB3", edgecolor=c.colour)
        axs[1].plot(x + off - 0.08, [d["total_geo"] for d in dec], "_", color="k", ms=18, mew=2)
        axs[1].set_xticks(x, [f"{p:g}%" for p in common])
        for i, ph in enumerate(common):
            for cc, ls in ((ref, "-"), (c, "--")):
                s = next(s for s in cc.snaps if s.phase == ph)
                h = np.sort(s.paths.GW_HI3.to_numpy())
                axs[2].plot(h, np.arange(1, len(h) + 1) / len(h), color=PHASE_COLOURS[i % 5], ls=ls, lw=1.3)
    axs[0].axhline(1, color=MUTED, lw=1)
    axs[0].set_yscale("log"); axs[0].set_xlabel("Systolic phase [%]")
    axs[0].set_ylabel("Mean GW HI3 ratio (95% bootstrap CI)")
    axs[0].set_title("a) Size effect per phase (band = scaling prediction)")
    axs[0].legend(fontsize=8)
    axs[1].axhline(0, color=MUTED, lw=0.8)
    axs[1].set_ylabel("log₁₀ ratio (blue: transit, pink: stress)")
    axs[1].set_title("b) Why: transit vs stress contribution")
    from matplotlib.lines import Line2D
    axs[2].legend(handles=[Line2D([], [], color=PHASE_COLOURS[i % 5], label=f"{p:g}%")
                           for i, p in enumerate(ref.table.phase_pct)], title="phase", fontsize=8)
    axs[2].set_xscale("log"); axs[2].set_xlabel("Per-path GW HI3 [%]"); axs[2].set_ylabel("ECDF")
    axs[2].set_title(f"c) Distributions (solid {ref.label}, dashed comparison)")
    note(fig, "Ratios compare matched phase fractions. Bootstrap resamples paths within each case (sampling spread, "
              "not CFD/model uncertainty).")
    fig.savefig(out / "07_size_comparison.png")
    plt.close(fig)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs=3, action="append", metavar=("LABEL", "FOLDER", "BPM"), required=True,
                    help="first case is the reference (e.g. 50 cc)")
    ap.add_argument("--systolic-fraction", type=float, default=0.5,
                    help="fraction of the cycle that is systole (default 0.5 -> 250 ms at 120 bpm)")
    ap.add_argument("--completed-only", action="store_true",
                    help="drop paths that do not end on the outlet plane from every analysis")
    ap.add_argument("--outdir", type=Path, default=None)
    a = ap.parse_args(argv)
    out = a.outdir or HERE / "outputs" / datetime.now().strftime("%Y%m%d_%H%M%S")
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    for i, (lab, folder, bpm) in enumerate(a.case):
        print(f"Case {lab}: {folder} @ {bpm} bpm")
        cases.append(Case(lab, folder, float(bpm), a.systolic_fraction, CASE_COLOURS[i % len(CASE_COLOURS)],
                          completed_only=a.completed_only))

    pd.concat([c.table for c in cases]).to_csv(out / "phase_summary.csv", index=False)
    for c in cases:
        for s in c.snaps:
            s.paths.assign(case=c.label, phase_pct=s.phase).to_csv(
                out / f"paths_{c.label.replace(' ', '')}_{s.phase:g}.csv", index=False)
    print("figures ...")
    hot, aud = [], []
    for c in cases:
        tag = "" if len(cases) == 1 else f"_{c.label.replace(' ', '')}"
        sub = out if not tag else out / c.label.replace(" ", "")
        sub.mkdir(exist_ok=True)
        fig_phase_mechanism(c, sub)
        fig_damage_timing(c, sub)
        hot.append(fig_damage_map(c, sub, plane=("x", "z")))
        fig_damage_map(c, sub, plane=("x", "y"))
        aud.append(fig_numerical_audit(c, sub))
    fig_stress_time_map(cases, out)
    fig_scaling(cases[0], out, cases).to_csv(out / "scaling_ratios.csv", index=False)
    pd.concat(hot).to_csv(out / "hotspots.csv", index=False)
    pd.concat(aud).to_csv(out / "numerical_audit.csv", index=False)
    if len(cases) > 1:
        fig_comparison(cases, out).to_csv(out / "size_comparison.csv", index=False)
    (out / "run.json").write_text(json.dumps({
        "created": datetime.now().isoformat(timespec="seconds"),
        "cases": [{"label": c.label, "folder": str(c.folder), "bpm": c.bpm, "systole_ms": c.Ts * 1e3,
                   "files": [s.name for s in c.snaps]} for c in cases],
        "systolic_fraction": a.systolic_fraction, "completed_only": a.completed_only}, indent=2))
    print(f"done -> {out}")


if __name__ == "__main__":
    main()
