"""Report-style 50 cc vs 35 cc comparison figures and tables (same style as report_figures.py).

    python compare_figures.py --case "50 cc 120 BPM" ..\\z_inputs 120 1.0 \\
                              --case "35 cc 120 BPM" "..\\35cc_inputs\\120BPM" 120 0.888 \\
                              --case "35 cc 171 BPM" "..\\35cc_inputs\\171BPM" 171 0.888 \\
                              --monitor "..\\35cc_inputs\\Shear Strain Rate Convergence Oscillation 35cc 83% 120BPM\\variableoutput.csv"

The first case is the reference. Writes to report_figs/ (override with --outdir):
  size_comparison           (main text) HI3 per phase, 35/50 ratios vs the scaling
                            bands, and measured transit-time / effective-stress ratios
                            vs their predictions
  size_comparison_completed same, outlet-reaching paths only
  size_dose_maps            where the dose is deposited, every case at 16.67/50/83.33 %
  size_lysis_sa             lysis margins and stress accumulation per case
  convergence_monitor       steady-solver oscillation of the valve shear-rate monitor
  size_mechanism            F2-style: HI3 against inflow and valve angle, every case
  size_lysis_envelopes      F7-style: stress-duration envelopes per phase, every case
  size_summary.csv          every number quoted in the comparison
  --only C1 C2 ...          rebuild a subset
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

import empirical_envelope as ee
import explore_core as ec
import explore_hemolysis as eh
import report_style as rs
import scaling

rs.apply()
HERE = Path(__file__).resolve().parent
W = rs.WIDTH
INK, MUTED = rs.INK, rs.MUTED
FORMATS = list(rs.FORMATS)
GW = ec.CONSTANTS["GW"]


def save(fig, out, name):
    rs.save(fig, out, name, FORMATS)


def style(label):
    return rs.CASE_STYLE.get(label, dict(color=MUTED, marker="o", ls="-"))


def logfmt(ax, axis="y"):
    (ax.yaxis if axis == "y" else ax.xaxis).set_minor_formatter(mt.NullFormatter())


class Run:
    """One case: snapshots (all or outlet-reaching paths) + metadata."""

    def __init__(self, label, folder, bpm, scale, completed_only):
        self.label, self.bpm, self.scale = label, float(bpm), float(scale)
        self.case = eh.Case(label, folder, bpm, 0.5, style(label)["color"], completed_only=completed_only)
        self.snaps = {round(s.phase, 2): s for s in self.case.snaps}
        self.all_snaps = {round(s.phase, 2): s for s in self.case.all_snaps}
        self.phases = sorted(self.snaps)


def boot_ratio(a, b, n=2000, seed=1):
    rng = np.random.default_rng(seed)
    ia = rng.integers(0, len(a), (n, len(a)))
    ib = rng.integers(0, len(b), (n, len(b)))
    r = b[ib].mean(1) / a[ia].mean(1)
    return np.percentile(r, [2.5, 97.5])


def geo(x):
    return float(np.exp(np.mean(np.log(x))))


def summary_table(runs):
    ref = runs[0]
    rows = []
    for r in runs:
        sc = scaling.ratios(35, 50, r.bpm, ref.bpm) if r is not ref else None
        for ph in r.phases:
            p = r.snaps[ph].paths
            a = ref.snaps[ph].paths
            hi, ha = p.GW_HI3.to_numpy(), a.GW_HI3.to_numpy()
            lo_ci, hi_ci = eh.boot_mean_ci(hi)
            row = {"case": r.label, "bpm": r.bpm, "length_scale": r.scale, "phase_pct": ph,
                   "n_paths": len(p), "GW_HI3_mean": hi.mean(), "GW_HI3_ci_lo": lo_ci, "GW_HI3_ci_hi": hi_ci,
                   "GW_HI3_median": np.median(hi), "HO_HI3_mean": p.HO_HI3.mean(), "TZ_HI3_mean": p.TZ_HI3.mean(),
                   "GW_HI2_mean": p.GW_HI2.mean(),
                   "GW_NIH_mg_100L": hi.mean() * (1 - ec.HCT) * ec.HB * 1000,
                   "transit_ms_median": p.transit_ms.median(), "transit_ms_geomean": geo(p.transit_ms),
                   "tau_eff_median_Pa": p.GW_tau_eff.median(), "tau_eff_geomean_Pa": geo(p.GW_tau_eff),
                   "tau_peak_p99_Pa": np.percentile(p.tau_peak, 99), "speed_p99_m_s": np.percentile(r.snaps[ph].vel, 99),
                   "SA_median": p.SA.median(), "SA_p95": np.percentile(p.SA, 95)}
            pa = r.all_snaps[ph].paths
            row["frac_paths_incomplete"] = float((~pa.completed).mean())
            row["GW_HI3_share_incomplete"] = float(pa.GW_HI3[~pa.completed].sum() / pa.GW_HI3.sum())
            if r is not ref:
                lo, hi2 = boot_ratio(ha, hi)
                d = eh.log_decomposition(a, p)
                row.update({"ratio_GW_HI3": hi.mean() / ha.mean(), "ratio_ci_lo": lo, "ratio_ci_hi": hi2,
                            "ratio_HO_HI3": p.HO_HI3.mean() / a.HO_HI3.mean(),
                            "ratio_TZ_HI3": p.TZ_HI3.mean() / a.TZ_HI3.mean(),
                            "ratio_GW_HI3_geomean": 10 ** d["total_geo"],
                            "log10_time_term": d["time"], "log10_stress_term": d["stress"],
                            "ratio_transit_geomean": geo(p.transit_ms) / geo(a.transit_ms),
                            "ratio_tau_eff_geomean": geo(p.GW_tau_eff) / geo(a.GW_tau_eff),
                            "ratio_SA_median": p.SA.median() / a.SA.median(),
                            "pred_T_ratio": sc["T_ratio"],
                            "pred_tau_ratio_viscous": sc["tau_ratio|viscous (τ∝μU/L)"],
                            "pred_tau_ratio_inertial": sc["tau_ratio|inertial (τ∝ρU²)"],
                            "pred_HI3_ratio_GW_lo": min(sc[f"HI3_ratio|GW|{k}"] for k in scaling.REGIMES),
                            "pred_HI3_ratio_GW_hi": max(sc[f"HI3_ratio|GW|{k}"] for k in scaling.REGIMES)})
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------

def c1_size_comparison(runs, tab, out, name="size_comparison"):
    ref = runs[0]
    fig, axs = plt.subplots(2, 2, figsize=(W, 5.0), constrained_layout=True)
    ax = axs[0, 0]
    for r in runs:
        t = tab[tab.case == r.label]
        st = style(r.label)
        y = t.GW_HI3_mean.to_numpy()
        ax.errorbar(t.phase_pct, y, yerr=[y - t.GW_HI3_ci_lo, t.GW_HI3_ci_hi - y], color=st["color"],
                    marker=st["marker"], ls=st["ls"], capsize=2.5, ms=5, label=r.label)
    ax.set_yscale("log"); logfmt(ax)
    ax.set_yticks([3e-5, 1e-4, 3e-4], ["3e-5", "1e-4", "3e-4"])
    ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel("Mean per-path HI3 [%] (GW)")
    ax.set_title("a) Predicted damage per pass")
    ax.legend(loc="upper center", fontsize=6.5)

    ax = axs[0, 1]
    bands = {}
    for k, r in enumerate(runs[1:]):
        t = tab[tab.case == r.label]
        st = style(r.label)
        lo, hi = t.pred_HI3_ratio_GW_lo.iloc[0], t.pred_HI3_ratio_GW_hi.iloc[0]
        ax.axhspan(lo, hi, color=st["color"], alpha=0.13, lw=0)
        if r.bpm > runs[0].bpm:
            ax.text(15.5, hi * 1.03, f"predicted, {r.bpm:g} BPM", color=st["color"], fontsize=6.3, va="bottom")
        else:
            ax.text(15.5, lo * 1.03, f"predicted, {r.bpm:g} BPM", color=st["color"], fontsize=6.3, va="bottom")
        off = (k - 0.5) * 1.2
        y = t.ratio_GW_HI3.to_numpy()
        ax.errorbar(t.phase_pct + off, y, yerr=[y - t.ratio_ci_lo, t.ratio_ci_hi - y], color=st["color"],
                    marker=st["marker"], ls="none", capsize=2.5, ms=5, label=f"{r.label} / 50 cc")
    ax.axhline(1, color=MUTED, lw=0.8)
    ax.set_yscale("log"); logfmt(ax)
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3], ["0.5", "0.75", "1", "1.5", "2", "3"])
    ax.set_ylim(0.45, 3.2)
    ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel("Mean HI3 ratio, 35/50 cc (95% CI)")
    ax.set_title("b) Size effect against the scaling prediction")
    ax.legend(loc="lower right", fontsize=6.3)

    for ax, col, pred_cols, ttl, ylab in (
            (axs[1, 0], "ratio_transit_geomean", ("pred_T_ratio",), "c) Transit time", "Transit-time ratio, 35/50 cc"),
            (axs[1, 1], "ratio_tau_eff_geomean", ("pred_tau_ratio_viscous", "pred_tau_ratio_inertial"),
             "d) Effective stress", "Effective-stress ratio, 35/50 cc")):
        for k, r in enumerate(runs[1:]):
            t = tab[tab.case == r.label]
            st = style(r.label)
            for j, pc in enumerate(pred_cols):
                ax.axhline(t[pc].iloc[0], color=st["color"], lw=1, ls=(":" if j else "--"))
            ax.plot(t.phase_pct, t[col], color=st["color"], marker=st["marker"], ls="none", ms=5.5,
                    label=r.label)
        ax.set_xlabel("Systolic phase [%]"); ax.set_ylabel(ylab); ax.set_title(ttl)
        ax.axhline(1, color=MUTED, lw=0.6)
    axs[1, 0].set_ylim(0.55, 1.2)
    axs[1, 0].text(0.02, 0.04, "dashed: predicted (T ∝ 1/BPM)", transform=axs[1, 0].transAxes,
                   fontsize=6.3, color=MUTED)
    axs[1, 1].set_ylim(0.6, 1.8)
    axs[1, 1].text(0.02, 0.04, "dashed: viscous, dotted: inertial prediction", transform=axs[1, 1].transAxes,
                   fontsize=6.3, color=MUTED)
    axs[1, 1].legend(loc="center right", fontsize=6.3)
    save(fig, out, name)


def c2_dose_maps(runs, out, phases=(16.67, 50.0, 83.33)):
    allxyz = np.vstack([s.xyz for r in runs for s in r.all_snaps.values()]) * 1e3
    lo, hi = allxyz.min(0) - 1, allxyz.max(0) + 1
    edges = [np.arange(lo[k], hi[k] + 1.0, 1.0) for k in range(3)]
    fig, axs = plt.subplots(len(runs), len(phases), figsize=(W, 1.55 * len(runs) + 0.5),
                            constrained_layout=True, sharex=True, sharey=True)
    cmap = plt.get_cmap("magma_r").copy(); cmap.set_under("#EDEDED")
    for i, r in enumerate(runs):
        for j, ph in enumerate(phases):
            s = r.all_snaps[ph]
            ax = axs[i, j]
            xyz = s.xyz * 1e3
            mid = xyz.copy(); mid[1:] = 0.5 * (xyz[1:] + xyz[:-1])
            inner = s.inner("GW")
            H, _, _ = np.histogram2d(mid[:, 0], mid[:, 2], bins=[edges[0], edges[2]], weights=inner)
            occ, _, _ = np.histogram2d(xyz[:, 0], xyz[:, 2], bins=[edges[0], edges[2]])
            H = np.where(occ > 0, np.maximum(H / inner.sum(), 1e-9), np.nan)
            m = ax.pcolormesh(edges[0], edges[2], H.T, cmap=cmap, norm=LogNorm(1e-5, 0.1), rasterized=True)
            H3, _ = np.histogramdd(mid, bins=edges, weights=inner)
            occ3, _ = np.histogramdd(xyz, bins=edges)
            vals = np.sort(H3[occ3 > 0])[::-1] / inner.sum()
            k50 = np.searchsorted(np.cumsum(vals), 0.5) + 1
            ax.text(0.03, 0.04, f"50% in {100 * k50 / len(vals):.2f}%", transform=ax.transAxes, fontsize=6,
                    bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
            ax.set_aspect("equal"); ax.grid(False)
            if i == 0:
                ax.set_title(f"{ph:g}% systole")
            if j == 0:
                ax.set_ylabel(f"{r.label}\nZ [mm]", fontsize=7)
            if i == len(runs) - 1:
                ax.set_xlabel("X [mm]")
    fig.colorbar(m, ax=axs.ravel().tolist(), shrink=0.7, pad=0.01, extend="min", label="Share of linear dose")
    save(fig, out, "size_dose_maps")


def c3_lysis_sa(runs, tab, out):
    curves = ee.lysis_curves()
    rows = []
    for r in runs:
        for ph in r.phases:
            env = ee.worst_envelope(r.all_snaps[ph])
            sf = ee.safety_factors(env, curves)
            rows.append({"case": r.label, "phase_pct": ph, "margin_RBC": sf["RBC"],
                         "margin_platelet": sf["Platelet"], "peak_stress_Pa": env.stress_Pa.max()})
    lys = pd.DataFrame(rows)
    lys.to_csv(out / "size_lysis_margins.csv", index=False)
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.3), constrained_layout=True)
    width = 0.26
    for k, r in enumerate(runs):
        st = style(r.label)
        l = lys[lys.case == r.label]
        x = np.arange(len(l)) + (k - 1) * width
        axs[0].bar(x, l.margin_RBC, width, color=st["color"], label=r.label)
        axs[1].bar(x, l.margin_platelet, width, color=st["color"])
        t = tab[tab.case == r.label]
        axs[2].plot(t.phase_pct + (k - 1) * 1.5, t.SA_median, color=st["color"], marker=st["marker"],
                    ls="none", ms=5)
        axs[2].plot(t.phase_pct + (k - 1) * 1.5, t.SA_p95, color=st["color"], marker=st["marker"],
                    ls="none", ms=5, mfc="white")
    ph = runs[0].phases
    for ax, ttl in ((axs[0], "a) Margin to RBC lysis"), (axs[1], "b) Margin to platelet lysis")):
        ax.set_xticks(range(len(ph)), [f"{p:g}" for p in ph], fontsize=6.5)
        ax.set_xlabel("Systolic phase [%]"); ax.set_title(ttl)
        ax.axhline(1, color=INK, lw=0.8, ls=":")
    axs[0].set_ylabel("Minimum safety factor")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="outside upper center", ncol=3, fontsize=6.8)
    axs[2].set_yscale("log"); logfmt(axs[2])
    axs[2].set_yticks([0.3, 1, 3], ["0.3", "1", "3"])
    axs[2].set_xlabel("Systolic phase [%]"); axs[2].set_ylabel("SA [dyn·s/cm²]")
    axs[2].set_title("c) SA median (●), P95 (○)")
    save(fig, out, "size_lysis_sa")
    return lys


def c4_monitor(monitor_csv, out, label="35 cc, 83.33%, 120 BPM"):
    t = open(monitor_csv).read().splitlines()
    i = next(k for k, l in enumerate(t) if l.startswith('"Accumulated'))
    d = pd.read_csv(io.StringIO("\n".join(t[i:]))).apply(pd.to_numeric, errors="coerce").dropna()
    x, y = d.iloc[:, 0].to_numpy(), d.iloc[:, 1].to_numpy()
    mean, sd = y.mean(), y.std()
    half = (y.max() - y.min()) / 2
    # dominant period from the autocorrelation
    yc = y - mean
    ac = np.correlate(yc, yc, "full")[len(yc) - 1:]
    ac /= ac[0]
    neg = np.argmax(ac < 0)
    period = neg + np.argmax(ac[neg:neg + 5000]) if neg else np.nan
    fig, ax = plt.subplots(figsize=(W, 2.0), constrained_layout=True)
    ax.plot(x, y, color=rs.DEVICE_COLOURS["35 cc"], lw=0.6)
    run = pd.Series(y).rolling(int(period) if np.isfinite(period) else 1000, min_periods=1).mean()
    ax.plot(x, run, color=INK, lw=1.2, label="running mean over one period")
    ax.axhline(mean, color=MUTED, ls="--", lw=0.8)
    ax.set_xlabel("Accumulated time step"); ax.set_ylabel("Area-averaged shear\nstrain rate at valve [1/s]")
    ax.set_title(f"Valve shear-rate monitor, {label}")
    ax.text(0.99, 0.97, f"mean {mean:.0f} 1/s, ±{100 * half / mean:.1f}% (SD {100 * sd / mean:.2f}%), "
            f"period ≈ {period:.0f} steps", transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            color=INK, bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.set_ylim(y.min() - 6, y.max() + 12)
    ax.legend(loc="upper left", fontsize=6.5)
    save(fig, out, "convergence_monitor")
    return {"mean": mean, "sd_pct": 100 * sd / mean, "half_range_pct": 100 * half / mean, "period_steps": period,
            "steps": [float(x.min()), float(x.max())],
            "implied_HI3_half_range_pct_GW": 100 * ((1 + half / mean) ** GW["alpha"] - 1)}


def c5_mechanism(runs, tab, out):
    """F2-style: HI3 against inflow and against valve angle, all cases, phases joined in order."""
    import operating_point as op
    import report_figures as rf
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.7), constrained_layout=True)
    for r in runs:
        st = style(r.label)
        t = tab[tab.case == r.label].sort_values("phase_pct")
        v = np.array([r.snaps[p].paths.v_inlet.mean() for p in t.phase_pct])
        g = op.geometry(r.scale)
        xh = [0.067, 0.25, 0.5, 0.75, 0.933]
        q = np.array([vi * op.dome(g["S"] * x, g["L"])[0] * 1e6 for vi, x in zip(v, xh)])
        ang = np.array([op.operating_point(p)["valve_angle_deg"] for p in t.phase_pct])
        y = t.GW_HI3_mean.to_numpy()
        for ax, x in ((axs[0], q), (axs[1], ang)):
            rf.smooth_path(ax, x, y, color=st["color"], lw=1.0)
            ax.plot(x, y, ls="none", marker=st["marker"], color=st["color"], ms=5.5, label=r.label, zorder=3)
        for xi, yi, p in zip(q, y, t.phase_pct):
            if p == 50.0:
                axs[0].annotate("50%", (xi, yi), xytext=(4, -8), textcoords="offset points",
                                fontsize=6, color=st["color"])
    for ax in axs:
        ax.set_yscale("log"); logfmt(ax)
        ax.set_yticks([3e-5, 1e-4, 3e-4], ["3e-5", "1e-4", "3e-4"])
        ax.set_ylabel("Mean HI3 [%] (GW)")
    axs[0].set_xlabel("Instantaneous inflow Q [mL/s]")
    axs[0].set_title("a) Against flow: same flow, more damage when smaller")
    axs[1].set_xlabel("Outlet-valve opening angle θ [°]")
    axs[1].set_title("b) Against valve opening")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="outside upper center", ncol=3, fontsize=6.8)
    for ax in axs:
        ax.text(0.98, 0.97, "lines: 16.67 → 83.33%", transform=ax.transAxes, ha="right", va="top",
                fontsize=6.3, color=MUTED)
    save(fig, out, "size_mechanism")


def c6_lysis_envelopes(runs, out):
    """F7-style: worst-case stress-duration envelope, one panel per phase, all cases."""
    curves = ee.lysis_curves()
    phases = runs[0].phases
    fig, axs = plt.subplots(1, len(phases), figsize=(W, 2.35), constrained_layout=True, sharex=True, sharey=True)
    for j, ph in enumerate(phases):
        ax = axs[j]
        for name, ls in (("RBC", "--"), ("Platelet", ":")):
            c = curves[name]
            ax.plot(c.time_s * 1e3, c.stress_Pa, color=INK, ls=ls, lw=1.0,
                    label=("red-cell lysis" if name == "RBC" else "platelet lysis") if j == 0 else None)
        for r in runs:
            st = style(r.label)
            env = ee.worst_envelope(r.all_snaps[ph]).sort_values(["duration_s", "stress_Pa"], ascending=[True, False])
            ax.plot(env.duration_s * 1e3, env.stress_Pa, color=st["color"], ls=st["ls"], lw=1.2,
                    label=r.label if j == 0 else None)
        ax.set_xscale("log"); ax.set_yscale("log"); logfmt(ax); logfmt(ax, "x")
        ax.set_title(f"{ph:g}%")
        ax.set_xlim(3e-3, 5e3); ax.set_ylim(5e-3, 2e3)
        ax.set_xticks([1e-2, 1, 1e2], ["0.01", "1", "100"])
        ax.set_xlabel("Exposure [ms]")
    axs[0].set_ylabel("Shear stress [Pa]")
    fig.legend(loc="outside upper center", ncol=5, fontsize=6.5)
    save(fig, out, "size_lysis_envelopes")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs=4, action="append", metavar=("LABEL", "FOLDER", "BPM", "SCALE"), required=True)
    ap.add_argument("--monitor", type=Path)
    ap.add_argument("--outdir", type=Path, default=HERE / "report_figs")
    ap.add_argument("--only", nargs="+")
    ap.add_argument("--formats", nargs="+", default=list(rs.FORMATS))
    a = ap.parse_args(argv)
    FORMATS[:] = a.formats
    a.outdir.mkdir(parents=True, exist_ok=True)
    want = set(x.upper() for x in (a.only or ["C1", "C2", "C3", "C4", "C5", "C6"]))
    runs = [Run(*c, completed_only=False) for c in a.case]
    tab = summary_table(runs)
    tab.to_csv(a.outdir / "size_summary.csv", index=False)
    if "C1" in want:
        c1_size_comparison(runs, tab, a.outdir)
        runs_c = [Run(*c, completed_only=True) for c in a.case]
        tab_c = summary_table(runs_c)
        tab_c.to_csv(a.outdir / "size_summary_completed_only.csv", index=False)
        c1_size_comparison(runs_c, tab_c, a.outdir, "size_comparison_completed")
    if "C2" in want:
        c2_dose_maps(runs, a.outdir)
    if "C3" in want:
        c3_lysis_sa(runs, tab, a.outdir)
    if "C5" in want:
        c5_mechanism(runs, tab, a.outdir)
    if "C6" in want:
        c6_lysis_envelopes(runs, a.outdir)
    if "C4" in want and a.monitor:
        m = c4_monitor(a.monitor, a.outdir)
        pd.Series(m).to_json(a.outdir / "convergence_monitor.json", indent=1)
        print("monitor", m)
    print(f"done -> {a.outdir}")


if __name__ == "__main__":
    main()
