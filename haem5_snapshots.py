"""Compare haem5 frozen-flow snapshots; originals and damage calculations are reused.

python haem5_snapshots.py z_inputs --empirical-dir outputs_empirical
Omit --empirical-dir to recompute empirical screening from the current exports.
"""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

import haem5
import haem5_empiricalthresholds as empirical

COLORS = ["#568FB5", "#8D85BB", "#B56FA0", "#D78D9E", "#CB9465"]
STYLES = ["-", "--", "-.", ":", (0, (5, 1))]
NOTE = ("Frozen-flow snapshots; exposure time is time along a path, not time within systole.\n"
        "All reconstructed paths; equal-path distribution statistics; bands show population spread.")

def natural_key(path):
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", path.name)]

def phase_from_name(path):
    match = re.search(r"_(\d+(?:\.\d+)?)%$", path.stem)
    if not match:
        raise ValueError("Cannot infer phase from " + path.name + "; supply --phases.")
    value = float(match.group(1))
    return value / 100 if value > 100 else value

def resample_paths(paths, grid, normalized=False):
    """Linear interpolation, no extrapolation; duplicate times use the last node."""
    values = np.full((len(paths), len(grid)), np.nan)
    for i, path in enumerate(paths):
        t, stress = path[:, 0], path[:, 1]
        if (not np.isfinite(path).all() or np.any(np.diff(t) < 0)
                or np.any(stress < 0) or t[-1] <= 0):
            raise ValueError("Paths require finite, nonnegative stress and increasing exposure time.")
        keep = np.r_[np.diff(t) > 0, True]
        t, stress = t[keep], stress[keep]
        if normalized:
            t = t / t[-1]
        values[i] = np.interp(grid, t, stress, left=np.nan, right=np.nan)
    return values

def population_stats(values, grid):
    count = np.isfinite(values).sum(axis=0)
    stats = np.full((5, len(grid)), np.nan)
    active = count > 0
    stats[:4, active] = np.nanpercentile(values[:, active], [10, 50, 90, 99], axis=0)
    stats[4, active] = np.nanmean(values[:, active], axis=0)
    return pd.DataFrame(dict(exposure=grid, n_present=count,
                             fraction_present=count / len(values),
                             p10_Pa=stats[0], median_Pa=stats[1],
                             p90_Pa=stats[2], p99_Pa=stats[3], mean_Pa=stats[4]))

def save(fig, output, name):
    fig.savefig(output / (name + ".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

def plot_sa(cases, output):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    positives = np.concatenate([c["summary"].SA_dyne_s_cm2.to_numpy() for c in cases])
    positives = positives[positives > 0]
    if not len(positives):
        raise ValueError("SA distribution needs at least one positive value.")
    lo, hi = np.log10(positives.min()), np.log10(positives.max())
    if hi - lo < .1:
        lo, hi = lo - .1, hi + .1
    grid = np.linspace(lo - .35, hi + .35, 500)
    for c in cases:
        sa = c["summary"].SA_dyne_s_cm2.to_numpy()
        positive = sa[sa > 0]
        if len(positive):
            # Density per log10 unit, scaled by positive fraction; zero mass stated in table.
            density = haem5._gaussian_kde_numpy(np.log10(positive), grid) * len(positive) / len(sa)
            axes[0].plot(10**grid, density, color=c["color"], ls=c["style"], label=c["label"])
        x = np.sort(sa)
        axes[1].step(x, np.arange(len(x), 0, -1) / len(x) * 100, where="pre",
                     color=c["color"], ls=c["style"], label=c["label"])
    axes[0].set(xscale="log", xlabel="Total path SA [dyne s/cm²]",
                ylabel="Density per log10(SA)", title="Stress accumulation distributions")
    axes[1].set(xscale="symlog", xlim=(0, max(positives) * 1.15),
                xlabel="Total path SA [dyne s/cm²]", ylabel="Paths with SA ≥ x [%]",
                title="Exposure tail (empirical survival curves)")
    for ax in axes:
        ax.axvline(haem5.HELLUMS_THRESHOLD_DYNE_S_CM2, color=".35", ls=":", lw=1)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Equal-path weighting • dotted reference: existing SA platelet criterion")
    save(fig, output, "sa_comparison")

def plot_empirical(cases, output):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.1), layout="constrained")
    for ax in axes:
        for curve, label, ls in ((empirical.RBC_THRESHOLD, "Red-cell lysis reference", "--"),
                                 (empirical.PLATELET_THRESHOLD, "Platelet lysis reference", ":")):
            ax.loglog(curve.time_s * 1000, curve.stress_Pa, color=".3", ls=ls, lw=1.5, label=label)
        for c in cases:
            env = c["envelope"].sort_values(["duration_s", "stress_Pa"], ascending=[True, False])
            ax.loglog(env.duration_s * 1000, env.stress_Pa,
                      color=c["color"], ls=c["style"], lw=2, label=c["label"])
        ax.set(xlabel="Longest contiguous exposure duration [ms]", ylabel="Shear stress [Pa]")
        ax.grid(which="major", alpha=.2)
    envs = pd.concat([c["envelope"] for c in cases])
    axes[0].set(xlim=(1e-3, 1e5), ylim=(min(1, envs.stress_Pa.min()/2), 2e5),
                title="Full empirical reference range")
    axes[1].set(xlim=(envs.duration_s.min()*1000/1.5, envs.duration_s.max()*1000*1.5),
                ylim=(envs.stress_Pa.min()/1.5, envs.stress_Pa.max()*1.5),
                title="Snapshot exposure detail")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Five-phase empirical exposure comparison" if len(cases) == 5 else "Empirical exposure comparison")
    fig.supxlabel("Approximate NIH traces; lysis screening only. Each envelope is a worst case across paths.", fontsize=9)
    save(fig, output, "empirical_stress_duration_all_snapshots")

def plot_histories(cases, output, points):
    max_t = max(p[-1, 0] for c in cases for p in c["paths"])
    # Log-spaced sampling resolves early exposures while retaining the complete residence-time tail.
    positive_t = min(p[p[:, 0] > 0, 0].min() for c in cases for p in c["paths"])
    grid = np.r_[0., np.geomspace(positive_t, max_t, points - 1)]
    all_stress = np.concatenate([p[:, 1] for c in cases for p in c["paths"]])
    positive_stress = all_stress[all_stress > 0]
    floor = positive_stress.min()/2 if len(positive_stress) else 1e-6
    ceiling = max(all_stress.max()*1.2, floor*10)
    bins = np.geomspace(floor, ceiling, 100)
    fig, axes = plt.subplots(len(cases), 1, figsize=(11, 2.65*len(cases)), sharex=True, sharey=True,
                             layout="constrained", squeeze=False)
    cloud, cloud_axes = plt.subplots(len(cases), 1, figsize=(11, 2.4*len(cases)), sharex=True, sharey=True,
                                     layout="constrained", squeeze=False)
    overlay, oa = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout="constrained",
                               gridspec_kw={"height_ratios": [3, 1]})
    normalized, na = plt.subplots(figsize=(10, 4.8), layout="constrained")
    norm = LogNorm(vmin=1/max(len(c["paths"]) for c in cases), vmax=1)
    for c, ax, ca in zip(cases, axes[:, 0], cloud_axes[:, 0]):
        matrix = resample_paths(c["paths"], grid)
        stats = population_stats(matrix, grid)
        stats.rename(columns={"exposure": "exposure_s"}).to_csv(output / (c["stem"]+"_shear_statistics.csv"), index=False)
        density = np.zeros((len(bins)-1, len(grid)))
        for j in range(len(grid)):
            finite = matrix[:, j][np.isfinite(matrix[:, j])]
            # Zero stress is displayed at the bottom log bin, rather than discarded.
            density[:, j] = np.histogram(np.maximum(finite, floor), bins=bins)[0] / len(matrix)
        time_edges = np.r_[0, (grid[:-1] + grid[1:])/2, grid[-1] + (grid[-1]-grid[-2])/2]
        mesh = ax.pcolormesh(time_edges*1000, bins, np.ma.masked_equal(density, 0),
                             cmap="Purples", norm=norm, shading="auto", rasterized=True)
        ax.fill_between(grid*1000, np.maximum(stats.p10_Pa, floor), np.maximum(stats.p90_Pa, floor),
                        color=c["color"], alpha=.15, label="P10–P90")
        for key, ls, label in (("median_Pa", "-", "Median"), ("p99_Pa", "--", "P99"), ("mean_Pa", ":", "Mean")):
            ax.plot(grid*1000, np.maximum(stats[key], floor), color=c["color"], ls=ls, lw=1.3, label=label)
        ax.set_title(c["label"] + f" | {len(matrix):,} paths", loc="left", fontsize=10)
        ax.set(yscale="log", ylabel="Shear [Pa]", ylim=(floor, ceiling))
        ax.legend(loc="upper right", fontsize=7, ncol=4)
        segments = [np.column_stack((p[:, 0]*1000, np.maximum(p[:, 1], floor))) for p in c["paths"]]
        ca.add_collection(LineCollection(segments, colors=c["color"], linewidths=.35, alpha=.025, rasterized=True))
        ca.plot(grid*1000, np.maximum(stats.median_Pa, floor), color=".15", lw=1.2, label="Median")
        ca.set(yscale="log", ylim=(floor, ceiling), ylabel="Shear [Pa]")
        ca.set_title(c["label"] + " | all paths, opacity 0.025", loc="left", fontsize=10)
        ca.legend(fontsize=8)
        oa[0].plot(grid*1000, stats.median_Pa, color=c["color"], ls=c["style"], label=c["label"]+" median")
        oa[0].plot(grid*1000, stats.mean_Pa, color=c["color"], lw=.8, alpha=.65)
        oa[1].plot(grid*1000, stats.fraction_present*100, color=c["color"], ls=c["style"], label=c["label"])
        ng = np.linspace(0, 1, points)
        ns = population_stats(resample_paths(c["paths"], ng, normalized=True), ng)
        ns.rename(columns={"exposure": "transit_fraction"}).to_csv(output/(c["stem"]+"_normalized_shear.csv"), index=False)
        na.plot(ng*100, ns.median_Pa, color=c["color"], ls=c["style"], label=c["label"])
    for ax in [*axes[:, 0], *cloud_axes[:, 0], *oa]:
        ax.set_xscale("symlog", linthresh=positive_t*1000)
        ax.set_xlim(0, max_t*1000)
        ax.grid(alpha=.15)
    axes[-1, 0].set_xlabel("Elapsed exposure time along path [ms] (linear near zero, logarithmic thereafter)")
    cloud_axes[-1, 0].set_xlabel("Elapsed exposure time along path [ms]")
    fig.colorbar(mesh, ax=list(axes[:, 0]), label="Fraction of all seeded paths per shear bin", shrink=.7)
    fig.suptitle("Shear history populations\n"+NOTE, fontsize=11)
    fig.supxlabel("Zero shear is shown at the lowest log bin. Interpolation is for visualization only.", fontsize=9)
    cloud.suptitle("All reconstructed shear histories\n"+NOTE, fontsize=11)
    oa[0].set_yscale("symlog", linthresh=.01)
    oa[0].set_ylim(bottom=0)
    oa[0].set(ylabel="Shear stress [Pa]", title="Medians (labelled styles) and means (thin solid lines)")
    oa[0].legend(fontsize=8)
    oa[1].set(ylabel="Paths present [%]", xlabel="Elapsed exposure time along path [ms]", ylim=(0, 105))
    overlay.suptitle("Population comparison • ended paths excluded; late-time statistics describe remaining paths")
    na.set_yscale("symlog", linthresh=.01)
    na.set_ylim(bottom=0)
    na.set(xlabel="Progress through each path's transit [%]", ylabel="Median shear [Pa]",
           title="Normalized transit comparison • absolute exposure duration is removed")
    na.legend(fontsize=8)
    na.grid(alpha=.2)
    for f, name in ((fig, "shear_density_panels"), (cloud, "shear_all_paths"),
                    (overlay, "shear_population_overlay"), (normalized, "shear_normalized_transit")):
        save(f, output, name)

def plot_phase(table, output):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), layout="constrained")
    for model, color, marker in (("HI2", COLORS[0], "o"), ("HI3", COLORS[2], "s")):
        axes[0].plot(table.phase_percent, table[model+"_percent"], marker=marker, color=color, label=model)
    axes[0].set(ylabel="Inlet-speed-weighted HI [%]", title="Predicted haemolysis")
    axes[0].legend()
    for key, label, style in (("sa_median", "Median", "-"), ("sa_p95", "P95", "--"), ("sa_p99", "P99", ":")):
        axes[1].plot(table.phase_percent, table[key], marker="o", ls=style, label=label, color=COLORS[1])
    axes[1].set(ylabel="Path SA [dyne s/cm²]", title="Stress accumulation (equal-path)")
    axes[1].legend()
    for key, label, style in (("peak_median", "Median", "-"), ("peak_p99", "P99", "--")):
        axes[2].plot(table.phase_percent, table[key], marker="o", ls=style, label=label, color=COLORS[2])
    axes[2].set(ylabel="Per-path peak stress [Pa]", title="Peak stresses (equal-path)")
    axes[2].legend()
    for ax in axes:
        ax.set_xlabel("Systolic phase [%]")
        ax.set_xticks(table.phase_percent)
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(alpha=.2)
    fig.suptitle("Snapshot trends • individual path damage calculated before averaging")
    save(fig, output, "snapshot_phase_trends")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=Path, default=Path("z_inputs"))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--phases", nargs="+", type=float, help="Percent phases in natural filename order; otherwise infer 1667%% as 16.67%%.")
    parser.add_argument("--empirical-dir", type=Path, help="Reuse newest saved envelope for each exact input stem. User must ensure exports match that run.")
    parser.add_argument("--points", type=int, default=400, help="Visualization samples per path; damage is computed on original nodes.")
    args = parser.parse_args()
    files = sorted(args.folder.glob("*.csv"), key=natural_key)
    if len(files) < 2:
        parser.error("Need at least two snapshot CSVs.")
    if args.points < 20:
        parser.error("--points must be at least 20.")
    phases = args.phases if args.phases is not None else [phase_from_name(p) for p in files]
    if len(phases) != len(files) or not np.isfinite(phases).all() or any(p < 0 or p > 100 for p in phases) or np.any(np.diff(phases) <= 0):
        parser.error("Supply one strictly increasing phase in [0, 100] per file.")
    output = args.output_dir or Path("outputs_snapshots") / datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
    output.mkdir(parents=True, exist_ok=False)
    record = {"status": "running", "files": [], "phases_percent": phases,
              "constants": haem5.ACTIVE_CONSTANTS, "points": args.points,
              "population_weighting": "equal path", "HI_weighting": "inlet speed proxy"}
    cases, rows = [], []
    try:
        for i, (path, phase) in enumerate(zip(files, phases)):
            print(f"[{i+1}/{len(files)}] {path.name} ({phase:g}% systole)", flush=True)
            result = haem5.run_pipeline(path)
            summary = result["summary"]
            summary.to_csv(output/(path.stem+"_streamlines.csv"), index=False)
            source = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            if args.empirical_dir:
                candidates = sorted(p for p in args.empirical_dir.glob("empirical_exposure_envelope_*.csv")
                                    if p.name.endswith("_"+path.stem+".csv"))
                if not candidates:
                    raise FileNotFoundError("No saved empirical envelope for "+path.name)
                source["empirical_envelope"] = str(candidates[-1].resolve())
                source["empirical_sha256"] = hashlib.sha256(candidates[-1].read_bytes()).hexdigest()
                source["empirical_input_match"] = "User-supplied cache; original input hash unavailable"
                envelope = pd.read_csv(candidates[-1])
                print("  Reusing "+candidates[-1].name, flush=True)
            else:
                print("  Computing empirical screening...", flush=True)
                lines, device, envelope = empirical.analyze_empirical_thresholds(result)
                lines.to_csv(output/(path.stem+"_empirical_paths.csv"), index=False)
                device.to_csv(output/(path.stem+"_empirical_device.csv"), index=False)
            if envelope.empty or not np.isfinite(envelope[["duration_s", "stress_Pa"]]).all().all() or (envelope[["duration_s", "stress_Pa"]] <= 0).any().any():
                raise ValueError("Empirical envelope must contain finite positive durations and stresses.")
            envelope.to_csv(output/(path.stem+"_empirical_envelope.csv"), index=False)
            record["files"].append(source)
            paths = [df[["t_cumulative", "shear"]].to_numpy() for df in result["streamlines"]]
            cases.append(dict(stem=path.stem, label=f"{phase:g}% systole", color=COLORS[i % len(COLORS)],
                              style=STYLES[i % len(STYLES)], summary=summary, paths=paths, envelope=envelope))
            rows.append(dict(file=path.name, phase_percent=phase, n_paths=len(paths),
                             HI2_percent=result["device_HI2_percent"], HI3_percent=result["device_HI3_percent"],
                             sa_median=summary.SA_dyne_s_cm2.median(), sa_p95=summary.SA_dyne_s_cm2.quantile(.95),
                             sa_p99=summary.SA_dyne_s_cm2.quantile(.99), sa_zero_fraction=float((summary.SA_dyne_s_cm2 == 0).mean()),
                             peak_median=summary.max_stress_Pa.median(), peak_p99=summary.max_stress_Pa.quantile(.99)))
            del result
        table = pd.DataFrame(rows)
        table.to_csv(output/"snapshot_summary.csv", index=False)
        print("Rendering combined figures...", flush=True)
        plot_sa(cases, output)
        plot_empirical(cases, output)
        plot_histories(cases, output, args.points)
        plot_phase(table, output)
        record["status"] = "complete"
        (output/"README.md").write_text(
            "# Snapshot comparison\n\n"+NOTE+"\n\n"
            "Seven PNGs compare SA distributions/tails, empirical envelopes (full and detail), "
            "shear density, all raw path histories, mean/median and population support, "
            "normalized transit medians, and phase trends. CSVs retain the numerical summaries.\n\n"
            "Phases inferred from filename suffixes (1667% = 16.67%); override with --phases. "
            "No full-cycle integration or transient particle tracking is performed.\n\n"
            "HI uses the unchanged haem5 model and inlet-speed weighting (equal-area inlet seeding assumption). "
            "All reconstructed paths are included; outlet completion is not checked. "
            "Shear interpolation does not change the damage calculations. Ended paths are excluded. "
            "Percentiles describe the population, not confidence intervals.\n\n"
            "SA density is per log10 unit; zero-SA fraction is in snapshot_summary.csv. "
            "Density panels show fractions of the ORIGINAL path population per shared log shear bin, "
            "so colour intensity also falls as paths end. Zero shear is placed at the lowest log bin.\n\n"
            "Empirical envelopes use approximate NIH lysis references, not activation thresholds or a "
            "validated variable-stress damage model. --empirical-dir trusts the supplied saved results; "
            "run.json records their exact filenames and hashes. Omit it to recalculate.\n",
            encoding="utf-8")
    except Exception as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        (output/"run.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print("Saved: "+str(output.resolve()), flush=True)

if __name__ == "__main__":
    main()
