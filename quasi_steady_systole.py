"""Combine frozen-flow CFD-Post snapshots into a quasi-steady systolic estimate.

Examples:
    python quasi_steady_systole.py snapshots
    python quasi_steady_systole.py snapshots --sampling midpoints
    python quasi_steady_systole.py snapshots --flow-rates-ml-s 10 120 200 120 10

Default: 250 ms systole, snapshots including both endpoints, natural filename
order. Endpoint quadrature gives the first/last snapshots half the interior
weight in time. Midpoints instead represent equal-duration bins. Use --files
to override order and --times-ms for explicit endpoint-inclusive sample times.

Without supplied flow rates, mean inlet speed is a relative Q proxy, assuming
constant inlet area across snapshots. --inlet-area-m2 gives absolute Q instead.
All modes retain haem5's equal-area, forward, normal inlet seeding assumption
for weighting paths within each snapshot. Seed count does not set snapshot Q.

These are weighted frozen-flow path estimates, NOT transient particle histories
or total haemoglobin mass released per beat. Path endpoints are not checked for
outlet arrival. haem5's velocity floor and damage-model assumptions still apply.
"""

import argparse
from datetime import datetime
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import haem5

METRICS = ("HI2_percent", "HI3_percent", "transit_ms", "SA_dyne_s_cm2", "max_stress_Pa")


def natural_key(path):
    return [float(p) if re.fullmatch(r"\d+(?:\.\d+)?", p) else p.casefold()
            for p in re.split(r"(\d+(?:\.\d+)?)", path.name)]


def snapshot_timing(count, systole_ms=250.0, sampling="endpoints", times_ms=None):
    """Return sample times and quadrature durations, both in milliseconds."""
    if count < 2:
        raise ValueError("Provide at least two snapshots from the same systole.")
    if not np.isfinite(systole_ms) or systole_ms <= 0:
        raise ValueError("Systole duration must be finite and positive.")
    if sampling not in ("endpoints", "midpoints"):
        raise ValueError("Unknown sampling convention.")
    if times_ms is not None:
        if sampling != "endpoints":
            raise ValueError("--times-ms requires endpoint sampling.")
        times = np.asarray(times_ms, dtype=float)
        if (times.shape != (count,) or not np.all(np.isfinite(times))
                or np.any(np.diff(times) <= 0)):
            raise ValueError("Supply one finite, strictly increasing time per snapshot.")
        if not np.isclose(times[0], 0) or not np.isclose(times[-1], systole_ms):
            raise ValueError("Explicit times must include 0 and the systole duration.")
        times = times.copy()
        times[0], times[-1] = 0.0, systole_ms
    elif sampling == "midpoints":
        times = (np.arange(count) + 0.5) * systole_ms / count
    else:
        times = np.linspace(0, systole_ms, count)
    edges = np.r_[0, (times[:-1] + times[1:]) / 2, systole_ms]
    return times, np.diff(edges)


def summarize_paths(summary, systole_ms):
    """Apply inlet-speed weights consistently to every reported path metric."""
    speeds = summary["mean_inlet_velocity"].to_numpy(dtype=float)
    values = summary[list(METRICS)].to_numpy(dtype=float)
    if (len(speeds) == 0 or not np.all(np.isfinite(speeds))
            or np.any(speeds < 0)):
        raise ValueError("Snapshot needs finite nonnegative inlet speeds.")
    if speeds.sum() == 0:
        row = {key: float("nan") for key in METRICS}
        row.update(SA_fraction_above_threshold=float("nan"),
                   fraction_transit_longer_than_systole=float("nan"),
                   mean_inlet_speed_m_s=0.0, n_streamlines=len(summary))
        return row, np.zeros(len(speeds))
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Snapshot contains invalid path metrics.")
    weights = speeds / speeds.sum()
    row = dict(zip(METRICS, weights @ values))
    row["SA_fraction_above_threshold"] = float(weights @ (
        summary["SA_dyne_s_cm2"].to_numpy() > haem5.HELLUMS_THRESHOLD_DYNE_S_CM2))
    row["fraction_transit_longer_than_systole"] = float(weights @ (
        summary["transit_ms"].to_numpy() > systole_ms))
    row["mean_inlet_speed_m_s"] = float(speeds.mean())
    row["n_streamlines"] = len(summary)
    return row, weights


def combine_snapshots(snapshots, flow):
    """Normalize Q_n * duration_n; HI is averaged after each path's damage calculation."""
    flow = np.asarray(flow, dtype=float)
    durations = snapshots["represented_ms"].to_numpy(dtype=float)
    if (flow.shape != durations.shape or not np.all(np.isfinite(flow))
            or np.any(flow < 0) or not np.all(np.isfinite(durations))
            or np.any(durations <= 0)):
        raise ValueError("Each snapshot needs a finite nonnegative flow and positive duration.")
    volume = flow * durations / 1000.0
    if volume.sum() <= 0:
        raise ValueError("Total represented flow volume must be positive.")
    weights = volume / volume.sum()
    columns = list(METRICS) + ["SA_fraction_above_threshold",
                              "fraction_transit_longer_than_systole"]
    active = weights > 0
    values = snapshots.loc[active, columns].to_numpy(dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("A positive-flow snapshot has undefined metrics (check inlet seed speeds).")
    result = dict(zip(columns, (weights[active] @ values).tolist()))
    for model in ("HI2", "HI3"):
        result[model + "_NIH_mg_per_100L"] = (
            result[model + "_percent"] * (1 - haem5.HCT) * haem5.HB * 1000)
    return weights, result


def plot_snapshots(snapshots, combined, output):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    times = snapshots["time_ms"]
    for col, color, marker in (("HI2_percent", "#9678B5", "o"),
                               ("HI3_percent", "#C7789B", "s")):
        axes[0, 0].plot(times, snapshots[col], marker=marker, color=color, label=col[:3])
        axes[0, 0].axhline(combined[col], color=color, linestyle="--", alpha=0.8)
    axes[0, 0].set_ylabel("Flow-weighted HI (%)")
    axes[0, 0].legend()
    for ax, col, label, color in (
            (axes[0, 1], "transit_ms", "Flow-weighted path transit (ms)", "#729EC2"),
            (axes[1, 0], "SA_dyne_s_cm2", "Flow-weighted SA (dyne s/cm²)", "#9678B5")):
        ax.plot(times, snapshots[col], "o-", color=color)
        ax.axhline(combined[col], color=color, linestyle="--")
        ax.set_ylabel(label)
    axes[1, 1].bar(times, snapshots["systole_weight"],
                   width=np.min(np.diff(times)) * 0.65, color="#C7789B", edgecolor="#75556A")
    axes[1, 1].set_ylabel("Snapshot share of represented volume")
    for ax in axes.flat:
        ax.set_xlabel("Time within systole (ms)")
        ax.grid(alpha=0.2)
    fig.suptitle("Quasi-steady systolic estimates\nDashed lines: combined volume-weighted means")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", type=Path, help="Folder containing only snapshots of one systole")
    parser.add_argument("--pattern", default="*.csv")
    parser.add_argument("--files", nargs="+", help="Snapshot filenames in chronological order, relative to folder")
    parser.add_argument("--systole-ms", type=float, default=250.0)
    parser.add_argument("--sampling", choices=["endpoints", "midpoints"], default="endpoints")
    parser.add_argument("--times-ms", nargs="+", type=float)
    flow_group = parser.add_mutually_exclusive_group()
    flow_group.add_argument("--flow-rates-ml-s", nargs="+", type=float,
                            help="Actual forward flow per snapshot, in the displayed file order")
    flow_group.add_argument("--inlet-area-m2", type=float, help="Same inlet area for every snapshot")
    parser.add_argument("--constants", choices=list(haem5.POWER_LAW_CONSTANTS), default="GW")
    parser.add_argument("--outdir", type=Path, help="New output directory (existing directories refused)")
    args = parser.parse_args(argv)
    try:
        if not args.folder.is_dir():
            raise ValueError(f"Snapshot folder does not exist: {args.folder}")
        files = ([args.folder / name for name in args.files] if args.files else
                 sorted((p for p in args.folder.glob(args.pattern) if p.is_file()), key=natural_key))
        if len(set(p.resolve() for p in files)) != len(files):
            raise ValueError("Duplicate snapshot files are not allowed.")
        for path in files:
            if not path.is_file():
                raise ValueError(f"Missing snapshot: {path}")
        times, durations = snapshot_timing(len(files), args.systole_ms, args.sampling, args.times_ms)
        if args.inlet_area_m2 is not None and (
                not np.isfinite(args.inlet_area_m2) or args.inlet_area_m2 <= 0):
            raise ValueError("Inlet area must be finite and positive.")
        if args.flow_rates_ml_s is not None:
            q = np.asarray(args.flow_rates_ml_s)
            if (len(q) != len(files) or not np.all(np.isfinite(q))
                    or np.any(q < 0) or q.sum() <= 0):
                raise ValueError("Supply one nonnegative finite flow per file, with positive total flow.")
        outdir = args.outdir or Path("outputs") / ("quasi_steady_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f"))
        if outdir.exists():
            raise ValueError(f"Output directory already exists: {outdir}")
    except ValueError as exc:
        parser.error(str(exc))

    print("Snapshot order and represented durations:")
    for path, time, duration in zip(files, times, durations):
        print(f"  {path.name}: t={time:g} ms, represented={duration:g} ms")
    print("Assumes equal-area inlet seeds and normal forward flow within each snapshot.")
    constants = haem5.POWER_LAW_CONSTANTS[args.constants]
    rows, path_tables = [], []
    for index, (path, time, duration) in enumerate(zip(files, times, durations), 1):
        print(f"Processing {index}/{len(files)}: {path.name}")
        try:
            with np.errstate(divide="ignore", invalid="ignore"):
                result = haem5.run_pipeline(path, **{k: constants[k] for k in ("A", "alpha", "beta")})
            row, weights = summarize_paths(result["summary"], args.systole_ms)
        except (ValueError, RuntimeError, KeyError, IndexError) as exc:
            parser.error(f"{path.name}: {exc}")
        rows.append({"snapshot": index, "file": path.name, "time_ms": time,
                     "represented_ms": duration, **row})
        table = result["summary"].copy()
        table.insert(0, "snapshot", index)
        table["within_snapshot_weight"] = weights
        path_tables.append(table)
        del result
    snapshots = pd.DataFrame(rows)
    if args.flow_rates_ml_s is not None:
        flow = np.asarray(args.flow_rates_ml_s)
        basis = "supplied snapshot forward flow (mL/s)"
    elif args.inlet_area_m2 is not None:
        flow = snapshots["mean_inlet_speed_m_s"].to_numpy() * args.inlet_area_m2 * 1e6
        basis = "inlet area times mean inlet speed (mL/s)"
    else:
        flow = snapshots["mean_inlet_speed_m_s"].to_numpy()
        basis = "relative Q from mean inlet speed; assumes constant inlet area"
    try:
        weights, combined = combine_snapshots(snapshots, flow)
    except ValueError as exc:
        parser.error(str(exc))
    absolute_flow = args.flow_rates_ml_s is not None or args.inlet_area_m2 is not None
    snapshots["flow_ml_s" if absolute_flow else "relative_flow_proxy_m_s"] = flow
    snapshots["systole_weight"] = weights
    combined.update(systole_ms=args.systole_ms, n_snapshots=len(files), constants=args.constants,
                    flow_basis=basis, sampling=args.sampling)
    if absolute_flow:
        snapshots["represented_volume_ml"] = flow * durations / 1000.0
        combined["represented_systolic_volume_ml"] = float(snapshots["represented_volume_ml"].sum())
    pooled = pd.concat(path_tables, ignore_index=True)
    pooled["systole_weight"] = pooled["within_snapshot_weight"] * pooled["snapshot"].map(
        dict(enumerate(weights, 1)))
    outdir.mkdir(parents=True, exist_ok=False)
    snapshots.to_csv(outdir / "snapshot_summary.csv", index=False)
    pd.DataFrame([combined]).to_csv(outdir / "systole_summary.csv", index=False)
    pooled.to_csv(outdir / "weighted_streamlines.csv", index=False)
    metadata = {"inputs": [str(p.resolve()) for p in files], "flow_basis": basis,
                "sampling": args.sampling, "times_ms": times.tolist(),
                "represented_ms": durations.tolist(), "constants": constants,
                "HCT": haem5.HCT, "HB_g_L": haem5.HB,
                "SA_threshold_dyne_s_cm2": haem5.HELLUMS_THRESHOLD_DYNE_S_CM2,
                "interpretation": "Volume-weighted frozen-flow path estimates, not transient histories or total damage per beat.",
                "assumptions": ["Equal-area, forward normal inlet seeding; first path node is inlet.",
                                "All exported paths included; outlet arrival is not verified.",
                                "haem5 minimum-velocity floor and damage models apply.",
                                "SA threshold exceedance is a model exposure statistic, not a clinical probability."]}
    (outdir / "run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    plot_snapshots(snapshots, combined, outdir / "systole_summary.png")
    print(f"\nFlow basis: {basis}")
    for metric in list(METRICS) + ["HI2_NIH_mg_per_100L", "HI3_NIH_mg_per_100L"]:
        print(f"  {metric}: {combined[metric]:.6g}")
    print(f"Volume-weighted fraction with transit > systole: {combined['fraction_transit_longer_than_systole']:.2%}")
    print("Quasi-steady estimates only: weighting does not reconstruct changing particle histories.")
    print(f"Saved results: {outdir.resolve()}")
    return outdir


if __name__ == "__main__":
    main()
