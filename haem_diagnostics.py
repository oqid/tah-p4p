"""Audit CFD-Post streamline exports using haem5's power-law conventions.

Run `python haem_diagnostics.py --help`. See HAEM_DIAGNOSTICS.md for scope,
manifest format, outlet classification, weighting and research interpretation.
Requires numpy/pandas and haem5's dependencies. No CFD solver is run.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from haem5 import POWER_LAW_CONSTANTS, HCT, HB


def read_export(path):
    """Strict numeric Line export parser; preserve connectivity order and IDs."""
    rows, pairs = [], []
    section = None
    header = None
    with open(path, encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith("["):
                section = line.strip("[]")
                if section == "Faces":
                    raise ValueError("Export plain Lines, not ribbons/tubes ([Faces]); orientation is ambiguous.")
                continue
            if section == "Data":
                if header is None:
                    header = [s.strip() for s in line.split(",")]
                    continue
                rows.append([float(v) for v in line.replace(",", " ").split()])
            elif section == "Lines":
                values = [float(v) for v in line.replace(",", " ").split()]
                if len(values) != 2 or any(v != int(v) for v in values):
                    raise ValueError("Expected two integer node IDs per [Lines] row.")
                pairs.append(tuple(map(int, values)))
    expected = ["Node Number", "X [ m ]", "Y [ m ]", "Z [ m ]",
                "varshear [ Pa ]", "Velocity [ m s^-1 ]"]
    norm = lambda s: "".join(s.lower().split())
    if not header or list(map(norm, header[:6])) != list(map(norm, expected)):
        raise ValueError("Expected haem5 export columns: " + ", ".join(expected))
    if not rows or not pairs:
        raise ValueError("Empty data or connectivity.")
    data = np.asarray(rows, dtype=float)
    if data.shape[1] != len(header) or not np.isfinite(data).all():
        raise ValueError("Inconsistent columns or nonfinite data.")
    ids = data[:, 0].astype(np.int64)
    if not np.array_equal(ids, data[:, 0]) or len(set(ids)) != len(ids):
        raise ValueError("Node IDs must be unique integers.")
    if (data[:, 4:6] < 0).any():
        raise ValueError("Negative stress or speed.")
    nodes = {int(i): r[1:6] for i, r in zip(ids, data)}
    paths = []
    current = []
    for a, b in pairs:
        if a not in nodes or b not in nodes:
            raise ValueError("Connectivity references an absent node.")
        if current and a == current[-1]:
            current.append(b)
        else:
            if current:
                paths.append(current)
            current = [a, b]
    paths.append(current)
    return nodes, paths


def segments(points, floor=True, midpoint=False):
    length = np.linalg.norm(np.diff(points[:, :3], axis=0), axis=1)
    speeds = points[:, 4]
    speed = (speeds[:-1] + speeds[1:]) / 2
    positive = speeds[speeds > 0]
    limit = max((positive.mean() if len(positive) else .1) * .01, .001)
    clipped = (speed < limit) & (length > 0)
    if not floor and np.any((speed == 0) & (length > 0)):
        return None, None, clipped
    denominator = np.maximum(speed, limit) if floor else speed
    dt = np.divide(length, denominator, out=np.zeros_like(length), where=denominator > 0)
    tau = (points[:-1, 3] + points[1:, 3]) / 2 if midpoint else points[1:, 3]
    return dt, tau, clipped


def damage(dt, tau, constants):
    c, a, b = (constants[k] for k in ("A", "alpha", "beta"))
    t = np.cumsum(dt)
    factor = np.zeros_like(t)
    np.power(t, b - 1, out=factor, where=t > 0)
    hi2 = np.sum(b * c * factor * tau ** a * dt)
    exact = np.sum(c * tau ** a * np.diff(np.r_[0., t ** b]))
    source = c ** (1 / b) * tau ** (a / b) * dt
    return float(hi2), float(source.sum() ** b), float(exact), source


def describe(values):
    return dict(zip(("min", "p50", "p90", "p95", "p99", "max"),
                    map(float, np.quantile(values, [0, .5, .9, .95, .99, 1]))))


def weighted(values, weights):
    # Do not silently discard undefined results from the selected population.
    values, weights = np.asarray(values), np.asarray(weights)
    active = weights > 0
    if not active.any() or not np.isfinite(values[active]).all():
        return float("nan")
    return float(np.average(values[active], weights=weights[active]))


def analyze(path, config, destination):
    nodes, ids = read_export(path)
    n = len(ids)
    weights = None
    if config.get("weights"):
        table = pd.read_csv(config["weights"])
        if not {"streamline_id", "flux_weight"}.issubset(table.columns):
            raise ValueError("Weights need streamline_id,flux_weight columns.")
        if table.streamline_id.duplicated().any() or set(table.streamline_id) != set(range(n)):
            raise ValueError("Weights must match every reconstructed zero-based streamline ID exactly.")
        weights = table.set_index("streamline_id").loc[range(n), "flux_weight"].to_numpy(float)
        if not np.isfinite(weights).all() or (weights < 0).any() or weights.sum() <= 0:
            raise ValueError("Flux weights must be finite, nonnegative and have a positive sum.")
    outlet = config.get("outlet_plane")
    if outlet:
        plane = np.asarray(outlet, dtype=float)
        if plane.shape != (6,) or not np.isfinite(plane).all() or np.linalg.norm(plane[3:]) == 0:
            raise ValueError("Outlet plane requires x,y,z,nx,ny,nz with nonzero outward normal.")
        origin, normal = plane[:3], plane[3:] / np.linalg.norm(plane[3:])
    records, source_bins, hot = [], np.zeros(4), []
    all_points = np.array(list(nodes.values()))
    floor_segments = 0
    for sid, path_ids in enumerate(ids):
        p = np.array([nodes[i] for i in path_ids])
        dt, tau, clipped = segments(p)
        floor_segments += int(clipped.sum())
        if dt.sum() <= 0:
            raise ValueError(f"Streamline {sid} has zero total travel time.")
        weight = float(weights[sid] if weights is not None else p[0, 4])
        status = "unverified"
        if outlet:
            distances = (p[:, :3] - origin) @ normal
            near_end = abs(distances[-1]) <= config.get("outlet_tolerance", 1e-4)
            forward = distances[-1] > distances[-2]
            status = "outlet_plane_reached" if near_end and forward and distances[0] < 0 else "incomplete_or_other_exit"
        record = {"streamline_id": sid, "start_node_id": path_ids[0], "end_node_id": path_ids[-1],
                  "n_points": len(p), "weight": weight, "status": status,
                  "transit_s": dt.sum(), "length_m": np.linalg.norm(np.diff(p[:, :3], axis=0), axis=1).sum(),
                  "floor_segments": int(clipped.sum()), "min_speed_m_s": p[:, 4].min(),
                  "max_stress_Pa": tau.max(), "SA_dyne_s_cm2": 10 * np.sum(tau * dt),
                  "repeated_node_ids": len(path_ids) - len(set(path_ids))}
        for name, point in (("start", p[0]), ("end", p[-1])):
            record.update({f"{name}_{axis}_m": float(v) for axis, v in zip("xyz", point[:3])})
        for name, constants in POWER_LAW_CONSTANTS.items():
            hi2, hi3, exact, source = damage(dt, tau, constants)
            record.update({f"{name}_HI2_percent": hi2, f"{name}_HI3_percent": hi3,
                           f"{name}_HI2_exact_percent": exact,
                           f"{name}_time_fraction_tau_above_fit": np.sum(dt[tau > constants['tau_max_Pa']]) / dt.sum(),
                           f"{name}_transit_above_fit": dt.sum() > constants['t_max_s']})
            if name == "TZ":
                record["TZ_time_fraction_tau_below_50Pa"] = np.sum(dt[tau < 50]) / dt.sum()
        gw = POWER_LAW_CONSTANTS["GW"]
        raw_dt, raw_tau, _ = segments(p, floor=False)
        record["GW_HI3_no_floor_percent"] = damage(raw_dt, raw_tau, gw)[1] if raw_dt is not None else np.nan
        mid_dt, mid_tau, _ = segments(p, midpoint=True)
        record["GW_HI3_midpoint_percent"] = damage(mid_dt, mid_tau, gw)[1]
        record["GW_HI2_reversed_percent"] = damage(dt[::-1], tau[::-1], gw)[0]
        for stride in (2, 4):
            selected = np.unique(np.r_[np.arange(0, len(p), stride), len(p) - 1])
            coarse_dt, coarse_tau, _ = segments(p[selected])
            record[f"GW_HI3_stride{stride}_percent"] = damage(coarse_dt, coarse_tau, gw)[1]
        _, _, _, source = damage(dt, tau, gw)
        if status != "incomplete_or_other_exit":
            source_bins += weight * np.histogram(np.cumsum(dt) / dt.sum(), [0, .25, .5, .75, 1], weights=source)[0]
            # Save highest source segments per trajectory for locating hotspots.
            for k in np.argsort(source)[-3:]:
                hot.append({"streamline_id": sid, "segment_index": int(k),
                            **dict(zip(("x_m", "y_m", "z_m"), (p[k, :3] + p[k+1, :3]) / 2)),
                            "tau_Pa": tau[k], "dt_s": dt[k], "time_fraction": np.cumsum(dt)[k] / dt.sum(),
                            "weighted_linearized_source": weight * source[k]})
        records.append(record)
    frame = pd.DataFrame(records)
    if frame.weight.sum() <= 0:
        raise ValueError("Zero inlet-speed weights; supply physical flux weights.")
    selected = frame[frame.status == "outlet_plane_reached"] if outlet else frame
    warnings = ["Viscous stress only; SST Reynolds stresses are not included in exported varshear.",
                "Steady snapshots do not reconstruct transient trajectories or full-cycle haemolysis.",
                "Mesh convergence, boundary-layer resolution, mass balance and physiological validation require solver data.",
                "Coefficient-range diagnostics are extrapolation flags, not evidence of biological validity."]
    if weights is None:
        warnings.append("Weights are inlet SPEED proxies; equal represented area and normal flow are unverified.")
    if not outlet:
        warnings.append("No outlet plane: all-export averages include potentially incomplete paths; not verified device HI.")
    else:
        warnings.append("Outlet test is geometric only: infinite plane, no aperture or velocity-vector verification; paths must end at the intended physical outlet.")
    expected = config.get("expected_seeds")
    if expected is not None and expected != n:
        warnings.append(f"Expected {expected} seeds; reconstructed {n} paths. Missing release flux is unknown; count difference is not a recirculation measurement.")
    if floor_segments:
        warnings.append(f"haem5 velocity floor affects {floor_segments} segments; baseline exposure is artificially shortened there.")
    if len(selected) == 0:
        warnings.append("No completed paths: completed-pass HI is undefined.")
    report = {"input": str(path.resolve()), "metadata": config,
              "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "n_nodes": len(nodes), "n_streamlines": n,
              "weighting": "supplied_flux" if weights is not None else "inlet_speed_proxy",
              "population": "outlet_plane_reached" if outlet else "all_exported_unverified",
              "selected_count": len(selected),
              "selected_weight_fraction_of_exported": float(selected.weight.sum() / frame.weight.sum()),
              "missing_seed_count": expected - n if expected is not None else None,
              "stress_Pa_point_sampled": describe(all_points[:, 3]),
              "speed_m_s_point_sampled": describe(all_points[:, 4]),
              "transit_s_all_paths": describe(frame.transit_s),
              "coordinate_extents_m": dict(zip("xyz", np.ptp(all_points[:, :3], axis=0).tolist())),
              "floor_segments": floor_segments, "warnings": warnings}
    metrics = {}
    for col in frame.columns:
        if "_percent" in col or "_time_fraction_" in col or "_transit_above_fit" in col:
            metrics[col] = weighted(selected[col], selected.weight)
    for name in POWER_LAW_CONSTANTS:
        for model in ("HI2", "HI3"):
            metrics[f"{name}_{model}_NIH_mg_per_100L"] = metrics[f"{name}_{model}_percent"] * (1-HCT)*HB*1000
    metrics["GW_HI3_all_exported_percent"] = weighted(frame.GW_HI3_percent, frame.weight)
    metrics["SA_selected_median_dyne_s_cm2"] = float(selected.SA_dyne_s_cm2.median())
    report["metrics"] = metrics
    report["GW_linearized_source_time_quartile_fractions"] = (source_bins/source_bins.sum()).tolist() if source_bins.sum() else None
    bpm = config.get("bpm", 120)
    report["transit_fraction_of_half_cycle_all_paths"] = describe(frame.transit_s / (30/bpm))
    # Subsamples of an existing export: representativeness diagnostic, not new-seed convergence.
    rng = np.random.default_rng(2026)
    subsamples = []
    for fraction in (.25, .5):
        count = max(1, int(len(selected) * fraction))
        if len(selected):
            estimates = []
            for _ in range(30):
                sample = selected.iloc[rng.choice(len(selected), count, replace=False)]
                estimates.append(weighted(sample.GW_HI3_percent, sample.weight))
            subsamples.append({"fraction": fraction, "count": count,
                               "HI3_percent_quantiles": describe(estimates)})
    report["existing_export_subsample_sensitivity_NOT_convergence"] = subsamples
    destination.mkdir(parents=True, exist_ok=False)
    frame.to_csv(destination / "streamlines.csv", index=False)
    pd.DataFrame(hot).to_csv(destination / "hotspot_segments.csv", index=False)
    write_json(destination / "diagnostics.json", report)
    lines = [f"# {config.get('label', path.stem)}", "", f"Input: {path}",
             f"Population: {report['population']}; weighting: {report['weighting']}",
             f"Paths: {n}; selected: {len(selected)}; selected/exported weight: {report['selected_weight_fraction_of_exported']:.6f}",
             "", "## Metrics", ""]
    lines.extend(f"- {key}: {value:.8g}" for key, value in metrics.items())
    lines.extend(["", "## Interpretation limits", ""] + [f"- {w}" for w in warnings])
    (destination / "report.md").write_text("\n".join(lines), encoding="utf-8")
    return {"label": config.get("label", path.stem), "input": str(path),
            "population": report["population"], "weighting": report["weighting"],
            "n_streamlines": n, "coverage_of_exported_weight": report["selected_weight_fraction_of_exported"],
            **{k: config.get(k) for k in ("size_cc", "phase", "geometry", "mesh_cells", "flow_L_min")}, **metrics}


def write_json(path, value):
    def clean(v):
        if isinstance(v, dict):
            return {k: clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [clean(x) for x in v]
        if isinstance(v, (float, np.floating)) and not np.isfinite(v):
            return None
        return v
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False), encoding="utf-8")


def self_test():
    from haem5 import compute_streamline_hemolysis
    c = POWER_LAW_CONSTANTS["GW"]
    analytic = c["A"] * 8 ** c["alpha"] * .2 ** c["beta"]
    errors = []
    for n in (4, 40, 400):
        dt, tau = np.full(n, .2/n), np.full(n, 8.)
        h2, h3, exact, _ = damage(dt, tau, c)
        np.testing.assert_allclose([h3, exact], analytic, rtol=1e-12)
        errors.append(abs(h2-analytic))
    assert errors[2] < errors[1] < errors[0]
    dt, tau = np.array([.1, .01, .02]), np.array([.2, 12., 4.])
    np.testing.assert_allclose(damage(dt, tau, c)[1], damage(dt[::-1], tau[::-1], c)[1])
    p = np.array([[0,0,0,.2,.1], [.01,0,0,12,.2], [.02,0,0,4,.001]], float)
    nodes = pd.DataFrame(p, columns=["x","y","z","shear","vel"])
    _, reference2, reference3 = compute_streamline_hemolysis(nodes, list(nodes.index))
    d, t, _ = segments(p)
    np.testing.assert_allclose(damage(d,t,c)[:2], [reference2, reference3], rtol=1e-12)
    p[:,4] = 0
    assert segments(p, floor=False)[0] is None
    print("PASS: uniform stress, refinement, ordering, haem5 agreement, zero-speed handling")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exports", nargs="*", type=Path)
    parser.add_argument("--manifest", type=Path, help="JSON with cases array; per-case settings override CLI defaults")
    parser.add_argument("--outdir", type=Path, default=Path("outputs/diagnostics"))
    parser.add_argument("--expected-seeds", type=int)
    parser.add_argument("--outlet-plane", nargs=6, type=float, metavar=("X","Y","Z","NX","NY","NZ"))
    parser.add_argument("--outlet-tolerance", type=float, default=1e-4, help="metres; default 0.1 mm")
    parser.add_argument("--bpm", type=float, default=120.)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        if not args.exports and not args.manifest:
            return
    cases = [{"path": str(p)} for p in args.exports]
    if args.manifest:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8-sig"))
        for case in manifest["cases"]:
            for key in ("path", "weights"):
                if case.get(key):
                    case[key] = str((args.manifest.parent / case[key]).resolve())
            cases.append(case)
    if not cases:
        parser.error("Provide exports or --manifest, or use --self-test.")
    defaults = {"expected_seeds": args.expected_seeds, "outlet_plane": args.outlet_plane,
                "outlet_tolerance": args.outlet_tolerance, "bpm": args.bpm,
                "flow_model": "SST", "stress_model": "molecular viscosity * scalar shear strain rate"}
    output = []
    for index, case in enumerate(cases):
        config = {**defaults, **case}
        if config["bpm"] <= 0 or config["outlet_tolerance"] < 0:
            raise ValueError("BPM must be positive and outlet tolerance nonnegative.")
        path = Path(config["path"])
        target = args.outdir / f"{index+1:02d}_{path.stem}"
        print(f"Analyzing {path} -> {target}", flush=True)
        output.append(analyze(path, config, target))
    pd.DataFrame(output).to_csv(args.outdir / "comparison.csv", index=False)
    print(f"Saved {args.outdir / 'comparison.csv'}; see per-case report.md and diagnostics.json.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
