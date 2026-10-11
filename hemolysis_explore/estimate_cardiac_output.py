"""Estimate outlet stroke volume and cycle-mean flow from quasi-steady snapshots.

Run from the repository root:
    python hemolysis_explore/estimate_cardiac_output.py
    python hemolysis_explore/estimate_cardiac_output.py --json-output outlet_co.json

Phase is percent of SYSTOLE, not percent of the complete cardiac cycle.
Defaults follow README.md / operating_point.py, including the final phase of 83%.
Outlet flows are supplied by the user, in ascending phase order. Assumed zero
flow at 0%/100% and throughout diastole; no diastolic leakage data are available.
The PCHIP/trapezoid difference is interpolation sensitivity, NOT an error bound.
Changing timing with fixed snapshot flows does not simulate a new operating point.
Requires numpy (already used by the surrounding analysis).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

DEFAULT_PHASES = (16.67, 33.33, 50.0, 66.67, 83.0)
DEFAULT_FLOWS = (9.32064, 16.35774, 20.18532, 19.28142, 12.1371)


def pchip_integral(x, y):
    """Exact integral of shape-preserving piecewise cubic Hermite interpolation.

    Interior slopes use a weighted harmonic mean; endpoint slopes use limited
    one-sided three-point estimates. These prevent overshoot on monotone spans.
    x must be strictly increasing and contain at least two points.
    """
    h = np.diff(x)
    delta = np.diff(y) / h
    slopes = np.zeros(len(x), dtype=float)
    if len(x) == 2:
        slopes[:] = delta[0]
    else:
        for k in range(1, len(x) - 1):
            left, right = delta[k - 1], delta[k]
            if left * right > 0:
                w1 = 2 * h[k] + h[k - 1]
                w2 = h[k] + 2 * h[k - 1]
                slopes[k] = (w1 + w2) / (w1 / left + w2 / right)

        def endpoint(h0, h1, d0, d1):
            slope = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
            if np.sign(slope) != np.sign(d0):
                return 0.0
            if np.sign(d0) != np.sign(d1) and abs(slope) > 3 * abs(d0):
                return 3 * d0
            return slope

        slopes[0] = endpoint(h[0], h[1], delta[0], delta[1])
        slopes[-1] = endpoint(h[-1], h[-2], delta[-1], delta[-2])
    # Integral of the Hermite basis on each interval, without grid approximation.
    return float(np.sum(h * (y[:-1] + y[1:]) / 2
                        + h**2 * (slopes[:-1] - slopes[1:]) / 12))


def estimate(phases=DEFAULT_PHASES, flows=DEFAULT_FLOWS, bpm=120.0,
             systolic_fraction=0.5, systole_seconds=None):
    """Return both integration estimates plus the explicit timing assumptions.

    For phase fraction f, I = integral_0^1 Q(f) df [L/min].
    SV [mL/beat] = I * Ts [s] * 1000/60.
    CO [L/min] = SV * bpm/1000, assuming zero diastolic outlet flow.
    """
    phases = np.asarray(phases, dtype=float)
    flows = np.asarray(flows, dtype=float)
    if (phases.ndim != 1 or flows.ndim != 1 or len(phases) != len(flows)
            or len(phases) < 2):
        raise ValueError("Supply at least two paired phases and outlet flows.")
    if not np.all(np.isfinite(phases)) or not np.all(np.isfinite(flows)):
        raise ValueError("All phases and flows must be finite.")
    if np.any(np.diff(phases) <= 0) or phases[0] <= 0 or phases[-1] >= 100:
        raise ValueError("Snapshot phases must increase strictly within (0, 100).")
    if not np.isfinite(bpm) or bpm <= 0:
        raise ValueError("Heart rate must be positive and finite.")
    cycle_seconds = 60.0 / bpm
    if systole_seconds is None:
        if not np.isfinite(systolic_fraction) or not 0 < systolic_fraction <= 1:
            raise ValueError("Systolic fraction must lie in (0, 1].")
        systole_seconds = cycle_seconds * systolic_fraction
    if not np.isfinite(systole_seconds) or not 0 < systole_seconds <= cycle_seconds:
        raise ValueError("Systole duration must be positive and no longer than a cycle.")
    x = np.r_[0.0, phases / 100.0, 1.0]
    y = np.r_[0.0, flows, 0.0]
    linear = float(np.sum(np.diff(x) * (y[:-1] + y[1:]) / 2))
    areas = {"pchip": pchip_integral(x, y), "trapezoid": linear}
    methods = {}
    for name, area in areas.items():
        sv = area * systole_seconds * 1000.0 / 60.0
        methods[name] = {"mean_systolic_flow_l_min": area,
                         "stroke_volume_ml": sv,
                         "cardiac_output_l_min": sv * bpm / 1000.0}
    return {"bpm": float(bpm), "systole_seconds": float(systole_seconds),
            "systolic_fraction": float(systole_seconds / cycle_seconds),
            "snapshot_phases_pct": phases.tolist(),
            "snapshot_outlet_flows_l_min": flows.tolist(),
            "assumed_endpoint_flows_l_min": [0.0, 0.0],
            "assumed_diastolic_flow_l_min": 0.0,
            "methods": methods,
            "limitations": "Quasi-steady reconstruction; interpolation spread is not "
                           "a confidence interval. Unmeasured transients and diastolic "
                           "regurgitation are not included."}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--phases", type=float, nargs="+", default=DEFAULT_PHASES,
                        help="Snapshot phases in percent of systole, ascending")
    parser.add_argument("--flows", type=float, nargs="+", default=DEFAULT_FLOWS,
                        help="Signed outlet flows in L/min; outward positive")
    parser.add_argument("--bpm", type=float, default=120.0)
    timing = parser.add_mutually_exclusive_group()
    timing.add_argument("--systolic-fraction", type=float, default=0.5)
    timing.add_argument("--systole-seconds", type=float)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()
    try:
        result = estimate(args.phases, args.flows, args.bpm,
                          args.systolic_fraction, args.systole_seconds)
    except ValueError as exc:
        parser.error(str(exc))
    print(f"Heart rate: {result['bpm']:g} bpm; systole: "
          f"{result['systole_seconds'] * 1000:g} ms "
          f"({result['systolic_fraction']:.1%} of cycle)")
    print("Assumptions: zero outlet flow at 0%/100% systole and during diastole.")
    print("Snapshot phases [%]:", ", ".join(f"{p:g}" for p in args.phases))
    for name, values in result["methods"].items():
        print(f"{name:10s}: stroke volume {values['stroke_volume_ml']:.3f} mL/beat; "
              f"cardiac output {values['cardiac_output_l_min']:.4f} L/min")
    print(result["limitations"])
    print("Timing changes reuse fixed flows; they do not predict a new CFD operating point.")
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"Saved {args.json_output}")


if __name__ == "__main__":
    main()
