"""
Independent empirical stress-exposure screening for haem5.py streamline results.

Purpose
-------
This companion analysis does NOT alter the Taskin/Giersiepen hemolysis model.
It adds a separate screening check against the empirical stress-exposure damage
boundaries reproduced as Fig. 26 in US Patent 5,924,975 from NIH/NHLBI guidance.

Why separate?
-------------
Taskin et al. use a continuous power-law model (HI = C t^alpha tau^beta) and
Lagrangian accumulation variants. Those equations encode stress and exposure-time
dependence but do not impose a hard zero-damage threshold. Fig. 26 is a different
kind of criterion: a stress-duration boundary. Keeping it separate avoids mixing
two distinct empirical models.

Method for time-varying CFD streamlines
---------------------------------------
The Fig. 26 boundary is for stress magnitude versus duration. For a streamline
whose stress varies with time, comparing local stress to *cumulative time since
inlet* is not physically equivalent to a sustained exposure. Instead, this script
constructs a conservative stress-duration exposure envelope:

    For each stress level S, find the longest CONTIGUOUS duration for which
    local shear stress >= S.

This produces a per-streamline exposure envelope (duration, stress). The envelope
is compared with the empirical RBC and platelet damage curves. A utilization
ratio > 1 means the exposure envelope crosses the corresponding empirical damage
boundary.

Outputs
-------
  empirical_threshold_summary_<run>.csv
      One row per streamline: exposure time, peak shear, RBC/platelet utilization,
      pass/fail, and inlet-velocity flow weight.

  empirical_threshold_device_<run>.csv
      Compact report-ready device summary.

  empirical_stress_duration_<run>.png
      Report-ready log-log plot: empirical damage boundaries + worst-case device
      exposure envelope.

Important source caveat
-----------------------
The patent reproduces Fig. 26 from NIH Publication 85-2185 and gives several
numerical anchors in its text, but it does not tabulate the full curves. The
curve arrays below are therefore APPROXIMATE digitizations/anchors intended for
engineering screening, not a replacement for the original NIH dataset. Keep this
qualification in the report.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import haem5


# ---------------------------------------------------------------------------
# 1. Empirical threshold curves from Fig. 26 / accompanying patent text
# ---------------------------------------------------------------------------
# Units used internally here: seconds and Pa.
# 1 dyne/cm^2 = 0.1 Pa.
#
# RBC anchors: approximate digitization of the solid "RED CELLS" curve, with
# the long-duration floor kept at ~150 Pa to match the patent's stated
# ~1500 dyne/cm^2 threshold (the patent later also describes ~2000 dyne/cm^2
# as a long-duration critical value; the figure/text are approximate).
RBC_THRESHOLD = pd.DataFrame({
    "time_s": np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4,
                         1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1,
                         1.0, 3.0, 10.0, 100.0]),
    "stress_Pa": np.array([6000, 4200, 2600, 1600, 950, 600,
                            380, 270, 210, 175, 155, 150,
                            150, 150, 150, 150], dtype=float),
})

# Platelet anchors: approximate digitization of the dashed "PLATELETS" curve.
# The 3 s, 35 Pa point is explicitly stated in the patent text
# (350 dyne/cm^2 at 3 seconds).
PLATELET_THRESHOLD = pd.DataFrame({
    "time_s": np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4,
                         1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1,
                         1.0, 3.0, 10.0, 30.0, 100.0]),
    "stress_Pa": np.array([100000, 60000, 30000, 17000, 9000, 5200,
                            3000, 1700, 900, 500, 220, 110,
                            55, 35, 24, 19, 17], dtype=float),
})

CURVES = {
    "RBC": RBC_THRESHOLD,
    "Platelet": PLATELET_THRESHOLD,
}


def _log_interp_threshold(time_s: np.ndarray | float, curve: pd.DataFrame) -> np.ndarray:
    """Log-log interpolate critical stress at exposure duration time_s."""
    t = np.asarray(time_s, dtype=float)
    t_safe = np.clip(t, curve["time_s"].iloc[0], curve["time_s"].iloc[-1])
    log_tau = np.interp(
        np.log10(t_safe),
        np.log10(curve["time_s"].to_numpy()),
        np.log10(curve["stress_Pa"].to_numpy()),
    )
    return 10.0 ** log_tau


def _longest_contiguous_duration_above(stress: np.ndarray, dt: np.ndarray, level: float) -> float:
    """Longest contiguous time for which segment stress is >= level."""
    active = stress >= level
    longest = 0.0
    current = 0.0
    # haem5 stores dt[0] = 0; rows 1: correspond to path segments.
    for is_active, dti in zip(active[1:], dt[1:]):
        if is_active:
            current += float(max(dti, 0.0))
            longest = max(longest, current)
        else:
            current = 0.0
    return longest


def build_exposure_envelope(df: pd.DataFrame, n_levels: int = 220) -> pd.DataFrame:
    """
    Build a stress-duration exposure envelope for one time-varying streamline.

    At each stress level S, duration = longest contiguous time with tau >= S.
    """
    stress = np.asarray(df["shear"], dtype=float)
    dt = np.asarray(df["dt"], dtype=float)
    finite = stress[np.isfinite(stress) & (stress > 0)]
    if len(finite) == 0:
        return pd.DataFrame(columns=["stress_Pa", "duration_s"])

    lo = max(float(np.min(finite)), 1e-6)
    hi = float(np.max(finite))
    if hi <= lo:
        levels = np.array([hi])
    else:
        levels = np.logspace(np.log10(lo), np.log10(hi), n_levels)

    durations = np.array([
        _longest_contiguous_duration_above(stress, dt, level)
        for level in levels
    ])
    out = pd.DataFrame({"stress_Pa": levels, "duration_s": durations})
    return out[out["duration_s"] > 0].reset_index(drop=True)


def criterion_utilization(envelope: pd.DataFrame, curve: pd.DataFrame) -> tuple[float, float, float]:
    """
    Return (max utilization, stress at max, duration at max).

    utilization = actual stress / empirical critical stress at that duration.
    >1 means the empirical boundary is crossed.
    """
    if envelope.empty:
        return 0.0, np.nan, np.nan
    crit = _log_interp_threshold(envelope["duration_s"].to_numpy(), curve)
    util = envelope["stress_Pa"].to_numpy() / crit
    i = int(np.nanargmax(util))
    return float(util[i]), float(envelope["stress_Pa"].iloc[i]), float(envelope["duration_s"].iloc[i])


def analyze_empirical_thresholds(result: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Analyze all haem5 streamlines and return per-line, device, worst-envelope tables."""
    rows = []
    envelopes = []

    # Same weighting convention used in haem5: inlet velocity as flux proxy.
    raw_weights = np.array([
        float(df["vel"].iloc[0]) if len(df) else 0.0
        for df in result["streamlines"]
    ])
    if np.sum(raw_weights) > 0:
        weights = raw_weights / np.sum(raw_weights)
    else:
        weights = np.full(len(raw_weights), 1.0 / max(len(raw_weights), 1))

    for i, (df, weight) in enumerate(zip(result["streamlines"], weights)):
        env = build_exposure_envelope(df)
        if not env.empty:
            env = env.copy()
            env["streamline"] = i
            envelopes.append(env)

        rbc_u, rbc_s, rbc_t = criterion_utilization(env, RBC_THRESHOLD)
        plt_u, plt_s, plt_t = criterion_utilization(env, PLATELET_THRESHOLD)

        rows.append({
            "streamline": i,
            "flow_weight": float(weight),
            "residence_time_s": float(df["t_cumulative"].iloc[-1]) if len(df) else 0.0,
            "peak_shear_Pa": float(df["shear"].max()) if len(df) else 0.0,
            "RBC_max_utilization": rbc_u,
            "RBC_margin_to_threshold": (1.0 / rbc_u) if rbc_u > 0 else np.inf,
            "RBC_exceeds": bool(rbc_u >= 1.0),
            "RBC_worst_stress_Pa": rbc_s,
            "RBC_worst_duration_s": rbc_t,
            "Platelet_max_utilization": plt_u,
            "Platelet_margin_to_threshold": (1.0 / plt_u) if plt_u > 0 else np.inf,
            "Platelet_exceeds": bool(plt_u >= 1.0),
            "Platelet_worst_stress_Pa": plt_s,
            "Platelet_worst_duration_s": plt_t,
        })

    per_line = pd.DataFrame(rows)

    device_rows = []
    for label, prefix in [("Red blood cells", "RBC"), ("Platelets", "Platelet")]:
        u = per_line[f"{prefix}_max_utilization"].to_numpy(dtype=float)
        exceed = per_line[f"{prefix}_exceeds"].to_numpy(dtype=bool)
        max_idx = int(np.nanargmax(u)) if len(u) else 0
        device_rows.append({
            "criterion": label,
            "max_utilization": float(np.nanmax(u)) if len(u) else np.nan,
            "minimum_safety_factor": float(1.0 / np.nanmax(u)) if len(u) and np.nanmax(u) > 0 else np.inf,
            "streamlines_exceeding_percent": float(100.0 * np.mean(exceed)) if len(exceed) else np.nan,
            "flow_weighted_exceedance_percent": float(100.0 * np.sum(per_line.loc[exceed, "flow_weight"])) if len(exceed) else np.nan,
            "worst_streamline": int(per_line.iloc[max_idx]["streamline"]) if len(per_line) else np.nan,
            "worst_stress_Pa": float(per_line.iloc[max_idx][f"{prefix}_worst_stress_Pa"]) if len(per_line) else np.nan,
            "worst_contiguous_duration_s": float(per_line.iloc[max_idx][f"{prefix}_worst_duration_s"]) if len(per_line) else np.nan,
            "status": "EXCEEDS" if np.any(exceed) else "BELOW THRESHOLD",
        })
    device = pd.DataFrame(device_rows)

    # Build a device worst-case exposure envelope: at each stress level, the
    # largest contiguous exposure duration found on any streamline.
    if envelopes:
        all_env = pd.concat(envelopes, ignore_index=True)
        stress_min = max(all_env["stress_Pa"].min(), 1e-6)
        stress_max = all_env["stress_Pa"].max()
        grid = np.logspace(np.log10(stress_min), np.log10(stress_max), 300)
        worst_d = np.zeros_like(grid)
        for df in result["streamlines"]:
            stress = np.asarray(df["shear"], dtype=float)
            dt = np.asarray(df["dt"], dtype=float)
            d = np.array([_longest_contiguous_duration_above(stress, dt, s) for s in grid])
            worst_d = np.maximum(worst_d, d)
        worst_envelope = pd.DataFrame({"stress_Pa": grid, "duration_s": worst_d})
        worst_envelope = worst_envelope[worst_envelope["duration_s"] > 0].reset_index(drop=True)
    else:
        worst_envelope = pd.DataFrame(columns=["stress_Pa", "duration_s"])

    return per_line, device, worst_envelope


def plot_empirical_thresholds(worst_envelope: pd.DataFrame, save_path: Path | None = None):
    """Report-ready stress-duration threshold plot."""
    fig, ax = plt.subplots(figsize=(8.2, 5.6))

    ax.loglog(RBC_THRESHOLD["time_s"], RBC_THRESHOLD["stress_Pa"],
              linewidth=2.2, label="RBC damage boundary (Fig. 26 approx.)")
    ax.loglog(PLATELET_THRESHOLD["time_s"], PLATELET_THRESHOLD["stress_Pa"],
              linewidth=2.2, linestyle="--", label="Platelet damage boundary (Fig. 26 approx.)")

    if not worst_envelope.empty:
        # Envelope is parameterized by stress; sort by duration for clean plotting.
        p = worst_envelope.sort_values("duration_s")
        ax.loglog(p["duration_s"], p["stress_Pa"], linewidth=2.4,
                  label="Worst-case device exposure envelope")

    ax.set_xlabel("Longest contiguous exposure duration [s]")
    ax.set_ylabel("Scalar shear stress [Pa]")
    ax.set_title("Empirical stress–exposure screening")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def print_device_summary(device: pd.DataFrame):
    print("\n" + "=" * 82)
    print("Empirical stress-exposure threshold screening (independent of Taskin HI)")
    print("=" * 82)
    cols = [
        "criterion", "max_utilization", "minimum_safety_factor",
        "streamlines_exceeding_percent", "flow_weighted_exceedance_percent",
        "worst_stress_Pa", "worst_contiguous_duration_s", "status",
    ]
    print(device[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nInterpretation: utilization = actual stress / empirical critical stress")
    print("  utilization < 1  -> below empirical boundary")
    print("  utilization >= 1 -> boundary crossed")
    print("NOTE: Fig. 26 curves are approximate digitizations/anchors for screening.")
    print("=" * 82)


def main():
    parser = argparse.ArgumentParser(
        description="Run haem5 plus independent empirical Fig. 26 stress-duration screening."
    )
    parser.add_argument("filepath", type=Path, help="CFD-Post Generic export CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_empirical"))
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG generation")
    args = parser.parse_args()

    result = haem5.run_pipeline(args.filepath)
    per_line, device, worst_envelope = analyze_empirical_thresholds(result)
    print_device_summary(device)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_id = f"{datetime.now():%Y-%m-%d_%H-%M-%S}_{Path(args.filepath).stem}"

    p_stream = args.output_dir / f"empirical_threshold_summary_{run_id}.csv"
    p_device = args.output_dir / f"empirical_threshold_device_{run_id}.csv"
    p_env = args.output_dir / f"empirical_exposure_envelope_{run_id}.csv"
    per_line.to_csv(p_stream, index=False)
    device.to_csv(p_device, index=False)
    worst_envelope.to_csv(p_env, index=False)

    print("Saved:")
    print(f"  {p_stream}")
    print(f"  {p_device}")
    print(f"  {p_env}")

    if not args.no_plot:
        p_plot = args.output_dir / f"empirical_stress_duration_{run_id}.png"
        plot_empirical_thresholds(worst_envelope, p_plot)
        print(f"  {p_plot}")


if __name__ == "__main__":
    main()
