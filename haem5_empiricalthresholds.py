"""
Independent empirical stress-exposure screening for haem5.py streamline results.

Purpose
-------
This companion analysis does NOT alter the Taskin/Giersiepen hemolysis model.
It adds a separate screening check against the empirical stress-exposure damage
boundaries in Fig. 3.6 of the NIH/NHLBI Guidelines for Blood-Material Interactions
(printed p. 78 / PDF p. 87), also reproduced as patent Fig. 26.

Why separate?
-------------
Taskin et al. use a continuous power-law model (HI = C tau^alpha t^beta) and
Lagrangian accumulation variants. Those equations encode stress and exposure-time
dependence but do not impose a hard zero-damage threshold. Fig. 26 is a different
kind of criterion: a stress-duration boundary. Keeping it separate avoids mixing
two distinct empirical models.

Method for time-varying CFD streamlines
---------------------------------------
The Fig. 26 boundary is for stress magnitude versus duration. For a streamline
whose stress varies with time, comparing local stress to *cumulative time since
inlet* is not physically equivalent to a sustained exposure. Instead, this script
constructs a longest-contiguous stress-duration exposure envelope:

    For each stress level S, find the longest CONTIGUOUS duration for which
    local shear stress >= S.

This produces a per-streamline exposure envelope (duration, stress). The envelope
is compared with the empirical RBC and platelet damage curves. A utilization
ratio > 1 means the exposure envelope crosses the corresponding empirical damage
boundary. This is a screening construction, not a validated variable-stress
damage model: it treats qualifying intervals as sustained stress S and ignores
separated repeated exposures. It is not guaranteed to bound biological damage.

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
The curves are approximate manual traces of the supplied NIH scan, not original
experimental measurements. See papers/figure_3_6_digitization.md. The Guidelines
describe significant lysis and also note activation/sublytic effects below these
boundaries. Comparisons outside the traced time range are left unassessed.
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
# 1. Approximate lysis curves traced from NIH Fig. 3.6
# ---------------------------------------------------------------------------
# Units used internally here: seconds and Pa.
# 1 dyne/cm^2 = 0.1 Pa.
#
# Use the primary-source trace rather than historical patent anchors.
# Calibration and limitations are recorded beside the source PDF.
_TRACE = pd.read_csv(Path(__file__).resolve().parent / "papers" / "figure_3_6_digitization.csv")
_TRACE["time_s"] = 10.0 ** (-6 + 8 * (_TRACE["x_px"] - 320) / 488)
_TRACE["stress_Pa"] = 0.1 * 10.0 ** (6 - 4 * (_TRACE["y_px"] - 268) / 471)
RBC_THRESHOLD = _TRACE.loc[_TRACE["cell_type"] == "RBC", ["time_s", "stress_Pa"]].reset_index(drop=True)
PLATELET_THRESHOLD = _TRACE.loc[_TRACE["cell_type"] == "Platelet", ["time_s", "stress_Pa"]].reset_index(drop=True)

CURVES = {
    "RBC": RBC_THRESHOLD,
    "Platelet": PLATELET_THRESHOLD,
}


def _log_interp_threshold(time_s: np.ndarray | float, curve: pd.DataFrame) -> np.ndarray:
    """Log-log interpolate critical stress at exposure duration time_s."""
    t = np.asarray(time_s, dtype=float)
    valid = (t >= curve["time_s"].iloc[0]) & (t <= curve["time_s"].iloc[-1])
    t_safe = np.clip(t, curve["time_s"].iloc[0], curve["time_s"].iloc[-1])
    log_tau = np.interp(
        np.log10(t_safe),
        np.log10(curve["time_s"].to_numpy()),
        np.log10(curve["stress_Pa"].to_numpy()),
    )
    return np.where(valid, 10.0 ** log_tau, np.nan)


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
    if not np.isfinite(util).any():
        return np.nan, np.nan, np.nan
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
        assessed = np.isfinite(u)
        max_idx = int(np.nanargmax(u)) if np.any(assessed) else 0
        maximum = float(np.nanmax(u)) if np.any(assessed) else np.nan
        device_rows.append({
            "criterion": label,
            "max_utilization": maximum,
            "minimum_safety_factor": float(1.0 / maximum) if maximum > 0 else np.nan,
            "streamlines_assessed": int(assessed.sum()),
            "streamlines_unassessed": int((~assessed).sum()),
            "streamlines_exceeding_percent": float(100.0 * np.mean(exceed)) if len(exceed) else np.nan,
            "flow_weighted_exceedance_percent": float(100.0 * np.sum(per_line.loc[exceed, "flow_weight"])) if len(exceed) else np.nan,
            "worst_streamline": int(per_line.iloc[max_idx]["streamline"]) if len(per_line) else np.nan,
            "worst_stress_Pa": float(per_line.iloc[max_idx][f"{prefix}_worst_stress_Pa"]) if len(per_line) else np.nan,
            "worst_contiguous_duration_s": float(per_line.iloc[max_idx][f"{prefix}_worst_duration_s"]) if len(per_line) else np.nan,
            "status": ("EXCEEDS LYSIS REFERENCE" if np.any(exceed) else
                       "INCOMPLETE COVERAGE" if not np.all(assessed) or not len(u) else
                       "BELOW LYSIS REFERENCE"),
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


def plot_empirical_thresholds(worst_envelope: pd.DataFrame,
                              save_path: Path | None = None):
    """Full-range report plot with exposure duration displayed in milliseconds."""
    from matplotlib.ticker import LogLocator, LogFormatterMathtext

    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    p = worst_envelope.replace([np.inf, -np.inf], np.nan).dropna()
    p = p[(p["duration_s"] > 0) & (p["stress_Pa"] > 0)]
    p = p.sort_values(["duration_s", "stress_Pa"], ascending=[True, False])
    ax.loglog(RBC_THRESHOLD["time_s"] * 1000, RBC_THRESHOLD["stress_Pa"],
              color="#bf5548", lw=1.6, label="Red-cell lysis")
    ax.loglog(PLATELET_THRESHOLD["time_s"] * 1000, PLATELET_THRESHOLD["stress_Pa"],
              color="#7757a5", lw=1.6, ls="--", label="Platelet lysis")
    if not p.empty:
        ax.loglog(p["duration_s"] * 1000, p["stress_Pa"], color="#087f8c", lw=1.8,
                  label="CFD exposure envelope", zorder=3)
    ax.set_xlabel("Exposure duration [ms]", fontsize=10)
    ax.set_ylabel("Shear stress [Pa]", fontsize=10)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=9)
    ax.grid(which="major", color="#e1e1e1", lw=0.5)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=7))
    ax.xaxis.set_major_formatter(LogFormatterMathtext())
    ax.yaxis.set_major_formatter(LogFormatterMathtext())
    ax.set_xlim(1e-3, 1e5)
    lowest = min(1., p["stress_Pa"].min() / 2) if not p.empty else 1.
    ax.set_ylim(10 ** np.floor(np.log10(lowest)), 2e5)
    ax.legend(loc="upper right", fontsize=8, frameon=True,
              facecolor="white", edgecolor="none", framealpha=0.95)
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pad_inches=0.03)
    return fig


def save_empirical_figures(worst_envelope: pd.DataFrame, base_path: Path):
    """Save the full-range figure as a PNG report asset."""
    png = base_path.with_suffix(".png")
    fig = plot_empirical_thresholds(worst_envelope, png)
    plt.close(fig)
    return [png]


def print_device_summary(device: pd.DataFrame):
    print("\n" + "=" * 82)
    print("Empirical stress-exposure threshold screening (independent of Taskin HI)")
    print("=" * 82)
    cols = [
        "criterion", "max_utilization", "streamlines_unassessed",
        "streamlines_exceeding_percent", "flow_weighted_exceedance_percent",
        "worst_stress_Pa", "worst_contiguous_duration_s", "status",
    ]
    print(device[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nInterpretation: utilization = actual stress / empirical critical stress")
    print("  utilization < 1  -> below empirical boundary")
    print("  utilization >= 1 -> boundary crossed")
    print("NOTE: Approximate NIH Fig. 3.6 traces; screening for lysis, not activation.")
    print("      Repeated separated exposures are not accumulated by this envelope.")
    print("=" * 82)


def main():
    parser = argparse.ArgumentParser(
        description="Run haem5 plus independent NIH Fig. 3.6 stress-duration lysis screening."
    )
    parser.add_argument("filepath", type=Path, help="CFD-Post Generic export CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_empirical"))
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG generation")
    args = parser.parse_args()
    # This command saves figures; it does not need an interactive GUI backend.
    plt.switch_backend("Agg")

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
        for path in save_empirical_figures(worst_envelope, p_plot):
            print(f"  {path}")


if __name__ == "__main__":
    main()
