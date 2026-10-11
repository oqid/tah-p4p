"""Vectorised worst-case stress–duration envelope (same definition as
haem5_empiricalthresholds.analyze_empirical_thresholds):

  for each stress level S, the longest CONTIGUOUS time any path spends with
  local stress >= S (segment i carries node stress tau_i and duration dt_i;
  the path's first node is excluded, as in haem5).

Compared against the approximate NIH Fig. 3.6 lysis traces (RBC, platelet),
read from ../papers/figure_3_6_digitization.csv (falls back to a local copy).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _trace_file():
    for p in (HERE.parent / "papers" / "figure_3_6_digitization.csv", HERE / "figure_3_6_digitization.csv"):
        if p.exists():
            return p
    raise FileNotFoundError("figure_3_6_digitization.csv not found")


def lysis_curves():
    t = pd.read_csv(_trace_file())
    t["time_s"] = 10.0 ** (-6 + 8 * (t["x_px"] - 320) / 488)          # same calibration as
    t["stress_Pa"] = 0.1 * 10.0 ** (6 - 4 * (t["y_px"] - 268) / 471)   # haem5_empiricalthresholds
    return {k: g[["time_s", "stress_Pa"]].reset_index(drop=True)
            for k, g in t.groupby("cell_type")}


def critical_stress(time_s, curve):
    """Log-log interpolated lysis stress at a duration (NaN outside the trace)."""
    t = np.asarray(time_s, float)
    ok = (t >= curve.time_s.iloc[0]) & (t <= curve.time_s.iloc[-1])
    lt = np.interp(np.log10(np.clip(t, curve.time_s.iloc[0], curve.time_s.iloc[-1])),
                   np.log10(curve.time_s), np.log10(curve.stress_Pa))
    return np.where(ok, 10 ** lt, np.nan)


def worst_envelope(snap, n_levels=300):
    tau, dt = snap.tau, snap.dt
    first = np.zeros(len(tau), bool)
    first[snap.offsets[:-1]] = True
    pos = tau[np.isfinite(tau) & (tau > 0)]
    seg = tau[~first & (dt > 0) & np.isfinite(tau)]          # stresses that carry duration
    levels = np.logspace(np.log10(max(pos.min(), 1e-6)), np.log10(seg.max()), n_levels)
    idx = np.arange(len(tau))
    worst = np.zeros(n_levels)
    for k, S in enumerate(levels):
        active = (tau >= S) & ~first
        cs = np.cumsum(np.where(active, dt, 0.0))
        brk = np.where(~active, idx, 0)
        last = np.maximum.accumulate(brk)
        worst[k] = (cs - cs[last]).max()
    env = pd.DataFrame({"stress_Pa": levels, "duration_s": worst})
    return env[env.duration_s > 0].reset_index(drop=True)


def safety_factors(env, curves):
    """Minimum margin (critical / actual stress) along the envelope, per curve."""
    out = {}
    for name, c in curves.items():
        u = env.stress_Pa.to_numpy() / critical_stress(env.duration_s.to_numpy(), c)
        out[name] = float(1 / np.nanmax(u)) if np.isfinite(u).any() else np.nan
    return out
