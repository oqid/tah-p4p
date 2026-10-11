"""Similarity-scaling hypotheses for a downsized (35 cc) chamber.

Before the 35 cc CFD exists, the 50 cc results plus dimensional analysis give a
falsifiable prediction for HI3. Assumptions (state them in the report):
  * geometric similarity: every length (incl. diaphragm stroke S and valve
    opening) scales by s = (V35/V50)^(1/3)
  * the same systolic fraction of the cycle and the same sinusoidal waveform,
    compared at the same normalised systolic phase
  * the same valve opening state at the same phase
  * per-path HI3 = A T^beta tau_eff^alpha  (exact identity, see explore_core)

Velocity:  U ~ S / Ts ~ s * BPM          (diaphragm-driven)
Transit:   T ~ L / U  ~ 1 / BPM
Stress (three bounding regimes, molecular viscous stress mu*gamma_dot):
  viscous  : tau ~ mu U / L                 -> ~ BPM            (Re-independent)
  bl       : tau ~ sqrt(rho mu) U^1.5 L^-0.5 -> ~ s * BPM^1.5   (laminar boundary layer / jet shear layer)
  inertial : tau ~ rho U^2                  -> ~ s^2 * BPM^2    (upper bound, stress follows dynamic pressure)
Each regime gives HI3 ratio = (T35/T50)^beta * (tau35/tau50)^alpha.
"""
from __future__ import annotations

import numpy as np

from explore_core import CONSTANTS

REGIMES = {
    # exponents (a, b) in tau ~ U^a * L^b
    "viscous (τ∝μU/L)": (1.0, -1.0),
    "boundary layer (τ∝U^1.5 L^-0.5)": (1.5, -0.5),
    "inertial (τ∝ρU²)": (2.0, 0.0),
}


def ratios(v_small=35.0, v_ref=50.0, bpm_small=120.0, bpm_ref=120.0):
    s = (v_small / v_ref) ** (1 / 3)
    f = bpm_small / bpm_ref
    U = s * f
    L = s
    T = L / U
    out = {"length_scale": s, "U_ratio": U, "T_ratio": T, "Re_ratio": U * L,
           "Q_ratio": U * L ** 2, "CO_ratio": (v_small / v_ref) * f,
           "Wo_ratio": L * np.sqrt(f)}
    for name, (a, b) in REGIMES.items():
        tau = U ** a * L ** b
        out[f"tau_ratio|{name}"] = tau
        out[f"SA_ratio|{name}"] = tau * T
        for k, c in CONSTANTS.items():
            out[f"HI3_ratio|{k}|{name}"] = T ** c["beta"] * tau ** c["alpha"]
    return out


def scenarios(bpm_ref=120.0, v_small=35.0, v_ref=50.0):
    matched_bpm = bpm_ref * v_ref / v_small
    return {
        f"Same rate ({bpm_ref:.0f} bpm)": ratios(v_small, v_ref, bpm_ref, bpm_ref),
        f"Matched output ({matched_bpm:.0f} bpm)": ratios(v_small, v_ref, matched_bpm, bpm_ref),
    }


if __name__ == "__main__":
    import pandas as pd
    pd.set_option("display.width", 200)
    print(pd.DataFrame(scenarios()).T.round(3).T.to_string())
