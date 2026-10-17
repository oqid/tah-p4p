"""Operating point of each steady systolic snapshot from the CAD/BC definitions.

Diaphragm (velocity inlet), Inventor parametric sketch revolved about its centre:
    x(t) = 2 L t,  y(t) = Amplitude * Pump * sin(2 pi t),  t in [0, 0.5]
    -> half-sine dome, radius R = L/2, height h = Amplitude * Pump, Pump = phase fraction
    L = 56.211 mm, Amplitude = 20 mm (50 cc)
Inlet speed (normal to the dome): v = pi S / (2 Ts) sin(pi phase), S = 20 mm
Q = v * A_dome(phase)  (velocity applied over the whole dome surface)
Outlet / valve orifice area A_out = 4.64768 cm^2.
Tilting-disc valve, angle theta from closed (0 = disc normal to the outlet axis):
    geometric open area ~ A_out (1 - cos theta)   (orifice minus projected disc)
    valve-jet speed proxy U_jet = Q / (A_out (1 - cos theta))
For 35 cc every length is scaled by s = 0.888 (areas s^2), same valve angles.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

RHO, MU = 1060.0, 0.0035
BASE = dict(L=56.211e-3, S=0.02, A_out=4.64768e-4)
VALVE_ANGLE_DEG = {16.67: 24.9, 33.33: 49.8, 50.0: 75.0, 66.67: 50.1, 83.33: 25.2}


def dome(h, L):
    """Surface area and volume of the revolved half-sine dome of height h."""
    R = L / 2
    r = np.linspace(0, R, 4001)
    y = h * np.cos(np.pi * r / L)
    dy = -h * np.pi / L * np.sin(np.pi * r / L)
    A = np.trapezoid(2 * np.pi * r * np.sqrt(1 + dy ** 2), r)
    V = np.trapezoid(2 * np.pi * r * y, r)
    return A, V


def geometry(scale=1.0):
    return {k: v * (scale ** 2 if k == "A_out" else scale) for k, v in BASE.items()}


def operating_point(phase_pct, bpm=120.0, scale=1.0, systolic_fraction=0.5, angle_deg=None):
    g = geometry(scale)
    Ts = 60.0 / bpm * systolic_fraction
    f = phase_pct / 100.0
    v = np.pi * g["S"] / (2 * Ts) * np.sin(np.pi * f)
    A_d, _ = dome(g["S"] * f, g["L"])
    Q = v * A_d
    if angle_deg is None:
        key = min(VALVE_ANGLE_DEG, key=lambda p: abs(p - phase_pct))
        ang = VALVE_ANGLE_DEG[key] if abs(key - phase_pct) <= 0.5 else None
    else:
        ang = angle_deg
    open_frac = 1 - np.cos(np.radians(ang)) if ang is not None else np.nan
    D = np.sqrt(4 * g["A_out"] / np.pi)
    U_out = Q / g["A_out"]
    U_jet = Q / (g["A_out"] * open_frac)
    return {"phase_pct": phase_pct, "bpm": bpm, "scale": scale, "systole_ms": Ts * 1e3,
            "diaphragm_speed_m_s": v, "dome_area_mm2": A_d * 1e6, "Q_ml_s": Q * 1e6,
            "Q_l_min": Q * 6e4, "valve_angle_deg": ang, "valve_open_frac": open_frac,
            "U_outlet_mean_m_s": U_out, "U_jet_proxy_m_s": U_jet,
            "Re_outlet": RHO * U_out * D / MU, "Re_jet_proxy": RHO * U_jet * D / MU}


def stroke(bpm=120.0, scale=1.0, systolic_fraction=0.5, n=2001):
    """Volume delivered by the inlet BC over systole, and the CAD dome sweep."""
    g = geometry(scale)
    Ts = 60.0 / bpm * systolic_fraction
    t = np.linspace(0, Ts, n)
    f = t / Ts
    v = np.pi * g["S"] / (2 * Ts) * np.sin(np.pi * f)
    A = np.array([dome(g["S"] * x, g["L"])[0] for x in f])
    sv = np.trapezoid(v * A, t)
    return {"bc_stroke_volume_ml": sv * 1e6, "bc_cardiac_output_l_min": sv * 1e3 * bpm,
            "cad_dome_swept_volume_ml": dome(g["S"], g["L"])[1] * 1e6,
            "t": t, "Q_ml_s": v * A * 1e6}


def table(phases=(16.67, 33.33, 50.0, 66.67, 83.33), **kw):
    return pd.DataFrame([operating_point(p, **kw) for p in phases])


if __name__ == "__main__":
    pd.set_option("display.width", 220)
    for lab, kw in (("50 cc 120 bpm", dict()), ("35 cc 120 bpm", dict(scale=0.888)),
                    ("35 cc 171 bpm", dict(scale=0.888, bpm=171.4))):
        st = stroke(**kw)
        print(f"\n{lab}: BC stroke volume {st['bc_stroke_volume_ml']:.1f} mL, CO "
              f"{st['bc_cardiac_output_l_min']:.2f} L/min, CAD dome sweep {st['cad_dome_swept_volume_ml']:.1f} mL")
        print(table(**kw).round(3).to_string(index=False))
