"""Measured flow through each steady snapshot: CFX massFlow() at the diaphragm
velocity inlet and at the outlet opening (source: TAH-P4P Mass Flow.xlsx, 17 Oct 2026).
Inlet and outlet agree to within 0.01 % in every case, so the inlet value is used.
The 35 cc 16.67 % values are from the corrected (scaled-geometry) re-runs.

    import massflow; massflow.q_l_min(scale, bpm, phase_pct)
"""
from __future__ import annotations

RHO = 1060.0
PHASES = (16.67, 33.33, 50.0, 66.67, 83.33)
# (scale, bpm): (inlet kg/s per phase, outlet kg/s magnitude per phase)
DATA = {
    (1.0, 120): ([0.162946, 0.288379, 0.356307, 0.339742, 0.212291],
                 [0.162947, 0.288374, 0.356307, 0.339738, 0.212293]),
    (0.888, 120): ([0.117592, 0.201893, 0.250064, 0.237796, 0.150072],
                   [0.117592, 0.201892, 0.250064, 0.237796, 0.150072]),
    (0.888, 171): ([0.162993, 0.288400, 0.356337, 0.339757, 0.212350],
                   [0.162992, 0.288402, 0.356337, 0.339736, 0.212351]),
}


def _key(scale, bpm):
    for (s, b) in DATA:
        if abs(s - scale) < 0.01 and abs(b - bpm) < 2:
            return (s, b)
    return None


def available(scale, bpm):
    return _key(scale, bpm) is not None


def _idx(phase_pct):
    i = min(range(5), key=lambda k: abs(PHASES[k] - phase_pct))
    if abs(PHASES[i] - phase_pct) > 0.5:
        raise KeyError(phase_pct)
    return i


def mdot(scale, bpm, phase_pct, where="inlet"):
    k = _key(scale, bpm)
    if k is None:
        raise KeyError((scale, bpm))
    return DATA[k][0 if where == "inlet" else 1][_idx(phase_pct)]


def q_ml_s(scale, bpm, phase_pct):
    return mdot(scale, bpm, phase_pct) / RHO * 1e6


def q_l_min(scale, bpm, phase_pct):
    return mdot(scale, bpm, phase_pct) / RHO * 6e4


def imbalance_pct(scale, bpm, phase_pct):
    a, b = mdot(scale, bpm, phase_pct, "inlet"), mdot(scale, bpm, phase_pct, "outlet")
    return 100 * abs(a - b) / a


if __name__ == "__main__":
    for (s, b) in DATA:
        qs = [q_l_min(s, b, p) for p in PHASES]
        im = max(imbalance_pct(s, b, p) for p in PHASES)
        print(f"scale {s} {b} BPM: Q = {[round(q, 2) for q in qs]} L/min, max in/out imbalance {im:.4f} %")
