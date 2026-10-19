"""Flow and cardiac output for every case, from the measured CFX massFlow (massflow.py).

Q(phase) = massFlow()@diaphragm inlet / rho; inlet and outlet agree within 0.01 %.
For reference the table also lists the analytic estimate (mean seed speed x dome
surface area, h = S * x_hat, lengths scaled by s) and the Sim Tracker outlet
integrals, whose 35 cc entries were mis-recorded and are not used.
Integration over systole uses estimate_cardiac_output.estimate (PCHIP and trapezoid).

    python cardiac_output_cases.py OUTDIR [FOLDER_50 FOLDER_35_120 FOLDER_35_171]
(defaults: ..\\z_inputs ..\\35cc_inputs\\120BPM ..\\35cc_inputs\\171BPM)
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import operating_point as op
import massflow as mf
from estimate_cardiac_output import estimate

PH = [16.67, 33.33, 50.0, 66.67, 83.33]
X_HAT = [0.067, 0.25, 0.5, 0.75, 0.933]          # normalised diaphragm position (Sim Tracker)
CASES = {  # label: (folder, bpm, scale, Sim Tracker outlet flow L/min)
    "50 cc 120 BPM": ("../z_inputs", 120, 1.0, [9.32064, 16.35774, 20.18532, 19.28142, 12.1371]),
    "35 cc 120 BPM": ("../35cc_inputs/120BPM", 120, 0.888, [13.69152, 13.65792, 15.62916, 16.04826, 13.59762]),
    "35 cc 171 BPM": ("../35cc_inputs/171BPM", 171, 0.888, [19.27374, 19.68138, 22.62198, 22.95432, 19.61862]),
}


def seed_speed(path):
    lines = open(path).read().splitlines()
    li = next(i for i, l in enumerate(lines) if l.strip() == "[Lines]")
    n = int(lines[li + 1].split(",")[1])
    di = next(i for i, l in enumerate(lines) if l.strip() == "[Data]")
    return float(np.loadtxt(lines[di + 2:di + 2 + n], delimiter=",")[:, 5].mean())


def main(out, folders=None):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    if folders:
        for lab, f in zip(CASES, folders):
            CASES[lab] = (f,) + CASES[lab][1:]
    import explore_core as ec
    rows, summary = [], {}
    for lab, (folder, bpm, s, q_tracker) in CASES.items():
        files = sorted(Path(folder).glob("*.csv"), key=lambda p: ec.phase_from_name(p))
        g = op.geometry(s)
        q_in = []
        for f, p, xh, qt in zip(files, PH, X_HAT, q_tracker):
            v = seed_speed(f)
            A, _ = op.dome(g["S"] * xh, g["L"])
            q_an = v * A * 6e4
            q = mf.q_l_min(s, bpm, p)
            q_in.append(q)
            rows.append({"case": lab, "phase_pct": p, "bpm": bpm, "length_scale": s,
                         "diaphragm_speed_cm_s": v * 100, "dome_area_mm2": A * 1e6,
                         "Q_massflow_l_min": q, "inlet_outlet_imbalance_pct": mf.imbalance_pct(s, bpm, p),
                         "Q_analytic_l_min": q_an, "analytic_over_massflow": q_an / q,
                         "Q_tracker_outlet_l_min": qt})
        e_in = estimate(PH, q_in, bpm)
        summary[lab] = {"massflow_based": e_in["methods"], "systole_ms": e_in["systole_seconds"] * 1e3}
    pd.DataFrame(rows).to_csv(out / "flows_by_phase.csv", index=False)
    (out / "cardiac_output.json").write_text(json.dumps(summary, indent=2))
    for lab, d in summary.items():
        i = d["massflow_based"]
        print(f"{lab}: massFlow SV {i['trapezoid']['stroke_volume_ml']:.1f}-{i['pchip']['stroke_volume_ml']:.1f} mL, "
              f"CO {i['trapezoid']['cardiac_output_l_min']:.2f}-{i['pchip']['cardiac_output_l_min']:.2f} L/min")
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "outputs/cardiac_output", sys.argv[2:5] or None)
