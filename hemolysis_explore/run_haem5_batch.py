"""Run haem5.run_pipeline (GW, HO, TZ) on each export and save a per-file summary,
plus verify explore_core reproduces haem5 per path. Usage: python run_haem5_batch.py OUT FILE..."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import haem5, explore_core as ec

out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
for f in map(Path, sys.argv[2:]):
    row = {"file": f.name, "phase_pct": ec.phase_from_name(f)}
    for k, c in haem5.POWER_LAW_CONSTANTS.items():
        r = haem5.run_pipeline(str(f), c["A"], c["alpha"], c["beta"])
        s = r["summary"]
        row.update({f"{k}_HI3_mean": s.HI3_percent.mean(), f"{k}_HI2_mean": s.HI2_percent.mean(),
                    f"{k}_HI3_vweighted": r["device_HI3_percent"],
                    f"{k}_NIH_HI3_vweighted_mg100L": r["device_HI3_NIH_mg_per_100L"]})
        if k == "GW":
            row.update(n_paths=len(s), SA_median=r["device_SA_median_dyne_s_cm2"],
                       transit_ms_median=s.transit_ms.median(), peak_stress_max=s.max_stress_Pa.max(),
                       inlet_speed_mean=s.mean_inlet_velocity.mean())
            s.to_csv(out / f"{f.stem}_haem5_GW_paths.csv", index=False)
            p = ec.build_snapshot(f).paths
            row["explore_core_max_rel_err_HI3"] = float(np.max(np.abs(p.GW_HI3.values - s.HI3_percent.values) / s.HI3_percent.values))
    (out / f"{f.stem}_haem5_summary.json").write_text(json.dumps(row, indent=1, default=float))
    print(f.name, "done", row["GW_HI3_mean"], row["explore_core_max_rel_err_HI3"], flush=True)
