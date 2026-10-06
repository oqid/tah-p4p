# tah-p4p

## Main analysis entry point

Use `run_haem.py` for routine investigations. It groups selected analyses in one
new folder under `outputs/runs/`, with a README, configuration/status record,
per-case console logs and combined case tables. Existing run folders are refused.

```powershell
# Default: printed core summary plus detailed diagnostics; no plots
python run_haem.py 1000p_100_alotpa_3djawn.csv

# Select extra analyses; plots are opt-in
python run_haem.py export.csv --analyses summary diagnostics empirical plots

# Named investigation with per-case geometry, phase, BPM, outlet and weights
python run_haem.py --manifest cases.json --outdir outputs/runs/50cc_vs_35cc

# Quick core summary and GW/HO/TZ comparison
python run_haem.py export.csv --analyses summary --compare-constants
```

`--manifest` uses the format documented in [HAEM_DIAGNOSTICS.md](HAEM_DIAGNOSTICS.md).
You can also pass several export paths. `--bpm`, `--expected-seeds`,
`--outlet-plane` and `--outlet-tolerance` provide defaults; manifest settings
override them. BPM is metadata here, not a change to the CFD field.

The core summary now prints the mean, P95, P99 and maximum of per-streamline
peak stresses, and saves each line's peak and transit time in `streamlines.csv`.
These are sampled stress peaks, not a sensitization prediction.

Volume flow rate and cardiac output are available in `haem5.py` and the runner:

```powershell
python run_haem.py export.csv --analyses summary --inlet-area-m2 0.0003 --stroke-volume-ml 50 --bpm 120
python haem5.py export.csv --mean-flow-rate-l-min 6
```

`--inlet-area-m2` estimates Q = area × mean inlet seed speed, reported in m³/s,
mL/s and L/min. This assumes seeds uniformly represent the entire inlet area,
start at the inlet and follow flow normal to that surface. Scalar speeds cannot
resolve reverse or oblique flow. A snapshot estimate is not cycle-mean cardiac
output. Cardiac output uses net delivered `--stroke-volume-ml` × `--bpm` / 1000,
or a supplied cycle-mean `--mean-flow-rate-l-min`. Chamber capacity alone is not
net stroke volume. Without inlet area, either cardiac output input also supplies
cycle-mean Q. Missing quantities are shown as unavailable. Manifest fields are
`inlet_area_m2`, `stroke_volume_ml`, and `mean_flow_rate_l_min`; the runner saves
these metrics in `summary.csv`, and haem5 saves them in `device_summary_*.csv`.

Core summary, empirical screening and plots use all reconstructed paths and
inlet-speed proxy weighting. Diagnostics can select completed paths and use
supplied flux weights; its report states the population and weighting. Use
diagnostics results for those controlled comparisons. `summary.csv` and
`comparison.csv` deliberately distinguish the two populations.

Use the scripts directly for specialist options: `haem5.py` for CAD overlays,
`haem_diagnostics.py` for numerical audits, `haem5_empiricalthresholds.py` for
standalone lysis screening and extra figure formats, and `sweeps/` for new-seed
convergence. `old/` and `haem5_recirculationtest.py` are not part of this runner.
The runner creates PNGs only when `plots` is selected; it does not duplicate PDF
variants. No existing scripts, inputs or results are moved.


| Script | Role |
|---|---|
| `haem5.py` | Core haemolysis/SA calculation, summary and plots |
| `haem_diagnostics.py` | Detailed numerical checks, path selection, weighting and comparison tables |
| `haem5_empiricalthresholds.py` | Separate empirical lysis screening |
| **`run_haem.py`** | Calls your selected analyses and groups their outputs and logs |

## Helpers

### Sweep

Sweep scripts, input exports, results and plots live in `sweeps/`. From the
project root, check pending inputs, run missing cases, then plot without TZ:

```powershell
python sweeps/run_sweep.py --resume --constants GW HO --dry-run
python sweeps/run_sweep.py --resume --constants GW HO
python sweeps/plot_sweep.py --constants GW HO
```

These commands also work from inside `sweeps/` using `python run_sweep.py`
and `python plot_sweep.py`. Default paths are relative to the scripts, rather
than the working directory. `sweep_N100000.csv` is included automatically;
`--resume` preserves existing results and skips completed constant sets.


### .STP Visualizer for Plots

To overlay a fluid-domain model on the 3D streamline figure:

```powershell
python -m pip install -r requirements-domain.txt
python haem5.py export.csv --fluid-domain fluid_domain.stp
```

Both `.stp` and `.step` are supported. The domain appears as a faint surface
with CAD edge outlines, receives the same rotation as the streamlines, and is
included in automatic plot framing. STEP units are converted to metres; the
model must use the same origin and coordinate axes as the CFD export.
`--domain-scale` applies an additional scale factor if needed.
Use `--domain-opacity 0.1` to adjust surface opacity or `--no-domain-border`
to hide the outlines. Omitting `--fluid-domain` keeps the usual streamline plot
and requires no CAD dependency. The overlay does not change the calculations.
