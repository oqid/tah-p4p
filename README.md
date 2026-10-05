# tah-p4p

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
