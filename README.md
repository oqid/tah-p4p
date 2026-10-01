# tah-p4p

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
