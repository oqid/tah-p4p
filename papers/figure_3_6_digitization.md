# Figure 3.6 curve provenance

Source: *Guidelines for Blood-Material Interactions*, NIH/NHLBI,
Fig. 3.6, printed p. 78 / PDF p. 87. Discussion: printed p. 77.
The caption attributes the figure to a modification of reference 44:
Hellums JD and Hardwick PA (1981), *Response of Platelets to Shear Stress —
A Review*, in *The Rheology of Blood, Blood Vessels and Associated Tissue*, p. 160.

`figure_3_6_digitization.csv` records manually selected curve centre points,
not original experimental measurements. The page was rendered using PyMuPDF
at a scale of 1.7, giving a 1027 × 1325 image. Calibration:

- x = 320 at 10^-6 s; x = 808 at 10^2 s.
- y = 268 at 10^6 dyn/cm²; y = 739 at 10^2 dyn/cm².
- time_s = 10^(-6 + 8 (x_px - 320) / 488).
- stress_Pa = 0.1 × 10^(6 - 4 (y_px - 268) / 471).

The scan has thick lines, a slightly tilted frame, and overlapping labels.
Sparse tracing and log-log interpolation introduce additional uncertainty.
No precision or uncertainty bound is implied by the floating-point conversion.
There is no fitted analytical equation in the figure. Points outside the traced
range are excluded from numerical comparisons rather than extrapolated.

Only red-cell and platelet curves are traced. PMN leukocytes and the
surface-induced-damage region are not modeled here.

The figure concerns significant cell lysis under shear exposure. The text on
printed p. 77 also discusses activation and sublytic effects; being below a
lysis curve does not establish absence of those effects. The longest-contiguous
CFD envelope is a screening construction in this code, not a model specified
or validated by the Guidelines, and ignores separated repeated exposures.
