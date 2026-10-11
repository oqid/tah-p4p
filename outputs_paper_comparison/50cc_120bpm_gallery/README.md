# Paper-inspired comparison

Open gallery.html to browse all twelve figures. discussion_metrics.csv contains
the values used for quantitative discussion. Per-case paths.csv and plot_data.npz
allow --replot without running CFD-export processing again.

## What the papers motivate

Sources: [Taskin et al.](../../papers/taskin-etal.pdf) and
[Marom et al.](../../papers/marom-etal.pdf).

Taskin et al. 2012, Figs. 4 and 7: actual path stress histories, shown here with
objective P50/P95/P99 HI3-rank path selection and accumulated HI3.
Figs. 6 and 9: model and coefficient sensitivity. We compare HI2 and HI3 with
GW, HO and TZ on identical trajectories; there is no Eulerian model or experimental
validation in these exports. Fig. 3 motivates the subset check, which is not
a substitute for new-seed or CFD discretization convergence.

Marom et al. 2014, Fig. 5: SA distribution shape, median and tail. The new SA
density uses linear axes and a shared bandwidth, with a separate tail panel.
The exceedance figure states the percentage of paths above each SA cutoff.
Figs. 4 and 6 motivate the spatial stress display. Our display projects vertices
along frozen-flow paths, not transient particle locations or a solved field contour.
Marom's transient chamber-seeded particles, flow models and operating conditions
differ from these inlet-seeded frozen-flow paths; published medians are not
validation targets for this pipeline.

## HI weighting

The primary estimate is arithmetic mean of individual path HI: sum(HI_i)/N.
The sensitivity figure additionally reports sum(v_inlet_i*HI_i)/sum(v_inlet_i).
Weighting does not change any path's damage calculation. Speed already enters
the transit-time calculation; changing operating flow can change both stress and
exposure. Taskin's flow-versus-HI trend does not establish a universal inverse
local-velocity law or determine the appropriate sampling weights.

An equal-path mean describes the sampled path population. It represents blood
throughput only if the seeding already represents equal transported amounts.
For equal-area inlet seeding, normal flux times area provides physical throughput
weights; scalar inlet speed is only a proxy under forward, normal-flow assumptions.
Double-weighting flux-proportional particle releases would be inappropriate.
All paths are included here; outlet arrival has not been checked.

HI3 is [sum(A^(1/beta)*tau^(alpha/beta)*dt)]^beta in this code's notation.
The right-hand dose-band plot partitions that additive inner dose, not final HI3.
The existing velocity floor and scalar-shear model are inherited from haem5.
Coefficients have limited calibration ranges; sensitivity is not validation.
In particular, TZ was fitted at 50-320 Pa; these paths include much lower stresses.

## 35 cc versus 50 cc discussion

First compare at 120 BPM in both geometries. Then compare a separate matched-
cardiac-output condition with the corresponding size-specific BPMs.
Compare like systolic phases and use identical seeding, stress definitions,
outlet criteria, mesh/segment resolution and model coefficients where possible.
Geometry-only causation cannot be inferred if flow/boundary conditions also differ.
Changing --bpm is metadata only: new CFD exports must embody the changed condition.
Cycle-mean cardiac output is metadata supplied by the user; nominal chamber
capacity alone does not establish delivered stroke volume or flow.

Questions for the figures:
- Is the bulk SA distribution shifted, or only its high-exposure tail?
- Does a longer transit or stronger stress explain the shift?
- Are a few paths dominating mean HI3?
- Do geometry rankings persist under HI2/HI3 and all three coefficient sets?
- Does seeding/weighting change the conclusion?
- Where are sampled high-stress regions, and do those regions move with phase?

Spatial axes retain physical millimetres. Different device coordinate systems
must be aligned upstream for anatomical location comparisons.
Rank correlations are descriptive, not causal tests. Subset bands describe the
fixed exported paths, not biological or CFD-model uncertainty.
