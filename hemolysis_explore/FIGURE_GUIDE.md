# Figure guide: which figures go in the report, and what each one argues

Your three favourites are the backbone of the set:
* `04_damage_map_xz` and `07_spatial_path_stress` → combined into **F3**
* `01_phase_mechanism` → rebuilt as **F2**, now using the valve angles

All report figures are 16 cm wide (page text width), with PNG at 300 dpi plus a vector PDF, in `report_figs/`. Regenerate them with:

```powershell
cd hemolysis_explore
python report_figures.py ..\z_inputs                 # all six
python report_figures.py ..\z_inputs --only F2       # just one, while editing
```

Where to edit:
* **Look** (fonts, sizes, grid, dpi): `p4p_report.mplstyle`.
* **Colours, markers, figure width:** `report_style.py`.
* **Content** (titles, labels, axis limits, label offsets, which phases F3 shows): the figure's own function in `report_figures.py`: `f1_model_states`, `f2_mechanism`, `f3_spatial`, `f4_timing`, `f5_completion`, `f6_scaling`.
* **Speed:** the first run parses the CSVs (about 1 min); later runs use `cache/` and take a few seconds.

The suggested order follows the argument:

**operating point → mechanism → location → method check → 35 cc prediction**

---

## F1 — The model through systole (Methods)

**Claim:** each snapshot is a different steady model. The diaphragm dome rises and its inlet speed follows a sine, while the outlet-valve disc opens and then closes again.

**Suggested caption:**
> The five steady systolic snapshots (50 cc, 120 bpm), drawn to scale from the CAD parameters.
> Top: diaphragm (velocity inlet) profile, the half-sine dome y = A·φ·cos(πr/L) with L = 56.2 mm and A = 20 mm. Grey lines show the other phases. Arrows show the prescribed normal inlet speed v = πS/(2T_s)·sin(πφ), drawn with length proportional to speed.
> Bottom: outlet-valve disc angle θ (0° = closed) in the 18 mm orifice. The disc pivot is drawn schematically.

**What to say:**
* Flow peaks at mid-systole, but that is also when the valve is most open (75°).
* Early and late in systole the flow is about half the mid-systole value, but it passes a disc only about 25° open.
* The integrated inflow is 55.5 mL per systole, or 6.7 L/min at 120 bpm (`operating_points.csv`).

## F2 — Mechanism (the centrepiece of the Results)

**Claim:** neither flow rate nor residence time explains the U-shaped HI3 through systole. It follows how far the outlet valve is open, through the stress it creates.

**Suggested caption:**
> (a) Mean per-path HI3 (Giersiepen–Wurzinger constants, 95% bootstrap CI) for all exported paths and for outlet-reaching paths only.
> (b) HI3 against instantaneous inflow. (c) HI3 against valve opening angle. In (b) and (c), a line joins the phases in order, from 16.67% to 83%; the arrowhead marks the direction.
> (d) Exact log decomposition of the change in geometric-mean HI3 relative to 50% systole. Since HI3 = A·T^β·τ_eff^α on every path, the change splits into a transit-time term β·Δ⟨ln T⟩ and a stress term α·Δ⟨ln τ_eff⟩.

**What to say:**
* As systole progresses, flow rises and then falls while HI3 falls and then rises (b). HI3 is highest at the lowest flows.
* Against valve angle (c), the two branches line up: the valve is narrowest at 16.67% and 83% and HI3 is highest there.
* The closing branch (66.67%, 83%) sits about 1.3–1.6× above the opening branch at the same angle. The angles match closely, so the difference comes from what changes alongside: the diaphragm is further expanded, and the inflow speed and chamber shape differ. It is worth a sentence.
* The stress term accounts for ×4.6 at 16.67% and ×12 at 83%. Transit time adds at most ×1.9 (d).
* Supporting check from the CFD itself: the 99th-percentile path speed ranks the phases the same way (83% 3.2 m/s > 16.67% 2.5 > 66.67% 2.1 > 33.33% 1.8 > 50% 1.6 m/s).
* HO and TZ constants rank the phases identically (`outputs/.../01_phase_mechanism`).

## F3 — Where the stress and the damage are (your favourite)

**Claim:** high stress follows the paths into the outflow tract, and nearly all of the damage is deposited in a small region at the valve.

**Suggested caption:**
> Top: local scalar stress along 200 randomly sampled paths at early, mid and late systole (X–Z projection). Bottom: linear HI3 dose, A^(1/β)·τ^(α/β)·dt, deposited per 1 mm column, as a share of each phase's total. Grey = region visited by paths.

**What to say:**
* Half of the dose sits in 0.03–0.45% of the visited volume, at the valve ring and in the outflow tract. The chamber body contributes essentially nothing.
* At 25° opening the dose forms a ring at the valve. At 75° it spreads into the outflow jet.
* The X–Y version and all five phases are in `outputs/.../04_damage_map_*.png` (Appendix).
* `hotspots.csv` gives voxel coordinates if you ever want a CFD-Post probe there. It isn't needed for the argument.

## F4 — Damage timing → why the quasi-steady approach is defensible (Discussion)

**Claim:** paths are slow to cross the chamber, but damage accumulates only in the last few tens of milliseconds. So violating the frozen-flow assumption over the long, slow part of a path barely changes HI3.

**Suggested caption:**
> Outlet-reaching paths only.
> (a) Cumulative share of linear dose against time remaining before the path exits.
> (b) Mean HI3 counting only the dose deposited within a final window, relative to the full-path value.
> Markers show the time left in systole after each snapshot, T_s(1 − φ).

**What to say:**
* 59–100% of paths take longer to exit than the time remaining in systole, so a frozen-field path is not physically realisable.
* Yet ≥ 90% of each phase's mean HI3 is accumulated within that remaining time.
* This is the honest way to state the limitation and bound its effect.

## F5 — Path completion (Methods verification / Limitations)

**Claim:** a small number of streamlines stop at the valve ring. They carry a disproportionate share of HI3, so results are reported with and without them.

**Suggested caption:**
> (a) Path end points. Grey points end on the outlet plane. Coloured markers end more than 1 mm from it, sized by path HI3.
> (b) Share of paths that stop early, and their share of summed HI3.

**What to say:**
* 1.4–9.6% of paths end on a ring about 16–20 mm across at the valve seat or disc edge. There the flow is slow (about 0.3 m/s) and the stress high (median about 49 Pa).
* They carry 18–45% of summed HI3.
* Excluding them lowers mean HI3 by 10–45%, but leaves the phase trend unchanged (F2a).
* Cause, confirmed: CFD-Post terminates a streamline when it makes hard contact with a wall. These paths hit the disc or seat while passing through the narrow valve gap.
  - Their exposure up to that point is real.
  - Their downstream exposure is missing, so their HI3 is a lower bound for those paths.
  - The high-stress near-wall contact also inflates their last segments.
* This is why the main figures report both versions.

## F6 — 35 cc hypothesis (end of Results or Discussion; replace with measured data later)

**Claim:** similarity scaling makes a falsifiable prediction for 35 cc. At the same rate, damage per pass stays the same or falls, but output drops. Restoring output roughly doubles damage per pass.

**Numbers for the text:**
* Every length is scaled by 0.888.
* At 120 bpm the inlet BC gives 38.9 mL per systole, or 4.66 L/min (×0.70).
* Matching 6.66 L/min needs 171 bpm.
* Predicted GW HI3 ratio: 0.56–1.0 at 120 bpm, and 1.8–2.4 at 171 bpm (1.5–2.0 with HO/TZ).
* For context, Marom et al. found a 1.56× higher SA median for 35 cc than 50 cc.

When the 35 cc runs exist, `explore_hemolysis.py` with two `--case` arguments produces `07_size_comparison` with these bands overlaid.

---

## F7 — Stress–duration screening against empirical lysis curves (Results, supporting, or Appendix)

**Claim:** the HI3 power law is not the only test. Measured against independent empirical lysis limits, every phase stays below both curves with margin.

This is the restyled version of `outputs_snapshots\systole_comparison\empirical_stress_duration_all_snapshots.png`. It uses the same calculation as `haem5_empiricalthresholds.py`, now vectorised in `empirical_envelope.py`; a check against the original loop gave identical durations. It writes `F7_envelope_<phase>.csv` and `F7_lysis_safety_factors.csv`.

**Suggested caption:**
> Worst-case stress–duration exposure envelope for each snapshot: for each stress level, the longest contiguous time any path spends at or above it. Shown against approximate traces of the NIH red-cell and platelet lysis boundaries (Guidelines for Blood–Material Interactions, Fig. 3.6).
> (a) Full reference range. (b) Detail.
> SF = minimum safety factor, the critical stress at that duration divided by the envelope stress.

**What to say:**
* Minimum safety factors are 3.3–5.1 against red-cell lysis and 7.5–12.3 against platelet lysis. The tightest is at 50% systole for red cells and 66.67% for platelets.
* The envelopes are similar across phases. The highest stresses (about 100–200 Pa) last only about 0.01–0.1 ms; long exposures (around 1 s) occur only below about 1 Pa.
* This complements F2: the phase trend in HI3 comes from how much of the flow passes through the high-stress region, not from any phase approaching a lysis limit.
* Caveats:
  - The traces are manual digitisations.
  - The method treats each contiguous interval as a sustained stress and ignores repeated exposures.
  - These are lysis boundaries; platelet activation can occur below them.
  - It uses all paths, including those that stop at the valve, as a worst case.

## Appendix candidates (already generated)

* `02_stress_time_map`: per-path τ_eff vs transit with iso-HI3 lines. A good bridge between F2 and F3 if you have room.
* `05_numerical_audit` c–d: stress-sampling and segment-length checks.
* The full five-phase `04_damage_map_xz` / `_xy`.
* From `haem5_paper_comparison`: SA density (01), Taskin model sensitivity (03) and subset stability (11).

## Methods notes from your answers (worth one sentence each in the report)

* **Diaphragm:** an Inventor parametric half-sine, x = 2Lt and y = A·Pump·sin(2πt) for t ∈ [0, 0.5], revolved, with L = 56.211 mm and A = 20 mm.
  - Pump equals the systolic fraction, so the dome *height* is linear in phase.
  - The inlet *speed* follows the sinusoidal displacement law. Mention both.
* **The inlet BC, not the CAD sweep, sets the flow.** Over the full stroke the dome's own volume change is only 23.0 mL. The velocity BC, applied normal to the whole dome surface, delivers 55.5 mL per systole.
  - That is fine for steady snapshots, but say so explicitly.
  - It is also why the delivered output (6.7 L/min) is close to the nominal 50 cc × 120 bpm.
* **Solver:** steady SST, physical timescale 5×10⁻⁴ s, converged after 1,234 accumulated pseudo-timesteps, no-slip walls.
  - Stress = molecular viscosity × strain rate. Reynolds stresses are not included.
  - Note that this likely *under*-estimates the stress experienced in the turbulent jet. That is acceptable for a relative 50/35 comparison, but it matters for absolute NIH.
* **Systole** is assumed to be 50% of the cycle (250 ms at 120 bpm).
