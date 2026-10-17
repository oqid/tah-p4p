# hemolysis_explore — extra analyses for the 50 cc → 35 cc hemolysis discussion

This folder adds analyses that the main scripts and `haem5_paper_comparison.py` don't produce yet. It uses a copy of `haem5.py` (unchanged) and its own segment-level loader, `explore_core.py`. That loader keeps every segment of every path, not only per-path totals, so it can show **where** and **when** along a path the damage builds up.

`test_explore_core.py` confirms that the loader gives the same per-path HI2, HI3, SA and transit time as `haem5.run_pipeline`. On `haemline_50cc_5000%.csv` the largest relative difference is 0 for HI3 and 1e-12 for HI2.

```powershell
cd hemolysis_explore
python explore_hemolysis.py --case "50 cc" ..\z_inputs 120
python explore_hemolysis.py --case "50 cc" ..\z_inputs 120 --completed-only
# when the 35 cc exports exist (first case = reference):
python explore_hemolysis.py --case "50 cc" ..\z_inputs 120 --case "35 cc" ..\z_inputs_35cc 120
python explore_hemolysis.py --case "50 cc" ..\z_inputs 120 --case "35 cc matched CO" ..\z_inputs_35cc_171 171
python scaling.py              # prints the similarity-scaling ratio table
python test_explore_core.py    # agreement check against haem5
```

The script needs only numpy, pandas and matplotlib (already in `.venv`). Parsed exports are cached in `cache/`, keyed on each file's hash, so a rerun takes about 20 s. File names must contain the phase, for example `_1667%` (16.67%). The systolic fraction defaults to 0.5, which gives a 250 ms systole at 120 bpm; change it with `--systolic-fraction`.

The results below come from `outputs/50cc_120bpm_all_paths/` (every exported path) and `outputs/50cc_120bpm_completed_only/` (only paths that reach the outlet). All figures use GW constants unless stated otherwise.

---

## 01 Phase mechanism: why HI3 is U-shaped through systole

| Phase | Mean GW HI3 [%] (95% CI) | Outlet-reaching paths only | Median transit [ms] | Median τ_eff [Pa] |
|---|---|---|---|---|
| 16.67 | 1.67e-4 (1.51–1.84e-4) | 1.51e-4 | 326 | 0.96 |
| 33.33 | 7.4e-5 (6.1–8.8e-5) | 4.2e-5 | 187 | 0.61 |
| 50 | 4.0e-5 (3.3–5.0e-5) | 2.3e-5 | 144 | 0.49 |
| 66.67 | 9.4e-5 (8.0–10.9e-5) | 5.9e-5 | 127 | 0.76 |
| 83 | 2.70e-4 (2.45–2.97e-4) | 2.33e-4 | 161 | 1.45 |

* **a)** HI3 is lowest at mid-systole and highest at the start and end of systole. The 83% phase is about 7× the 50% phase. GW, HO and TZ all rank the phases in the same order. The ejection-weighted value over the five phases is 1.08e-4% (8.1e-5% using only outlet-reaching paths). For this weighting, flow at each phase is taken as proportional to the mean diaphragm-inlet speed, which follows the sine profile: 0.063, 0.109, 0.126, 0.110 and 0.064 m/s.
* **b)** Across the five phases, HI3 falls as instantaneous flow rises (log–log slope about −2). This runs against the "more flow, more damage" expectation in your early notes.
* **c)** Each path obeys the exact identity **HI3 = A·T^β·τ_eff^α**. Taking logs splits the change from mid-systole exactly into a transit-time term and a stress term:
  - At 83%, the stress term is about ×12 and the transit term only about ×1.1.
  - At 16.67%, the stress term is about ×4.6 and the transit term about ×1.9.
  - **So the phase trend is mainly a stress effect, not a residence-time effect.**

  Stress is higher at early and late systole even though flow is lower there. That suggests the valve configuration at each phase (a partly open disc means a narrower gap and a stronger jet) matters more than the diaphragm speed. Can you confirm the valve angle or gap in each phase's CAD? (See the questions at the end.)

## 02 Stress–time map (τ_eff vs transit, with iso-HI3 lines)

This differs from figure 05 in `haem5_paper_comparison`, which uses time-averaged stress. Here the y-axis is τ_eff, the stress that, combined with the path's transit time, reproduces that path's HI3 exactly. That means the dashed iso-HI3 lines can be read directly off the plot.

* For a typical path, τ_eff is about 0.5–1.5 Pa.
* The damage-weighted medians (hollow markers) sit at τ_eff of about 6–10 Pa with similar or only slightly longer transit times. **Predicted damage comes from a high-stress minority of paths, not from long-residence paths.**
* The vertical blue line marks the time left in systole after each snapshot.
* The arrows show where the 50% median would move in a 35 cc device under the scaling hypotheses (see 06).

## 03 Damage timing, and why the frozen-flow assumption is less damaging than it looks

Many paths take longer to exit than the time left in systole:

| Phase | Paths taking longer than the remaining systole |
|---|---|
| 16.67% | 84% |
| 33.33% | 59% |
| 50% | 62% |
| 66.67% | 92% |
| 83% (42 ms left) | 100% |

So most paths are not physically possible within one systole. However:

* **a)** For all paths, **67–98% of the linear dose is deposited in the final 25 ms before exit**, in the valve and outflow tract. For outlet-reaching paths only, the figure is 26–97%; their dose builds up over the final ~50–100 ms instead.
* **b)** If you count only the dose deposited within the time remaining in systole, mean HI3 keeps **≥ 90% of its full value at every phase**, and ≥ 99% at 16.67% and 83%. The slow, unphysical early part of each path happens at low stress, so it contributes almost nothing.

This is a good defence of the quasi-steady method in the report. The transit time is not physically consistent with a single systole, but the predicted damage hardly depends on that inconsistency.

## 04 Damage map: where the dose is deposited

These are 1 mm projected maps of where the linear HI3 dose is deposited, in the X–Z view (comparable to `07_spatial_path_stress`) and the X–Y view. They come with a dose-concentration curve and `hotspots.csv`, which lists the 20 highest-dose 3D voxels per phase so you can place probe points in CFD-Post.

* **50% of the dose sits in 0.03–0.45% of the volume the paths visit**, or 0.25–0.77% using only outlet-reaching paths.
* The chamber body contributes essentially nothing. The dose is concentrated around the outlet valve and the outflow tract, roughly X 0–20 mm, Z −3 to 13 mm.

This supports a design discussion focused on the valve. Mean stress over the whole chamber is the wrong metric here.

## 05 Numerical audit: important for the report

* **a–b) Some paths never reach the outlet, and they carry a large share of HI3.** The outlet is found as a plane fitted to the path end points; paths that reach it stop exactly on it. Paths that end more than 1 mm off that plane stopped inside the domain:

  | Phase | Paths that stop inside the domain | Share of summed HI3 they carry |
  |---|---|---|
  | 16.67% | 9.6% | 18% |
  | 33.33% | 3.4% | 45% |
  | 50% | 1.4% | 45% |
  | 66.67% | 3.5% | 39% |
  | 83% | 9.0% | 21% |

  Their end points form a ring about 16 mm across, which looks like the valve seat or disc edge. At their end points they are slow (median about 0.3 m/s, compared with 1.3 m/s at the outlet) and at high stress (median about 49 Pa, compared with 0.4 Pa). These are most likely streamlines that run into the wall or the gap and stop.

  How to handle them in the report:
  - Present both versions (the `all_paths` and `completed_only` runs). Removing them lowers mean HI3 by 10–45%, but **the U-shape and the phase ranking remain**.
  - Say that their HI3 is uncertain: the exposure is probably real (blood near the valve), but the path is cut short.
* **Velocity floor:** never active (under 3e-6 of total time). The 1% speed floor does not affect any result here.
* **c) Stress sampling:** 2–9% of the dose sits in segments where τ changes by more than 2× between consecutive points, with the most at 50%. That is acceptable, but it is the next thing to refine (more export segments near the valve).
* **d)** *Corrected:* an earlier version of this panel plotted densities on log-spaced bins, which exaggerated the short segments. Weighted by dose, the median segment length is about 1 mm at most phases (0.1 mm at 50%), similar to the overall median of 1.2–1.45 mm. The export is therefore no finer where the damage accumulates, and refining the segment count near the valve is the next numerical check.
* **Sampling spread:** the bootstrap 95% CI on mean HI3 is about ±10–25%, because the distribution is heavy-tailed. A test run comparing 50 cc with itself gave ratio CIs of about 0.7–1.4. **A 35 cc vs 50 cc difference smaller than about ±40% at a single phase cannot be distinguished from path-sampling noise.** Report the geometric mean or median alongside the mean.

## 06 Scaling hypothesis for 35 cc (testable before the CFD exists)

Assumptions:
* Uniform geometric similarity, with every length scaled by (35/50)^(1/3) = 0.888, including the diaphragm stroke and the valve gap.
* The same systolic fraction, and the same valve state at the same phase.
* Stress is bounded between viscous scaling (τ ∝ μU/L) and inertial scaling (τ ∝ ρU²).

| Scenario | U | T | Re | Predicted HI3 ratio, GW | Predicted HI3 ratio, HO/TZ |
|---|---|---|---|---|---|
| Same rate, 120 bpm (CO ×0.70) | ×0.89 | ×1.0 | ×0.79 | **0.56–1.0** | 0.62–1.0 |
| Matched CO, 171 bpm | ×1.27 | ×0.70 | ×1.13 | **1.8–2.4** | 1.5–2.0 |

For comparison, Marom et al. 2014 report a 35/50 median-SA ratio of 1.56, though at 164 and 140 bpm respectively.

The framing for the report: at the same rate, downsizing should not increase damage per pass, but it reduces output by 30%. Recovering that output means a faster rate, which the scaling predicts roughly doubles HI3 per pass. When the 35 cc results arrive, figure 07 overlays these bands on the measured ratios. If the measured ratios fall outside the bands, that is itself a finding: the geometry or valve state did not scale uniformly.

## 07 Size comparison (generated automatically once a second case is supplied)

For each phase, this figure shows:
* the 35/50 ratio of mean HI3, with a bootstrap 95% CI and the scaling band behind it;
* the same log split as figure 01c, attributing the size effect to transit time or to stress;
* overlaid ECDFs of per-path HI3.

`size_comparison.csv` contains the numbers.

## Output files

`phase_summary.csv` has one row per phase with every number quoted above, including the GW/HO/TZ mean, CI, geometric mean, P50/P95/P99, NIH conversion, completeness and timing metrics. The other files:

* `paths_*.csv`: per-path T, HI3 (all constant sets), τ_eff, SA, peak stress, end distance from the outlet plane, and the completed flag
* `numerical_audit.csv`
* `hotspots.csv`
* `scaling_ratios.csv`
* `run.json`

## Report figures and operating points

* `report_figures.py` (styled by `p4p_report.mplstyle` + `report_style.py`; `--only F2` rebuilds one figure) writes the curated report set to `report_figs/` (F1–F6, PNG + PDF), and `operating_points.csv` (inflow Q = v·A_dome, valve angle, CFD P99 speed and HI3 for each phase).
* `FIGURE_GUIDE.md` gives captions, the claim each figure supports, and the methods notes.
* `operating_point.py` also contains a simple valve-gap "jet proxy". It was dropped from the figures: the real gap geometry (an 18 mm orifice, with flow on both sides of the disc) makes it too crude.

Key result once the valve angles were known: HI3 falls as inflow rises, but tracks the valve opening angle. The narrowest openings (≈25°) give the highest HI3, and the CFD's own peak path speeds rank the phases the same way.

## Notes

* CFD-Post ends streamlines on hard wall contact. Paths that stop at the valve ring therefore have a truncated (lower-bound) exposure, so results are shown with and without them.

## 35 cc comparison (17 Oct 2026)

```powershell
cd hemolysis_explore
python run_haem5_batch.py outputs\35cc_haem5 ..\35cc_inputs\120BPM\*.csv ..\35cc_inputs\171BPM\*.csv
python explore_hemolysis.py --case "50 cc" ..\z_inputs 120 --case "35 cc 120 BPM" ..\35cc_inputs\120BPM 120 --case "35 cc 171 BPM" ..\35cc_inputs\171BPM 171 --outdir outputs\35cc_comparison_all_paths
python cardiac_output_cases.py outputs\cardiac_output        # optional: folder_50 folder_35_120 folder_35_171
python compare_figures.py --case "50 cc 120 BPM" ..\z_inputs 120 1.0 --case "35 cc 120 BPM" ..\35cc_inputs\120BPM 120 0.888 --case "35 cc 171 BPM" ..\35cc_inputs\171BPM 171 0.888 --monitor "..\35cc_inputs\Shear Strain Rate Convergence Oscillation 35cc 83% 120BPM\variableoutput.csv"
```

* `outputs/35cc_RESULTS.md` collects every number from this run.
* The loader now reads phases from names like `35cc 120BPM 16.67%.csv`. It snaps the 50 cc `_8300%` file to 83.33%.
* Flow and cardiac output come from the inlet (seed speed × dome area). The Sim Tracker's 35 cc outlet integrals violate mass conservation, so they are not used.

### Per-case report figures and the three-case figures

```powershell
# full F1–F7 + A1–A5 set for each 35 cc case (--scale sets the geometry in F1 and the valve diameter)
python report_figures.py ..\35cc_inputs\120BPM --label "35 cc 120 BPM" --bpm 120 --scale 0.888 --outdir report_figs\35cc_120BPM
python report_figures.py ..\35cc_inputs\171BPM --label "35 cc 171 BPM" --bpm 171 --scale 0.888 --outdir report_figs\35cc_171BPM
```

`compare_figures.py` writes into `report_figs`:

| Figure | File | What it shows |
|---|---|---|
| C1 | `size_comparison`, `size_comparison_completed` | 35/50 HI3 ratio per phase with CIs and the similarity bands |
| C2 | `size_dose_maps` | dose maps of all three cases side by side |
| C3 | `size_lysis_sa` | SA distributions and lysis safety factors |
| C4 | `convergence_monitor` | the 35 cc solver oscillation |
| C5 | `size_mechanism` | the F2 b/c mechanism plot for all three cases: HI3 against inflow and against valve angle |
| C6 | `size_lysis_envelopes` | the F7 empirical-lysis plot with one panel per phase and all three cases |

Use `--only C5` (for example) to rebuild just one of them.
