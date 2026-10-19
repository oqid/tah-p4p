> **SUPERSEDED (17 Oct 2026):** the 16.67 % numbers here come from the unscaled-geometry runs, and flows are analytic. Use report_figs/size_summary*.csv, size_lysis_margins.csv and outputs/cardiac_output (massFlow-based) instead.

# 35 cc vs 50 cc results (run 17 Oct 2026)

**Cases**
* 50 cc at 120 BPM (`z_inputs`).
* 35 cc at 120 BPM and at 171 BPM (`35cc_inputs/120BPM`, `171BPM`). Every length is scaled by 0.888, and the valve angles are the same at each phase.
* About 4,850–4,920 paths per phase.

**Constants and averaging:** GW constants unless stated. Averages are taken over all exported paths unless stated.

**How to regenerate**
* `run_haem5_batch.py`
* `explore_hemolysis.py` (three cases)
* `cardiac_output_cases.py`
* `compare_figures.py`

The command lines are in the README.

## Operating point (inflow = diaphragm speed × dome area; mass conservation)

| Case | Q at each phase (L/min) | Stroke volume (mL) | Cardiac output (L/min) |
|---|---|---|---|
| 50 cc 120 BPM | 9.4, 16.7, 20.6, 19.9, 12.5 | 54.9–56.4 | 6.59–6.76 |
| 35 cc 120 BPM | 6.8, 11.7, 14.5, 13.9, 8.9 | 38.7–39.8 | 4.65–4.77 (×0.70) |
| 35 cc 171 BPM | 9.4, 16.7, 20.7, 19.9, 12.5 | 38.5–39.6 | 6.59–6.77 (×1.00) |

Ranges span trapezoid and PCHIP interpolation.

The **tracker's outlet integrals for 35 cc are not used.**
* They exceed this inflow by 8–106%, which is impossible in a converged steady incompressible solution.
* For 50 cc the tracker agrees with the inflow to within 3%.
* Using the tracker values would give 6.05–6.33 L/min at 120 BPM and 8.7–9.1 L/min at 171 BPM.

## Predicted haemolysis per pass

Each ratio is 35 cc divided by 50 cc at the same phase, with a bootstrap 95% CI in brackets.

| Phase | 35 cc 120 BPM | 35 cc 171 BPM |
|---|---|---|
| 16.67 | 0.99 (0.87–1.14) | 1.83 (1.59–2.13) |
| 33.33 | 0.72 (0.58–0.90) | 1.43 (1.12–1.80) |
| 50.00 | 0.94 (0.73–1.23) | 1.87 (1.43–2.45) |
| 66.67 | 1.01 (0.81–1.26) | 1.94 (1.56–2.44) |
| 83.33 | 1.01 (0.87–1.16) | 1.86 (1.60–2.14) |
| Prediction (GW) | 0.56–1.0 | 1.79–2.38 |

* **Outlet-reaching paths only:** 0.84–1.05 at 120 BPM and 1.62–2.14 at 171 BPM.
* **Other constant sets at 171 BPM:** HO gives 1.31–1.63 (predicted 1.55–1.96). TZ gives 1.38–1.72 (predicted 1.61–2.04).
* **H1 is supported:** at the same rate, downsizing did not increase damage per pass. Four of the five phases are indistinguishable from 1. 33.33% is lower, at 0.72; its CI excludes 1 for all paths but only just for outlet-reaching paths.
* **H2 is supported:** at matched output, damage per pass roughly doubled. The geometric-mean ratio is 1.58–1.83, and the arithmetic means are 1.43–1.94.

## Mechanism (exact log decomposition, geometric means)

| | Transit-time ratio | Effective-stress ratio |
|---|---|---|
| 35 cc 120 BPM, measured | 1.00 at 33.33–83.33% (1.09 at 16.67%) | 1.01–1.03 at 33.33–83.33% (0.94 at 16.67%) |
| 35 cc 120 BPM, predicted | 1.00 | 1.00 (viscous) or 0.79 (inertial) |
| 35 cc 171 BPM, measured | 0.70 at 33.33–83.33% (0.79 at 16.67%) | 1.42–1.44 at 33.33–83.33% (1.31 at 16.67%) |
| 35 cc 171 BPM, predicted | 0.70 | 1.43 (viscous) or 1.61 (inertial) |

* **The flow scaled exactly as viscous similarity predicts.** The stress increase drives the rise at 171 BPM: the stress term is +0.37–0.39 in log10, against −0.12 for transit time.
* **16.67% departs from this somewhat in both 35 cc cases.** The 35 cc inlet speed there was 3% above exact scaling (the tracker used f = 0.172), but that alone doesn't explain the gap.

## Where the damage occurs, and the paths that stop early

* **Location is unchanged:** dose maps (`size_dose_maps`) show the same valve ring and outflow tract. Half the dose sits in 0.37–0.49% of the visited volume at 16.67% and 83.33%, and 0.05% at 50%.
* **Early-stopping paths:** 1.4–9.6% of paths end at the valve in 50 cc and 1.4–8.9% in 35 cc (both 35 cc cases). They carry 17–45% of summed HI3, which is similar across cases.

## Platelets and lysis

* **SA medians are unchanged by size or rate**, at 0.39–0.73 dyn·s/cm² (SA median ratio 1.00–1.06).
* **Minimum margin to RBC lysis:** 3.3–5.1 (50 cc), 3.4–6.2 (35 cc 120 BPM), 2.5–4.0 (35 cc 171 BPM). The lowest is 2.5 at 83.33% for 171 BPM, where the peak path stress is 348 Pa.
* **Minimum margin to platelet lysis:** 7.5–12.3, 7.6–16.3 and 6.1–12.6 for the same three cases.

## Combined over systole (flow-weighted, Eq. hi-qs)

| | GW NIH (mg/100 L) | HO | TZ |
|---|---|---|---|
| 50 cc | 11.8 (unweighted 13.5) | 0.19 | 2.09 |
| 35 cc 120 BPM | 11.4 (unweighted 13.1) | 0.18 | 2.04 |
| 35 cc 171 BPM | 21.4 (unweighted 24.6) | 0.29 | 3.37 |

* **Haemoglobin released per minute** is proportional to HI × output. Relative to 50 cc it is about 0.68 at 120 BPM and about 1.81 at 171 BPM.

## Solver steadiness (35 cc, 83.33%, 120 BPM)

* The area-averaged shear strain rate at the valve oscillates by ±0.7% (SD 0.28%) about 5,033 1/s, with a period of about 1,330 steps.
* Over the last 20,000 steps the running mean is flat.
* With GW constants, ±0.7% in stress gives about ±1.7% in HI3, far below the ±40% path-sampling resolution.
