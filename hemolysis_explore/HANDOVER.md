# Handover: writing the 35 cc results and discussion (as of 17 Oct 2026)

Read this first, then `claude/35cc-results.md` (every number) and `claude/report-plan-notes.md` (to-do list and the supervisor's request). This note covers what exists only in the previous chat: the state of Overleaf, the decisions made, inconsistencies still to fix, table-ready numbers and the writing plan.

---

## 1. Where everything lives

**On the user's computer** (`C:\Users\kingk\Documents\tah-p4p`, reached through the linked-computer tools):

* `hemolysis_explore\` is the analysis package. It contains:
  * scripts: `explore_core.py`, `report_figures.py`, `compare_figures.py`, `explore_hemolysis.py`, `empirical_envelope.py`, `operating_point.py`, `cardiac_output_cases.py`, `run_haem5_batch.py` and the user's own `estimate_cardiac_output.py`
  * the style files `report_style.py` and `p4p_report.mplstyle`
  * `README.md`, which has the commands
  * `FIGURE_GUIDE.md`, which has captions and talking points for every figure, including C1–C6
* `hemolysis_explore\report_figs\` holds the 50 cc figures F1–F7 and A1–A5, plus the three-case figures:
  * `size_comparison`
  * `size_comparison_completed`
  * `size_dose_maps`
  * `size_lysis_sa`
  * `convergence_monitor`
  * `size_mechanism`
  * `size_lysis_envelopes`
* `report_figs\35cc_120BPM\` and `report_figs\35cc_171BPM\` hold full per-case F1–F7 and A1–A5 sets.
* `report_figs\size_summary.csv`, `size_summary_completed_only.csv` and `size_lysis_margins.csv` are the source data for the tables below.
* `hemolysis_explore\outputs\cardiac_output\flows_by_phase.csv` and `cardiac_output.json` hold the flows and cardiac output.
* `hemolysis_explore\overleaf_upload\figures\` holds every PDF/PNG the report needs. **The Overleaf connection is text-only**, so the user uploads figure files to the Overleaf `figures/` folder by hand. Tell them which ones are new.
* Inputs:
  * `z_inputs\` (50 cc)
  * `35cc_inputs\120BPM`, `35cc_inputs\171BPM`
  * the 35 cc solver monitor `35cc_inputs\Shear Strain Rate Convergence Oscillation 35cc 83% 120BPM\variableoutput.csv`
* House rules:
  * Deleting files in the user's folders is off. Move unwanted files to `tah-p4p\_to_delete\` instead.
  * The connection refuses files over 30 MB, so larger transfers go as split zips that are extracted on the computer.
  * Never print the git token in `overleaf_config.json`.

**Overleaf project "Final Year Report"** (through the `overleaf__*` tools). The relevant files are listed below.

| File | State |
|---|---|
| `sections/methods/methods.tex` | Written. Includes the scaling-hypothesis subsection (`subsec:scaling-hypothesis`), the post-processing sensitivity (`subsubsec:post-processing-sensitivity`), the geometry (`subsubsec:geometry-35cc`) and the equations `eq:taueff`, `eq:hi-qs`, `eq:nih`, `eq:motion`. |
| `sections/results/00-results-overview.tex` | Written. Already points to `subsec:res-35cc`. |
| `sections/results/01-haemodynamics.tex` | 50 cc flow, written. Table `tab:haemodynamics`. **Ends with a commented-out note proposing CFD-Post contours. This is the natural home for the supervisor's screenshots (§6).** |
| `sections/results/02-haemolysis-50cc.tex` | 50 cc, fully written. Tables `tab:hi-phase` and `tab:lysis`; figures `fig:F2` and `fig:F3`. |
| `sections/results/03-size-comparison.tex` | **Skeleton only: `\tbd{}` placeholders, the table `tab:size-comparison` and the figure `fig:size-comparison` (`size_comparison.pdf`).** The header comment and the "jet speed" wording are out of date; rewrite the section. |
| `sections/05-discussion.tex` | Written for 50 cc. `subsec:hemolysis-trends` is `\provisional` with `\tbd` placeholders and a commented structure to follow. `subsec:viability` has a provisional 35 cc paragraph. |
| `sections/06-conclusion.tex` | Item 5 and the closing sentence are `\tbd`. |
| `appendix/appendixA.tex` | Appendix A: figures A.1–A.10 with labels `fig:app-scaling`, `fig:pp-sweeps`, `fig:app-resolution`, `fig:app-model`, `fig:app-stress-time`, `fig:app-dose`, `fig:F4`, `fig:F5`, `fig:app-sa`, `fig:F7`. Appendix B covers anatomical fit and 35 cc output. |
| `sections/abstract.tex` | Check whether it needs the 35 cc result once the results section is written. |

Macros in use: `\reportfig{file}{width}`, `\tbd{...}` (red placeholder), `\provisional`, `\HIthree`, `\HItwo`, `\taueff`. Citation keys already used include `Marom2014`, `Brockhaus202234`, `Taskin2012`, `Zaman2026`, `NHLBI1985`, `Hellums1994`, `Fukamachi2018383` and `Nascimbene2024661`.

Style of the existing text:
* British spelling (haemolysis), past tense for results, short paragraphs.
* Numbers are given with their CIs, and figures are referenced as `Figure~\ref{fig:F2}(a)`.
* Edit section by section with `overleaf__update_section` or `edit_file`, then show the user.

---

## 2. Decisions made (and why)

* **Flow is inlet-based:** seed speed × diaphragm dome area, from `operating_point.py`. The dome is y = h·cos(πr/L), with L = 56.211 mm × scale and S = 20 mm × scale.
  * The Sim Tracker's 35 cc outlet integrals (areaInt of u·n at the Opening) exceed this inflow by 8–106%. That is impossible in a converged incompressible solution, so they are not used. The user agreed and doesn't know the cause.
  * For 50 cc the tracker agrees with the inflow to within 3%.
* **35 cc geometry:** every length is ×0.888 and the valve angles are unchanged.
* **Cases:** 35 cc at 120 BPM (same rate) and at 171 BPM (matched output, chosen as 120/0.70).
* **Phases:** 16.67, 33.33, 50, 66.67 and 83.33% of systole. The last 50 cc file is named `_8300%` and is treated as 83.33%.
* **Main measure:** mean GW HI3 over all paths, with a 95% bootstrap CI. Outlet-reaching-only results are reported as the second population, as in the 50 cc section. HO and TZ ratios give the sensitivity to the constants.
* **Mechanism:** the exact log decomposition HI3 = A·T^β·τ_eff^α, splitting each change into a transit-time term and a stress term.
* **No jet-speed proxy.** The user rejected it. The u99 (99th-percentile path speed) is fine, since it's already used in `tab:haemodynamics`.
* **Each comparison is per pass.** Haemoglobin release per minute (HI × output) is a secondary point.

---

## 3. Inconsistencies to resolve before or while writing (ask the user where marked)

1. **50 cc flows in the report use the tracker outlet values; 35 cc must use the inlet.** The written 50 cc numbers are:
   * Q = 9.32, 16.36, 20.19, 19.28, 12.14 L/min
   * stroke volume 53.5–55.0 mL
   * cardiac output 6.4–6.6 L/min

   The inlet-based 50 cc values are:
   * Q = 9.36, 16.66, 20.64, 19.88, 12.51 L/min
   * stroke volume 54.9–56.4 mL
   * cardiac output 6.59–6.76 L/min

   These figures appear in `01-haemodynamics`, `tab:size-comparison`, Discussion `subsec:viability` and Conclusion item 1. **Ask the user** whether to switch everything to inlet-based (recommended, for consistency; the change is under 3%) or keep the 50 cc outlet values and footnote the 35 cc method.
2. **Appendix B "Cardiac output"** predicts 2.5–3.4 L/min for 35 cc at 120 BPM from a scaled real-device fill (Nascimbene). The simulation delivers 4.65–4.77 L/min because the inlet boundary condition prescribes a full stroke. Add a sentence reconciling the two: the simulated value is the boundary condition's full-stroke output, not a filling prediction.
3. **The `03-size-comparison` opening paragraph asks about "differences in mesh size"** between 35 cc and 50 cc. That isn't known. **Ask the user** how the 35 cc mesh was made (scaled copy or re-meshed, and the element count).
4. **The 16.67% phase departs from exact scaling** in both 35 cc cases (transit ratio 1.09 and 0.79; stress ratio 0.94 and 1.31). The 35 cc inlet speed there is 3% above exact scaling (the tracker used f = 0.172), but that alone doesn't explain it. Report it honestly as a minor deviation.
5. **Discussion text still says** "a same-rate ratio between 0.56 and 1.0 may not be resolvable" and "~±40% resolution". The measured ratios are about 1.0, and only 33.33% is distinguishable. Update to match.

---

## 4. Table-ready numbers

GW HI3 values are in units of 10⁻⁵ % (same as `tab:hi-phase`), with 95% bootstrap CIs. Ratios are 35 cc divided by 50 cc at the same phase.

| Phase (% systole) | 16.67 | 33.33 | 50.00 | 66.67 | 83.33 |
|---|---|---|---|---|---|
| Q inlet, 50 cc (L/min) | 9.36 | 16.66 | 20.64 | 19.88 | 12.51 |
| Q inlet, 35 cc 120 BPM | 6.76 | 11.67 | 14.51 | 13.94 | 8.86 |
| Q inlet, 35 cc 171 BPM | 9.38 | 16.66 | 20.66 | 19.89 | 12.51 |
| Paths, 35 cc (both cases) | 4,854 | 4,878 | 4,914 | 4,899 | 4,856 |
| Mean GW HI3, 50 cc | 16.7 | 7.40 | 4.03 | 9.38 | 27.0 |
| Mean GW HI3, 35 cc 120 | 16.5 (14.8–18.2) | 5.35 (4.68–6.08) | 3.79 (3.17–4.51) | 9.45 (8.12–11.1) | 27.3 (24.7–30.0) |
| Mean GW HI3, 35 cc 171 | 30.5 (27.1–34.0) | 10.6 (9.02–12.3) | 7.54 (6.27–8.91) | 18.2 (15.7–21.3) | 50.2 (45.2–55.5) |
| Ratio, 120 BPM (CI) | 0.99 (0.87–1.14) | 0.72 (0.58–0.90) | 0.94 (0.73–1.23) | 1.01 (0.81–1.26) | 1.01 (0.87–1.16) |
| Ratio, 171 BPM (CI) | 1.83 (1.59–2.13) | 1.43 (1.12–1.80) | 1.87 (1.43–2.45) | 1.94 (1.56–2.44) | 1.86 (1.60–2.14) |
| Ratio, outlet-reaching only, 120 | 0.96 (0.84–1.10) | 0.84 (0.71–0.99) | 1.01 (0.81–1.26) | 1.00 (0.84–1.19) | 1.05 (0.91–1.21) |
| Ratio, outlet-reaching only, 171 | 1.82 (1.56–2.12) | 1.62 (1.36–1.94) | 2.14 (1.65–2.77) | 1.94 (1.62–2.32) | 1.85 (1.60–2.14) |
| HO / TZ ratio, 120 | 1.00 / 0.98 | 0.80 / 0.82 | 0.97 / 0.97 | 1.01 / 1.02 | 1.03 / 1.02 |
| HO / TZ ratio, 171 | 1.56 / 1.60 | 1.31 / 1.38 | 1.58 / 1.67 | 1.63 / 1.72 | 1.60 / 1.65 |
| GW NIH, 35 cc 120 / 171 (mg/100 L) | 17.3 / 32.0 | 5.6 / 11.1 | 4.0 / 7.9 | 9.9 / 19.1 | 28.7 / 52.7 |
| Median transit, 50 / 120 / 171 (ms) | 326 / 355 / 257 | 187 / 185 / 130 | 144 / 143 / 100 | 127 / 127 / 89 | 161 / 159 / 113 |
| Median τ_eff, 50 / 120 / 171 (Pa) | 0.96 / 0.89 / 1.24 | 0.61 / 0.63 / 0.86 | 0.49 / 0.53 / 0.72 | 0.76 / 0.78 / 1.08 | 1.45 / 1.48 / 2.10 |
| u99 path speed, 50 / 120 / 171 (m/s) | 2.51 / 2.31 / 3.20 | 1.82 / 1.64 / 2.33 | 1.64 / 1.47 / 2.08 | 2.12 / 1.91 / 2.71 | 3.22 / 2.91 / 4.09 |
| Median SA, 50 / 120 / 171 (dyn·s/cm²) | 0.66 / 0.66 / 0.66 | 0.43 / 0.44 / 0.43 | 0.38 / 0.40 / 0.39 | 0.46 / 0.47 / 0.47 | 0.72 / 0.72 / 0.73 |
| Factor below RBC threshold, 50 / 120 / 171 | 5.1 / 5.5 / 3.9 | 4.1 / 6.2 / 4.0 | 3.3 / 4.2 / 3.8 | 4.7 / 4.6 / 3.2 | 4.2 / 3.4 / 2.5 |
| Factor below platelet threshold, 50 / 120 / 171 | 9.9 / 8.4 / 6.1 | 9.2 / 13.2 / 7.8 | 12.3 / 16.3 / 12.6 | 7.5 / 7.6 / 6.3 | 9.6 / 8.5 / 7.1 |
| Peak path stress, 50 / 120 / 171 (Pa) | 146 / 141 / 212 | 87 / 79 / 115 | 107 / 96 / 131 | 115 / 125 / 202 | 195 / 242 / 348 |
| Paths ending early, 35 cc (%) | 7.7 | 3.1 | 1.4–1.6 | 2.9–3.0 | 8.9 |
| Their share of summed HI3, 120 / 171 | 18 / 17% | 36 / 37% | 41 / 37% | 39 / 39% | 19 / 22% |

**Cardiac output** (trapezoid and PCHIP interpolation):
* 50 cc: 6.59–6.76 L/min
* 35 cc 120 BPM: 4.65–4.77 L/min (×0.70), stroke volume 38.7–39.8 mL
* 35 cc 171 BPM: 6.59–6.77 L/min (×1.00)

**Combined flow-weighted GW NIH** (Eq. hi-qs), unweighted in brackets:

| Case | GW | HO | TZ |
|---|---|---|---|
| 50 cc | 11.8 (13.5) | 0.19 | 2.09 |
| 35 cc 120 BPM | 11.4 (13.1) | 0.18 | 2.04 |
| 35 cc 171 BPM | 21.4 (24.6) | 0.29 | 3.37 |

**Haemoglobin released per minute** (HI × output) relative to 50 cc: about 0.68 at 120 BPM and about 1.81 at 171 BPM.

**Mechanism.** Geometric-mean ratios over 33.33–83.33%, with 16.67% in brackets:

| Case | Transit time | τ_eff | Prediction |
|---|---|---|---|
| 120 BPM | 1.00 (1.09) | 1.01–1.03 (0.94) | T 1.00; τ 1.00 viscous or 0.79 inertial |
| 171 BPM | 0.70 (0.79) | 1.42–1.44 (1.31) | T 0.70; τ 1.43 viscous or 1.61 inertial |

* The flow scales as **viscous similarity** predicts.
* At 171 BPM the stress term contributes +0.37–0.39 in log10, against −0.12 for transit time.

**Predicted bands** (methods): 0.56–1.0 at 120 BPM and 1.79–2.38 at 171 BPM for GW; for HO/TZ at 171 BPM, 1.55–1.96 / 1.61–2.04.

**Dose location** is unchanged (`size_dose_maps`). Half of the dose sits in:
* 0.37–0.49% of the visited volume at 16.67% and 83.33%
* 0.05% at 50%

**Solver steadiness** (35 cc, 83.33%, 120 BPM):
* The valve shear strain rate oscillates by ±0.7% (SD 0.28%) about 5,033 1/s, with a period of about 1,330 steps.
* That is about ±1.7% in HI3, far below the sampling resolution.

**Literature contrast:**
* Marom 2014 found the median SA rose 1.56× (50→35 cc at 140→164 BPM). Here the SA median ratio is 1.00–1.06.
* The SA median is insensitive in this model because it is dominated by low-stress bulk transit.
* Brockhaus 2022 found +22% haemolysis for a 20% smaller ReinHeart at a higher rate.

---

## 5. Figure plan for the 35 cc section

Captions and talking points for every figure are in `FIGURE_GUIDE.md`.

**Main text** (`03-size-comparison`):
* **`size_comparison.pdf`** (C1) replaces the placeholder `fig:size-comparison`. It shows the ratio per phase with CIs against the prediction bands.
* **`size_mechanism.pdf`** (C5) is the three-case version of F2 b/c: HI3 against inflow and against valve angle.
  * 120 BPM shifts left at the same height.
  * 171 BPM has the same Q as 50 cc but sits about 1.8× higher.
  * The U shape is kept in all three.
* **`size_lysis_envelopes.pdf`** (C6) is the three-case version of F7. The 171 BPM envelope comes closest to the RBC curve (factor 2.5 at 83.33%), but every case stays below both curves.
  * It could go in the appendix if space is tight.

**Appendix A**, added after the existing entries:
* `size_comparison_completed`
* `size_dose_maps` (C2)
* `size_lysis_sa` (C3)
* `convergence_monitor` (C4)

The per-case F-sets in `report_figs\35cc_*` are for checking only; don't add them unless asked.

After writing, list for the user which new PDFs to upload from `overleaf_upload\figures`. `size_mechanism.pdf` and `size_lysis_envelopes.pdf` are new. Copy any other figure you use there first.

---

## 6. Supervisor request (still to do; needs the user to export from CFD-Post)

The supervisor wants one or two flow/contour screenshots to discuss stagnation and recirculation (washout). Suggestions:
* A velocity contour with in-plane vectors or streamlines on the XZ mid-plane at 50% systole.
* The same view near closure (16.67% or 83.33%), or a low-velocity iso-volume (|u| < 0.05 m/s).
* If possible, the same plane and colour scale for 50 cc and 35 cc 171 BPM.

Where they go: the end of `01-haemodynamics.tex`, where the commented note already sits. The discussion link is that haemolysis dose sits at the valve and outflow tract, while slow regions are a thrombosis (washout) concern. Platelet SA is low everywhere, but stagnation is a separate risk that the Lagrangian model doesn't capture.

---

## 7. Suggested order of work

1. Ask the user about §3 items 1 and 3. Fix item 1 everywhere once decided.
2. Write `03-size-comparison.tex`:
   * opening paragraph: cases, paths and mesh
   * same rate: operating point, ratios, mechanism
   * matched output: the same three
   * the table, from the §4 numbers
   * the figures C1, C5 and optionally C6
   * the short combined-estimate subsection
3. Discussion `subsec:hemolysis-trends`:
   * H1 is supported (four of five phases about 1; 33.33% at 0.72).
   * H2 is supported (about 1.8–1.9×; inside the GW band; HO/TZ slightly below their bands).
   * The mechanism is viscous similarity.
   * The valve still dominates.
   * Contrast with Marom's SA and Brockhaus.
   * The trade-off: same rate is fine per pass but gives 4.7 L/min; matched output gives about 2× damage per pass and about 1.8× haemoglobin per minute, with the RBC margin down to 2.5.
4. Discussion `subsec:viability`: replace the provisional paragraph with the simulated 4.65–4.77 L/min and assess its adequacy against the 2–9 L/min requirement.
5. Update the uncertainty paragraph (§3 item 5) and Appendix B (§3 item 2).
6. Conclusion item 5 and the closing sentence. Then check the abstract.
7. Add the appendix entries and give the user the upload list.
8. Later: the CFD screenshots (§6), once the user exports them.
