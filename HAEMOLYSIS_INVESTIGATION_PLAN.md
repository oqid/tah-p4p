# SynCardia haemolysis investigation plan

Prepared 6 October 2026 from the project scripts and collected literature.

## Aim and scope

The project already has most of the haemolysis calculation machinery. The main missing element is a structured comparison pipeline that explains **why haemolysis changes** and establishes confidence in those changes.

Organise the investigation around three questions:

1. How does predicted haemolysis vary through systole in the 50 cc model?
2. What changes when the device is reduced to 35 cc, and which changes explain the haemolysis difference?
3. How does lowering BPM alter haemolysis, transit time and delivered flow in the 35 cc model?

The strongest defensible outcome is quantitative model predictions and controlled comparisons, with absolute haemolysis validation wherever sufficiently matched experimental data exist. Taskin found substantial errors in absolute haemolysis predictions, so using its preferred formulation does not itself validate the SynCardia predictions. See [Taskin et al.](https://researchportal.bath.ac.uk/en/publications/evaluation-of-eulerian-and-lagrangian-models-for-hemolysis-estima/).

## 1. Separate geometry from operation

A 50-to-35 cc comparison can change chamber geometry, valve dimensions, flow rate, velocity, transit time and BPM simultaneously. Without separating those factors, it can demonstrate a difference without explaining it.

| Comparison | What to hold consistent | What it establishes |
|---|---|---|
| 50 cc across systole | Geometry, BPM, numerical settings | Phase dependence |
| 50 vs 35 cc at intended operation | Corresponding phase and clearly defined operating conditions | Combined size-and-operation effect |
| 50 vs 35 cc at matched instantaneous flow | Flow, corresponding chamber/valve configuration, fluid properties | Geometry effect at that flow |
| 35 cc across BPM | Geometry, waveform construction and systolic fraction | Rate-related operating effect |
| 35 cc baseline vs modified design | Operating conditions and numerical settings | Effect of the specific design change |

Matched-flow cases need not represent intended clinical operation; they are useful controlled investigations.

Lowering BPM does not preserve cardiac output automatically. Ideally:

```text
Q_mean [L/min] = SV [mL] * BPM / 1000
```

For full effective stroke volumes, 50 cc at 120 BPM gives 6 L/min, whereas 35 cc gives 4.2 L/min. Matching 6 L/min with 35 cc would require approximately 171 BPM. Lowering BPM further reduces output unless another parameter compensates.

This does not invalidate the investigation. It means conclusions should describe the **haemolysis–output trade-off**, rather than treating reduced haemolysis alone as an improvement.

## 2. Decompose stress and transit-time effects

Mean velocity and mean stress alone will not adequately explain the haemolysis results. HI3 weights high stresses nonlinearly.

For each trajectory, define:

```text
T = sum_i(dt_i)

tau_eff = [sum_i(tau_i^(alpha/beta) * dt_i) / T]^(beta/alpha)

HI3 = A * T^beta * tau_eff^alpha
```

The final identity holds exactly for the implemented HI3 calculation.

This explains competing effects:

- Longer transit increases predicted haemolysis at fixed effective stress.
- Higher effective stress increases predicted haemolysis at fixed transit time.
- Faster flow may shorten transit while increasing stress; the net effect depends on both.

For GW, the stress exponent is 2.416 and the time exponent is 0.785. This makes high-stress exposure particularly influential within this model.

**Code addition:** export total transit time, time-weighted mean stress and effective stress for every trajectory and coefficient set.

**Graph:** transit time versus effective stress, coloured by trajectory HI3, comparing phases or device sizes. This separates long-exposure paths from high-stress paths.

Apply this identity per trajectory. Device-averaged haemolysis cannot generally be reconstructed from device-averaged transit and stress.

## 3. Main graphs to add

| Graph | Question it answers | Required addition |
|---|---|---|
| HI3/NIH versus normalised systolic phase, with flow alongside | When is predicted haemolysis highest? | Numeric phase metadata and case plotting |
| 35/50 cc haemolysis ratio versus phase | Does downsizing affect every phase equally? | Matched-case comparison |
| Transit and effective stress versus phase | What explains the phase or size difference? | Effective-stress calculation |
| Flux-weighted HI3 and transit-time cumulative distributions | Does the average hide highly exposed populations? | Distribution plotting; statistics mostly exist |
| Cumulative share of predicted haemolysis versus cumulative flow fraction, sorting paths by HI3 | Does a small flow fraction dominate damage? | Contribution analysis |
| Spatial maps of linearised HI3 source | Where should design changes be investigated? | More complete segment-source export |
| NIH, transit time and mean output versus BPM | What improves or worsens when rate is lowered? | Operating-condition sweep |
| Haemolysis versus output, labelled by BPM/design | Which cases offer a useful trade-off? | Combined case table |

Stress-coloured streamlines show where stress is high. Source-coloured streamlines additionally account for exposure duration and nonlinear stress weighting.

The diagnostics currently exports only the three largest source segments per selected path. That identifies candidate hotspots but is insufficient for complete spatial or regional accounting.

## 4. Investigate valves and leakage explicitly

Marom identifies high stresses around valves and regurgitant flow, including differences between leakage through a closed valve and regurgitation during valve closure. These are relevant mechanisms to investigate. See [Marom et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4245394/).

Distinguish at least:

- Chamber bulk.
- Open outlet valve and downstream jet.
- Closed inlet-valve leakage regions.
- Connectors or abrupt area changes.

**Code addition:** assign segment midpoints to geometry-specific regions and sum their linearised HI3 source.

Report these as fractions of linearised source, rather than additive fractions of final haemolysis, because the final exponent is nonlinear.

Check trajectory coverage: inlet-seeded paths may poorly represent chamber blood entering a leakage route or remaining trapped in the chamber. Do not infer negligible damage merely because few exported paths visit that region. Additional releases and separate exit populations may be needed.

## 5. Separate verification and validation

### Calculation verification

The diagnostics already contains analytical checks and integration sensitivities. Retain these.

Correct the opening HI3 equation in `haem5.py`: at the time of this review, its docstring contradicts the implemented calculation, which sums the linearised source first and exponentiates afterwards.

### Numerical verification

Demonstrate convergence using new CFD solutions and trajectory exports. The particle-count sweep is useful but does not replace:

- Mesh refinement, particularly near valve gaps and jets.
- Trajectory integration/export-spacing refinement.
- Release-location and release-count sensitivity.
- Velocity-floor sensitivity in slower-flow cases.

Monitor HI3, transit-time tails and hotspot locations alongside pressure and velocity. Numerical variation should be smaller than the size or BPM effect being claimed.

### Physical validation

Compare matched flow, pressure and velocity measurements where available, then haemolysis where the experimental setup is sufficiently comparable. Agreement with another numerical paper is a consistency check, not experimental validation.

Audit two existing summary items:

- Establish the source and operating conditions behind `SYNCARDIA_50CC_NIH = 37.15`.
- Remove or justify the wording “ASTM acceptable limit” for the hard-coded 100 mg/100 L. An ASTM testing method does not automatically establish a universal clinical acceptance threshold. FDA also describes the lack of broadly established bench safety criteria: [FDA discussion](https://www.fda.gov/science-research/advancing-regulatory-science/biokinetic-modeling-hemolysis-and-renal-injury-cardiopulmonary-bypass-patients-establishing-safety).

The collected Realheart paper reports SynCardia haemolysis from a complete systemic-and-pulmonary circulation setup. Its NIH is contextual evidence, not a directly matched validation target for one chamber at one systolic snapshot. See [Zaman et al.](https://www.nature.com/articles/s41598-026-68178-2).

Audit the HI-to-NIH conversion units and blood-property definitions before making that comparison.

## 6. Make the quasi-steady assumption measurable

The example input analysed on 6 October 2026 had an unweighted mean transit time of approximately 66 ms. Compare transit with how quickly the actual flow and valve configuration change.

Add:

- Actual systolic duration.
- Snapshot time and normalised systolic phase.
- Transit-time mean, P95 and P99 relative to systolic duration.
- A comparison with the timescale over which prescribed flow/configuration changes appreciably.

A half-cycle is only the systolic duration when systole occupies 50% of the cycle.

Early and late systole may deserve special attention because valve motion or rapidly changing flow can undermine frozen-field trajectories. Sampling more steady states improves phase resolution; it does not turn them into transient particle histories.

If calculating a combined snapshot result, use ejected-volume weighting:

```text
HI_QS = sum_k(HI_k * Q_k * delta_t_k) / sum_k(Q_k * delta_t_k)
```

Label it an **ejection-weighted quasi-steady estimate**, not validated whole-cycle haemolysis.

## 7. Prioritised code work

Build on `haem_diagnostics.py` rather than adding all investigation logic to `haem5.py`.

1. **Expand the case manifest and comparison CSV.** Carry through BPM, phase fraction, systolic fraction, effective stroke volume, instantaneous flow, pressure head, valve configuration and geometry description. At review time, BPM is used in diagnostics but omitted from its top-level comparison CSV.
2. **Add explanatory trajectory metrics.** Transit in ms, time-weighted stress, effective stress and each path’s contribution to weighted device HI.
3. **Add complete segment-source and region analysis.** Keep the existing small hotspot export as a lightweight option.
4. **Create a case-comparison plotting script.** Generate phase, size, BPM, distribution and trade-off graphs from one combined table.
5. **Add solver-summary inputs.** Flow, pressure, mass balance and volume-based low-velocity measures require solver exports; scalar streamline data alone cannot establish them.

Changing BPM metadata alone will not change calculated haemolysis. Each BPM case needs a corresponding CFD field generated from its revised operating waveform.

## 8. Manageable sequence of new simulations

Start with:

1. Five representative systolic phases for the 50 cc baseline.
2. The corresponding five phases for the 35 cc baseline.
3. Two additional BPM settings for the 35 cc model at three informative phases—initially early, peak and late systole.
4. A few matched-flow size comparisons to isolate geometry effects.
5. One targeted design change, chosen from source maps, at the most informative operating conditions.

This gives **16 primary cases before matched-flow, convergence and design studies**. Refine phase sampling where results change rapidly.

Keep platelet SA and empirical lysis screening as supporting analyses. The report priorities are reliable HI3 comparisons, stress–time explanations, valve-region attribution and the haemolysis–output trade-off. These additions explain mechanisms behind the results rather than simply reporting a series of NIH values.
