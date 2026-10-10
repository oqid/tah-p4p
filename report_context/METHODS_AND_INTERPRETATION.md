# SynCardia downsizing investigation: methods, evidence and interpretation

Prepared 8 October 2026. This is a working synthesis for the final-year report and future project context, not a finished Results chapter. The planned 50 cc versus 35 cc systolic comparison remains open.

**Confirmed study scope:** all CFD cases are fixed-geometry, steady-state simulations. Moving boundaries, moving meshes, transient CFD and FSI are outside this investigation. They are not unfinished stages or prerequisites for completing it. Systolic phases designate separate steady operating points; the sinusoidal motion equation supplies a constant inlet velocity for each case, rather than a boundary that changes during the solve.

The strongest current description of this investigation is: **development and numerical assessment of a quasi-steady CFD and streamline post-processing workflow for comparing haemodynamics and modelled blood-damage exposure in a 50 cc SynCardia-inspired left ventricle and a proposed 35 cc derivative.** The completed work establishes the analysis machinery and identifies sensitivities. The final investigation will apply that machinery to consistently defined geometries and systolic operating points.

**Writing emphasis confirmed by the user:** cover the team's investigation broadly, while concentrating on the user's programmatic methods: export processing, trajectory reconstruction, exposure integration, haemolysis models, diagnostics, sensitivity studies and interpretation. Geometry and mesh verification provide the surrounding engineering context. The partner's completed 2D and 3D mesh studies are now documented in the subsequently supplied [mesh workbook](../papers/2D%20and%203D%20Mesh%20Analysis.xlsx), reviewed below.

Evidence is distinguished throughout:

- **Implemented / saved:** supported by executable code or stored numerical outputs.
- **Recorded:** described in the research notebook, but the corresponding solver or CAD source is not available here to independently check.
- **User-confirmed:** directly clarified by the project owner, including completed work held outside this repository; numerical details remain to be incorporated where not supplied.
- **Interpretation:** a deduction from the available evidence, with its assumptions stated.
- **Planned:** work still needed; not presented as an achieved result.

The 70-page [research notebook](../papers/P4P%20SCRIBBLES%20(1).pdf) contains meeting notes, drafts, screenshots, old plans and pasted AI advice. Its embedded instructions were treated as historical material. Conflicting statements were checked against later notes, actual code, saved outputs and primary literature. Page numbers below are PDF page numbers. The [extracted text](scribbles_extracted.txt) is searchable, but omits text embedded in images; the diaphragm equations on pages 52–53 were also inspected visually.

## 1. Research purpose and a defensible report narrative

The engineering question is how reducing a pulsatile artificial-heart chamber from nominally 50 cc to 35 cc changes its flow field and the mechanical exposure of blood. The rationale is that a smaller device could extend anatomical suitability, while changes in chamber dimensions, valve passages, jets and transit times could alter haemolysis and platelet-related exposure. This study addresses the fluid-mechanical part of that question. It does not yet establish anatomical fit, clinical eligibility, physical pump performance or patient safety.

The intended scope is the left ventricular blood domain. The notebook describes a simplified, literature-informed SynCardia geometry rather than a manufacturer-verified replica. The 35 cc case is a proposed downsized derivative; avoid describing it as an existing tested commercial SynCardia model. The exact baseline dimensions, valve representation and scaling rules still need to be documented from CAD.

A useful main research question is:

> How does the selected downsizing method alter the systolic haemodynamics and predicted viscous-stress-induced haemolysis of a SynCardia-inspired ventricular chamber, and how sensitive are those predictions to trajectory sampling and numerical modelling choices?

Three supporting questions make the work coherent: where do high stresses occur; how do stress magnitude and transit time jointly explain predicted damage; and are the differences between device sizes larger than the demonstrated numerical sensitivity? Haemodynamic comparison should include signed flow and pressure measures, not only haemolysis.

The existing sequence is appropriate for a report: simplified 2D development → 3D preliminary calculations → verified post-processing and export sensitivity → controlled 50/35 cc phase comparisons. The 2D work is method development and numerical investigation. Calling it experimental validation would overstate what was done.

## 2. What has been done, and what remains open

| Component | Evidence and current status |
|---|---|
| Initial post-processing proof of concept | Notebook pp. 4–8 records an aerofoil export used to exercise the parser, exposure calculation and plots. Its haemolysis values are not heart results. |
| Simplified 2D ventricular simulations | Notebook pp. 18–33 records boundary-condition repairs, blood-property changes, geometry smoothing, mesh refinement and haemolysis calculations. A retained `export_18pa.csv` and diagnostics support the post-processing results. |
| Preliminary 3D analysis | Notebook pp. 34–35 and saved exports show 3D streamline calculations. These are development cases; complete solver metadata is missing. |
| Haemolysis and stress accumulation | Implemented in `haem5.py`, with GW/HO/TZ coefficient choices and both HI2/HI3 accumulators. |
| Numerical diagnostics | Implemented in `haem_diagnostics.py`, including analytic checks, alternative quadratures, weighting, endpoint classification and distributions. |
| Release-count study | Saved results extend to 300,000 requested paths, with a denser 1,000–10,000 study. |
| Segment and tolerance studies | Actual exports and numerical results are present, including segment settings 25–5,000 and grid-relative tolerances 0.02–0.0025. These studies are no longer merely proposed. |
| Empirical stress-duration screening | Implemented, with digitised reference-curve provenance and saved screening outputs. |
| Diaphragm and valve CAD | User confirms several CAD models exist. A sinusoidal diaphragm profile is revolved and Boolean-subtracted from the fluid domain. A separate SynHall-inspired valve is simplified to a disk with translation/rotation used to represent its configuration. CAD assemblies and detailed dimensions were not reviewed here; configurations are fixed during each steady solve. |
| Partner's mesh studies | Supplied workbook confirms 2D and 3D studies recording velocity, pressure, wall shear and y+. Richardson extrapolation uses maximum velocity; 3D gives p = 2.0185 and an extrapolated maximum speed of 2.56353 m/s. Haemolysis convergence is not tabulated. See detailed assessment below. |
| Reproducible case runner | `run_haem.py` groups results, logs, configurations and input hashes in separate run folders. |
| Controlled 50 cc versus 35 cc phase matrix | Planned. No adequately labelled completed comparison was identified in this repository. |
| Moving-boundary / transient CFD / FSI | Explicitly outside the user-confirmed study scope. Earlier notebook suggestions do not describe the current investigation or its remaining work. |
| Experimental validation | Not established by the supplied material. Early notebook proposals for physical loops must not be written as completed methods. |

No CFX solver project/results or Inventor assemblies were found in the file inventory used for this review. The user confirms multiple CAD models exist separately. The mesh workbook now supplies mesh parameters and tabulated solver outputs, although the underlying solver files were not inspected. The STEP file provides a plotting overlay; it does not establish all simulation settings. Accordingly, recorded workbook/notebook settings are distinguished from independently recovered solver configuration.

## 3. Geometry development and downsizing

The user confirms that the current CAD work includes several chamber configurations. The diaphragm shape is constructed by revolving a sinusoidal profile, then using a Boolean subtraction to remove that solid from the fluid domain. This constructs the stationary fluid-domain geometry for a chosen diaphragm configuration. The spatial sinusoidal profile used to make the CAD surface and the temporal sinusoidal displacement prescription used to set inlet speed are distinct parts of the method; record both equations and their parameters rather than assuming they are the same function.

A separate valve model is based on the SynHall valve and simplified to a translating/rotating disk. Here, translation and rotation describe how the disk configuration is represented or positioned in CAD. The disk remains fixed in each steady CFD calculation. Its diameter, thickness, position, angle, clearances and any omitted hinge detail should be documented for the final geometries. Do not carry the early 2D rectangular-valve simplification forward as a description of this later disk model.

The 2D geometry used simplified valve blockages and an idealised diaphragm/chamber. Early notes identify an outlet accidentally assigned as a wall, difficulty resolving narrow valve regions and poor mesh quality. The outlet was changed to an opening; later revisions removed a valve gap and added common 0.2 mm fillets, with 0.4 mm diaphragm-corner fillets (pp. 19–29). These changes enabled a roughly 110,000-element calculation.

Interpret these edits carefully. Filleting can improve mesh quality, but also changes local acceleration, separation and stress. Removing a gap removes a potential leakage jet. Differences between such revisions cannot be attributed purely to numerical refinement. Report these as geometry simplifications and retain a diagram showing the final features. A closed valve modelled as an impermeable wall excludes closed-valve leakage by construction; it does not demonstrate that leakage is negligible in the real device.

For geometrically similar downsizing, volume and length obey

\[
s_L=(35/50)^{1/3}=0.887904,\qquad s_A=s_L^2\approx0.788374,\qquad s_V=0.7.
\]

Thus, a uniform 11.21% reduction in lengths gives a 30% reduction in volume. Scaling every length by 0.7 would instead give 17.15 cc from 50 cc. This is a geometric calculation, not evidence that the current 35 cc CAD actually uses uniform scaling.

Specify which features scale: chamber height/diameter, diaphragm stroke, valve diameter/thickness, hinge gaps, connector diameters, and fillet radii. A design retaining valve dimensions while shrinking the chamber is a distinct downsizing strategy. Record end-diastolic and end-systolic fluid volumes, because nominal chamber capacity is not the effective ejected stroke volume.

**Interpretive hypothesis:** smaller passages at matched flow may increase velocity gradients and stress, while transit time may decrease. At the same BPM with a smaller delivered stroke volume, flow also changes. Therefore a size effect is conditional on the operating constraint; “smaller necessarily causes greater haemolysis” is not a result yet.

## 4. CFD assumptions and systolic operating conditions

The notebook records ANSYS CFX, blood density changed to 1,060 kg/m³, dynamic viscosity changed to 0.0035 Pa·s, and an incompressible Newtonian approximation (pp. 20–21). The analysis metadata describes an SST turbulence model, but this is a configurable/default label in Python, not proof of the solver configuration. Confirm the CFX model and version from the solver setup. SST should not be described loosely as “SST k-epsilon or something.”

The implied continuum calculation solves conservation of mass and momentum for the prescribed geometry and boundaries. Newtonian modelling removes shear-dependent viscosity; incompressibility removes density changes driven by pressure. The post-processor follows passive trajectories in the solved velocity field and does not resolve individual blood cells, cell deformation, cell–wall collisions or feedback from blood damage to the flow.

The early representative systolic case used a diaphragm-boundary speed near 0.126 m/s and a recorded outlet pressure around 100 mmHg. These are historical development settings. Confirm pressure reference, opening direction/backflow settings, wall treatment, inlet turbulence, residual targets and actual converged balances for each final case. A prescribed pressure level is not the same quantity as the pressure difference required to drive flow.

### Quasi-steady diaphragm prescription

The latest notebook plan is to sample **20%, 40%, 60% and 80% of systole** (pp. 53–54), superseding an earlier 0/25/50/75/100% proposal (p. 49). The user confirmed these phases on 8 October 2026, together with the intention to investigate both same-BPM and matched-flow comparisons. These are planned cases, not saved completed results.

For BPM \(B\), cycle period \(T_c=60/B\), systolic fraction \(f_s\), and \(T_s=f_sT_c\), the recorded half-cosine displacement is

\[
x(t)=\frac S2[1-\cos(\pi t/T_s)],\qquad
v(t)=\frac{\pi S}{2T_s}\sin(\pi t/T_s),\qquad 0\leq t\leq T_s.
\]

With the notebook assumptions \(S=20\) mm, \(B=120\) and \(f_s=0.5\), systole is 0.25 s and peak displacement speed is 0.12566 m/s.

| Systolic fraction | Diaphragm displacement | Displacement speed |
|---:|---:|---:|
| 0.20 | 1.910 mm | 0.07386 m/s |
| 0.40 | 6.910 mm | 0.11951 m/s |
| 0.60 | 13.090 mm | 0.11951 m/s |
| 0.80 | 18.090 mm | 0.07386 m/s |

This prescription creates an interesting comparison: 20/80% and 40/60% have equal displacement speeds but different chamber configurations. Differences between those pairs could help distinguish geometry effects from speed effects, provided their imposed volume fluxes are also checked. Equal displacement speed does not ensure equal flow for a changing dome shape.

For fluid-domain volume \(V(x)\), ideal ejection with the other port sealed obeys

\[
Q(t)=-\frac{dV}{dt}=-\frac{dV}{dx}\frac{dx}{dt}.
\]

A translating piston of constant area gives \(Q=Av\); a deforming dome generally does not. Extract CAD volumes over the prescribed displacement, derive the swept-volume rate, and make the flow-driving boundary consistent with it. This is especially important after downsizing.

**User-confirmed implementation, 8 October 2026:** the diaphragm surface is set as an **inlet**, with prescribed velocity equal to the sinusoidal diaphragm displacement speed at the selected cycle position, using 20 mm total travel. The user confirms the resulting peak of approximately 12.6 cm/s. Thus the implemented prescription is the phase-specific speed above; deriving and matching the CAD swept-volume flux is a recommended consistency check, not a step already demonstrated.

Each simulation holds its geometry and prescribed inlet velocity fixed. Any separate CAD configurations used to represent different phases are independent stationary domains. The displacement equation defines the intended operating point; no diaphragm boundary moves during a CFD solve. Here, “quasi-steady” means a collection of these steady-state cases at representative systolic phases.

In a fixed-domain steady calculation, this inlet is an approximation to displacement-driven flow. It injects surrogate fluid into a frozen chamber; it is not a dynamically moving impermeable diaphragm. The exact inlet velocity direction/profile and phase-specific CAD-volume mapping still need recording. The imposed inlet speed can remain the stated method while its ability to reproduce the intended swept-volume flow is assessed explicitly.

Steady streamlines can represent trajectories in a frozen field. Independent phase snapshots do not follow the same blood through changing chamber and valve configurations. They cannot resolve valve closure impulses, transient leakage, filling history, repeated passages or FSI. Sampling more phases improves phase coverage without removing that limitation.

### Flow rate and cardiac output

Prefer the signed normal flux from CFD:

\[
Q_{out}=\int_{A_{out}}\mathbf{u}\cdot\mathbf{n}\,dA.
\]

The notebook correctly replaces area times mean speed with the normal-component integral (p. 47). A scalar speed average overestimates net flux when motion is oblique or partly reversed. The current Python area-times-seed-speed option is also conditional on normal flow and representative equal-area inlet seeding.

The recorded preliminary value \(0.000472695\,\mathrm{m^3/s}\) converts to **28.36 L/min instantaneously**. It was first attached to the speed-based calculation; the notebook says the normal-flux correction also gave about 28 L/min but does not retain a separate exact number. This should trigger a stroke/area/boundary-condition audit, not an immediate statement that cardiac output is 28 L/min.

For illustration only, a half-sine ejection flow with a constant-area piston and 50% systolic fraction has \(\bar Q_{cycle}=Q_{peak}/\pi\). A peak of 28.36 L/min would then imply approximately 9.03 L/min and 75.2 mL/beat at 120 BPM. A changing dome or a non-peak sample invalidates that shortcut. The useful conclusion is that the preliminary flow must be reconciled with actual swept volume before its results become the final baseline.

Use \(CO=SV_{net}B/1000\) in L/min for net stroke volume in mL. Full effective strokes of 50 and 35 mL at 120 BPM give 6 and 4.2 L/min, respectively; those are arithmetic examples, not demonstrated outputs. Matching 6 L/min with 35 mL would require about 171.4 BPM. Changing BPM metadata in Python does not change a CFD field.

## 5. Streamline export and exposure reconstruction

The export workflow is recorded on pp. 4–5 and implemented in [haem5.py](../haem5.py): define scalar `varshear`, seed forward streamlines, and export coordinates, stress, speed, node numbers and connectivity using CFD-Post Generic CSV. Use plain Line exports for final diagnostics. The core parser has a Faces reconstruction fallback, but diagnostics rejects ribbon/tube Faces data because physical direction may be lost.

Nodes are reconstructed into ordered paths from connectivity. For segment \(i\),

\[
\Delta s_i=\|\mathbf{x}_i-\mathbf{x}_{i-1}\|,\quad
\bar u_i=(u_i+u_{i-1})/2,\quad
\Delta t_i=\frac{\Delta s_i}{\max(\bar u_i,u_{floor})}.
\]

The code uses \(u_{floor}=\max(0.01\bar u_{positive,path},0.001\,\mathrm{m/s})\), with a fallback positive-speed mean of 0.1 m/s when no positive speeds exist. Cumulative transit is \(t_i=\sum_{k\leq i}\Delta t_k\). The first node has zero duration. Stress is sampled at the downstream endpoint of each segment.

The floor raises low speeds; it does not set particles to zero velocity. Consequently it can shorten estimated residence in very slow regions and reduce accumulated exposure. It was inactive in the audited preliminary examples, which is reassuring for those exports only. It does not establish adequate capture of stagnant or trapped blood.

The trajectories use variable exposure increments inferred from spatial distance and speed, not uniformly spaced physical timesteps. Their total time is an estimated passage time through the represented path, not automatically chamber washout time or blood age over multiple beats.

### Scalar stress

The exported expression is molecular dynamic viscosity multiplied by scalar shear strain rate. Under the standard incompressible Newtonian definitions, writing \(S=(\nabla u+\nabla u^T)/2\) gives \(\tau_{scalar}=\mu\sqrt{2S:S}\). With viscous tensor \(\sigma=2\mu S\) and zero trace, this equals

\[
\left\{\frac{(\sigma_{xx}-\sigma_{yy})^2+(\sigma_{yy}-\sigma_{zz})^2+(\sigma_{zz}-\sigma_{xx})^2}{6}
+\sigma_{xy}^2+\sigma_{yz}^2+\sigma_{zx}^2\right\}^{1/2}.
\]

This algebra explains the notebook's comparison to Taskin's scalar convention. Confirm that the actual export uses these definitions and molecular viscosity. The [CFX theory definitions](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/cfx_thry/i1298716.html) identify the tensor and scalar strain-rate quantities. This is a stress within the fluid, distinct from wall shear stress. SST affects the resolved mean flow, but the damage scalar excludes explicit Reynolds-stress contributions; the method therefore estimates damage from the resolved viscous exposure, not all turbulent blood trauma.

## 6. Haemolysis equations and model choices

Using repository notation, uniform-stress haemolysis is \(HI[\%]=A\tau^\alpha t^\beta\), with stress in Pa and time in seconds. Here alpha is the stress exponent and beta is the time exponent; Taskin's paper uses the symbols in the opposite roles. Match exponents by their physical variable when writing the report.

| Set | A, for percent output | Stress exponent alpha | Time exponent beta | Calibration range recorded in code |
|---|---:|---:|---:|---|
| GW | 3.62 × 10⁻⁵ | 2.416 | 0.785 | stress <255 Pa; duration <0.7 s |
| HO | 1.80 × 10⁻⁶ | 1.991 | 0.765 | stress <700 Pa; duration <0.7 s |
| TZ | 1.228 × 10⁻⁵ | 1.9918 | 0.6606 | 50–320 Pa; duration <1.5 s |

These are alternative empirical parameterisations, not repeated measurements or a statistical confidence interval. The repository defaults to GW. Retain HO as a material sensitivity comparison. TZ is poorly matched to the predominantly low-stress preliminary exposure: the 100-path diagnostic reports 99.990% of its weighted residence below 50 Pa. It can remain in an appendix rather than determining the main conclusion.

### HI2: cumulative-time derivative form

\[
HI2=\sum_i\beta A\,t_i^{\beta-1}\tau_i^\alpha\Delta t_i.
\]

Because beta is below one, identical stress-duration pairs receive smaller contributions later in the trajectory. This makes HI2 sensitive to where high stress occurs relative to the chosen starting time. The diagnostic exact alternative,

\[
HI2_{exact}=\sum_i A\tau_i^\alpha(t_i^\beta-t_{i-1}^\beta),
\]

integrates the cumulative-time factor analytically for piecewise-constant stress. It checks quadrature; it does not remove the model's ordering dependence.

### HI3: primary linearised accumulation

The executed implementation is

\[
HI3=\left[\sum_i A^{1/\beta}\tau_i^{\alpha/\beta}\Delta t_i\right]^\beta
=A\left[\sum_i\tau_i^{\alpha/\beta}\Delta t_i\right]^\beta.
\]

The sum is taken **before** the final exponent. The opening docstring of `haem5.py` still prints the contradictory expression that exponentiates each segment before summing. Use the executed equation above. The code itself already implements the correct order.

HI3 exactly recovers the constant-stress power law and is invariant to subdividing unchanged constant-stress exposure. This does not make CFD trajectory results independent of spatial resolution: refining a path can alter stress peaks, duration and geometry. HI3 also loses information about the ordering of stress-duration pairs; it is not a biological memory model.

A helpful explanatory identity, derived from this implementation, is

\[
T=\sum_i\Delta t_i,\qquad
\tau_{eff}=\left[\frac{\sum_i\tau_i^{\alpha/\beta}\Delta t_i}{T}\right]^{\beta/\alpha},
\qquad HI3=A\,T^\beta\tau_{eff}^{\alpha}.
\]

This separates transit and stress effects for each path. For GW, the linearised integrand weights stress to approximately the power 3.08. A small number of short high-stress intervals can therefore dominate damage even when mean stress is low. The identity cannot generally be applied to averaged transit and stress to reproduce the device average. Effective stress is a proposed additional explanatory export, not a field already present in the saved diagnostic CSV.

Taskin supports examining accumulator and coefficient sensitivity, but does not validate this TAH. Its Table 3 gives particularly different errors across coefficient sets: HI3's aggregate percentage error is 57% with HO and 5,805% with GW for its tested devices. These errors are not uncertainty bounds for your cases. Selecting GW simply because its NIH is closer to a desired literature number would be circular. See [Taskin paper](../papers/taskin-etal.pdf) and its [publication record](https://researchportal.bath.ac.uk/en/publications/evaluation-of-eulerian-and-lagrangian-models-for-hemolysis-estima/).

### Population averaging and NIH

The core pipeline averages per-path HI using normalised starting-speed weights. This approximates flux weighting only when seeds represent equal inlet areas and flow is normal and forward. Equally spaced points on a line do not automatically represent equal areas over a 3D inlet. Diagnostics can instead read supplied normal-flux/area or mass-flux weights and select outlet-completed paths.

Without an outlet plane, all reconstructed paths are included and labelled **completion unverified**. An endpoint-plane test additionally requires an upstream start and outward approach; it is a conservative geometric test, not knowledge of the solver's termination reason or outlet aperture. Missing releases have unknown flux. Their number cannot be interpreted as recirculation, nor can the surviving paths be assumed to cover all inlet flow.

For HI expressed numerically as percent, the repository conversion is

\[
NIH_{mg/100L}=1000\,HI[\%]\,(1-Hct)\,Hb,
\]

with Hct = 0.30 and Hb = 150 g/L: a multiplier of 105,000. For example, 0.00018326456% becomes 19.24278 mg/100 L. The factor of 100 for percent and the 100 L reporting volume cancel. This conversion is conditional on the adopted HI definition \(HI=100\Delta Hb_{plasma}/Hb\); if HI is instead defined as a mass fraction of total whole-blood haemoglobin released, the plasma-volume factor must be reconsidered. State the definition and confirm the Hb basis, especially when haemodilution is involved. The arithmetic reproduces the repository; it is not experimental NIH measurement.

## 7. Platelet-related exposure and empirical screening

The implemented stress accumulation is \(SA=10\sum_i\tau_i\Delta t_i\), reported in dyn·s/cm² because 1 Pa = 10 dyn/cm². It measures accumulated mechanical exposure. The core reports an unweighted path median and the count fraction above 35 dyn·s/cm²; diagnostics additionally provides explicitly weighted distributions. Do not mix those populations when comparing figures.

Marom's full-cycle TAH study reported median SA values of 1.12, 2.19 and 3.41 dyn·s/cm² for 70, 50 and 35 cc models, respectively. It provides a relevant qualitative size-comparison context, but its simulation scope differs from a single frozen systolic field. Your lower preliminary median is not proof that the simplified model is less thrombogenic. See [Marom et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4245394/) and [local paper](../papers/marom-etal.pdf).

Describe 35 dyn·s/cm² as the selected SA reference criterion. Exceedance is neither an activated-platelet fraction nor a clinical thrombosis probability. Stress waveform, repetition and sensitisation can matter beyond a single integrated SA; these are not modelled here. [Sheriff et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4854792/) provides experimental motivation for that limitation, without validating a numerical threshold for this device.

The separate empirical screen uses manually digitised RBC and platelet **lysis** curves from *Guidelines for Blood-Material Interactions*, Figure 3.6, printed p. 78 / PDF p. 87. For each stress level S, it finds the longest contiguous duration with stress at least S, then compares S with the interpolated critical stress at that duration. The maximum ratio is the utilisation. The code samples 220 stress levels per path and uses log-log reference interpolation; exposures outside the traced duration domain are unassessed.

This is an independent screening construction, not an extra term in HI3. It ignores separated repeated exposures and does not reproduce a validated variable-stress lysis model. Being below the curves does not establish zero haemolysis, absence of platelet activation or clinical safety. The reciprocal utilisation is named `minimum_safety_factor` in the saved CSV but should be called a reference-curve margin in the report. See [digitisation provenance](../papers/figure_3_6_digitization.md).

## 8. Existing results and what they mean

### Preliminary 2D and 3D calculations

| Dataset | Reconstructed paths | GW HI2 (%) | GW HI3 (%) | Code-converted HI3 NIH, mg/100 L | Unweighted median SA, dyn·s/cm² |
|---|---:|---:|---:|---:|---:|
| 2D `export_18pa` | 969 | 4.2083 × 10⁻⁶ | 7.3226 × 10⁻⁶ | 0.7689 | 0.3372 |
| 3D 100-path development export | 100 | 8.4979 × 10⁻⁵ | 1.8326 × 10⁻⁴ | 19.2428 | 0.6747 |
| 3D 1,000-request release export | 985 | 1.2867 × 10⁻⁴ | 2.2810 × 10⁻⁴ | 23.9506 | 0.6433 |

Sources: [2D diagnostics](../outputs/diagnostics_18pa/01_export_18pa/report.md), [3D diagnostics](../outputs/runs/runner_check_2026-10-06/01_1000p_100_alotpa_3djawn/diagnostics/report.md), and [release sweep](../sweeps/sweep_results.csv). These averages use inlet-speed proxies with unverified outlet completion. Geometry and operating changes prevent interpreting the rows as a controlled 2D/3D effect or a size comparison.

The notebook's “7% haemolysis” beside the 18 Pa 2D case is inconsistent with the retained calculated HI; use the explicit percentage above. The early expectation that a small NIH necessarily meant broken code was also too strong: different geometry, stress, boundary conditions and empirical coefficients can all produce a low value.

The two inputs `z_inputs/INPUT_TEST_3D.csv` and `z_inputs/1000p_100_alotpa_3djawn.csv` have identical SHA256 hashes. Their matching runner results are duplicate analysis of one export, not independent simulations or replicates. The filename containing “1000p” reconstructs 100 paths; do not infer release count or mesh size from that label.

In this 100-path example, point-sampled stress has mean 2.52 Pa and maximum 75.63 Pa. Unweighted transit has mean 65.98 ms, P95 108.47 ms, P99 121.14 ms and maximum 137.20 ms. With an assumed 250 ms systole, these are respectively about 26%, 43%, 48% and 55% of systolic duration. This makes the frozen-field assumption material: some paths last through a substantial fraction of systole. It is a timescale warning, not a measured quasi-steady error.

GW HI3 is about 2.16 times HI2. Replacing the HI2 time quadrature with its exact piecewise expression changes the device result only about 0.34%; reversing stress-duration pairs changes HI2 about 49%. This supports an accumulator-ordering explanation rather than blaming the HI2/HI3 difference entirely on time quadrature. For HI3, removing the speed floor changes nothing, mean-endpoint stress changes the result about −1.40%, and stride-2/stride-4 thinning changes it about +0.50%/+2.00%. Thinning existing data does not establish true convergence.

The 100-path GW/HO HI3 ratio is about 68.45, much larger than those local quadrature changes. Coefficient choice is therefore a major source of absolute prediction variation in this example. It should be disclosed alongside numerical sensitivity rather than hidden by reporting GW alone.

An additional calculation from the saved path table shows that the five highest-HI3 paths represent about 5.01% of the inlet-speed proxy weight but contribute **73.44%** of the weighted GW HI3. The top ten contribute 85.11%. This is preliminary evidence that a small high-exposure tail controls the mean; it helps explain why seeding matters and why the median can appear stable while haemolysis changes. It does not identify the anatomical hotspot without spatial inspection. Reproducible values are in [path contributions](preliminary_100_path_contributions.csv).

Empirical screening of the same 100 paths gives maximum RBC utilisation 0.1845 at approximately 39.93 Pa for 4.08 ms, and platelet-lysis utilisation 0.07072 at approximately 36.74 Pa for 5.07 ms. None crosses the selected traced curves. This supports only the statement “no assessed path crossed the selected empirical lysis reference.”

### Release-count sensitivity

| Requested paths | Reconstructed paths | GW HI3 NIH | HO HI3 NIH |
|---:|---:|---:|---:|
| 100 | 100 | 19.2428 | 0.28112 |
| 1,000 | 985 | 23.9506 | 0.36165 |
| 5,000 | 4,809 | 28.1824 | 0.41450 |
| 10,000 | 9,709 | 28.7538 | 0.42177 |
| 30,000 | 29,451 | 29.8096 | 0.43533 |
| 100,000 | 97,331 | 28.1020 | 0.41310 |
| 300,000 | 292,782 | 28.1119 | 0.41327 |

NIH units are mg/100 L under the repository conversion. Table assembled from the three release-study CSVs; it assumes their notebook-described common development field and settings. The underlying solver provenance should be retained to confirm that assumption.

The two largest GW estimates differ by approximately **0.035%**, and HO by **0.041%**. This supports practical stability of the mean against release count for this flow field. The 30,000-path GW value is still about 6.04% above the 300,000 reference, so the old interpretation of 30,000 as an adequate converged reference should be updated.

The 5,000-request result is approximately +0.25% from the 300,000 reference for GW HI3 and +0.30% for HO HI3. That makes it a useful provisional production setting, conditional on confirming the other metrics and difficult new cases. It does not mean every count above 5,000 is within 5%: the 7,000 and 9,000 GW results are about +6.11% and +7.04%. Such nonmonotonic behaviour is compatible with changing seed locations sampling concentrated high-damage regions; it does not alone prove the mechanism.

At 300,000 requested paths, 292,782 are reconstructed: 97.594% by count. The unrepresented 2.406% cannot be translated into missing flow or trapped blood. The saved SA exceedance fraction is 0.0001981, equivalent to about 0.0198% of paths (58 paths), despite a median near 0.6375. Thus “zero exceedance” in a small sample does not establish absence of a rare tail. The stored `P_SA_above_hellums` column is a fraction, not already a percentage.

### Configured segment count

At 5,000 requested paths, all segment settings reconstruct 4,809 paths, but their damage estimates differ strongly:

| Configured segments | GW HI3 NIH, mg/100 L | Difference from the 5,000-setting result |
|---:|---:|---:|
| 25 | 1.4688 | −94.79% |
| 50 | 4.2347 | −84.97% |
| 100 | 23.3021 | −17.32% |
| 200 | 27.8235 | −1.27% |
| 500 | 28.1795 | −0.0103% |
| 1,000 | 28.1824 | 0% |
| 2,000 | 28.1824 | 0% |
| 5,000 | 28.1824 | reference |

The S1000/S2000/S5000 files are **byte-identical**, checked by SHA256. This is evidence that raising the configured count no longer changes this export. It is not three independently finer integrations. The very low S25/S50 values may reflect truncation before high-stress regions or path completion, as well as resolution; equal path counts cannot distinguish these possibilities. Check endpoint positions, path lengths and transit times before describing the sweep as pure integration convergence.

A 1,000-segment setting is therefore a defensible candidate for avoiding the observed export cap in this development case, with completion checks in all final cases. The available data do not establish a universal 1,000-segment requirement.

### Tracking tolerance

| Grid-relative tolerance | Paths | GW HI3 NIH | GW change from 0.01 | HO HI3 NIH |
|---:|---:|---:|---:|---:|
| 0.0200 | 4,798 | 27.6979 | −1.72% | 0.40902 |
| 0.0100 | 4,809 | 28.1824 | reference | 0.41450 |
| 0.0050 | 4,810 | 26.2662 | −6.80% | 0.39492 |
| 0.0025 | 4,811 | 27.3980 | −2.78% | 0.40640 |

This is CFD-Post tracking tolerance, not a CFX nonlinear residual tolerance. Results vary nonmonotonically, with slightly different exported populations. Retaining the notebook's 0.01 setting may be a pragmatic choice, but “all differences <5%” is unsupported. Mean variation over settings is not a reliable error bound. If the final size effect is only a few percent, this observed sensitivity could be comparable to the effect and needs further checking.

The [sensitivity figure](export_sensitivity.png) and its [PDF version](export_sensitivity.pdf) summarise the three studies. Separate denominators are stated in the caption; the sweeps must not be added into a total uncertainty estimate.

### Mesh study recorded in the notebook

Page 32 lists nominal element sizes 0.175, 0.100 and 0.0571 mm; cell counts 110,008, 330,188 and 998,156; and a stress quantity of 18, 24 and 25 Pa. Nearby notes refer to maxima, so confirm whether this is maximum stress or another consistently defined quantity.

The notebook reports a fine-grid GCI near 2.4%, using assumed order p = 2 and safety factor 1.25. The arithmetic is reproducible: with refinement ratio near 1.75, \(1.25(1/25)/(1.75^2-1)\approx2.42\%\). However, the rounded sequence 18→24→25 itself implies an observed order near 3.20, not 2. An assumed second-order estimate and a measured-order estimate are different claims. Rounded values, uncertain metric definition and missing solver records prevent treating this as a verified uncertainty bound.

This study supports decreasing changes in one stress metric on the simplified 2D geometry. It does not establish mesh independence of 3D HI3, tail exposure or valve-gap flow. It is worth retaining as development evidence, with the calculation convention made explicit.

### Partner's supplied 2D and 3D mesh workbook

The subsequently supplied [2D and 3D Mesh Analysis.xlsx](../papers/2D%20and%203D%20Mesh%20Analysis.xlsx) establishes completed studies in **both dimensions**. Its four worksheets are `y+ calculator`, `Inflation height calculator` (hidden), `2D Base Ref` and `3D Mesh Analysis`. Cell values, stored formula results and embedded method illustrations were inspected. The underlying CFD was not rerun, and the workbook was not modified. Extracted cell values/formulas are preserved in [mesh_workbook_cells.json](mesh_workbook_cells.json).

The workbook records mean/maximum velocity, pressure and wall shear stress, y+ measures, mesh counts and some solve times. The Richardson calculation uses **maximum velocity**, not wall shear, scalar `varshear`, HI or NIH. The latter blood-damage quantities are absent from these tables. This is complementary to the user's trajectory/export studies: it investigates the spatial CFD mesh rather than the number or resolution of exported paths.

#### Near-wall mesh design

The y+ calculator uses density 1,060 kg/m³, viscosity 0.0035 Pa·s, a trial speed of 1.57 m/s and a characteristic span of 0.007537 m, yielding Re ≈ 3,584. Its target y+ is 1. It estimates an initial first-layer height of 64.92 μm, then divides by a measured maximum y+ of 6.27 to propose 10.35 μm. This is a recorded design iteration, not evidence that every final wall location meets y+ = 1.

The inflation calculator uses geometric layer growth, h_n = h_1 r^(n−1), and total thickness h_1(r^n−1)/(r−1). The 3D sheet specifies global/local first-layer heights of 30/18.4 μm, 17 maximum layers and growth rate 1.1. These are shown as common settings rather than separate values for each 3D mesh. The actual layer counts and any changes between meshes require solver/mesh records to confirm.

#### 3D results: strongest evidence for maximum-velocity convergence

Use the coherent table in `3D Mesh Analysis!C24:F34`. The earlier block at rows 14–23 has base-column values apparently pasted against mismatched labels (for example, 3,983,051 beside `yplusval`); it should not be used for physical interpretation. The later table correctly aligns those base values with their quantities.

| Quantity | Coarse | Base | Fine |
|---|---:|---:|---:|
| Global element size, mm | 1.000 | 0.500 | 0.250 |
| Local element size, mm | 0.7071 | 0.3536 | 0.1768 |
| Elements | 749,820 | 3,983,051 | 22,781,658 |
| Nodes | 232,893 | 1,186,948 | 6,098,820 |
| Maximum velocity, m/s | 2.45723 | 2.53729 | 2.55705 |
| Average velocity, m/s | 0.681712 | 0.698183 | 0.703475 |
| Average pressure, Pa | 14,828.0 | 14,804.8 | 14,810.0 |
| Maximum pressure, Pa | 16,384.8 | 16,203.7 | 16,217.8 |
| Average wall shear, Pa | 2.17233 | 1.91189 | 2.60200 |
| Maximum wall shear, Pa | 403.893 | 424.659 | 548.196 |
| Maximum y+ | 13.2491 | 7.69725 | 5.80777 |
| Average y+ | 0.255544 | 0.133912 | 0.179148 |
| Recorded solve time | — | 4 hours | 49 hours |

Using nominal spacing ratio r = 2, with phi_1 fine, phi_2 base and phi_3 coarse, the workbook evaluates:

\[
p=\frac{\ln[(\phi_3-\phi_2)/(\phi_2-\phi_1)]}{\ln r}=2.0184987,
\qquad \phi_{ext}=\phi_1+\frac{\phi_1-\phi_2}{r^p-1}=2.5635253\;\mathrm{m/s}.
\]

Independent arithmetic reproduced the workbook's p, extrapolation and signed relative differences from the extrapolated value: **4.1464% coarse, 1.0234% base and 0.2526% fine**. These are Richardson extrapolation differences. They are not formal GCI values: the displayed calculation contains no GCI safety factor. Richardson extrapolation is the basis of Roache-style GCI methods, but the names should not be interchanged. See [NASA's grid-convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html). The exact reference used by the partner remains to be recorded.

The base-to-fine maximum-velocity change is +0.779%, and the average-velocity change is +0.758%. The recorded run times increase from 4 to 49 hours, illustrating the computational trade-off, though hardware and solver settings are not given. This supports a practical argument for the base mesh **for the velocity quantities assessed**, conditional on the refinement assumptions. It does not establish that the base mesh has already been selected for all final cases.

Important qualifications are visible in the same table. Maximum wall shear rises **29.09%**, and average wall shear rises **36.10%**, from base to fine. Average wall shear is nonmonotonic across the three meshes. Thus velocity stability coexists with appreciable shear sensitivity. Wall shear is not the scalar viscous stress integrated along the blood trajectories, so these changes neither quantify HI3 error nor prove that HI3 changes by the same amount. They do show why velocity convergence alone cannot establish convergence of the damage calculation.

Average and maximum pressure change by only +0.035% and +0.087% from base to fine, but these percentages use the tabulated pressure levels as denominators. They do not establish equivalent relative accuracy in the smaller pressure difference across the device. The workbook does not state the averaging domains/weights or pressure reference; record those definitions for the report.

The fine mesh has average y+ below 1 but maximum y+ = 5.81. Therefore do not write that y+ < 1 everywhere. Nominal global/local spacings are halved, but element counts imply effective three-dimensional refinement ratios of approximately 1.745 and 1.788 rather than 2. This does not automatically invalidate a nominal-spacing study, especially with inflation/local controls, but the chosen definition of h and consistency of refinement should be documented. The inverse-node-count chart axis is a plotting coordinate, not a physical 3D mesh spacing.

#### 2D results: stable bulk quantities, less clear extrapolation

`2D Base Ref` contains two populated result blocks, rows 16–25 and 32–43. Their differing element counts and results should be preserved as separate series until their mesh revisions are identified. The upper block spans 7,501–334,579 elements. The lower block spans 9,977–695,470 elements and supplies the explicit maximum-velocity Richardson calculation. Both share a base result of 150,118 elements; this alone does not make the series interchangeable.

For the lower block, maximum velocity across the seven levels is 1.64661, 1.67414, 1.71726, 1.72267, 1.73045, 1.72846 and 1.73128 m/s. The final four values are close but nonmonotonic. From the 150,118-element base to the 695,470-element finest case, maximum velocity changes +0.048%, average velocity −0.165%, average pressure +0.042%, maximum pressure +2.93%, average wall shear +69.48% and maximum wall shear +25.80%.

The workbook's first extrapolation triplet uses lower-block maximum velocities from 16,977, 36,185 and 70,700 elements and nominal r = 1.5. It gives p = 5.1194 and extrapolated speed 1.723446 m/s. The next triplet gives p = −0.8960, and the next produces `#NUM!` because the successive differences change sign. Consequently the workbook's listed 2D “convergence errors” are differences from an extrapolation based on that earlier triplet, not a demonstrated reliable finest-grid asymptotic error bound. The high first p, negative next p and nonmonotonic sequence warrant explicit qualification. The workbook suggests mesh quality as a possible explanation; this is a hypothesis, not a demonstrated cause.

These completed studies substantiate a report narrative of **mesh refinement and quantity-specific convergence assessment**. The 3D maximum-velocity study has a consistent decreasing-increment trend under its stated nominal spacing ratio, while wall shear remains more sensitive and the 2D extrapolation is less robust. Neither table directly establishes haemolysis mesh convergence. The earlier notebook 18/24/25 Pa study remains separate historical evidence and should not be merged with these wall-shear series.

## 9. Literature comparisons and claims to correct

The value **37.15 ± 12.42 mg/100 L** is traceable to Zaman et al., *In Vitro Hemolytic Performance of the Novel Realheart Total Artificial Heart Versus SynCardia Total Artificial Heart With Human Blood*. The study uses a modified systemic circulation rig with targeted 5 ± 0.25 L/min, 100 mmHg and diluted human blood at Hct 30 ± 2%. It is an in vitro result, not a clinical patient measurement. The accessible abstract does not establish the 50 cc identification used in the code constant; confirm model size from the full methods before retaining that label. [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/42010686/).

The local `papers/realheart.pdf` is a **different** full-body circulation paper. Its reported SynCardia mgNIH is **226.0 ± 30.6 mg/100 L**, with systemic and pulmonary circulation at a targeted 6.2 ± 0.25 L/min. It is not the source of 37.15. [Full-body study](https://www.nature.com/articles/s41598-026-68178-2).

Numerical proximity to either experimental value is contextual consistency at most until geometry, operating waveform, circuit scope, blood properties and exposure definitions are matched. Do not tune model coefficients to force agreement and then call that independent validation.

The code still prints “ASTM acceptable limit: 100 mg/100 L.” No basis for treating that number as a universal pass/fail limit was established here. The published scope of [ASTM F1841-25](https://store.astm.org/f1841-25.html) concerns a blood-pump haemolysis assessment protocol. A testing protocol should not be equated to a universal clinical acceptance threshold without a specific supporting clause and applicable context.

| Historical wording or implication | More accurate report wording |
|---|---|
| “Validated because NIH is near the literature value” | “Compared with literature for context; matched physical validation remains outstanding.” |
| “Flow-weighted device HI” for current unverified runs | “Mean HI of reconstructed paths, weighted by inlet-speed proxies; completion unverified.” |
| “7% haemolysis” for the 18 Pa case | Explicit saved HI3 = 7.3226 × 10⁻⁶%. |
| “No activated platelets / safe below Hellums” | “No sampled paths exceeded the selected SA reference in this export.” |
| “Whole-cycle haemolysis from systolic snapshots” | “Ejection-weighted quasi-steady estimate over the represented systolic intervals.” |
| “Fixed time segments” | “Variable exposure increments reconstructed from distance and mean endpoint speed.” |
| “All tolerance effects are below 5%” | “Observed GW HI3 change reaches −6.80% relative to tolerance 0.01.” |
| “Independent simulations named INPUT_TEST_3D and 1000p…” | “Identical exports analysed through different runs.” |
| “First experimental evidence / physical loop validation” | Omit unless actual experiments and a substantiated novelty review are added. |

The old presentation claims about regulatory status, clinical eligibility dimensions, survival and treatment should be sourced separately if used in the Introduction; this methods synthesis does not verify them.

## 10. Completing the 50 cc / 35 cc investigation

The user has confirmed **two comparison families**, both at 20%, 40%, 60% and 80% of systole:

| Comparison | Prescribed constraint | Main question |
|---|---|---|
| Same BPM | Both sizes use the same beat rate and their own defined stroke/volume histories; keep systolic fraction consistent unless explicitly varied | How does downsizing change flow and exposure at the same operating rate? |
| Matched flow | Both sizes have the same instantaneous volume flux at corresponding phases, with the driving prescription documented | How does geometry change flow structure and exposure at the same imposed flow? |

The first family gives eight size–phase cases. If the four 50 cc baseline fields are reused as the matched-flow references, four additional 35 cc fields would give twelve unique CFD cases supporting both comparisons. Otherwise the number depends on the selected reference waveform. This is a planning count, not a claim that those cases exist. Distinguish matching instantaneous phase flow from matching cycle-mean cardiac output. If the smaller device's BPM or systolic duration is adjusted to achieve the matched waveform, document that change; if flux is imposed independently of a physically consistent stroke, describe the case as a controlled geometry comparison. A reduced-BPM study or design modification is an optional extension after these two families are reliable.

| Final-case item | Value to insert |
|---|---|
| CAD identification and dimensions | [50 cc model revision; 35 cc revision; actual phase volumes] |
| Downsizing method | [uniform/selective scaling; valve and connector exceptions] |
| Diaphragm displacement and shape | Confirmed S = 20 mm and diaphragm-as-inlet sinusoidal speed; [shape equation; V(x); 35 cc stroke/scaling; velocity direction/profile] |
| Operation | [BPM; systolic fraction; intended net stroke; phase fractions] |
| Boundary conditions | [imposed normal flux/speed; driving surface; opening pressure and reference; valve state] |
| Fluid and turbulence | [density; viscosity; Newtonian assumption; exact CFX turbulence/wall treatment] |
| Mesh | 2D/3D workbook reviewed; [confirm chosen production mesh, case/phase correspondence, topology, gap resolution, quality, averaging definitions, refinement convention and exact reference] |
| Solver convergence | [residual definition/target; mass imbalance; stable Q, pressure and stress monitors] |
| Trajectory export | [seeding surface and distribution; requested/recovered count; count cap; tolerance; exit termination] |
| Analysis population | [selected exits; coverage; supplied flux weights or documented proxy] |
| Validation | [matched experimental/benchmark quantities, or explicitly unavailable] |

Prioritise the remaining two weeks around the following dependencies:

1. **Reconcile volume, motion and flow.** Freeze the two CAD variants, confirm V(x), derive Q for each phase and select the operating constraint. Resolve the preliminary 28 L/min issue before using that field as a physiological baseline.
2. **Establish one traceable reference case.** Record solver conditions, signed fluxes, pressures and trajectory completion. Use the existing 5,000-path / 1,000-segment / 0.01-tolerance settings as candidates, with their measured limitations visible.
   Use the partner's reviewed 2D/3D workbook for the mesh evidence, with its maximum-velocity convergence separated from the still-unmeasured HI3 mesh sensitivity.
3. **Run the same-BPM matrix and the additional matched-flow cases consistently.** Keep blood properties, metric definitions, plotting scales and analysis populations aligned. Repeat release/tolerance/completion checks in the smaller geometry and the most difficult phases; existing convergence does not transfer automatically.
4. **Quantify mechanisms.** Compare stresses and passage-time distributions, map candidate hotspots onto the geometry, and relate HI3 differences to effective stress and transit. Add region/source totals only if full segment data are exported; the current hotspot CSV holds only three selected segments per path.
5. **Check whether claimed effects exceed observed numerical changes.** If they do not, report the ranking as unresolved at the present resolution. Do not assign significance tests to thousands of streamlines as though they were independent experimental replicates.

For each final case, report signed outlet flow, relevant chamber/outlet pressure statistics and pressure difference, peak jet speed and location, stress distributions, transit distributions, GW and HO HI3, HI2 sensitivity, SA, reconstructed/selected counts and represented flux coverage. Pressure, volume-weighted stagnation and mass balance require solver field/surface exports; they cannot be recovered reliably from the current scalar streamline files alone.

Candidate figures are: geometry and boundary schematics; matched-scale velocity and scalar-stress contours; HI3 versus phase with flow alongside; 35/50 ratios versus phase; transit versus effective stress; cumulative weighted damage contribution; and the existing numerical sensitivity figure. Streamline density is not a quantitative measure of recirculating volume.

### Combining phase snapshots

If positive ejection flux \(Q_k\), representative interval \(\Delta t_k\) and compatible per-case HI averages are available,

\[
HI_{QS}=\frac{\sum_k HI_kQ_k\Delta t_k}{\sum_kQ_k\Delta t_k}.
\]

The denominator is represented ejected volume. With only 20/40/60/80% samples, specify how the unsampled ends of systole and interval weights are handled; equal weights are not implied. The volume history can define each interval's ejected volume directly. Handle reverse flow and multiple exits explicitly rather than allowing signed weights to cancel incompatible populations. Label the result an ejection-weighted quasi-steady estimate; filling and changing-flow histories are absent.

## 11. Adaptable prose for the final report

**Team context and contribution — user-confirmed:**

> The wider investigation combines CAD development, mesh assessment, steady-state CFD and programmatic blood-damage analysis. Multiple CAD configurations were constructed using a revolved sinusoidal diaphragm profile subtracted from the fluid domain, together with a simplified SynHall-inspired disk valve positioned to represent the selected configuration. The project partner undertook two- and three-dimensional mesh studies monitoring velocity, pressure, wall shear stress and y+. Richardson extrapolation of the three-dimensional maximum velocity yielded an observed order of 2.02 using the nominal spacing ratio of two. The base-mesh velocity differed from the extrapolated value by 1.02%, while wall shear remained appreciably mesh-sensitive. The emphasis of the present account is the programmatic workflow linking CFD streamline exports to exposure histories, haemolysis predictions and numerical diagnostics. All chamber and valve configurations remain fixed within each steady-state simulation.

**Methods overview — supported now:**

> A computational workflow was developed to investigate the relationship between ventricular downsizing, haemodynamics and mechanical blood-damage exposure in a SynCardia-inspired artificial-heart chamber. Initial simplified two-dimensional calculations were used to develop the geometry, boundary conditions and streamline post-processing. The workflow was subsequently applied to preliminary three-dimensional steady flow fields. Diaphragm-driven ejection was approximated by prescribing an inlet velocity on the diaphragm surface, corresponding to the instantaneous speed of a sinusoidal displacement over a 20 mm stroke. Scalar viscous stress and velocity magnitude were exported along CFD-Post streamlines, and segment exposure durations were reconstructed from spatial distance and average endpoint speed. Haemolysis was evaluated using a linearised power-law accumulation method, with alternative accumulation and coefficient choices retained for sensitivity assessment. Stress accumulation and empirical stress-duration screening provided complementary measures of mechanical exposure.

**Verification paragraph — supported now:**

> Numerical assessment distinguished calculation verification from CFD and trajectory sensitivity. Analytical tests checked the uniform-stress solution, subdivision behaviour and agreement between implementations. Release-count studies demonstrated close agreement of mean HI3 between 100,000 and 300,000 requested paths for the development field. A lower release count of 5,000 produced a similar mean in that case, although intermediate counts varied nonmonotonically. Increasing the configured segment count beyond 1,000 produced identical exports. Tracking-tolerance sensitivity remained measurable, including a 6.8% change in GW HI3 relative to the selected 0.01 setting. These checks support a provisional analysis configuration but do not establish convergence for every geometry and phase.

**Planned comparison paragraph — convert to past tense only after completion:**

> The comparative study will evaluate the 50 cc baseline and 35 cc derivative using separate fixed-geometry, steady-state simulations representing four systolic fractions: 0.20, 0.40, 0.60 and 0.80. The diaphragm inlet velocity will remain constant within each solve and will be selected from the prescribed sinusoidal motion at the corresponding phase. One comparison will hold BPM constant while using each geometry's defined stroke prescription; a second will compare corresponding geometries at matched instantaneous flow. The phase-specific chamber configuration and driving prescription will be documented for each family, using [final verified mesh and trajectory settings]. Haemolysis differences will be interpreted alongside pressure, flow, stress and passage-time distributions.

**Results placeholder:**

> At [phase], downsizing changed [signed Q / pressure metric] from [value] to [value]. The [flux-weighted or explicitly proxy-weighted] GW HI3 changed by [percentage], while HO predicted [direction and percentage]. The associated transit-time and effective-stress changes were [values]. Spatial inspection located the principal exposure change in [region]. The observed size effect was [larger than/comparable to] the changes observed during [specific numerical checks].

**Limitations paragraph — applicable to the current scope:**

> The predictions represent passive trajectories within prescribed steady fields and do not reconstruct transient blood histories through filling, ejection and valve motion. Simplified valve geometry and excluded leakage routes may omit relevant exposure mechanisms. The stress measure includes molecular viscous stress derived from the resolved mean velocity field, while biological responses are represented by empirical exposure models. Absolute NIH depends strongly on the selected coefficients and conversion conventions. These calculations therefore support conditional engineering comparisons within the modelled conditions; validation of full-device haemocompatibility requires appropriately matched experimental evidence.

## 12. Reproducibility and evidence map

The current main entry point is [run_haem.py](../run_haem.py); [README](../README.md) and [diagnostics guide](../HAEM_DIAGNOSTICS.md) explain manifest and output conventions. Core summaries/plots/empirical screens use all reconstructed paths and speed-proxy weights. Diagnostics can use a selected completed population and supplied weights. Use the diagnostic comparison for controlled populations, and always state which result is plotted.

The [investigation plan](../HAEMOLYSIS_INVESTIGATION_PLAN.md) is useful prior planning, not an execution record. Some items it proposes now exist; other claims need the clarifications above. Historical scripts in `old/`, the recirculation diagnostic and the portable pipeline copy document development but should not silently replace the current analysis implementation.

This review reran the existing `test_haem_diagnostics` and `test_empirical_thresholds` suites: **8 tests passed**. This verifies the tested calculations and handling, not the underlying CFD or biology. No production solver cases or large haemolysis sweeps were rerun. Analysis code and existing results were left unchanged.

New evidence can be regenerated with:

```powershell
.venv/Scripts/python.exe report_context/build_evidence.py
```

The script derives tables and a sensitivity figure from saved CSVs, computes preliminary high-damage-path contributions, and records source/segment-export SHA256 hashes in [evidence provenance](evidence_provenance.json). It does not simulate flow or generate new blood-damage predictions. The colour palette uses restrained purple and blue with distinct markers.

| Source | Main use |
|---|---|
| Notebook pp. 4–8 | Initial export workflow and non-heart proof of concept |
| Notebook pp. 18–34 | 2D development, material settings, mesh work and early 3D transition |
| Notebook pp. 35–51 | Release/tolerance studies, empirical screening and flow concern |
| Notebook pp. 52–54 | Latest volume-motion formulation and four-phase plan |
| `papers/2D and 3D Mesh Analysis.xlsx` | Partner's completed mesh studies; y+/inflation design, bulk-flow and wall-shear results, maximum-velocity Richardson extrapolation |
| `haem5.py` | Executed equations, velocity floor, proxy averaging, SA and NIH conversion |
| `haem_diagnostics.py` | Verification variants, distributions, completion and weighting |
| `haem5_empiricalthresholds.py` | Longest-contiguous exposure screen |
| `sweeps/sweep_results.csv`, `300k_sweep_results.csv`, `a1to10k_sweep_results.csv` | Release-count evidence |
| `sweeps/segment_sweeps/segment_sweep_results.csv` | Segment-setting evidence |
| `sweeps/tolerance_sweep_results.csv` | Tracking-tolerance evidence |
| `papers/taskin-etal.pdf` | Power-law definitions, coefficient sets and method limitations |
| `papers/marom-etal.pdf`, `papers/nihms-769456.pdf` | TAH SA context and repeated-exposure limitations |

**Update boundary:** future results should replace the planned case descriptions and placeholders while preserving the distinction between observed data, numerical assumptions and biological interpretation. No conclusion on whether the 35 cc model outperforms or underperforms the 50 cc model is supported yet.
