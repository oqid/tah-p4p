# Streamline haemolysis diagnostics

`haem_diagnostics.py` implements the numerical conventions used by `haem5.py`:
GW/HO/TZ coefficients, endpoint stress, segment distance divided by average
endpoint speed, the original velocity floor, HI2 and correctly linearized HI3.
It leaves haem5 unchanged. It supports both 2D and 3D **plain Line** CFD-Post
exports with the same columns as export_18pa.csv (including node numbers,
varshear in Pa and speed in m/s). Export forward trajectories from the inlet.
Ribbon/Tube [Faces] exports are deliberately rejected because reconstruction
can lose physical direction. It requires haem5 and its installed dependencies.

## Run

```powershell
python haem_diagnostics.py --self-test
python haem_diagnostics.py export_18pa.csv --expected-seeds 1000 --outdir outputs/audit_2d
python haem_diagnostics.py case50.csv case35.csv --outdir outputs/audit_scaling
python haem_diagnostics.py --manifest diagnostic_cases.json --outdir outputs/audit_cases
```

Use a new output directory for each run; existing case directories are not
overwritten. CLI defaults apply to all cases; manifest fields override them.
Manifest paths are relative to the manifest file. Example (replace filenames):

```json
{
  "cases": [
    {
      "path": "case50.csv",
      "label": "50cc baseline, peak ejection",
      "size_cc": 50,
      "geometry": "baseline",
      "phase": "peak_ejection",
      "mesh_cells": 1000000,
      "bpm": 120,
      "expected_seeds": 4000
    },
    {
      "path": "case35.csv",
      "label": "35cc scaled, peak ejection",
      "size_cc": 35,
      "geometry": "scaled",
      "phase": "peak_ejection",
      "bpm": 120,
      "expected_seeds": 4000
    }
  ]
}
```

Optional manifest fields:

- `outlet_plane`: `[x,y,z,nx,ny,nz]`, point in metres and **outward** normal.
- `outlet_tolerance`: endpoint distance tolerance in metres (default 0.0001).
- `weights`: path to a CSV with `streamline_id,flux_weight` for **every**
  reconstructed path; IDs are zero-based connectivity order, as in streamlines.csv.
  Use represented inlet normal flux times area, or mass flux. Any consistent
  units work. Verify correspondence using start_node_id and start coordinates.
- `flow_L_min`, `flow_model`, `stress_model`, and additional metadata are retained.
  Default flow_model is SST and stress_model is molecular viscosity times shear
  strain rate; change these for other exports.

Equivalent common outlet settings can be supplied with
`--outlet-plane X Y Z NX NY NZ --outlet-tolerance 0.0001`.
Use different planes per case if scaling changes outlet position/orientation.

## Outputs and interpretation

- `comparison.csv`: case metadata, selected-population HI/NIH, sensitivity results.
- Per-case `report.md`: principal values and interpretation limits.
- `diagnostics.json`: input SHA256, settings, percentiles, timing, fit-range flags,
  source distribution and deterministic subsampling results. Undefined results
  are JSON null, not zero.
- `streamlines.csv`: endpoints, weights, status, residence times, stress peaks,
  SA, all three coefficient sets and numerical sensitivity diagnostics.
- `hotspot_segments.csv`: three largest GW linearized-source segments per
  selected trajectory, with coordinates for comparison against CFD geometry.
  This is a candidate hotspot sample, not the complete spatial damage field.

Without an outlet plane, the selected population is **all exported paths,
completion unverified**. With a plane, a path must end within tolerance, start
upstream and approach in the outward direction. Other paths are excluded from
completed-pass averages and remain in streamlines.csv. This conservative test
does not know the outlet aperture or true solver termination reason. Inspect
endpoints against your geometry, and test tolerance sensitivity. The export must
terminate at the intended outlet; a trajectory continued far beyond it fails
this endpoint test rather than accumulating downstream damage as device damage.

Without a weight file, starting speed is used solely to reproduce haem5's
**inlet-speed proxy weighting**. This is not verified flow weighting.
Coverage uses only exported weights: missing seeds' flux cannot be inferred
from a count. It is not total inlet coverage when releases are missing.

Numerical variants:

- HI2 exact: analytically integrates the time factor for each piecewise-constant
  endpoint stress. This tests HI2 quadrature, not biological validity.
- HI2 reversed: reverses stress-duration pairs to expose ordering sensitivity;
  it is not another physical simulation.
- HI3 no floor: removes haem5's minimum-speed clamp. Positive-distance segments
  with zero mean speed make this diagnostic undefined; they are never silently
  removed from a weighted average.
- HI3 midpoint: uses mean endpoint stress as a stress-quadrature sensitivity.
- HI3 stride2/stride4: recomputes paths after retaining every second/fourth point
  and the endpoints. Includes chord-length changes. Large differences indicate
  resolution sensitivity; small differences cannot establish solver/export
  convergence or recover unsampled stress peaks.
- Random 25%/50% subsets, 30 repeats with seed 2026: existing-export sampling
  sensitivity. These are not confidence intervals, new-seed convergence or
  independent numerical solutions.

Point stress percentiles are not volume-weighted. Transit and SA percentiles are
unweighted over paths. Fit stress flags use per-path residence-time fractions,
then the selected path weights; total transit exceeding a calibration duration
is only a screening flag for variable stress. TZ explicitly checks its 50 Pa
lower calibration bound; no undocumented GW/HO lower bound is assumed.

Time-quartile source fractions refer to the weighted **linearized** HI3 source,
not additive fractions of final haemolysis. Transit relative to a half-cycle
assumes sinusoidal motion; it diagnoses the frozen-flow approximation, not its
error. NIH uses Hct=0.30 and Hb=150 g/L from haem5. There is no safety threshold
or automatic clinical acceptability verdict.

## Report-ready investigation

1. Define the comparison: 50cc baseline versus 35cc scaled, then isolate geometry
   changes. Prespecify operating conditions and acceptable performance criteria.
   Distinguish matched flow/head comparisons from intended operation at different
   stroke volumes or rates. Keeping 120 bpm for both does not keep cardiac output
   equal. Scale all dimensions consistently only where the actual design does so.
2. At each chosen ejection phase, justify chamber shape, valve position, inlet
   volume flux and pressure. Compare phase trends without treating independent
   snapshots as actual transient histories. A full-cycle claim requires inlet,
   outlet and valve conditions appropriate to filling and ejection.
3. Demonstrate mesh, release-population and integration convergence with new
   CFD/trajectory exports. Monitor HI3 as well as hydrodynamics; numerical
   uncertainty must be smaller than the geometry effect being claimed.
4. Report viscous-only SST stress and verify its scalar convention. Use HI3 as
   the primary result, HI2 and coefficient sets as sensitivity cases; report
   coverage, hotspots, damage distributions and operating-point robustness.
5. Validate matched hydrodynamics where measurements exist. Use Marom for
   qualitative flow/SA context and Taskin for haemolysis methodology, not as
   absolute NIH validation. Absolute acceptable haemolysis needs a justified
   matched experimental/reference criterion; CFD alone supports conditional
   comparative predictions, not a clinical safety finding.

Do not average phase HI values and call that transient whole-cycle haemolysis.
If an approximate ejected-volume-weighted snapshot summary is used, define the
phase intervals and flux weighting and label it a quasi-steady approximation.
Evaluate each outlet separately if there is reverse flow or multiple exits.




---

Run a single export:

`python haem_diagnostics.py export_18pa.csv --expected-seeds 1000 --outdir outputs/my_audit`

Or compare several:

`python haem_diagnostics.py case50.csv case35.csv --outdir outputs/scaling_audit`

A documented JSON manifest supports different outlet planes, flux weights, geometries, sizes, and cycle phases for each case.

The script reports:

- HI2/HI3 and NIH using GW, HO, and TZ coefficients.
- Stress, velocity, residence-time and calibration-range diagnostics.
- Velocity-floor, stress-sampling and HI2 integration sensitivity.
- Export-resolution and trajectory-subsampling sensitivity.
- Trajectory endpoints, optional outlet-completion classification and flux weighting.
- Coordinates of candidate damage hotspots.
- Per-case reports and a combined comparison CSV.

Verification: analytical checks and four automated tests passed, including synthetic 2D/3D cases and incomplete-path weighting. Your existing export reproduces both original HI values exactly. Its generated report (outputs/diagnostics_18pa/01_export_18pa/report.md) also shows that removing the velocity floor changes nothing; using mean endpoint stress changes HI3 by approximately −0.46%.

Outlet completion remains unverified until you supply the outlet plane. Subsampling diagnostics are explicitly distinguished from actual mesh or new-seed convergence.
