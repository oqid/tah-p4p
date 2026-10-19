# Valve cutaway figures (35 cc, 120 BPM)

Flow at the systolic snapshots with the valve geometry in view, in the style of Dumont et al. (2007), Fig. 5.
Built only from data already on disk: the pathline exports (`35cc_inputs/120BPM`) and the per-phase STP fluid
domains (`E:\35cc Retry\3D 35cc STP Files`). No new ANSYS exports were needed.

Valve angles in titles are the design values from `hemolysis_explore/operating_point.py` (24.9/49.8/75/50.1/25.2°),
as in F3. Measured from the STP disc against the outlet-tube axis they come out 3–4° larger at low opening, because
the tube axis is not exactly the valve-housing axis.

## Figures (`figures/`, PNG + PDF)
**+Y (outflow) up, as in the CAD drawings**
| File | Content |
|---|---|
| `cutaway3d_Yup_stress` | Rear half of the domain cut on the mid-plane, viewed normal to the cut. 800 paths per phase coloured by local scalar stress on the F3 scale (magma, 0.01–100 Pa). The wake vortices on the minor-orifice side show clearly. |
| `cutaway3d_Yup_regions_front` | Outflow tract regions: reversed axial flow (blue isosurface, u·n̂ < −0.03 m/s), no path reached (lilac), path segments ≥ 15 Pa (pink) and grey context paths. |
| `cutaway3d_Yup_regions_oblique` | The same as above, from downstream and to the side, so the disc face is visible. |
| `F3_companion_section` | F3 phases (16.67/50/83.33%). Mid-plane section; top row local stress (F3 scale), bottom row share of linear GW HI3 dose per mm² of the ±1 mm slab (F3 dose colour map). |
| `F3_companion_section3d` | The same two fields drawn on the mid-plane inside the translucent half-shell, with the disc and pivots in 3D. |

**Earlier versions (original CFD orientation)**: `valve_section_zoom/full_35cc_120` (speed) and `valve_cutaway3d_35cc_120`.

## Method notes (for captions)
* The section plane is perpendicular to the disc hinge axis and passes through both ports, so it is the chamber's mid-plane.
* Section fields average every path node within ±1 mm of the plane onto a 0.5 mm grid (Gaussian, σ = 0.7 mm). The lines show in-plane direction only.
* Flow direction comes from consecutive points on each path, because the exports have no vector components.
* Grey or lilac means no inlet-seeded path came near. In a steady (frozen) field such paths cannot enter a closed recirculation zone, so a large void beside the disc is consistent with separated flow. Small voids may just be gaps between the ~4,900 paths.
* The 3D isosurfaces come from path nodes binned on a 0.6 mm voxel grid (σ = 0.8 mm), downstream of the valve only.
* Dose per node follows `explore_core.Snapshot.inner` (A^(1/β) τ^(α/β) dt, GW constants, the same velocity floor).
* Paths that end before the outlet (7.7% at 16.67%, 8.9% at 83.33%) stop about 0.13 mm from the disc or the valve ring, mostly in the disc–ring gap. They run into the wall; they are not trapped in the wake.

## Rerun (from this folder)
Needs `numpy scipy matplotlib gmsh pyvista`.
```
python load.py; python geom.py                   # CSV -> ph*.npz, STP -> geom*.pkl (CSV_DIR / STP_DIR env vars)
python fields.py; python fields_y.py; python volumes.py
python fig_section.py; python fig_f3companion.py
for k in 16 33 50 66 83: python render3d.py k; python render_v2.py k front_valve stress;
                         python render_vol.py k front_valve; python render_vol.py k oblique_down
for k in 16 50 83:       python render_plane.py k stress; python render_plane.py k dose
python fig_3d.py; python fig_film.py; python fig_f3companion3d.py
```
Camera presets live in `VIEWS` in `render_v2.py`. Disc face IDs are 47–51 and the pivot and ring IDs are in `common.py`.

## Round 3 additions
| File | Content |
|---|---|
| `cutaway3d_Yup_stress_34above` | 3/4 view (swung 30° and raised 22°). Stress-coloured paths, a faint ghosted backflow volume with a dashed outline, and an open circular arrow giving the sense of recirculation. |
| `cutaway3d_Yup_stress_backflow_front` | The original front view with the same backflow outline and recirculation arrow. |
| `cutaway3d_Yup_flow_features` | Labelled key at 16.67% and 50%: major/minor-orifice jets, shear layer, separation at the disc edge, recirculation zone (backflow), and the thin separated layer at the housing step. |
| `F3_companion_section3d` | Section planes now 72% opaque (colour bars match). |

* Backflow outline: the screen silhouette of the axial-velocity < −0.03 m/s isosurface from `volumes.py` (largest region, rear half), drawn dashed on top.
* Recirculation arrow: placed at the centroid of that region. Its sense is anticlockwise viewed from +H at every phase: the jet (major orifice, −E2 side) drives the near side downstream and the wall side returns upstream. A velocity-gradient check gave the same sign at four of five phases. At 83.33% the fit window took in the jet and returned the opposite sign.
* Feature points (`features.py`, `feats3d`): jet = mean of the top 2% of path-point speeds 2–12 mm past the disc centre (major = −E2 side); shear layer = midpoint of jet and recirculation centre (indicative); separation = disc tip on the cut plane. These label positions are approximate, for orientation; they are not measured quantities.
* Rerun: `python render_v4.py <k> <front|sw30e22|sw-25e-15>`, then `python fig_annot.py <tag> <name>` or `python fig_annot.py key sw30e22 <name>`.

## Round 4
* Backflow is now a soft multiply-blend shade (#8DB8E8 at 55%, mask blurred 2.5 px) from a second render pass of the backflow volume with the same camera (`render_mask.py`). It tints the paths without covering them. The dashed outlines and ghost surface are gone, and the rotation arrow has been removed. The recirculation turns anticlockwise in these views at every phase, so state that in the caption.
* `34below` was dropped. To rebuild: `python render_mask.py <k> <tag> full`, then `fig_annot.py` as before.
