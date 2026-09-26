# export_18pa

Input: export_18pa.csv
Population: all_exported_unverified; weighting: inlet_speed_proxy
Paths: 969; selected: 969; selected/exported weight: 1.000000

## Metrics

- GW_HI2_percent: 4.2083241e-06
- GW_HI3_percent: 7.3225609e-06
- GW_HI2_exact_percent: 4.2301599e-06
- GW_time_fraction_tau_above_fit: 0
- GW_transit_above_fit: 0
- HO_HI2_percent: 1.6076096e-07
- HO_HI3_percent: 2.5307717e-07
- HO_HI2_exact_percent: 1.6213873e-07
- HO_time_fraction_tau_above_fit: 0
- HO_transit_above_fit: 0
- TZ_HI2_percent: 1.2547198e-06
- TZ_HI3_percent: 2.5983781e-06
- TZ_HI2_exact_percent: 1.2771149e-06
- TZ_time_fraction_tau_above_fit: 0
- TZ_transit_above_fit: 0
- TZ_time_fraction_tau_below_50Pa: 1
- GW_HI3_no_floor_percent: 7.3225609e-06
- GW_HI3_midpoint_percent: 7.2891593e-06
- GW_HI2_reversed_percent: 5.399302e-06
- GW_HI3_stride2_percent: 7.2068939e-06
- GW_HI3_stride4_percent: 7.3953174e-06
- GW_HI2_NIH_mg_per_100L: 0.44187404
- GW_HI3_NIH_mg_per_100L: 0.76886889
- HO_HI2_NIH_mg_per_100L: 0.016879901
- HO_HI3_NIH_mg_per_100L: 0.026573103
- TZ_HI2_NIH_mg_per_100L: 0.13174558
- TZ_HI3_NIH_mg_per_100L: 0.2728297
- GW_HI3_all_exported_percent: 7.3225609e-06
- SA_selected_median_dyne_s_cm2: 0.33724743

## Distribution statistics

| Quantity (unit) | Population / weighting | Mean | Std | P50 | P90 | P95 | P99 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| scalar_shear_stress (Pa) | all_exported_nodes / point_sampled | 0.669271 | 0.963986 | 0.402598 | 1.35627 | 1.99478 | 5.19322 | 18.0609 |
| speed (m/s) | all_exported_nodes / point_sampled | 0.851994 | 0.448894 | 0.904642 | 1.39602 | 1.41899 | 1.4319 | 1.59468 |
| transit_time (s) | all_exported_paths / equal_path | 0.133468 | 0.0516876 | 0.099756 | 0.217442 | 0.23125 | 0.247528 | 0.387053 |
| transit_time (s) | all_exported_unverified / equal_path | 0.133468 | 0.0516876 | 0.099756 | 0.217442 | 0.23125 | 0.247528 | 0.387053 |
| transit_time (s) | all_exported_unverified / inlet_speed_proxy | 0.133388 | 0.0514387 | 0.099756 | 0.217359 | 0.23127 | 0.245002 | 0.387053 |
| path_maximum_stress (Pa) | all_exported_unverified / equal_path | 2.56115 | 2.76434 | 1.32632 | 6.12377 | 8.60422 | 12.3706 | 18.0609 |
| path_maximum_stress (Pa) | all_exported_unverified / inlet_speed_proxy | 2.55623 | 2.7511 | 1.32632 | 6.11774 | 8.61833 | 12.3441 | 18.0609 |
| stress_accumulation (dyne*s/cm2) | all_exported_unverified / equal_path | 0.58944 | 0.47701 | 0.337247 | 1.24759 | 1.58379 | 2.30544 | 3.42269 |
| stress_accumulation (dyne*s/cm2) | all_exported_unverified / inlet_speed_proxy | 0.588484 | 0.473846 | 0.337247 | 1.24658 | 1.58626 | 2.28462 | 3.42269 |
| GW_HI2 (%) | all_exported_unverified / equal_path | 4.22623e-06 | 8.3623e-06 | 6.61443e-07 | 1.26162e-05 | 2.23635e-05 | 4.49891e-05 | 5.44403e-05 |
| GW_HI2 (%) | all_exported_unverified / inlet_speed_proxy | 4.20832e-06 | 8.30734e-06 | 6.61443e-07 | 1.2611e-05 | 2.24407e-05 | 4.47033e-05 | 5.44403e-05 |
| GW_HI3 (%) | all_exported_unverified / equal_path | 7.36767e-06 | 1.61006e-05 | 1.13989e-06 | 1.98368e-05 | 3.88666e-05 | 8.70527e-05 | 0.00013242 |
| GW_HI3 (%) | all_exported_unverified / inlet_speed_proxy | 7.32256e-06 | 1.59198e-05 | 1.13989e-06 | 1.97989e-05 | 3.88914e-05 | 8.69464e-05 | 0.00013242 |

- Speed percentiles describe velocity magnitude, not volumetric flow rate.
- Point-sampled stresses/speeds depend on export sampling; they are not volume- or residence-time-weighted.
- Weighted path statistics use the same supplied flux or inlet-speed proxy as the HI average.
- Unweighted percentiles use linear interpolation; weighted percentiles use the inverse empirical CDF.
- Standard deviation describes distribution spread, not uncertainty in a device prediction; points/paths are not independent experimental replicates.
- Weighted count is the number of positive-weight paths, not an effective sample size.

## Interpretation limits

- Viscous stress only; SST Reynolds stresses are not included in exported varshear.
- Steady snapshots do not reconstruct transient trajectories or full-cycle haemolysis.
- Mesh convergence, boundary-layer resolution, mass balance and physiological validation require solver data.
- Coefficient-range diagnostics are extrapolation flags, not evidence of biological validity.
- Weights are inlet SPEED proxies; equal represented area and normal flow are unverified.
- No outlet plane: all-export averages include potentially incomplete paths; not verified device HI.
- Expected 1000 seeds; reconstructed 969 paths. Missing release flux is unknown; count difference is not a recirculation measurement.