# 1000p_100_alotpa_3djawn

Input: C:\Users\kingk\Documents\tah-p4p\1000p_100_alotpa_3djawn.csv
Population: all_exported_unverified; weighting: inlet_speed_proxy
Paths: 100; selected: 100; selected/exported weight: 1.000000

## Metrics

- GW_HI2_percent: 8.4978821e-05
- GW_HI3_percent: 0.00018326456
- GW_HI2_exact_percent: 8.5271547e-05
- GW_time_fraction_tau_above_fit: 0
- GW_transit_above_fit: 0
- HO_HI2_percent: 1.334049e-06
- HO_HI3_percent: 2.6773605e-06
- HO_HI2_exact_percent: 1.3452316e-06
- HO_time_fraction_tau_above_fit: 0
- HO_transit_above_fit: 0
- TZ_HI2_percent: 1.0812829e-05
- TZ_HI3_percent: 3.1189159e-05
- TZ_HI2_exact_percent: 1.1014839e-05
- TZ_time_fraction_tau_above_fit: 0
- TZ_transit_above_fit: 0
- TZ_time_fraction_tau_below_50Pa: 0.99990151
- GW_HI3_no_floor_percent: 0.00018326456
- GW_HI3_midpoint_percent: 0.00018070058
- GW_HI2_reversed_percent: 0.00012674628
- GW_HI3_stride2_percent: 0.00018417285
- GW_HI3_stride4_percent: 0.00018693477
- GW_HI2_NIH_mg_per_100L: 8.9227762
- GW_HI3_NIH_mg_per_100L: 19.242779
- HO_HI2_NIH_mg_per_100L: 0.14007514
- HO_HI3_NIH_mg_per_100L: 0.28112286
- TZ_HI2_NIH_mg_per_100L: 1.135347
- TZ_HI3_NIH_mg_per_100L: 3.2748617
- GW_HI3_all_exported_percent: 0.00018326456
- SA_selected_median_dyne_s_cm2: 0.67468414

## Distribution statistics

| Quantity (unit) | Population / weighting | Mean | Std | P50 | P90 | P95 | P99 | Max |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| scalar_shear_stress (Pa) | all_exported_nodes / point_sampled | 2.52449 | 4.71049 | 1.16179 | 5.94983 | 9.10862 | 24.1222 | 75.6325 |
| speed (m/s) | all_exported_nodes / point_sampled | 1.76874 | 0.857869 | 1.94463 | 2.79235 | 2.97249 | 3.14335 | 3.37991 |
| transit_time (s) | all_exported_paths / equal_path | 0.0659774 | 0.0227536 | 0.0598051 | 0.096018 | 0.108472 | 0.121141 | 0.137197 |
| transit_time (s) | all_exported_unverified / equal_path | 0.0659774 | 0.0227536 | 0.0598051 | 0.096018 | 0.108472 | 0.121141 | 0.137197 |
| transit_time (s) | all_exported_unverified / inlet_speed_proxy | 0.0659822 | 0.0227561 | 0.0598306 | 0.0967688 | 0.113231 | 0.137197 | 0.137197 |
| path_maximum_stress (Pa) | all_exported_unverified / equal_path | 7.85324 | 12.0577 | 3.72398 | 17.0075 | 28.462 | 55.8962 | 75.6325 |
| path_maximum_stress (Pa) | all_exported_unverified / inlet_speed_proxy | 7.85653 | 12.0624 | 3.74803 | 17.0165 | 35.4446 | 75.6325 | 75.6325 |
| stress_accumulation (dyne*s/cm2) | all_exported_unverified / equal_path | 0.856355 | 0.674205 | 0.674684 | 1.50692 | 2.19271 | 3.57961 | 4.22865 |
| stress_accumulation (dyne*s/cm2) | all_exported_unverified / inlet_speed_proxy | 0.856608 | 0.674463 | 0.67796 | 1.54893 | 2.2055 | 4.22865 | 4.22865 |
| GW_HI2 (%) | all_exported_unverified / equal_path | 8.49083e-05 | 0.000314925 | 7.64554e-06 | 0.000103597 | 0.000302333 | 0.00178498 | 0.00252837 |
| GW_HI2 (%) | all_exported_unverified / inlet_speed_proxy | 8.49788e-05 | 0.000315061 | 7.80332e-06 | 0.000127558 | 0.000379495 | 0.00252837 | 0.00252837 |
| GW_HI3 (%) | all_exported_unverified / equal_path | 0.000183098 | 0.000737932 | 1.19903e-05 | 0.000223209 | 0.000606678 | 0.00488694 | 0.00540425 |
| GW_HI3 (%) | all_exported_unverified / inlet_speed_proxy | 0.000183265 | 0.000738322 | 1.22579e-05 | 0.000317344 | 0.000734973 | 0.00540425 | 0.00540425 |

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