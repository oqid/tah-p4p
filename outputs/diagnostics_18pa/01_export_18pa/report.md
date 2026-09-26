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

## Interpretation limits

- Viscous stress only; SST Reynolds stresses are not included in exported varshear.
- Steady snapshots do not reconstruct transient trajectories or full-cycle haemolysis.
- Mesh convergence, boundary-layer resolution, mass balance and physiological validation require solver data.
- Coefficient-range diagnostics are extrapolation flags, not evidence of biological validity.
- Weights are inlet SPEED proxies; equal represented area and normal flow are unverified.
- No outlet plane: all-export averages include potentially incomplete paths; not verified device HI.
- Expected 1000 seeds; reconstructed 969 paths. Missing release flux is unknown; count difference is not a recirculation measurement.