# Snapshot comparison

Frozen-flow snapshots; exposure time is time along a path, not time within systole.
All reconstructed paths; equal-path distribution statistics; bands show population spread.

Seven PNGs compare SA distributions/tails, empirical envelopes (full and detail), shear density, all raw path histories, mean/median and population support, normalized transit medians, and phase trends. CSVs retain the numerical summaries.

Phases inferred from filename suffixes (1667% = 16.67%); override with --phases. No full-cycle integration or transient particle tracking is performed.

HI uses the unchanged haem5 model and inlet-speed weighting (equal-area inlet seeding assumption). All reconstructed paths are included; outlet completion is not checked. Shear interpolation does not change the damage calculations. Ended paths are excluded. Percentiles describe the population, not confidence intervals.

SA density is per log10 unit; zero-SA fraction is in snapshot_summary.csv. Density panels show fractions of the ORIGINAL path population per shared log shear bin, so colour intensity also falls as paths end. Zero shear is placed at the lowest log bin.

Empirical envelopes use approximate NIH lysis references, not activation thresholds or a validated variable-stress damage model. --empirical-dir trusts the supplied saved results; run.json records their exact filenames and hashes. Omit it to recalculate.
