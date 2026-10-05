"""Checks for contiguous exposure screening and primary-source curve units."""
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from haem5_empiricalthresholds import (
    _longest_contiguous_duration_above, _log_interp_threshold,
    RBC_THRESHOLD, PLATELET_THRESHOLD, criterion_utilization,
    analyze_empirical_thresholds, plot_empirical_thresholds,
)


class ExposureTests(unittest.TestCase):
    def test_separated_exposures_are_not_added(self):
        stress = np.array([0., 50., 50., 5., 50.])
        dt = np.array([0., .002, .003, .001, .004])
        self.assertAlmostEqual(_longest_contiguous_duration_above(stress, dt, 40), .005)
        self.assertAlmostEqual(_longest_contiguous_duration_above(stress, dt, 5), .010)

    def test_source_units_and_no_out_of_range_extrapolation(self):
        self.assertGreater(RBC_THRESHOLD["stress_Pa"].iloc[-1], 100)
        self.assertLess(RBC_THRESHOLD["stress_Pa"].iloc[-1], 200)
        self.assertLess(PLATELET_THRESHOLD["stress_Pa"].iloc[-1], 20)
        self.assertTrue(np.isnan(_log_interp_threshold(1e-9, RBC_THRESHOLD)))
        env = pd.DataFrame({"duration_s": [1e-9], "stress_Pa": [100.]})
        self.assertTrue(np.isnan(criterion_utilization(env, RBC_THRESHOLD)[0]))

    def test_device_frontier_can_combine_different_paths(self):
        def path(stress, duration):
            return pd.DataFrame({"shear": [stress, stress], "dt": [0., duration],
                                 "vel": [1., 1.], "t_cumulative": [0., duration]})
        _, device, envelope = analyze_empirical_thresholds(
            {"streamlines": [path(50., .002), path(10., .1)]})
        high = envelope[envelope["stress_Pa"] > 10]
        low = envelope[envelope["stress_Pa"] <= 10]
        np.testing.assert_allclose(high["duration_s"], .002)
        np.testing.assert_allclose(low["duration_s"], .1)
        self.assertTrue((device["streamlines_unassessed"] == 0).all())
        fig = plot_empirical_thresholds(envelope)
        fig.canvas.draw()
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual(fig.axes[0].get_xlabel(), "Exposure duration [ms]")
        np.testing.assert_allclose(fig.axes[0].lines[2].get_xdata(),
                                   envelope.sort_values(["duration_s", "stress_Pa"],
                                                        ascending=[True, False])["duration_s"] * 1000)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
