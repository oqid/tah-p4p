"""Regression checks for snapshot alignment and the shared empirical plot."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import haem5_snapshots as snapshots

class SnapshotTests(unittest.TestCase):
    def test_ended_paths_do_not_become_zero_stress(self):
        paths = [np.array([[0., 2.], [1., 4.]]),
                 np.array([[0., 10.], [2., 20.]])]
        grid = np.array([0., .5, 1., 1.5, 2., 3.])
        values = snapshots.resample_paths(paths, grid)
        self.assertTrue(np.isnan(values[0, 3:]).all())
        self.assertAlmostEqual(values[0, 1], 3.)
        stats = snapshots.population_stats(values, grid)
        self.assertEqual(stats.n_present.tolist(), [2, 2, 2, 1, 1, 0])
        self.assertAlmostEqual(stats.median_Pa.iloc[3], 17.5)
        self.assertTrue(np.isnan(stats.mean_Pa.iloc[-1]))
        self.assertAlmostEqual(stats.fraction_present.iloc[3], .5)

    def test_duplicate_times_and_normalized_alignment(self):
        path = np.array([[0., 1.], [0., 3.], [2., 7.]])
        values = snapshots.resample_paths([path], np.array([0., .5, 1.]), normalized=True)
        np.testing.assert_allclose(values, [[3., 5., 7.]])

    def test_time_reversal_is_rejected(self):
        with self.assertRaises(ValueError):
            snapshots.resample_paths([np.array([[0., 1.], [2., 2.], [1., 3.]])], np.array([0., 1.]))

    def test_filename_phase(self):
        self.assertEqual(snapshots.phase_from_name(Path("haemline_50cc_1667%.csv")), 16.67)
        self.assertEqual(snapshots.phase_from_name(Path("haemline_50cc_83.0%.csv")), 83.)
        with self.assertRaises(ValueError):
            snapshots.phase_from_name(Path("unknown.csv"))

    def test_density_and_empirical_render_with_zero_stress_and_different_durations(self):
        cases = []
        for i in range(2):
            cases.append(dict(stem=str(i), label=str(i), color=snapshots.COLORS[i], style="-",
                              paths=[np.array([[0., 0.], [.001*(i+1), 5.], [.002*(i+1), 1.]])],
                              envelope=pd.DataFrame(dict(duration_s=[.001, .003], stress_Pa=[10., 1.]))))
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            snapshots.plot_histories(cases, output, 25)
            snapshots.plot_empirical(cases, output)
            self.assertEqual(len(list(output.glob("*.png"))), 5)
            for path in output.glob("*.png"):
                self.assertGreater(path.stat().st_size, 1000)

if __name__ == "__main__":
    unittest.main()
