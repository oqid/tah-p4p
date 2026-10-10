"""Numerical and CLI checks for quasi-steady snapshot aggregation."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import quasi_steady_systole as qs


class QuasiSteadyTests(unittest.TestCase):
    def test_time_quadrature(self):
        times, durations = qs.snapshot_timing(5)
        np.testing.assert_allclose(times, [0, 62.5, 125, 187.5, 250])
        np.testing.assert_allclose(durations, [31.25, 62.5, 62.5, 62.5, 31.25])
        times, durations = qs.snapshot_timing(5, sampling="midpoints")
        np.testing.assert_allclose(times, [25, 75, 125, 175, 225])
        np.testing.assert_allclose(durations, [50] * 5)
        _, durations = qs.snapshot_timing(3, times_ms=[0, 100, 250])
        np.testing.assert_allclose(durations, [50, 125, 75])
        for times in ([0, 0, 250], [10, 100, 250], [0, 250]):
            with self.assertRaises(ValueError):
                qs.snapshot_timing(3, times_ms=times)

    def test_weighting_is_by_volume_not_seed_count(self):
        def paths(count, speed, value):
            return pd.DataFrame({**{key: [value] * count for key in qs.METRICS},
                                 "mean_inlet_velocity": [speed] * count})
        a, _ = qs.summarize_paths(paths(2, 1, 10), 250)
        b, _ = qs.summarize_paths(paths(20, 3, 30), 250)
        snapshots = pd.DataFrame([a, b]).assign(represented_ms=[125, 125])
        weights, combined = qs.combine_snapshots(snapshots, snapshots.mean_inlet_speed_m_s)
        np.testing.assert_allclose(weights, [0.25, 0.75])
        self.assertAlmostEqual(combined["HI3_percent"], 25)
        # Uniform scaling of Q leaves the average unchanged.
        _, scaled = qs.combine_snapshots(snapshots, [100, 300])
        self.assertAlmostEqual(scaled["HI3_percent"], 25)
        _, zero = qs.combine_snapshots(snapshots, [0, 3])
        self.assertAlmostEqual(zero["HI3_percent"], 30)
        with self.assertRaises(ValueError):
            qs.combine_snapshots(snapshots, [0, 0])
        with self.assertRaises(ValueError):
            qs.combine_snapshots(snapshots, [-1, 3])

    def test_zero_flow_endpoint(self):
        paths = pd.DataFrame({**{key: [10.] for key in qs.METRICS},
                              "mean_inlet_velocity": [0.]})
        zero, weights = qs.summarize_paths(paths, 250)
        self.assertEqual(weights.sum(), 0)
        paths["mean_inlet_velocity"] = 1.
        flowing, _ = qs.summarize_paths(paths, 250)
        snapshots = pd.DataFrame([zero, flowing]).assign(represented_ms=[125, 125])
        _, combined = qs.combine_snapshots(snapshots, [0, 1])
        self.assertAlmostEqual(combined["HI3_percent"], 10)
        with self.assertRaises(ValueError):
            qs.combine_snapshots(snapshots, [1, 1])

    def test_within_snapshot_flux(self):
        summary = pd.DataFrame({**{key: [10, 30] for key in qs.METRICS},
                                "mean_inlet_velocity": [1, 3]})
        row, weights = qs.summarize_paths(summary, 20)
        self.assertAlmostEqual(row["transit_ms"], 25)
        self.assertAlmostEqual(row["fraction_transit_longer_than_systole"], 0.75)
        np.testing.assert_allclose(weights, [0.25, 0.75])

    def test_cli_real_parser_and_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary) / "snapshots"
            folder.mkdir()
            # Constant shear and speed, one 0.1 m path: transit=100 ms.
            export = ("[Name]\ntest\n[Data]\n"
                      "Node Number, X [ m ], Y [ m ], Z [ m ], varshear [ Pa ], Velocity [ m s^-1 ], X [ m ], Y [ m ], Z [ m ]\n"
                      "0, 0, 0, 0, 10, 1, 0, 0, 0\n"
                      "1, 0.1, 0, 0, 10, 1, 0.1, 0, 0\n"
                      "[Lines]\n0, 1\n")
            for number in [10, 2, 1, 4, 3]:
                (folder / f"snapshot_{number}.csv").write_text(export)
            out = Path(temporary) / "results"
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                qs.main([str(folder), "--outdir", str(out), "--flow-rates-ml-s",
                         "0", "100", "100", "100", "0"])
            snapshots = pd.read_csv(out / "snapshot_summary.csv")
            total = pd.read_csv(out / "systole_summary.csv").iloc[0]
            pooled = pd.read_csv(out / "weighted_streamlines.csv")
            self.assertEqual(list(snapshots.file), [f"snapshot_{n}.csv" for n in [1, 2, 3, 4, 10]])
            self.assertAlmostEqual(total.transit_ms, 100)
            self.assertAlmostEqual(total.SA_dyne_s_cm2, 10)
            c = qs.haem5.POWER_LAW_CONSTANTS["GW"]
            expected_hi3 = c["A"] * 10 ** c["alpha"] * 0.1 ** c["beta"]
            self.assertAlmostEqual(total.HI3_percent, expected_hi3)
            self.assertAlmostEqual(total.represented_systolic_volume_ml, 18.75)
            self.assertAlmostEqual(pooled.systole_weight.sum(), 1)
            self.assertTrue((out / "systole_summary.png").stat().st_size > 1000)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                qs.main([str(folder), "--outdir", str(out)])


if __name__ == "__main__":
    unittest.main()
