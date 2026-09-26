"""Focused checks of the diagnostic population/weighting logic and 3D parsing."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import haem_diagnostics as audit


class DiagnosticsTests(unittest.TestCase):
    def export(self, root, three_d=False):
        path = root / "fixture.csv"
        z = .002 if three_d else 0
        path.write_text(
            "[Name]\nfixture\n[Data]\n"
            "Node Number, X [ m ], Y [ m ], Z [ m ], varshear [ Pa ], Velocity [ m s^-1 ]\n"
            f"0,0,0,0,2,1\n1,.01,0,{z},2,1\n"
            "2,0,.001,0,8,1\n3,.005,.001,0,8,1\n"
            "[Lines]\n0,1\n2,3\n", encoding="utf-8")
        return path

    def test_analytical_and_haem5(self):
        audit.self_test()

    def test_completion_and_flux_weights_in_2d_and_3d(self):
        for three_d in (False, True):
            with self.subTest(three_d=three_d), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                path = self.export(root, three_d)
                weights = root / "weights.csv"
                weights.write_text("streamline_id,flux_weight\n0,1\n1,3\n")
                config = {"weights": str(weights), "outlet_plane": [.01,0,0,1,0,0],
                          "outlet_tolerance": 1e-6, "expected_seeds": 3}
                result = audit.analyze(path, config, root / "out")
                self.assertEqual(result["coverage_of_exported_weight"], .25)
                c = audit.POWER_LAW_CONSTANTS["GW"]
                time = np.hypot(.01, .002 if three_d else 0)
                expected = c["A"] * 2 ** c["alpha"] * time ** c["beta"]
                self.assertAlmostEqual(result["GW_HI3_percent"], expected)
                report = json.loads((root / "out/diagnostics.json").read_text())
                self.assertEqual(report["missing_seed_count"], 1)

    def test_no_completed_paths_are_undefined(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = self.export(root)
            result = audit.analyze(path, {"outlet_plane": [1,0,0,1,0,0]}, root / "out")
            self.assertTrue(np.isnan(result["GW_HI3_percent"]))
            report = json.loads((root / "out/diagnostics.json").read_text())
            self.assertIsNone(report["metrics"]["GW_HI3_percent"])

    def test_bad_connectivity_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.export(Path(temp))
            path.write_text(path.read_text().replace("2,3\n", "2,99\n"))
            with self.assertRaisesRegex(ValueError, "absent node"):
                audit.read_export(path)


if __name__ == "__main__":
    unittest.main()
