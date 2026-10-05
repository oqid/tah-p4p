"""CLI regressions for the relocated sweep scripts."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

SWEEPS = Path(__file__).resolve().parent
ROOT = SWEEPS.parent


class SweepPathTests(unittest.TestCase):
    def run_cli(self, script, args, cwd):
        result = subprocess.run([sys.executable, str(SWEEPS / script), *args],
                                cwd=cwd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def test_default_discovery_and_resume_from_both_folders(self):
        for cwd in (ROOT, SWEEPS):
            output = self.run_cli("run_sweep.py", ["--resume", "--constants", "GW", "HO",
                                                  "--dry-run"], cwd)
            self.assertIn(str(SWEEPS / "sweep_results.csv"), output)
            self.assertIn("N=100000: sweep_N100000.csv", output)
            self.assertIn("N=10: sweep_N10.csv -> already done", output)

    def test_small_sweep_resume_and_plot_from_unrelated_folder(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            csv = folder / "results.csv"
            args = ["--min-n", "10", "--max-n", "10", "--constants", "GW", "HO",
                    "--out", str(csv)]
            self.run_cli("run_sweep.py", args, folder)
            before = csv.read_bytes()
            output = self.run_cli("run_sweep.py", [*args, "--resume"], folder)
            self.assertIn("already done, skipping", output)
            self.assertEqual(csv.read_bytes(), before)
            rows = pd.read_csv(csv)
            self.assertEqual(set(rows["constants"]), {"GW", "HO"})
            self.assertEqual(len(rows), 2)
            png = folder / "plots" / "sweep.png"
            self.run_cli("plot_sweep.py", [str(csv), "--constants", "GW", "HO",
                                           "--out", str(png)], folder)
            self.assertGreater(png.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
