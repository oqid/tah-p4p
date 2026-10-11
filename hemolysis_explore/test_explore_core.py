"""Check explore_core reproduces haem5 per path, on one export.
Run:  python test_explore_core.py ..\\z_inputs\\haemline_50cc_5000%.csv
"""
import sys
from pathlib import Path

import explore_core as ec

ec.check_constants()
f = Path(sys.argv[1]) if len(sys.argv) > 1 else sorted((Path(__file__).parent.parent / "z_inputs").glob("*.csv"))[0]
r = ec.verify_against_haem5(f)
print(r)
assert all(v < 1e-9 for k, v in r.items() if k.startswith("max_rel")), "mismatch with haem5"
print("OK: per-path HI2/HI3/SA/transit identical to haem5")
