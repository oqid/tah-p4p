"""
Run haem5.py's hemolysis pipeline over every sweep_N<count>.csv export and
collect HI2 / HI3 (and NIH) for each constant set into one results CSV.

It imports haem5 directly (no subprocess), so haem5's __main__ block never
runs: no plots and no output files are produced by haem5 itself, and you do
not need to touch GENERATE_OUTPUTS. The per-file text table from
haem5.run_pipeline_all_constants is still printed.

Usage:
    python sweeps/run_sweep.py --resume --constants GW HO
    python sweeps/run_sweep.py --dry-run --resume          # list pending inputs
    python sweeps/run_sweep.py <folder> --max-n 3000        # custom input folder

Inputs and default output live beside this script, independent of working directory.
Output: sweeps/sweep_results.csv (one row per file per constant set).
"""

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd

SWEEP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SWEEP_DIR))
import haem5


def n_from_name(path):
    m = re.search(r"N(\d+)", path.stem)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description="Run haem5 over a particle-count sweep.")
    ap.add_argument("folder", type=Path, nargs="?", default=SWEEP_DIR,
                    help="folder containing the sweep CSVs (default: beside this script)")
    ap.add_argument("--pattern", default="sweep_N*.csv", help="glob (default: sweep_N*.csv)")
    ap.add_argument("--constants", nargs="+", default=["GW", "HO", "TZ"],
                    choices=list(haem5.POWER_LAW_CONSTANTS), help="constant sets to run")
    ap.add_argument("--out", type=Path, default=SWEEP_DIR / "sweep_results.csv")
    ap.add_argument("--min-n", type=int, default=None)
    ap.add_argument("--max-n", type=int, default=None)
    ap.add_argument("--resume", action="store_true",
                    help="skip (N, constants) combinations already in --out")
    ap.add_argument("--dry-run", action="store_true",
                    help="list inputs and pending constant sets without running or writing")
    args = ap.parse_args()
    if not args.folder.is_dir():
        ap.error(f"Input folder does not exist: {args.folder}")

    files = []
    for p in args.folder.glob(args.pattern):
        n = n_from_name(p)
        if n is None:
            continue
        if args.min_n is not None and n < args.min_n:
            continue
        if args.max_n is not None and n > args.max_n:
            continue
        files.append((n, p))
    files.sort()  # smallest first, so a crash on the big ones keeps the small ones
    if not files:
        sys.exit(f"No files matching {args.pattern} in {args.folder}")

    rows = []
    done = set()
    if args.resume and args.out.exists():
        prev = pd.read_csv(args.out)
        rows = prev.to_dict("records")
        done = set(zip(prev["N_requested"], prev["constants"]))

    if args.dry_run:
        print(f"Inputs: {args.folder.resolve()}")
        print(f"Results: {args.out.resolve()}")
        for n, path in files:
            todo = [c for c in args.constants if (n, c) not in done]
            print(f"N={n}: {path.name} -> {', '.join(todo) if todo else 'already done'}")
        return

    args.out.parent.mkdir(parents=True, exist_ok=True)

    for n, path in files:
        todo = [c for c in args.constants if (n, c) not in done]
        if not todo:
            print(f"[N={n}] already done, skipping")
            continue

        print(f"\n##### N = {n}   ({path.name}) #####")
        t0 = time.time()
        try:
            all_results, comparison = haem5.run_pipeline_all_constants(path, tuple(todo))
        except Exception as exc:  # keep going if one file is bad
            print(f"  FAILED on {path.name}: {exc}")
            continue

        for name in todo:
            res = all_results[name]
            rows.append({
                "N_requested": n,
                "n_streamlines": len(res["summary"]),
                "constants": name,
                "HI2_percent": res["device_HI2_percent"],
                "HI3_percent": res["device_HI3_percent"],
                "HI2_NIH_mg_per_100L": res["device_HI2_NIH_mg_per_100L"],
                "HI3_NIH_mg_per_100L": res["device_HI3_NIH_mg_per_100L"],
                "SA_median_dyne_s_cm2": res["device_SA_median_dyne_s_cm2"],
                "P_SA_above_hellums": res["device_SA_prob_above_hellums"],
                "frac_tau_over_range": comparison.loc[name, "frac_segments_tau_over_range"],
                "frac_t_over_range": comparison.loc[name, "frac_segments_t_over_range"],
                "file": path.name,
            })
        del all_results, comparison  # free memory before the next (bigger) file

        # save after every file so nothing is lost if a later file crashes
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"  done in {time.time() - t0:.1f} s, saved {args.out}")

    if not rows:
        sys.exit("No results were produced; see the failures above.")
    out = pd.DataFrame(rows).sort_values(["constants", "N_requested"])
    print("\n" + "=" * 70)
    print(out[["N_requested", "n_streamlines", "constants", "HI2_percent",
               "HI3_percent", "HI2_NIH_mg_per_100L", "HI3_NIH_mg_per_100L"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("=" * 70)


if __name__ == "__main__":
    main()
