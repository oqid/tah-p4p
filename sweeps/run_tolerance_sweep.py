"""Process fixed-pathline exports generated at different grid-relative tolerances.

The sweep keeps the requested pathline count fixed while varying the CFD
integration tolerance. It processes existing exports; it does not rerun CFD.
"""

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd

SWEEP_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = SWEEP_DIR / "tolerance_sweep"
sys.path.insert(0, str(SWEEP_DIR.parent))
import haem5


def metadata_from_name(path):
    match = re.fullmatch(
        r"step_sweep_N(?P<pathlines>\d+)_T(?P<step>T|div\d+)?_GT(?P<tolerance>[\d.]+)",
        path.stem,
    )
    if not match:
        return None
    groups = match.groupdict()
    tolerance = float(groups["tolerance"])
    # The base tolerance T is 0.01 grid-relative. The plain T filename also
    # labels the 2T export, so derive labels from the actual exported value.
    ratio = tolerance / 0.01
    labels = {2.0: "2T", 1.0: "T", 0.5: "T/2", 0.25: "T/4"}
    label = labels.get(round(ratio, 8), f"{ratio:g}T")
    return {
        "pathlines": int(groups["pathlines"]),
        "tolerance_label": label,
        "tolerance_grid_relative": tolerance,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Process haemolysis over a fixed-pathline tolerance sweep."
    )
    ap.add_argument("folder", type=Path, nargs="?", default=DEFAULT_INPUT_DIR,
                    help="folder containing tolerance exports")
    ap.add_argument("--pathlines", type=int, default=5000,
                    help="requested pathline count (default: 5000)")
    ap.add_argument("--constants", nargs="+", default=["GW", "HO"],
                    choices=list(haem5.POWER_LAW_CONSTANTS),
                    help="constant sets to run (default: GW HO; TZ excluded)")
    ap.add_argument("--out", type=Path, default=SWEEP_DIR / "tolerance_sweep_results.csv")
    ap.add_argument("--resume", action="store_true",
                    help="skip tolerance and constant combinations already in --out")
    ap.add_argument("--dry-run", action="store_true",
                    help="list inputs and pending constant sets without running")
    args = ap.parse_args()
    if args.pathlines <= 0:
        ap.error("--pathlines must be positive")
    if not args.folder.is_dir():
        ap.error(f"Input folder does not exist: {args.folder}")

    files = []
    for path in args.folder.glob("step_sweep_N*_T*_GT*.csv"):
        metadata = metadata_from_name(path)
        if metadata and metadata["pathlines"] == args.pathlines:
            files.append((metadata["tolerance_grid_relative"], path, metadata))
    files.sort(key=lambda item: item[0], reverse=True)
    if not files:
        sys.exit(f"No tolerance exports for N={args.pathlines} in {args.folder}")

    rows = []
    done = set()
    if args.resume and args.out.exists():
        previous = pd.read_csv(args.out)
        rows = previous.to_dict("records")
        done = set(zip(previous["tolerance_grid_relative"], previous["constants"]))

    if args.dry_run:
        print(f"Inputs: {args.folder.resolve()}")
        print(f"Results: {args.out.resolve()}")
        for tolerance, path, metadata in files:
            todo = [c for c in args.constants if (tolerance, c) not in done]
            print(f"N={args.pathlines}, {metadata['tolerance_label']} (grid-relative {tolerance:g}): "
                  f"{path.name} -> {', '.join(todo) if todo else 'already done'}")
        return

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for tolerance, path, metadata in files:
        todo = [c for c in args.constants if (tolerance, c) not in done]
        if not todo:
            print(f"[grid-relative tolerance={tolerance:g}] already done, skipping")
            continue

        print(f"\n##### N={args.pathlines}, {metadata['tolerance_label']} "
              f"(grid-relative tolerance={tolerance:g}) ({path.name}) #####")
        started = time.time()
        try:
            all_results, comparison = haem5.run_pipeline_all_constants(path, tuple(todo))
        except Exception as exc:
            print(f"  FAILED on {path.name}: {exc}")
            continue

        for name in todo:
            result = all_results[name]
            rows.append({
                "N_requested": args.pathlines,
                "n_streamlines": len(result["summary"]),
                "tolerance_label": metadata["tolerance_label"],
                "tolerance_grid_relative": tolerance,
                "constants": name,
                "HI2_percent": result["device_HI2_percent"],
                "HI3_percent": result["device_HI3_percent"],
                "HI2_NIH_mg_per_100L": result["device_HI2_NIH_mg_per_100L"],
                "HI3_NIH_mg_per_100L": result["device_HI3_NIH_mg_per_100L"],
                "SA_median_dyne_s_cm2": result["device_SA_median_dyne_s_cm2"],
                "P_SA_above_hellums": result["device_SA_prob_above_hellums"],
                "frac_tau_over_range": comparison.loc[name, "frac_segments_tau_over_range"],
                "frac_t_over_range": comparison.loc[name, "frac_segments_t_over_range"],
                "file": path.name,
            })
        del all_results, comparison
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"  done in {time.time() - started:.1f} s, saved {args.out}")

    if not rows:
        sys.exit("No results were produced; see the failures above.")
    output = pd.DataFrame(rows).sort_values(["constants", "tolerance_grid_relative"],
                                            ascending=[True, False])
    print("\n" + "=" * 80)
    print(output[["N_requested", "n_streamlines", "tolerance_label",
                  "tolerance_grid_relative", "constants", "HI2_percent", "HI3_percent",
                  "HI2_NIH_mg_per_100L", "HI3_NIH_mg_per_100L"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("=" * 80)


if __name__ == "__main__":
    main()
