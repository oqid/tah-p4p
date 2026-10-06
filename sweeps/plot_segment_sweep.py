"""Plot haemolysis results against requested segment count."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

LINESTYLES = {"HI2": "--", "HI3": ":"}
MARKERS = {"GW": "o", "HO": "s", "TZ": "*"}
COLOURS = {"GW": "tab:blue", "HO": "tab:orange", "TZ": "tab:green"}
SWEEP_DIR = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path, help="segment_sweep_results.csv")
    ap.add_argument("--pathlines", type=int, default=None,
                    help="requested pathline count to plot (default: all in file)")
    ap.add_argument("--hi", nargs="+", default=["HI2", "HI3"], choices=["HI2", "HI3"])
    ap.add_argument("--constants", nargs="+", default=["GW", "HO", "TZ"],
                    choices=["GW", "HO", "TZ"])
    ap.add_argument("--metric", choices=["hi", "nih"], default="hi")
    ap.add_argument("--out", type=Path, default=SWEEP_DIR / "hi_vs_segments.png")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if not args.show:
        plt.switch_backend("Agg")
    if not args.csv.is_file():
        ap.error(f"Results file does not exist: {args.csv}")
    data = pd.read_csv(args.csv)
    if args.pathlines is not None:
        data = data[data["N_requested"] == args.pathlines]
    if data.empty:
        ap.error("No results match the selected pathline count")
    suffix = "_percent" if args.metric == "hi" else "_NIH_mg_per_100L"
    ylabel = "HI (%)" if args.metric == "hi" else "NIH (mg/100L)"
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for constant in args.constants:
        subset = data[data["constants"] == constant].sort_values("segments_requested")
        if subset.empty:
            print(f"(no rows for {constant}, skipping)")
            continue
        for model in args.hi:
            x = subset["segments_requested"]
            y = subset[f"{model}{suffix}"]
            ax.plot(x, y, linestyle=LINESTYLES[model], color="black",
                    marker=MARKERS[constant], markerfacecolor="none",
                    markersize=6, linewidth=1.3, label=f"{model}-{constant}")
            print(f"\n{model}-{constant}")
            print(pd.DataFrame({"segments": x, "value": y}).to_string(
                index=False, float_format=lambda v: f"{v:.4g}"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of segments per pathline", fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.grid(True, which="major", linestyle=":", alpha=0.5)
    ax.legend(fontsize=9, ncol=3, framealpha=1, edgecolor="black",
              loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"\nSaved {args.out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
