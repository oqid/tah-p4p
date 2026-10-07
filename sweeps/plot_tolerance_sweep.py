"""Plot haemolysis metrics against grid-relative tolerance."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter

SWEEP_DIR = Path(__file__).resolve().parent
TOLERANCE_ORDER = [0.02, 0.01, 0.005, 0.0025]
TOLERANCE_LABELS = {0.02: "2T\n0.02", 0.01: "T\n0.01",
                    0.005: "T/2\n0.005", 0.0025: "T/4\n0.0025"}
LINESTYLES = {"HI2": "--", "HI3": ":"}
MARKERS = {"GW": "o", "HO": "s"}
# Muted orange and pink from the lesbian pride flag palette.
COLOURS = {"GW": "#D98255", "HO": "#C56A91"}
PATHLINES_PLOT_YLIM = (3.1744528705483715e-7, 3.923862981256825e-4)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path, nargs="?",
                    default=SWEEP_DIR / "tolerance_sweep_results.csv")
    ap.add_argument("--metric", choices=["nih", "hi"], default="hi",
                    help="HI percent (default) or NIH mg/100L")
    ap.add_argument("--out", type=Path, default=SWEEP_DIR / "hi_vs_tolerance.png")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if not args.show:
        plt.switch_backend("Agg")
    if not args.csv.is_file():
        ap.error(f"Results file does not exist: {args.csv}")

    data = pd.read_csv(args.csv)
    data = data[data["constants"].isin(["GW", "HO"])]
    if data.empty:
        ap.error("No GW or HO results found in the selected file")
    observed = sorted(data["tolerance_grid_relative"].unique(), reverse=True)
    if len(observed) != 4:
        ap.error(f"Expected four tolerance values; found {observed}")

    suffix = "_NIH_mg_per_100L" if args.metric == "nih" else "_percent"
    ylabel = "NIH (mg/100L)" if args.metric == "nih" else "HI (%)"
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for constant in ["GW", "HO"]:
        subset = data[data["constants"] == constant].sort_values(
            "tolerance_grid_relative", ascending=False
        )
        if subset.empty:
            continue
        for metric in ["HI2", "HI3"]:
            ax.plot(subset["tolerance_grid_relative"], subset[f"{metric}{suffix}"],
                    linestyle=LINESTYLES[metric], color=COLOURS[constant],
                    marker=MARKERS[constant], markerfacecolor="none", markersize=6,
                    linewidth=1.5, label=f"{metric}-{constant}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(PATHLINES_PLOT_YLIM)
    ax.xaxis.set_major_locator(FixedLocator(TOLERANCE_ORDER))
    ax.xaxis.set_major_formatter(
        FixedFormatter([TOLERANCE_LABELS[value] for value in TOLERANCE_ORDER])
    )
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.022, 0.0022)
    ax.set_xlabel("grid-relative integration tolerance", fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.grid(True, which="major", linestyle=":", alpha=0.5)
    ax.legend(fontsize=9, ncol=2, framealpha=1, edgecolor="black",
              loc="upper center", bbox_to_anchor=(0.5, -0.20))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Saved {args.out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
