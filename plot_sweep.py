"""
Plot HI (%) vs number of path lines from run_sweep.py's sweep_results.csv,
in the style of the particle-count convergence figure (log-log, one curve per
model x constant set).

Usage:
    python plot_sweep.py                                  # HI2 + HI3, GW + HO + TZ
    python plot_sweep.py --hi HI2 HI3 --constants GW HO   # pick models / constants
    python plot_sweep.py --metric nih                     # plot NIH (mg/100L) instead of HI%
    python plot_sweep.py --x requested                    # x = requested N, not actual streamline count
    python plot_sweep.py --colour --show                  # colour by constant set, open a window
    python plot_sweep.py --linear-y                       # linear y axis

Also prints, per curve, the % change of each level relative to the finest
level and relative to the previous level, so you can state the N at which
the result is converged to your tolerance.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

LINESTYLES = {"HI2": "--", "HI3": ":"}          # paper: HI1 dotted, HI2 dashed
MARKERS = {"GW": "o", "HO": "s", "TZ": "*"}
COLOURS = {"GW": "tab:blue", "HO": "tab:orange", "TZ": "tab:green"}


def main():
    ap = argparse.ArgumentParser(description="Plot HI vs number of path lines.")
    ap.add_argument("csv", type=Path, nargs="?", default=Path("sweep_results.csv"))
    ap.add_argument("--hi", nargs="+", default=["HI2", "HI3"], choices=["HI2", "HI3"])
    ap.add_argument("--constants", nargs="+", default=["GW", "HO", "TZ"],
                    choices=["GW", "HO", "TZ"])
    ap.add_argument("--metric", choices=["hi", "nih"], default="hi",
                    help="hi = HI%% (default), nih = mg/100L")
    ap.add_argument("--x", choices=["streamlines", "requested"], default="streamlines",
                    help="x axis: actual streamlines exported (default) or requested N")
    ap.add_argument("--colour", action="store_true", help="colour by constant set")
    ap.add_argument("--linear-y", action="store_true")
    ap.add_argument("--out", type=Path, default=Path("hi_vs_pathlines.png"))
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    xcol = "n_streamlines" if args.x == "streamlines" else "N_requested"
    suffix = "_percent" if args.metric == "hi" else "_NIH_mg_per_100L"
    ylabel = "HI (%)" if args.metric == "hi" else "NIH (mg/100L)"

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for c in args.constants:
        sub = df[df["constants"] == c].sort_values(xcol)
        if sub.empty:
            print(f"(no rows for {c}, skipping)")
            continue
        for model in args.hi:
            y = sub[f"{model}{suffix}"]
            colour = COLOURS[c] if args.colour else "black"
            ax.plot(sub[xcol], y, linestyle=LINESTYLES[model], color=colour,
                    marker=MARKERS[c], markersize=6, markerfacecolor="none",
                    linewidth=1.3, label=f"{model}-{c}")

            ref = y.iloc[-1]
            chg_ref = (y / ref - 1) * 100
            chg_prev = y.pct_change() * 100
            print(f"\n{model}-{c}  (reference = finest level, {sub[xcol].iloc[-1]:g} lines)")
            print(pd.DataFrame({"lines": sub[xcol].values, "value": y.values,
                                "% vs finest": chg_ref.values,
                                "% vs previous": chg_prev.values})
                  .to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    ax.set_xscale("log")
    if not args.linear_y:
        ax.set_yscale("log")
    ax.set_xlabel("number of path lines", fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.grid(True, which="major", linestyle=":", alpha=0.5)
    ax.legend(fontsize=9, ncol=3, framealpha=1, edgecolor="black",
              loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"\nSaved {args.out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
