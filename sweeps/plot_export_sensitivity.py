"""Compact three-panel summary of the streamline export sensitivity studies
(body figure replacing the old Table 6), in the style of Taskin et al. 2012,
Fig. 3C-D: mean HI3 normalised by the reference setting of each study.

    python sweeps/plot_export_sensitivity.py

Reads sweep_results.csv + a1to10k_sweep_results.csv + 300k_sweep_results.csv
(streamline count), segment_sweeps/segment_sweep_results.csv and
tolerance_sweep_results.csv. Writes export_sensitivity.pdf/.png here and to
hemolysis_explore/overleaf_upload/figures/.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hemolysis_explore"))
import report_style as rs  # noqa: E402

COLOURS = {"GW": "#D98255", "HO": "#C56A91"}   # as hi_vs_pathlines / hi_vs_tolerance
MARKERS = {"GW": "o", "HO": "s"}
LS = {"GW": "-", "HO": "--"}
BAND = 0.05                                     # acceptance criterion
UPLOAD = HERE.parent / "hemolysis_explore" / "overleaf_upload" / "figures"


def streamline_data():
    parts = [pd.read_csv(HERE / f) for f in
             ("sweep_results.csv", "a1to10k_sweep_results.csv", "300k_sweep_results.csv")]
    d = pd.concat(parts).drop_duplicates(["N_requested", "constants"], keep="first")
    return d.rename(columns={"N_requested": "x"}), 300_000, 5_000


def segment_data():
    d = pd.read_csv(HERE / "segment_sweeps" / "segment_sweep_results.csv")
    return d.rename(columns={"segments_requested": "x"}), 1_000, 1_000


def tolerance_data():
    d = pd.read_csv(HERE / "tolerance_sweep_results.csv")
    return d.rename(columns={"tolerance_grid_relative": "x"}), 0.01, 0.01


def thousands(v, _):
    return f"{v:,.0f}" if v >= 1 else f"{v:g}"


def panel(ax, data, ref, used, metric="HI3"):
    for c in ("GW", "HO"):
        s = data[data.constants == c].sort_values("x")
        r = s.loc[s.x == ref, f"{metric}_percent"].iloc[0]
        y = s[f"{metric}_percent"] / r
        ax.plot(s.x, y, color=COLOURS[c], ls=LS[c], marker=MARKERS[c], ms=3.2,
                mfc="white", mew=1.0, lw=1.2, label=f"{c}", zorder=3)
        print(f"  {metric}-{c}: " + ", ".join(
            f"{x:g}:{(v - 1) * 100:+.2f}%" for x, v in zip(s.x, y)))
    ax.axhspan(1 - BAND, 1 + BAND, color="#E6E6E6", alpha=0.8, lw=0, zorder=0)
    ax.axhline(1, color=rs.MUTED, lw=0.6, zorder=1)
    ax.axvline(used, color=rs.INK, lw=0.8, ls=":", zorder=2.5)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.2)


def main():
    rs.apply()
    fig, axes = plt.subplots(1, 3, figsize=(rs.WIDTH, 2.2), sharey=True,
                             gridspec_kw=dict(wspace=0.12))
    specs = [
        (streamline_data, "Number of streamlines", "(a)"),
        (segment_data, "Maximum segments per streamline", "(b)"),
        (tolerance_data, "Tracking tolerance (grid-relative)", "(c)"),
    ]
    for ax, (fn, xlabel, tag) in zip(axes, specs):
        data, ref, used = fn()
        print(tag, xlabel, f"(ref {ref:g}, used {used:g})")
        panel(ax, data, ref, used)
        ax.set_xlabel(xlabel)
        ax.set_title(tag, loc="left")
        ax.xaxis.set_major_formatter(FuncFormatter(thousands))

    # (a) label reference and production points
    a = axes[0]
    a.set_xlim(6, 6e5)
    a.xaxis.set_major_locator(FixedLocator([10, 100, 1e3, 1e4, 1e5]))
    a.xaxis.set_minor_locator(NullLocator())
    a.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: {10: "10", 100: "100", 1e3: "1k", 1e4: "10k", 1e5: "100k"}.get(v, "")))

    b = axes[1]
    b.set_xlim(18, 7000)
    b.xaxis.set_major_locator(FixedLocator([25, 100, 1000, 5000]))
    b.xaxis.set_minor_locator(NullLocator())
    b.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: {25: "25", 100: "100", 1000: "1k", 5000: "5k"}.get(v, "")))

    c = axes[2]
    c.set_xlim(0.025, 0.002)                    # finer tolerance to the right
    c.xaxis.set_major_locator(FixedLocator([0.02, 0.01, 0.005, 0.0025]))
    c.xaxis.set_minor_locator(NullLocator())
    c.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

    axes[0].set_ylabel(r"Mean HI$_3$ / reference value")
    axes[0].legend(loc="lower right", handlelength=2.2, title=None,
                   bbox_to_anchor=(1.0, 0.12))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.21, top=0.9)
    UPLOAD.mkdir(parents=True, exist_ok=True)
    for out in (HERE, UPLOAD):
        for ext in ("pdf", "png"):
            fig.savefig(out / f"export_sensitivity.{ext}")
    print("saved", HERE / "export_sensitivity.pdf", "and", UPLOAD)

    # HI2 check for the caption
    print("HI2 check:")
    for fn, *_ in specs:
        data, ref, used = fn()
        fig2, ax2 = plt.subplots()
        panel(ax2, data, ref, used, metric="HI2")
        plt.close(fig2)


if __name__ == "__main__":
    main()
