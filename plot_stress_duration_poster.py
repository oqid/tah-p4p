"""Poster version of haem5_empiricalthresholds.plot_empirical_thresholds.

Run with .venv/Scripts/python plot_stress_duration_poster.py
Reads saved results; does not rerun or modify the haemolysis analysis.
Exports an exact 3:4 canvas as PNG, PDF and SVG beside the original plot.

Reference: Zaman et al., doi:10.1111/aor.70149 (mean +/- SD, n=6).
The reference's 50 cc designation is not verified, so it is not printed.
NIH uses the repository's conversion convention, Hct=0.30, Hb=150 g/L;
it is a model-derived comparison, not measured NIH or matched validation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, LogFormatterMathtext
import numpy as np
import pandas as pd

from haem5_empiricalthresholds import RBC_THRESHOLD, PLATELET_THRESHOLD


ROOT = Path(__file__).resolve().parent
DEFAULT_CASE = ROOT / "outputs/runs/2026-10-06_17-34-08_756471/01_INPUT_TEST_3D"
HCT, HB = 0.30, 150.0


def make_poster(case_dir: Path):
    p = pd.read_csv(case_dir / "empirical/exposure_envelope.csv")
    p = p.replace([np.inf, -np.inf], np.nan).dropna()
    p = p[(p.duration_s > 0) & (p.stress_Pa > 0)]
    p = p.sort_values(["duration_s", "stress_Pa"], ascending=[True, False])
    if p.empty:
        raise ValueError("No positive, finite exposure-envelope data.")
    summary = pd.read_csv(case_dir.parent / "summary.csv")
    case = case_dir.name.split("_", 1)[1]
    rows = summary.loc[summary["case"] == case]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one summary row for {case!r}.")
    hi = float(rows.iloc[0]["device_HI3_percent"])
    nih = hi * (1 - HCT) * HB * 1000
    if not np.isfinite(hi) or not np.isclose(
        nih, float(rows.iloc[0]["device_HI3_NIH_mg_per_100L"]), rtol=1e-8
    ):
        raise ValueError("Stored NIH does not match the adopted conversion.")

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 20,
                         "axes.linewidth": 1.8, "svg.fonttype": "none",
                         "pdf.fonttype": 42}):
        fig = plt.figure(figsize=(9, 12), facecolor="white")
        ax = fig.add_axes([0.16, 0.305, 0.81, 0.67])
        ax.loglog(RBC_THRESHOLD.time_s * 1000, RBC_THRESHOLD.stress_Pa,
                  color="#D989B5", lw=3.5, label="Red-cell lysis")
        ax.loglog(PLATELET_THRESHOLD.time_s * 1000, PLATELET_THRESHOLD.stress_Pa,
                  color="#9B8AC4", lw=3.5, ls="--", label="Platelet lysis")
        ax.loglog(p.duration_s * 1000, p.stress_Pa, color="#78A9D1", lw=4,
                  label="CFD exposure\nenvelope", zorder=3)
        ax.set_xlabel("Exposure duration [ms]", fontsize=25, labelpad=12)
        ax.set_ylabel("Shear stress [Pa]", fontsize=25, labelpad=8)
        ax.tick_params(which="both", direction="in", top=True, right=True,
                       labelsize=21, width=1.5)
        ax.tick_params(which="major", length=8, pad=8)
        ax.tick_params(which="minor", length=4)
        ax.grid(which="major", color="#e1e1e1", lw=0.9)
        ax.set_axisbelow(True)
        ax.xaxis.set_major_locator(FixedLocator(10.0 ** np.arange(-3, 6, 2)))
        ax.yaxis.set_major_locator(FixedLocator(10.0 ** np.arange(-2, 6, 2)))
        ax.xaxis.set_major_formatter(LogFormatterMathtext())
        ax.yaxis.set_major_formatter(LogFormatterMathtext())
        ax.set_xlim(1e-3, 1e5)
        ax.set_ylim(10 ** np.floor(np.log10(min(1., p.stress_Pa.min() / 2))), 2e5)
        ax.legend(loc="upper right", fontsize=19, handlelength=2.0,
                  labelspacing=0.6, frameon=True, facecolor="white",
                  edgecolor="none", framealpha=0.95)

        # Compact calculation and side-by-side comparison, with no title/caption.
        fig.add_artist(Line2D([0.08, 0.97], [0.215, 0.215],
                             transform=fig.transFigure, color="#cccccc", lw=1))
        fig.text(0.525, 0.185,
                 r"$\mathrm{NIH}=1000\,\mathrm{HI}_3[\%]\,(1-\mathrm{Hct})\,\mathrm{Hb}$",
                 ha="center", va="center", fontsize=21)
        fig.text(0.525, 0.151, r"$\mathrm{Hct}=0.30\quad\mathrm{Hb}=150\ \mathrm{g/L}$",
                 ha="center", va="center", fontsize=18)
        fig.text(0.29, 0.111, "This test (HI3, GW)", ha="center", fontsize=19)
        fig.text(0.75, 0.111, "SynCardia (in vitro)", ha="center", fontsize=19)
        fig.text(0.29, 0.075, f"{nih:.2f}", ha="center", fontsize=30, weight="bold")
        fig.text(0.75, 0.075, r"$37.15\pm12.42$", ha="center", fontsize=30)
        fig.text(0.29, 0.049, "mg/100 L", ha="center", fontsize=18)
        fig.text(0.75, 0.049, "mg/100 L", ha="center", fontsize=18)
        exponent = int(np.floor(np.log10(hi))) if hi > 0 else 0
        mantissa = hi / 10 ** exponent
        fig.text(0.29, 0.019,
                 rf"$\mathrm{{HI}}_3={mantissa:.3f}\times10^{{{exponent}}}\,\%$",
                 ha="center", fontsize=18)
        fig.text(0.75, 0.019, "Mean ± SD", ha="center", fontsize=18)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--output", type=Path, help="Output stem (without extension)")
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    fig = make_poster(args.case_dir.resolve())
    stem = args.output or args.case_dir / "empirical/stress_duration_poster"
    stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        path = stem.with_suffix("." + extension)
        # No tight crop: preserve the exact 3:4 aspect ratio in every format.
        fig.savefig(path, dpi=args.dpi, facecolor="white", bbox_inches=None)
        print(path.resolve())
    plt.close(fig)


if __name__ == "__main__":
    main()
