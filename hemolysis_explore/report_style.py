"""One place for the look of the report figures.

    import report_style as rs
    rs.apply()                      # loads p4p_report.mplstyle
    fig, ax = plt.subplots(figsize=(rs.WIDTH, 2.5))

Change fonts/sizes/grid in p4p_report.mplstyle; change colours/markers here.
Colours: phases use a lightness-ordered purple -> pink -> peach ramp (bi/trans
flag inspired, see AGENTS.md) and every phase also has its own marker, so
colour is never the only cue. Device colours are a validated blue/pink pair.
"""
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STYLE = HERE / "p4p_report.mplstyle"

WIDTH = 6.3                     # inches ≈ 16 cm text width; use WIDTH/2 for half-width figures

PHASE_COLOURS = ["#2D1E5C", "#5B3A8E", "#9B4F96", "#D0668F", "#EE9A6B"]   # 16.67 … 83 %
PHASE_MARKERS = ["o", "s", "D", "^", "v"]
DEVICE_COLOURS = {"50 cc": "#2E6DB4", "35 cc": "#C2477F"}
TIME_TERM, STRESS_TERM = "#8FB8DE", "#D98CB3"                            # F2d bars
INK, MUTED, GRID = "#1F1F1F", "#6B6B6B", "#E6E6E6"
FORMATS = ("png", "pdf")


def apply():
    plt.style.use(str(STYLE))


def save(fig, outdir, name, formats=FORMATS):
    for ext in formats:
        fig.savefig(Path(outdir) / f"{name}.{ext}")
    plt.close(fig)
