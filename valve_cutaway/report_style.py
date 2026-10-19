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
# three-case comparison (validated set: CVD dE >= 8.5; orange needs labels/legend)
CASE_STYLE = {
    "50 cc 120 BPM": dict(color="#2E6DB4", marker="o", ls="-"),
    "35 cc 120 BPM": dict(color="#C2477F", marker="s", ls="--"),
    "35 cc 171 BPM": dict(color="#D98A1E", marker="^", ls="-."),
}
TIME_TERM, STRESS_TERM = "#8FB8DE", "#D98CB3"                            # F2d bars
INK, MUTED, GRID = "#1F1F1F", "#6B6B6B", "#E6E6E6"
FORMATS = ("png", "pdf")


FONT_DIR = HERE / "fonts"          # Helvetica Neue LT Pro .otf files (the user's own copy)


def _register_fonts():
    """Register the bundled Helvetica Neue and make it the default sans-serif."""
    from matplotlib import font_manager
    names = []
    for f in sorted(FONT_DIR.glob("*.otf")) + sorted(FONT_DIR.glob("*.ttf")):
        font_manager.fontManager.addfont(str(f))
        names.append(font_manager.FontProperties(fname=str(f)).get_name())
    if names:
        fam = names[0]
        # family list = per-glyph fallback (symbols such as ∝, →, θ missing from Helvetica)
        plt.rcParams["font.family"] = [fam, "DejaVu Sans"]
    return names


def apply():
    plt.style.use(str(STYLE))
    _register_fonts()


def save(fig, outdir, name, formats=FORMATS):
    for ext in formats:
        fig.savefig(Path(outdir) / f"{name}.{ext}")
    plt.close(fig)
