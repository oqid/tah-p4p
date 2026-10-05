"""Generate a geometry-only comparison: python preview_domain_styles.py model.stp."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fluid_domain import (load_fluid_domain, camera_border_collection,
                          domain_surface_collection)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--output", type=Path, default=Path("outputs_stp/domain_styles.png"))
    args = parser.parse_args()
    mesh = load_fluid_domain(args.model)
    vertices = mesh.vertices[:, [0, 2, 1]].copy()
    vertices[:, 1] *= -1
    borders = mesh.border_segments[:, :, [0, 2, 1]].copy()
    borders[:, :, 1] *= -1
    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    pad = 0.02 * max(hi - lo)
    fig = plt.figure(figsize=(14, 10), facecolor="white")
    styles = [("Current transparent shell", "flat", 0.05),
              ("Soft matte lighting", "matte", 0.18),
              ("Subtle satin highlight", "glossy", 0.18),
              ("Shape study: stronger matte surface", "matte", 0.65)]
    for i, (title, shading, opacity) in enumerate(styles, 1):
        ax = fig.add_subplot(2, 2, i, projection="3d", computed_zorder=False)
        ax.add_collection3d(domain_surface_collection(vertices, mesh.triangles,
                            shading=shading, alpha=opacity, zorder=1), autolim=False)
        ax.add_collection3d(camera_border_collection(vertices, mesh.triangles, borders,
                            colors="#455a64", linewidth=0.7, alpha=0.45, zorder=2), autolim=False)
        ax.set_xlim(lo[0] - pad, hi[0] + pad)
        ax.set_ylim(lo[1] - pad, hi[1] + pad)
        ax.set_zlim(lo[2] - pad, hi[2] + pad)
        ax.set_box_aspect(hi - lo + 2 * pad)
        ax.set_proj_type("ortho")
        ax.view_init(elev=15, azim=-80)
        ax.set_axis_off()
        ax.set_title(f"{title}\nSurface opacity {opacity}; line opacity 0.45", fontsize=12)
    fig.suptitle("Your STEP model — same camera in every panel (geometry only)", fontsize=16)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.9, wspace=0, hspace=0.08)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
