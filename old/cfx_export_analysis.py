"""
quick_diag.py - fast sanity check on a CFD-Post "Generic" export.

Usage:
    python quick_diag.py path/to/export.csv

Answers, in order:
  1. Are we reading the columns we think we are? (units bug check)
  2. What does the raw node-level shear distribution look like?
  3. What does the raw node-level velocity distribution look like?
  4. Is shear correlated with velocity magnitude the way it should be?
  5. Is there a suspicious floor of near-zero shear nodes?

Run this BEFORE trusting any HI/SA number.
"""

import sys
import re
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Reuse the same parser as haem5.py so we're definitely looking at the same
# numbers. If you'd rather keep this script standalone, paste the parser in.
# ---------------------------------------------------------------------------
def parse_cfd_post_export(filepath):
    with open(filepath, "r") as f:
        text = f.read()

    sections = {}
    current = None
    buf = []
    for line in text.splitlines():
        m = re.match(r"^\s*\[(.+?)\]\s*$", line)
        if m:
            if current is not None:
                sections[current] = buf
            current = m.group(1)
            buf = []
        else:
            if current is not None:
                buf.append(line)
    if current is not None:
        sections[current] = buf

    data_lines = [l for l in sections["Data"] if l.strip() != ""]
    header = data_lines[0]

    # --- KEY: print the header so you can see what columns actually exist ---
    print("Raw header line from [Data]:")
    print(f"  {header!r}")

    n_cols = len(header.split("\t")) if "\t" in header else len(header.split())

    rows = []
    for l in data_lines[1:]:
        parts = l.replace(",", "\t").split()
        parts = [p for p in parts if p != ""]
        if len(parts) == n_cols:
            node_id = None
            vals = parts
        else:
            node_id = int(float(parts[0]))
            vals = parts[1:]
        rows.append((node_id, vals))

    parsed = []
    auto_id = 0
    for node_id, vals in rows:
        vals = [float(v) for v in vals]
        x, y, z, shear, vel = vals[0], vals[1], vals[2], vals[3], vals[4]
        if node_id is None:
            node_id = auto_id
        auto_id = max(auto_id, node_id) + 1
        parsed.append((node_id, x, y, z, shear, vel))

    nodes = pd.DataFrame(parsed, columns=["node_id", "x", "y", "z", "shear", "vel"])
    nodes = nodes.drop_duplicates(subset="node_id").set_index("node_id").sort_index()
    return nodes


def describe(name, arr, unit):
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    print(f"\n--- {name} [{unit}] ---")
    print(f"  n            = {len(arr)}")
    print(f"  min          = {arr.min():.6g}")
    print(f"  p01          = {np.percentile(arr, 1):.6g}")
    print(f"  p25          = {np.percentile(arr, 25):.6g}")
    print(f"  median       = {np.median(arr):.6g}")
    print(f"  p75          = {np.percentile(arr, 75):.6g}")
    print(f"  p95          = {np.percentile(arr, 95):.6g}")
    print(f"  p99          = {np.percentile(arr, 99):.6g}")
    print(f"  max          = {arr.max():.6g}")
    print(f"  mean         = {arr.mean():.6g}")
    print(f"  std          = {arr.std():.6g}")
    # "thin high end" check: how much of the mass is above 5x the median?
    med = np.median(arr)
    if med > 0:
        frac_above_5x = np.mean(arr > 5 * med)
        print(f"  frac > 5*median = {frac_above_5x*100:.2f} %")
    print(f"  frac == 0    = {np.mean(arr == 0)*100:.2f} %")
    print(f"  frac < 1e-6  = {np.mean(arr < 1e-6)*100:.2f} %")


def main(filepath):
    print("=" * 72)
    print(f"Diagnostic on: {filepath}")
    print("=" * 72)

    nodes = parse_cfd_post_export(filepath)

    # -----------------------------------------------------------------
    # 1. Column sanity - the most common "low shear" cause is a units
    #    mismatch or reading the wrong column. Print the first few rows
    #    so you can eyeball whether shear looks like Pa or dyne/cm^2.
    # -----------------------------------------------------------------
    print("\nFirst 5 rows (node_id, x, y, z, shear, vel):")
    print(nodes.head().to_string())

    # -----------------------------------------------------------------
    # 2. Distributions
    # -----------------------------------------------------------------
    describe("Shear", nodes["shear"].to_numpy(), "unknown - check header")
    describe("Velocity", nodes["vel"].to_numpy(), "unknown - check header")

    # -----------------------------------------------------------------
    # 3. Physical plausibility checks
    # -----------------------------------------------------------------
    shear = nodes["shear"].to_numpy()
    vel = nodes["vel"].to_numpy()

    print("\n--- Plausibility checks ---")

    # Shear magnitude: in a blood pump / TAH you expect peak tau of order
    # 10-500 Pa in jets and gaps, median of order 0.1-5 Pa in bulk flow.
    # If the median is <0.01 Pa or the max is <1 Pa, something is off.
    med_shear = np.median(shear)
    max_shear = np.max(shear)
    if med_shear < 0.01:
        print(f"  [!] median shear = {med_shear:.4g} - suspiciously low.")
        print("      If your header says [Pa] but this is ~1000x smaller than")
        print("      expected, you may be reading dyne/cm^2 or a strain-rate")
        print("      magnitude instead of a stress. Check the export header.")
    if max_shear < 1.0:
        print(f"  [!] max shear = {max_shear:.4g} - no high-shear region at all.")
        print("      This is almost certainly an export/seeding problem, not")
        print("      a physics result - every blood pump has a high-shear gap.")

    # Velocity magnitude: expect peak of order 1-5 m/s in a TAH.
    max_vel = np.max(vel)
    if max_vel < 0.1:
        print(f"  [!] max velocity = {max_vel:.4g} - seems too low for a TAH.")
        print("      Check units (m/s vs mm/s) and whether you exported the")
        print("      velocity MAGNITUDE or just one component.")

    # -----------------------------------------------------------------
    # 4. Correlation: shear should scale roughly with |grad u|, so in a
    #    shear-thinning or Newtonian flow the high-shear nodes should be
    #    a subset of the high-velocity nodes (plus boundary layers).
    #    If shear is uncorrelated with velocity, you may be reading a
    #    wall-normal component or a wrong variable entirely.
    # -----------------------------------------------------------------
    if len(shear) > 10:
        corr = np.corrcoef(shear, vel)[0, 1]
        print(f"\n  corr(shear, velocity) = {corr:.3f}")
        if abs(corr) < 0.1:
            print("  [!] shear and velocity are essentially uncorrelated.")
            print("      That is unusual for a Newtonian flow field and")
            print("      suggests one of the two columns is not what you")
            print("      think it is.")

    # -----------------------------------------------------------------
    # 5. Spatial check: if your high-shear nodes are all clustered at
    #    the inlet/outlet, you're probably only seeding the core flow.
    # -----------------------------------------------------------------
    print("\n--- Spatial distribution of high-shear nodes (top 5%) ---")
    p95 = np.percentile(shear, 95)
    hi = nodes[shear >= p95]
    print(f"  threshold p95 = {p95:.6g}")
    for ax in ("x", "y", "z"):
        print(f"  {ax}: [{hi[ax].min():.4g}, {hi[ax].max():.4g}]  "
              f"(all nodes: [{nodes[ax].min():.4g}, {nodes[ax].max():.4g}])")
    print("  If the high-shear bounds are much narrower than the full domain,")
    print("  your seeding is only capturing one region.")

    # -----------------------------------------------------------------
    # 6. Bounding box, to catch unit errors (mm vs m) in the coordinates
    # -----------------------------------------------------------------
    print("\n--- Domain bounding box ---")
    for ax in ("x", "y", "z"):
        span = nodes[ax].max() - nodes[ax].min()
        print(f"  {ax}: [{nodes[ax].min():.4g}, {nodes[ax].max():.4g}]  span = {span:.4g}")
    print("  A TAH is roughly 0.05-0.10 m across. If spans are ~100x that,")
    print("  coordinates are in mm. If ~1000x, in micrometres.")

    print("\n" + "=" * 72)
    print("Done. If shear median < 0.01 Pa or max < 1 Pa, do NOT trust")
    print("any HI/SA output from haem5.py until this is resolved.")
    print("=" * 72)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python quick_diag.py <export.csv>")
        sys.exit(1)
    main(sys.argv[1])