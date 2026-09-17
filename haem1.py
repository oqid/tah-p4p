"""
Hemolysis pipeline for CFD-Post streamline exports (CFX / CFD-Post "Generic" export format).

Implements the discretized cumulative Heuser-Opitz / Giersiepen power-law model
described in section 2.1 of the MECHENG 700 mid-year report:

    d(HI%)_i = beta * A * t_i^(beta-1) * tau_i^alpha * dt_i
    HI% = sum_i d(HI%)_i

then converts to NIH via:

    NIH = HI% * (1 - Hct) * Hb

Usage:
    python hemolysis_pipeline.py path/to/export.csv

Author: generated for Sevan Dalzell, MECHENG 700 capstone
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

# ---------------------------------------------------------------------------
# 1. Model constants
# ---------------------------------------------------------------------------

# Giersiepen et al. power-law coefficients (commonly cited values; confirm
# against the exact source you cited in the report before trusting results).
A_COEF = 3.62e-5      # units depend on source - check your reference!
ALPHA = 2.416
BETA = 0.785

# Blood properties used for NIH conversion (from your report, section 2.1)
HCT = 0.30             # hematocrit (volume fraction)
HB = 150.0             # g/L

# Literature comparison points
SYNCARDIA_50CC_NIH = 37.15   # mg/100L, +/- 12.42
ASTM_LIMIT_NIH = 100.0       # mg/100L  (0.1 g/100L expressed in mg/100L)


# ---------------------------------------------------------------------------
# 2. Parsing the CFD-Post "Generic" export format
# ---------------------------------------------------------------------------

def parse_cfd_post_export(filepath):
    """
    Parses a CFD-Post Generic CSV export with [Name]/[Data]/[Lines] sections.

    Returns
    -------
    nodes : pd.DataFrame indexed by node id, columns: x, y, z, shear, vel
    line_pairs : list of (start_node, end_node) tuples, in file order
    """
    with open(filepath, "r") as f:
        text = f.read()

    # Split into bracketed sections, e.g. [Name], [Data], [Lines]
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

    if "Data" not in sections or "Lines" not in sections:
        raise ValueError(
            "Expected [Data] and [Lines] sections in the export file. "
            "Make sure 'Node Numbers' and 'Line and Face Connectivity' "
            "were both ticked in the CFD-Post export dialog."
        )

    # --- Parse [Data] ---
    data_lines = [l for l in sections["Data"] if l.strip() != ""]
    header = data_lines[0]
    n_cols = len(header.split("\t")) if "\t" in header else len(header.split())

    rows = []
    for l in data_lines[1:]:
        parts = l.replace(",", "\t").split()
        parts = [p for p in parts if p != ""]
        if len(parts) == n_cols:
            # No leading node-id column -> this is a release-point row,
            # assign it a node id in file order later.
            node_id = None
            vals = parts
        else:
            # Leading node-id column present
            node_id = int(float(parts[0]))
            vals = parts[1:]
        rows.append((node_id, vals))

    # Expected value columns: X, Y, Z, shear, Velocity, X(dup), Y(dup), Z(dup)
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

    # --- Parse [Lines] ---
    line_pairs = []
    for l in sections["Lines"]:
        parts = l.replace(",", "\t").split()
        parts = [p for p in parts if p != ""]
        if len(parts) >= 2:
            a, b = int(float(parts[0])), int(float(parts[1]))
            line_pairs.append((a, b))

    return nodes, line_pairs


def reconstruct_streamlines(line_pairs):
    """
    Walks the [Lines] connectivity pairs and splits them into individual
    streamlines. A new streamline starts whenever the next pair's start
    node does not match the previous pair's end node.

    Returns
    -------
    list of lists, each inner list is the ordered sequence of node ids
    for one streamline.
    """
    streamlines = []
    current = None

    for start, end in line_pairs:
        if current is None:
            current = [start, end]
        elif start == current[-1]:
            current.append(end)
        else:
            streamlines.append(current)
            current = [start, end]

    if current is not None:
        streamlines.append(current)

    return streamlines


# ---------------------------------------------------------------------------
# 3. Hemolysis calculation
# ---------------------------------------------------------------------------

def compute_streamline_hemolysis(nodes, streamline_node_ids,
                                  A=A_COEF, alpha=ALPHA, beta=BETA):
    """
    Builds a per-streamline dataframe (position, shear, velocity, arc length,
    dt, cumulative time) and evaluates the discretized hemolysis sum.

    Returns
    -------
    df : pd.DataFrame for this streamline (one row per node along the path)
    HI_percent : float, total HI% accumulated along this streamline
    """
    pts = nodes.loc[streamline_node_ids].copy()
    pts = pts.reset_index()

    xyz = pts[["x", "y", "z"]].to_numpy()
    seg_len = np.zeros(len(pts))
    seg_len[1:] = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    pts["arc_length"] = seg_len

    # Average velocity across each segment (trapezoidal), guard against v=0
    vel = pts["vel"].to_numpy()
    seg_vel = np.zeros(len(pts))
    seg_vel[1:] = 0.5 * (vel[1:] + vel[:-1])
    seg_vel[seg_vel < 1e-6] = np.nan  # avoid divide-by-zero -> becomes NaN dt

    dt = np.zeros(len(pts))
    dt[1:] = seg_len[1:] / seg_vel[1:]
    dt = np.nan_to_num(dt, nan=0.0)
    pts["dt"] = dt
    pts["t_cumulative"] = np.cumsum(dt)

    # Discretized Giersiepen sum: dHI_i = beta * A * t_i^(beta-1) * tau_i^alpha * dt_i
    t = pts["t_cumulative"].to_numpy()
    tau = pts["shear"].to_numpy()
    dt_arr = pts["dt"].to_numpy()

    with np.errstate(divide="ignore", invalid="ignore"):
        t_pow = np.where(t > 0, t ** (beta - 1), 0.0)
    dHI = beta * A * t_pow * (tau ** alpha) * dt_arr
    dHI = np.nan_to_num(dHI, nan=0.0, posinf=0.0, neginf=0.0)

    pts["dHI_percent"] = dHI
    pts["HI_percent_cumulative"] = np.cumsum(dHI)

    HI_percent_total = pts["HI_percent_cumulative"].iloc[-1] if len(pts) else 0.0
    return pts, HI_percent_total


def run_pipeline(filepath, A=A_COEF, alpha=ALPHA, beta=BETA):
    """
    Full pipeline: parse -> reconstruct -> compute per-streamline hemolysis
    -> flow-weighted device average -> NIH conversion.
    """
    nodes, line_pairs = parse_cfd_post_export(filepath)
    streamline_ids = reconstruct_streamlines(line_pairs)

    results = []
    per_streamline_dfs = []
    for sid_list in streamline_ids:
        # keep only node ids that actually exist in the parsed data table
        sid_list = [n for n in sid_list if n in nodes.index]
        if len(sid_list) < 2:
            continue
        df, HI_total = compute_streamline_hemolysis(nodes, sid_list, A, alpha, beta)
        per_streamline_dfs.append(df)
        results.append({
            "n_points": len(sid_list),
            "HI_percent": HI_total,
            "mean_inlet_velocity": df["vel"].iloc[0],
        })

    summary = pd.DataFrame(results)

    if len(summary) == 0:
        raise RuntimeError("No valid streamlines were reconstructed - check the export file.")

    # Flow-weighted average: weight each streamline by its inlet velocity
    # (proxy for the flow rate it represents, since inlet sample points are
    # evenly spaced -> equal area weighting, so velocity approximates flux share)
    weights = summary["mean_inlet_velocity"].to_numpy()
    weights = weights / weights.sum()
    device_HI_percent = float(np.sum(summary["HI_percent"].to_numpy() * weights))

    device_NIH = device_HI_percent * (1 - HCT) * HB  # g/100L (per report eq. 3 units)
    device_NIH_mg = device_NIH * 1000.0  # convert g/100L -> mg/100L for comparison

    return {
        "nodes": nodes,
        "streamlines": per_streamline_dfs,
        "summary": summary,
        "device_HI_percent": device_HI_percent,
        "device_NIH_mg_per_100L": device_NIH_mg,
    }


# ---------------------------------------------------------------------------
# 4. Visualization
# ---------------------------------------------------------------------------

def plot_streamlines_3d(streamline_dfs, save_path=None, max_lines=None):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")

    all_shear = np.concatenate([df["shear"].to_numpy() for df in streamline_dfs])
    vmin, vmax = np.nanmin(all_shear), np.nanmax(all_shear)
    cmap = plt.get_cmap("turbo")

    dfs = streamline_dfs if max_lines is None else streamline_dfs[:max_lines]
    for df in dfs:
        pts = df[["x", "y", "z"]].to_numpy()
        shear = df["shear"].to_numpy()
        norm_shear = (shear - vmin) / (vmax - vmin + 1e-30)
        for i in range(len(pts) - 1):
            ax.plot(pts[i:i+2, 0], pts[i:i+2, 1], pts[i:i+2, 2],
                     color=cmap(norm_shear[i]), linewidth=1.2)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("Local shear stress [Pa]")

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("Reconstructed streamlines coloured by local shear stress")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    return fig


def plot_shear_vs_time(streamline_dfs, save_path=None, n_lines=8):
    fig, ax = plt.subplots(figsize=(9, 5))
    for df in streamline_dfs[:n_lines]:
        ax.plot(df["t_cumulative"], df["shear"], marker="o", markersize=2, alpha=0.8)
    ax.set_xlabel("Cumulative exposure time, t [s]")
    ax.set_ylabel("Local shear stress, τ [Pa]")
    ax.set_title(f"Shear stress vs exposure time (first {n_lines} streamlines)")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    return fig


def plot_hi_histogram(summary, save_path=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(summary["HI_percent"], bins=30, color="#c0392b", edgecolor="black", alpha=0.85)
    ax.set_xlabel("Streamline HI% (unweighted)")
    ax.set_ylabel("Number of streamlines")
    ax.set_title("Distribution of per-streamline hemolysis index")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    return fig


def print_summary(result):
    print("=" * 60)
    print(f"Streamlines processed:        {len(result['summary'])}")
    print(f"Device HI% (flow-weighted):   {result['device_HI_percent']:.6e} %")
    print(f"Device NIH:                   {result['device_NIH_mg_per_100L']:.4f} mg/100L")
    print(f"SynCardia 50cc literature:    {SYNCARDIA_50CC_NIH} +/- 12.42 mg/100L")
    print(f"ASTM acceptable limit:        {ASTM_LIMIT_NIH} mg/100L")
    print("=" * 60)


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python hemolysis_pipeline.py <export.csv>")
        sys.exit(1)

    filepath = sys.argv[1]
    result = run_pipeline(filepath)
    print_summary(result)

    plot_streamlines_3d(result["streamlines"], save_path="streamlines_3d.png")
    plot_shear_vs_time(result["streamlines"], save_path="shear_vs_time.png")
    plot_hi_histogram(result["summary"], save_path="hi_histogram.png")

    result["summary"].to_csv("streamline_summary.csv", index=False)
    print("Saved: streamlines_3d.png, shear_vs_time.png, hi_histogram.png, streamline_summary.csv")
