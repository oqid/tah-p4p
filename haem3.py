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

Author: Claude and Sevan Dalzell, MECHENG 700 
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

    if "Data" not in sections or not ({"Lines", "Faces"} & sections.keys()):
        raise ValueError(
            "Expected [Data] plus a [Lines] or [Faces] section in the export "
            "file. Make sure 'Node Numbers' and 'Line and Face Connectivity' "
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

    # --- Parse [Lines], or fall back to [Faces] ---
    if "Lines" in sections:
        line_pairs = []
        for l in sections["Lines"]:
            parts = l.replace(",", "\t").split()
            parts = [p for p in parts if p != ""]
            if len(parts) >= 2:
                a, b = int(float(parts[0])), int(float(parts[1]))
                line_pairs.append((a, b))
        return nodes, line_pairs

    nodes, line_pairs = _reconstruct_positions_from_faces(nodes, sections["Faces"])
    return nodes, line_pairs


def _reconstruct_positions_from_faces(nodes, faces_lines):
    """
    CFD-Post exports a streamline plotted as a Ribbon/Tube (rather than a plain
    Line) as a [Faces] section instead of [Lines]: every point along the path
    is duplicated across two parallel rails (e.g. a tiny z-offset), and each
    quad face stitches an adjacent pair of rail-points together.

    Within any one face, the two SHORT edges are "rungs" joining the same
    physical point across its two duplicate node ids; the two LONG edges join
    that point to its neighbour along the streamline. Classifying edges by
    relative length (not an absolute distance threshold) lets this collapse
    the duplicates back into single points and rebuild the ordered line pairs
    regardless of the export's rail spacing.
    """
    faces = []
    for l in faces_lines:
        parts = [p for p in l.replace(",", "\t").split() if p != ""]
        if len(parts) >= 4:
            faces.append(tuple(int(float(p)) for p in parts[:4]))

    if not faces:
        return nodes, []

    coords = nodes[["x", "y", "z"]]

    def edge_len(a, b):
        return np.linalg.norm(coords.loc[a].to_numpy() - coords.loc[b].to_numpy())

    parent = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    long_edges = set()
    for n0, n1, n2, n3 in faces:
        edges = [(n0, n1), (n1, n2), (n2, n3), (n3, n0)]
        lens = [edge_len(a, b) for a, b in edges]
        order = np.argsort(lens)
        for i in order[:2]:      # two shortest -> same-point "rung"
            union(*edges[i])
        for i in order[2:]:      # two longest -> streamline-direction edge
            long_edges.add(edges[i])

    pos_of = {n: find(n) for n in nodes.index}
    collapsed = nodes.copy()
    collapsed["_pos"] = collapsed.index.map(pos_of)
    pos_nodes = collapsed.groupby("_pos")[["x", "y", "z", "shear", "vel"]].mean()
    pos_nodes.index.name = "node_id"

    adj = {}
    for a, b in long_edges:
        pa, pb = pos_of[a], pos_of[b]
        if pa == pb:
            continue
        adj.setdefault(pa, set()).add(pb)
        adj.setdefault(pb, set()).add(pa)

    line_pairs = []
    seen = set()
    for n in adj:
        if n in seen:
            continue
        comp = set()
        stack = [n]
        while stack:
            m = stack.pop()
            if m in comp:
                continue
            comp.add(m)
            stack.extend(adj.get(m, ()))
        seen |= comp

        leaves = [m for m in comp if len(adj[m]) == 1]
        start = leaves[0] if leaves else next(iter(comp))

        path = [start]
        visited = {start}
        cur = start
        while True:
            nxts = [m for m in adj.get(cur, ()) if m not in visited]
            if not nxts:
                break
            nxt = nxts[0]
            path.append(nxt)
            visited.add(nxt)
            cur = nxt

        for a, b in zip(path[:-1], path[1:]):
            line_pairs.append((a, b))

    return pos_nodes, line_pairs


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
    """
    Builds a per-streamline dataframe and evaluates the discretized hemolysis sum.
    """
    pts = nodes.loc[streamline_node_ids].copy()
    pts = pts.reset_index()

    xyz = pts[["x", "y", "z"]].to_numpy()
    seg_len = np.zeros(len(pts))
    seg_len[1:] = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    pts["arc_length"] = seg_len

    # Average velocity across each segment
    vel = pts["vel"].to_numpy()
    seg_vel = np.zeros(len(pts))
    seg_vel[1:] = 0.5 * (vel[1:] + vel[:-1])
    
    # --- FIX: Better velocity handling ---
    # Calculate a reasonable minimum velocity (1% of mean velocity in this streamline)
    mean_vel = np.mean(vel[vel > 0]) if np.any(vel > 0) else 0.1
    min_vel = max(mean_vel * 0.01, 0.001)  # 1% of mean, but at least 0.001 m/s
    
    # Cap velocities at the minimum to avoid infinite exposure time
    seg_vel = np.maximum(seg_vel, min_vel)
    
    dt = np.zeros(len(pts))
    dt[1:] = seg_len[1:] / seg_vel[1:]
    pts["dt"] = dt
    pts["t_cumulative"] = np.cumsum(dt)

    # Discretized Giersiepen sum
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

def _axis_range(streamline_dfs, col):
    lo = min(df[col].min() for df in streamline_dfs)
    hi = max(df[col].max() for df in streamline_dfs)
    return lo, hi, hi - lo


def plot_streamlines_3d(streamline_dfs, save_path=None, max_lines=None,
                         flat_tolerance=1e-3):
    """
    Plots reconstructed streamlines coloured by local shear stress.

    Automatically detects near-planar (2D CFD case) data - where one axis
    has negligible range compared to the other two - and switches to a
    clean 2D plot instead of a distorted 3D one. For genuinely 3D data,
    keeps proper axis proportions (equal aspect box) so the geometry
    isn't stretched.
    """
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    x_lo, x_hi, x_range = _axis_range(streamline_dfs, "x")
    y_lo, y_hi, y_range = _axis_range(streamline_dfs, "y")
    z_lo, z_hi, z_range = _axis_range(streamline_dfs, "z")

    # all_shear = np.concatenate([df["shear"].to_numpy() for df in streamline_dfs])
    # vmin, vmax = np.nanmin(all_shear), np.nanmax(all_shear)
    # cmap = plt.get_cmap("turbo")
    # norm = plt.Normalize(vmin=vmin, vmax=vmax)

    all_shear = np.concatenate([df["shear"].to_numpy() for df in streamline_dfs])
    from matplotlib.colors import LogNorm
    vmin, vmax = np.nanmin(all_shear), np.nanmax(all_shear)
    cmap = plt.get_cmap("turbo")
    norm = LogNorm(vmin=vmin, vmax=vmax)  # ← use LogNorm instead of Normalize

    dfs = streamline_dfs if max_lines is None else streamline_dfs[:max_lines]
    in_plane_scale = max(x_range, y_range, 1e-12)

    # A "2D case" is one where the out-of-plane axis barely varies compared
    # to the other two - e.g. a single-cell-thick CFX slab extrusion.
    is_flat_z = z_range < flat_tolerance * in_plane_scale
    is_flat_y = y_range < flat_tolerance * max(x_range, z_range, 1e-12)

    if is_flat_z or is_flat_y:
        # --- 2D rendering ---
        if is_flat_z:
            xa, ya, xlabel, ylabel = "x", "y", "X [m]", "Y [m]"
        else:
            xa, ya, xlabel, ylabel = "x", "z", "X [m]", "Z [m]"

        fig, ax = plt.subplots(figsize=(11, 6))
        for df in dfs:
            pts = df[[xa, ya]].to_numpy()
            shear = df["shear"].to_numpy()
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            seg_colors = cmap(norm(0.5 * (shear[:-1] + shear[1:])))
            lc = LineCollection(segs, colors=seg_colors, linewidth=1.4)
            ax.add_collection(lc)

        ax.set_xlim(min(df[xa].min() for df in dfs), max(df[xa].max() for df in dfs))
        ax.set_ylim(min(df[ya].min() for df in dfs), max(df[ya].max() for df in dfs))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title("Reconstructed streamlines coloured by local shear stress\n"
                      "(rendered in 2D - out-of-plane axis range was negligible)")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.85, pad=0.02)
        cbar.set_label("Local shear stress [Pa]")

    else:
        # --- genuine 3D rendering, with real proportions preserved ---
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection="3d")

        for df in dfs:
            pts = df[["x", "y", "z"]].to_numpy()
            shear = df["shear"].to_numpy()
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            seg_colors = cmap(norm(0.5 * (shear[:-1] + shear[1:])))
            lc = Line3DCollection(segs, colors=seg_colors, linewidth=1.2)
            ax.add_collection3d(lc)

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_zlim(z_lo, z_hi)
        # Preserve true relative proportions instead of stretching to a cube -
        # this is what prevents the "distorted box" artifact.
        ax.set_box_aspect((max(x_range, 1e-9), max(y_range, 1e-9), max(z_range, 1e-9)))

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
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

    # After run_pipeline in your main script:
    print("\n=== Streamline Quality Check ===")
    print(f"Total streamlines: {len(result['streamlines'])}")
    print(f"Average points per streamline: {np.mean([len(df) for df in result['streamlines']])}")
    print(f"Min/Max velocity: {result['nodes']['vel'].min():.3f} / {result['nodes']['vel'].max():.3f} m/s")
    print(f"Min/Max shear: {result['nodes']['shear'].min():.1f} / {result['nodes']['shear'].max():.1f} Pa")

    # Check for zero-velocity issues
    zero_vel = result['nodes'][result['nodes']['vel'] < 0.001]
    if len(zero_vel) > 0:
        print(f"WARNING: {len(zero_vel)} nodes with velocity < 0.001 m/s")
        print(f"  These are in regions that might be underestimated")
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