"""


This is a variation of haem5.py where im checking for if too many of the streamlines
get floored by the min velocity, or recirculate too many times and get ignored



Hemolysis pipeline for CFD-Post streamline exports (CFX / CFD-Post "Generic" export format).

Implements two Lagrangian power-law damage-accumulation models, following the
naming/derivation in Taskin et al., "Evaluation of Eulerian and Lagrangian
Models for Hemolysis Estimation," ASAIO J 2012;58:363-372:

HI2 (their Eq. 4) - discretized temporal derivative of the base Giersiepen
power law, described in section 2.1 of the MECHENG 700 mid-year report:

    d(HI%)_i = beta * A * t_i^(beta-1) * tau_i^alpha * dt_i
    HI% = sum_i d(HI%)_i

Note: t_i is *cumulative* time since the start of the streamline, so because
beta < 1 (beta-1 is negative), a given (tau, dt) pair contributes LESS damage
the later it occurs along the path. Taskin et al. identify this as HI2's key
weakness - it systematically under-predicts hemolysis when the high-shear
region occurs mid- or late-path rather than at the inlet, and it had the
lowest correlation coefficients of the Lagrangian methods they tested.

HI3 (their Eq. 5, after Garon & Farinas 2004) - sum of purely local segment
damage, with no dependence on cumulative path time:

    HI%_i = [ A^(1/beta) * tau_i^(alpha/beta) * dt_i ]^beta
    HI%   = sum_i HI%_i

HI3 does not have the late-path suppression artifact of HI2, is independent
of time-step size, and gave the smallest relative error against experiment
of the Lagrangian methods in Taskin et al. (Table 3). Both are computed here
so they can be compared directly on the same streamline data.

Both convert to NIH via:

    NIH = HI% * (1 - Hct) * Hb

Usage:
    python hemolysis_pipeline.py path/to/export.csv

Useful comparison macro:
    python haem5.py export.csv --compare-constants

Author: Claude and Sevan Dalzell, MECHENG 700 
"""

import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)
from tqdm import tqdm

# ---------------------------------------------------------------------------
# 1. Model constants
# ---------------------------------------------------------------------------

POWER_LAW_CONSTANTS = {
    "GW": {"A": 3.62e-5, "alpha": 2.416, "beta": 0.785,
           "tau_max_Pa": 255.0, "t_max_s": 0.700},
    "HO": {"A": 1.80e-6, "alpha": 1.991, "beta": 0.765,
           "tau_max_Pa": 700.0, "t_max_s": 0.700},
    "TZ": {"A": 1.228e-5, "alpha": 1.9918, "beta": 0.6606,
           "tau_max_Pa": 320.0, "t_max_s": 1.500},
}

ACTIVE_CONSTANTS = "GW"
A_COEF = POWER_LAW_CONSTANTS[ACTIVE_CONSTANTS]["A"]
ALPHA = POWER_LAW_CONSTANTS[ACTIVE_CONSTANTS]["alpha"]
BETA = POWER_LAW_CONSTANTS[ACTIVE_CONSTANTS]["beta"]

HCT = 0.30             # hematocrit (volume fraction)
HB = 150.0             # g/L

SYNCARDIA_50CC_NIH = 37.15   # mg/100L, +/- 12.42
ASTM_LIMIT_NIH = 100.0       # mg/100L  (0.1 g/100L expressed in mg/100L)

PA_S_TO_DYNE_S_CM2 = 10.0

HELLUMS_THRESHOLD_DYNE_S_CM2 = 35.0

SA_LITERATURE_MEDIANS = {
    "70cc_TAH": 1.12,
    "50cc_TAH": 2.19,
    "35cc_TAH": 3.41,
}

GENERATE_OUTPUTS = True
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

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

    data_lines = [l for l in sections["Data"] if l.strip() != ""]
    header = data_lines[0]
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
        for i in order[:2]:
            union(*edges[i])
        for i in order[2:]:
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
    dt, cumulative time) and evaluates two discretized Lagrangian hemolysis
    sums along it - HI2 and HI3, per Taskin et al. 2012 (Eqs. 4 and 5).
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
    mean_vel = np.mean(vel[vel > 0]) if np.any(vel > 0) else 0.1
    min_vel = max(mean_vel * 0.01, 0.001)  # 1% of mean, but at least 0.001 m/s

    # --- NEW: per-segment truncation diagnostic (see diagnose_streamline_coverage) ---
    truncated_mask = (vel > 0) & (vel < min_vel)
    if np.any(truncated_mask):
        num_truncated = np.sum(truncated_mask)
        print(f"[TRUNCATION] Streamline: {num_truncated} values truncated to min_vel={min_vel:.6f} "
              f"(mean_vel={mean_vel:.6f}, original min={np.min(vel[truncated_mask]):.6f})")
        frac = num_truncated / np.sum(vel > 0)
        if frac > 0.05:  # only warn if >5% truncated
            print(f"[TRUNCATION WARNING] {frac*100:.1f}% of values truncated!")

    # Cap velocities at the minimum to avoid infinite exposure time.
    # --- NEW: record which segments were floored so coverage diagnostics can see it ---
    seg_vel_floored = seg_vel < min_vel
    seg_vel = np.maximum(seg_vel, min_vel)
    pts["seg_vel_floored"] = seg_vel_floored

    dt = np.zeros(len(pts))
    dt[1:] = seg_len[1:] / seg_vel[1:]
    pts["dt"] = dt
    pts["t_cumulative"] = np.cumsum(dt)

    t = pts["t_cumulative"].to_numpy()
    tau = pts["shear"].to_numpy()
    dt_arr = pts["dt"].to_numpy()

    # --- HI2: discretized temporal derivative (Taskin et al. Eq. 4) ---
    with np.errstate(divide="ignore", invalid="ignore"):
        t_pow = np.where(t > 0, t ** (beta - 1), 0.0)

    dHI2 = beta * A * t_pow * (tau ** alpha) * dt_arr
    dHI2 = np.nan_to_num(dHI2, nan=0.0, posinf=0.0, neginf=0.0)

    pts["dHI2_percent"] = dHI2
    pts["HI2_percent_cumulative"] = np.cumsum(dHI2)

    # --- HI3: local/linearized sum (Taskin et al. Eq. 5, after Garon & Farinas) ---
    inner_term = A ** (1.0 / beta) * (tau ** (alpha / beta)) * dt_arr
    inner_term = np.nan_to_num(inner_term, nan=0.0, posinf=0.0, neginf=0.0)
    inner_cumsum = np.cumsum(inner_term)

    pts["HI3_inner_term"] = inner_term
    pts["HI3_percent_cumulative"] = inner_cumsum ** beta

    # Backwards-compatible aliases (old column/variable names -> HI2)
    pts["dHI_percent"] = pts["dHI2_percent"]
    pts["HI_percent_cumulative"] = pts["HI2_percent_cumulative"]

    HI2_percent_total = pts["HI2_percent_cumulative"].iloc[-1] if len(pts) else 0.0
    HI3_percent_total = pts["HI3_percent_cumulative"].iloc[-1] if len(pts) else 0.0
    return pts, HI2_percent_total, HI3_percent_total


def compute_streamline_SA(pts):
    """
    Linear stress-accumulation (SA) model, as used in Marom et al. 2014.
    """
    tau_pa = pts["shear"].to_numpy()
    dt_arr = pts["dt"].to_numpy()

    dSA_pa_s = tau_pa * dt_arr
    dSA_dyne_s_cm2 = dSA_pa_s * PA_S_TO_DYNE_S_CM2

    pts = pts.copy()
    pts["dSA_dyne_s_cm2"] = dSA_dyne_s_cm2
    pts["SA_cumulative_dyne_s_cm2"] = np.cumsum(dSA_dyne_s_cm2)

    SA_total_dyne_s_cm2 = pts["SA_cumulative_dyne_s_cm2"].iloc[-1] if len(pts) else 0.0
    return pts, SA_total_dyne_s_cm2


# --- NEW: per-streamline coverage diagnostic ---------------------------------
def diagnose_streamline_coverage(pts, sid_index=None):
    """
    Per-streamline diagnostic: how much of the path was rescued by the
    velocity floor, and how long the path actually is in time.

    If most streamlines have (a) similar short t_cumulative AND (b) a
    substantial fraction of floored segments, the problem is upstream:
    your CFD-Post seeding isn't capturing the slow / stagnation regions,
    so the hemolysis model is being asked to extrapolate into territory
    it never sees. No amount of min_vel tuning fixes that - you need to
    re-seed streamlines (denser near walls, valve hinges, recirculation
    zones) and re-export.

    Parameters
    ----------
    pts : pd.DataFrame
        Per-streamline dataframe from compute_streamline_hemolysis, must
        have "seg_vel_floored" (bool) and "t_cumulative" columns.
    sid_index : int or str, optional
        Identifier for reporting.

    Returns
    -------
    dict with per-streamline stats, or None if the streamline is too short.
    """
    n_seg = len(pts)
    if n_seg < 2:
        return None

    floored = pts["seg_vel_floored"].to_numpy()
    # segment 0 has no meaningful velocity (it's the seed point), so
    # exclude it from the floored fraction
    floored_segs = floored[1:]
    frac_floored = float(np.mean(floored_segs)) if len(floored_segs) else 0.0

    t_total = float(pts["t_cumulative"].iloc[-1])
    t_median_seg = float(np.median(pts["dt"].to_numpy()[1:])) if n_seg > 1 else 0.0
    t_max_seg = float(np.max(pts["dt"].to_numpy()[1:])) if n_seg > 1 else 0.0

    return {
        "streamline_id": sid_index,
        "n_segments": n_seg - 1,
        "frac_segments_floored": frac_floored,
        "t_total_s": t_total,
        "t_median_segment_s": t_median_seg,
        "t_max_segment_s": t_max_seg,
        "arc_length_m": float(pts["arc_length"].sum()),
    }


# --- NEW: aggregate coverage report -----------------------------------------
def _report_coverage(coverage_df):
    """
    Prints the aggregate coverage diagnostic.

    Interpretation logic:
      - frac_segments_floored is small (<~1%) and t_total is well-distributed
        across streamlines -> the velocity floor is a rare safety net; HI/SA
        numbers are being driven by real flow data. Fine.
      - frac_segments_floored is large (>10% on many streamlines) OR t_total
        is tightly clustered near zero -> the floor is doing real work, and
        the seed points are not resolving low-velocity regions. The fix is
        upstream in CFD-Post seeding, not in the hemolysis math.
    """
    if coverage_df is None or len(coverage_df) == 0:
        print("\n[COVERAGE] No streamlines to diagnose.")
        return

    frac = coverage_df["frac_segments_floored"].to_numpy()
    t_tot = coverage_df["t_total_s"].to_numpy()

    print("\n" + "=" * 70)
    print("Streamline coverage diagnostic")
    print("=" * 70)
    print(f"  Streamlines analysed:              {len(coverage_df)}")
    print(f"  Segments floored (all lines):      {frac.mean()*100:6.2f} %  "
          f"(median per-line: {np.median(frac)*100:.2f} %)")
    print(f"  Streamlines with >10% floored:     "
          f"{np.mean(frac > 0.10)*100:6.2f} %")
    print(f"  Streamlines with >50% floored:     "
          f"{np.mean(frac > 0.50)*100:6.2f} %")
    print("-" * 70)
    print(f"  Total transit time t_total [s]:")
    print(f"    min / median / max:              "
          f"{t_tot.min():.4e} / {np.median(t_tot):.4e} / {t_tot.max():.4e}")
    print(f"    std / mean:                      {t_tot.std():.4e} / {t_tot.mean():.4e}")
    print(f"    coefficient of variation:        "
          f"{t_tot.std()/max(t_tot.mean(),1e-12):.3f}")
    print("-" * 70)

    heavy_floor = np.mean(frac > 0.10) > 0.20
    tight_time = (t_tot.std() / max(t_tot.mean(), 1e-12)) < 0.25

    if heavy_floor and tight_time:
        print("  VERDICT: Floor is doing real work AND transit times are")
        print("           tightly clustered -> you are NOT capturing the")
        print("           slow/recirculation regions. Fix upstream in")
        print("           CFD-Post seeding, not in the hemolysis math.")
    elif heavy_floor:
        print("  VERDICT: Floor is doing real work on a meaningful fraction")
        print("           of segments. Some slow regions exist but are")
        print("           under-resolved. Consider denser seeding there.")
    elif tight_time:
        print("  VERDICT: Floor is rarely hit, but transit times are very")
        print("           uniform - suspicious. Check whether the export")
        print("           contains more than one distinct path family.")
    else:
        print("  VERDICT: Coverage looks healthy - floor is a rare safety")
        print("           net and transit times span a plausible range.")
    print("=" * 70)


def run_pipeline(filepath, A=A_COEF, alpha=ALPHA, beta=BETA):
    """
    Full pipeline: parse -> reconstruct -> compute per-streamline hemolysis
    -> flow-weighted device average -> NIH conversion.
    """
    nodes, line_pairs = parse_cfd_post_export(filepath)
    streamline_ids = reconstruct_streamlines(line_pairs)

    results = []
    per_streamline_dfs = []
    coverage_diags = []          # --- NEW
    for i, sid_list in enumerate(tqdm(streamline_ids, desc="Hemolysis", unit="line")):
        # keep only node ids that actually exist in the parsed data table
        sid_list = [n for n in sid_list if n in nodes.index]
        if len(sid_list) < 2:
            continue
        df, HI2_total, HI3_total = compute_streamline_hemolysis(nodes, sid_list, A, alpha, beta)
        df, SA_total = compute_streamline_SA(df)
        per_streamline_dfs.append(df)

        # --- NEW: collect coverage diagnostic for this streamline
        diag = diagnose_streamline_coverage(df, sid_index=i)
        if diag is not None:
            coverage_diags.append(diag)

        results.append({
            "n_points": len(sid_list),
            "HI2_percent": HI2_total,
            "HI3_percent": HI3_total,
            "HI_percent": HI2_total,
            "SA_dyne_s_cm2": SA_total,
            "mean_inlet_velocity": df["vel"].iloc[0],
        })

    summary = pd.DataFrame(results)

    if len(summary) == 0:
        raise RuntimeError("No valid streamlines were reconstructed - check the export file.")

    weights = summary["mean_inlet_velocity"].to_numpy()
    weights = weights / weights.sum()
    device_HI2_percent = float(np.sum(summary["HI2_percent"].to_numpy() * weights))
    device_HI3_percent = float(np.sum(summary["HI3_percent"].to_numpy() * weights))
    device_HI_percent = device_HI2_percent

    def _to_nih_mg(hi_percent):
        nih_g = hi_percent * (1 - HCT) * HB  # g/100L (per report eq. 3 units)
        return nih_g * 1000.0  # convert g/100L -> mg/100L for comparison

    device_NIH_mg = _to_nih_mg(device_HI2_percent)
    device_HI2_NIH_mg = _to_nih_mg(device_HI2_percent)
    device_HI3_NIH_mg = _to_nih_mg(device_HI3_percent)

    sa_values = summary["SA_dyne_s_cm2"].to_numpy()
    device_SA_median = float(np.median(sa_values))
    device_SA_prob_above_threshold = float(
        np.mean(sa_values > HELLUMS_THRESHOLD_DYNE_S_CM2)
    )

    # --- NEW: aggregate coverage diagnostic + report
    coverage_df = pd.DataFrame(coverage_diags)
    _report_coverage(coverage_df)

    return {
        "nodes": nodes,
        "streamlines": per_streamline_dfs,
        "summary": summary,
        "device_HI_percent": device_HI_percent,
        "device_NIH_mg_per_100L": device_NIH_mg,
        "device_HI2_percent": device_HI2_percent,
        "device_HI3_percent": device_HI3_percent,
        "device_HI2_NIH_mg_per_100L": device_HI2_NIH_mg,
        "device_HI3_NIH_mg_per_100L": device_HI3_NIH_mg,
        "device_SA_median_dyne_s_cm2": device_SA_median,
        "device_SA_prob_above_hellums": device_SA_prob_above_threshold,
        "coverage_diagnostics": coverage_df,   # --- NEW
    }


def check_constant_range_coverage(result, constants_name):
    """
    Reports what fraction of streamline data falls outside the tau/t range
    that a given power-law constant set was actually fit over.
    """
    limits = POWER_LAW_CONSTANTS[constants_name]
    all_tau = np.concatenate([df["shear"].to_numpy()[1:] for df in result["streamlines"]])
    all_t = np.concatenate([df["t_cumulative"].to_numpy()[1:] for df in result["streamlines"]])

    frac_tau_over = float(np.mean(all_tau > limits["tau_max_Pa"])) if len(all_tau) else 0.0
    frac_t_over = float(np.mean(all_t > limits["t_max_s"])) if len(all_t) else 0.0

    print(f"[{constants_name}] fitted range: tau < {limits['tau_max_Pa']:.0f} Pa, "
          f"t < {limits['t_max_s']*1000:.0f} ms")
    print(f"  segments with tau over range: {frac_tau_over*100:.1f}%")
    print(f"  segments with t over range:   {frac_t_over*100:.1f}%")
    return frac_tau_over, frac_t_over


def run_pipeline_all_constants(filepath, constant_names=("GW", "HO", "TZ")):
    """
    Runs the full pipeline once per named constant set.
    """
    all_results = {}
    rows = []
    for name in constant_names:
        c = POWER_LAW_CONSTANTS[name]
        res = run_pipeline(filepath, A=c["A"], alpha=c["alpha"], beta=c["beta"])
        all_results[name] = res
        frac_tau_over, frac_t_over = check_constant_range_coverage(res, name)
        rows.append({
            "constants": name,
            "HI2_percent": res["device_HI2_percent"],
            "HI3_percent": res["device_HI3_percent"],
            "HI2_NIH_mg_per_100L": res["device_HI2_NIH_mg_per_100L"],
            "HI3_NIH_mg_per_100L": res["device_HI3_NIH_mg_per_100L"],
            "frac_segments_tau_over_range": frac_tau_over,
            "frac_segments_t_over_range": frac_t_over,
        })

    comparison = pd.DataFrame(rows).set_index("constants")
    print("\n" + "=" * 70)
    print("Constant-set comparison (device-level, flow-weighted):")
    print(comparison.to_string(float_format=lambda x: f"{x:.4g}"))
    print(f"SynCardia 50cc literature NIH: {SYNCARDIA_50CC_NIH} +/- 12.42 mg/100L")
    print(f"ASTM acceptable limit:         {ASTM_LIMIT_NIH} mg/100L")
    print("=" * 70)

    return all_results, comparison


# ---------------------------------------------------------------------------
# 4. Visualization
# ---------------------------------------------------------------------------

def _axis_range(streamline_dfs, col):
    lo = min(df[col].min() for df in streamline_dfs)
    hi = max(df[col].max() for df in streamline_dfs)
    return lo, hi, hi - lo


def plot_streamlines_3d(streamline_dfs, save_path=None, max_lines=None,
                         flat_tolerance=1e-3):
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    x_lo, x_hi, x_range = _axis_range(streamline_dfs, "x")
    y_lo, y_hi, y_range = _axis_range(streamline_dfs, "y")
    z_lo, z_hi, z_range = _axis_range(streamline_dfs, "z")

    all_shear = np.concatenate([df["shear"].to_numpy() for df in streamline_dfs])
    from matplotlib.colors import LogNorm
    vmin, vmax = np.nanmin(all_shear), np.nanmax(all_shear)
    cmap = plt.get_cmap("turbo")
    norm = LogNorm(vmin=vmin, vmax=vmax)

    dfs = streamline_dfs if max_lines is None else streamline_dfs[:max_lines]
    in_plane_scale = max(x_range, y_range, 1e-12)

    is_flat_z = z_range < flat_tolerance * in_plane_scale
    is_flat_y = y_range < flat_tolerance * max(x_range, z_range, 1e-12)

    if is_flat_z or is_flat_y:
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
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    ax, ax2 = axes
    ax.hist(summary["HI2_percent"], bins=30, color="#c0392b", edgecolor="black", alpha=0.85)
    ax.set_xlabel("Streamline HI2% (unweighted)")
    ax.set_ylabel("Number of streamlines")
    ax.set_title("HI2 (temporal-derivative form)")
    ax.grid(alpha=0.3)

    ax2.hist(summary["HI3_percent"], bins=30, color="#2980b9", edgecolor="black", alpha=0.85)
    ax2.set_xlabel("Streamline HI3% (unweighted)")
    ax2.set_title("HI3 (local/linearized form)")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    return fig


def _gaussian_kde_numpy(data, x_grid):
    data = np.asarray(data, dtype=float)
    n = len(data)
    std = np.std(data, ddof=1) if n > 1 else 1.0
    if std <= 0:
        std = 1.0
    bandwidth = 1.06 * std * n ** (-1 / 5)
    bandwidth = max(bandwidth, 1e-6)

    diffs = (x_grid[:, None] - data[None, :]) / bandwidth
    kernel = np.exp(-0.5 * diffs ** 2) / np.sqrt(2 * np.pi)
    pdf = kernel.sum(axis=1) / (n * bandwidth)
    return pdf


def plot_sa_pdf(summary, save_path=None, label=None, threshold=HELLUMS_THRESHOLD_DYNE_S_CM2,
                 ax=None, color=None):
    try:
        from scipy.stats import gaussian_kde
        use_scipy = True
    except ImportError:
        use_scipy = False

    if isinstance(summary, pd.DataFrame):
        summaries = [summary]
        labels = [label if label else "SA distribution"]
    else:
        summaries = summary
        labels = label if label else [f"Case {i+1}" for i in range(len(summaries))]

    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=(8, 5.5))
    else:
        fig = ax.figure

    colors = plt.get_cmap("tab10").colors
    for i, (s, lab) in enumerate(zip(summaries, labels)):
        sa = s["SA_dyne_s_cm2"].to_numpy()
        sa = sa[np.isfinite(sa)]
        if len(sa) < 2:
            continue
        kde = gaussian_kde(sa) if use_scipy else None
        x_grid = np.linspace(0, max(sa.max() * 1.1, threshold * 1.2), 400)
        pdf = kde(x_grid) if use_scipy else _gaussian_kde_numpy(sa, x_grid)
        c = color if (color and len(summaries) == 1) else colors[i % len(colors)]
        ax.plot(x_grid, pdf, color=c, linewidth=2, label=lab)

        median_sa = np.median(sa)
        ax.axvline(median_sa, color=c, linestyle="--", linewidth=1.2, alpha=0.8)
        ax.text(median_sa, ax.get_ylim()[1] if not created_fig else 0, f"{median_sa:.2f}",
                color=c, fontsize=9, ha="center", va="bottom")

    ax.axvline(threshold, color="black", linestyle=":", linewidth=1.5,
               label=f"Hellums threshold ({threshold:.0f} dyne\u00b7s/cm\u00b2)")

    ax.set_xlabel("Stress accumulation, SA [dyne\u00b7s/cm\u00b2]")
    ax.set_ylabel("Probability density function")
    ax.set_title("Stress-accumulation distribution (thrombogenic footprint)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    return fig


def print_summary(result):
    print("=" * 60)
    print(f"Streamlines processed:        {len(result['summary'])}")
    print(f"Constants used:                A={A_COEF:.4g}  alpha={ALPHA:.4g}  beta={BETA:.4g} "
          f"({ACTIVE_CONSTANTS})")
    print("-" * 60)
    print("HI2 (temporal-derivative form, Taskin et al. Eq. 4):")
    print(f"  Device HI2% (flow-weighted): {result['device_HI2_percent']:.6e} %")
    print(f"  Device NIH (from HI2):       {result['device_HI2_NIH_mg_per_100L']:.4f} mg/100L")
    print("HI3 (local/linearized form, Taskin et al. Eq. 5):")
    print(f"  Device HI3% (flow-weighted): {result['device_HI3_percent']:.6e} %")
    print(f"  Device NIH (from HI3):       {result['device_HI3_NIH_mg_per_100L']:.4f} mg/100L")
    print("-" * 60)
    print(f"SynCardia 50cc literature:    {SYNCARDIA_50CC_NIH} +/- 12.42 mg/100L")
    print(f"ASTM acceptable limit:        {ASTM_LIMIT_NIH} mg/100L")
    print("-" * 60)
    print(f"Device SA (median):           {result['device_SA_median_dyne_s_cm2']:.4f} dyne*s/cm^2")
    print(f"Hellums threshold:            {HELLUMS_THRESHOLD_DYNE_S_CM2} dyne*s/cm^2")
    print(f"P(SA > Hellums threshold):    {result['device_SA_prob_above_hellums']*100:.3f} %")
    print(f"Marom et al. medians (ref):   {SA_LITERATURE_MEDIANS}")
    print("=" * 60)

    print("\n=== Streamline Quality Check ===")
    print(f"Total streamlines: {len(result['streamlines'])}")
    print(f"Average points per streamline: {np.mean([len(df) for df in result['streamlines']])}")
    print(f"Min/Max velocity: {result['nodes']['vel'].min():.3f} / {result['nodes']['vel'].max():.3f} m/s")
    print(f"Min/Max shear: {result['nodes']['shear'].min():.1f} / {result['nodes']['shear'].max():.1f} Pa")

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
        print("Usage: python hemolysis_pipeline.py <export.csv> [--compare-constants]")
        sys.exit(1)

    filepath = sys.argv[1]

    if "--compare-constants" in sys.argv:
        all_results, comparison = run_pipeline_all_constants(filepath)
        result = all_results[ACTIVE_CONSTANTS]
    else:
        result = run_pipeline(filepath)

    print_summary(result)

    if GENERATE_OUTPUTS:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        try:
            script_name = Path(__file__).stem
        except NameError:
            script_name = "haemolysis_pipeline_anyversion"

        run_id = (
            f"{datetime.now():%Y-%m-%d_%H-%M-%S}_"
            f"{script_name}_"
            f"{Path(filepath).stem}"
        )

        output_paths = {
            "streamlines": OUTPUT_DIR / f"streamlines_3d_{run_id}.png",
            "shear":       OUTPUT_DIR / f"shear_vs_time_{run_id}.png",
            "histogram":   OUTPUT_DIR / f"hi_histogram_{run_id}.png",
            "sa_pdf":      OUTPUT_DIR / f"sa_pdf_{run_id}.png",
            "summary":     OUTPUT_DIR / f"streamline_summary_{run_id}.csv",
            "coverage":    OUTPUT_DIR / f"coverage_diagnostics_{run_id}.csv",   # --- NEW
        }
        plot_streamlines_3d(result["streamlines"], save_path=output_paths["streamlines"])
        plot_shear_vs_time(result["streamlines"], save_path=output_paths["shear"])
        plot_hi_histogram(result["summary"], save_path=output_paths["histogram"])
        plot_sa_pdf(result["summary"], save_path=output_paths["sa_pdf"], label="This device")

        result["summary"].to_csv(output_paths["summary"], index=False)
        # --- NEW: persist coverage diagnostics
        if len(result["coverage_diagnostics"]) > 0:
            result["coverage_diagnostics"].to_csv(output_paths["coverage"], index=False)

        print("Saved outputs:")
        for output_path in output_paths.values():
            print(f"  {output_path}")
    else:
        print("Output generation disabled (GENERATE_OUTPUTS=False).")