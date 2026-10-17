"""Fast, vectorised loader for CFD-Post streamline exports (haem5-compatible).

Re-implements haem5.compute_streamline_hemolysis / compute_streamline_SA with
numpy arrays so every *segment* (not just every path) is available for
spatial, temporal and decomposition analyses. Numerics are identical to haem5:

  dt_i  = |x_i - x_{i-1}| / max(0.5 (v_i + v_{i-1}), v_floor)
  v_floor = max(0.01 * mean(v > 0 on that path), 0.001 m/s)
  stress sampled at the segment END node (tau_i), as in haem5
  HI3 = [ sum A^(1/b) tau_i^(a/b) dt_i ]^b      (Taskin Eq. 5, Garon-Farinas)
  HI2 = sum b A t_i^(b-1) tau_i^a dt_i            (Taskin Eq. 4)
  SA  = 10 * sum tau_i dt_i   [dyne s / cm^2]

verify_against_haem5() checks per-path agreement on a file.
Only numpy / pandas are required (no scipy).
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Same constants as haem5.POWER_LAW_CONSTANTS (copied so this file has no
# side-effect imports; check_constants() asserts they still match).
CONSTANTS = {
    "GW": dict(A=3.62e-5, alpha=2.416, beta=0.785),
    "HO": dict(A=1.80e-6, alpha=1.991, beta=0.765),
    "TZ": dict(A=1.228e-5, alpha=1.9918, beta=0.6606),
}
PA_S_TO_DYNE = 10.0
HCT, HB = 0.30, 150.0


def check_constants():
    import haem5
    for k, c in CONSTANTS.items():
        h = haem5.POWER_LAW_CONSTANTS[k]
        assert np.isclose(h["A"], c["A"]) and np.isclose(h["alpha"], c["alpha"]) \
            and np.isclose(h["beta"], c["beta"]), k


SNAP_PHASES = (16.67, 33.33, 50.0, 66.67, 83.33)   # the five simulated sixths of systole


def phase_from_name(path: Path) -> float | None:
    """Phase in % systole from the file name.

    'haemline_50cc_1667%.csv' -> 16.67   (digits without a point = hundredths)
    '35cc 120BPM 16.67%.csv'  -> 16.67   (explicit decimal point)
    '..._8300%.csv'           -> 83.33   (snapped: within 0.5 of a simulated sixth)
    """
    m = re.search(r"(\d+(?:\.\d+)?)%", path.name)
    if not m:
        return None
    txt = m.group(1)
    val = float(txt) if "." in txt else int(txt) / 100.0
    near = min(SNAP_PHASES, key=lambda p: abs(p - val))
    return near if abs(near - val) <= 0.5 else val


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _read_export(path: Path):
    with open(path, "r") as f:
        lines = f.readlines()
    idx = {}
    for i, l in enumerate(lines):
        s = l.strip()
        if s.startswith("[") and s.endswith("]"):
            idx[s[1:-1]] = i
    if "Data" not in idx or "Lines" not in idx:
        raise ValueError(f"{path.name}: need [Data] and [Lines] sections "
                         "(export streamlines as plain Lines with node numbers)")
    d0, l0 = idx["Data"] + 2, idx["Lines"] + 1  # skip [Data] header row
    data_rows = [r for r in lines[d0:idx["Lines"]] if r.strip()]
    arr = np.loadtxt(data_rows, delimiter=",", ndmin=2)
    end_lines = min([v for v in idx.values() if v > idx["Lines"]], default=len(lines))
    pair_rows = [r for r in lines[l0:end_lines] if r.strip()]
    pairs = np.loadtxt(pair_rows, delimiter=",", dtype=np.int64, ndmin=2)
    if arr.shape[1] == 9:       # node id + 8 values
        nid = arr[:, 0].astype(np.int64)
        vals = arr[:, 1:6]
    else:                       # no node ids -> file order (haem5 behaviour)
        nid = np.arange(len(arr))
        vals = arr[:, 0:5]
    # keep first occurrence of each node id (haem5 drop_duplicates)
    _, first = np.unique(nid, return_index=True)
    first.sort()
    nid, vals = nid[first], vals[first]
    return nid, vals, pairs


def _reconstruct(pairs: np.ndarray):
    """Split consecutive (start, end) pairs into node sequences (haem5 logic)."""
    a, b = pairs[:, 0], pairs[:, 1]
    breaks = np.r_[True, a[1:] != b[:-1]]
    starts = np.flatnonzero(breaks)
    ends = np.r_[starts[1:], len(pairs)]
    seqs = [np.r_[a[s], b[s:e]] for s, e in zip(starts, ends)]
    return seqs


# ---------------------------------------------------------------------------
# Case container
# ---------------------------------------------------------------------------

@dataclass
class Snapshot:
    name: str
    phase: float | None
    # node-level arrays, all paths concatenated
    pid: np.ndarray        # path index of each node
    xyz: np.ndarray        # (n, 3) [m]
    tau: np.ndarray        # scalar stress [Pa]
    vel: np.ndarray        # speed [m/s]
    dt: np.ndarray         # segment duration ending at this node [s] (0 at path start)
    t: np.ndarray          # cumulative time [s]
    floored: np.ndarray    # bool: segment speed was clipped to the floor
    offsets: np.ndarray    # start index of each path (len n_paths + 1)
    paths: pd.DataFrame = field(default=None)

    @property
    def n_paths(self):
        return len(self.offsets) - 1

    def inner(self, const="GW"):
        """Linearised HI3 dose per segment: A^(1/b) tau^(a/b) dt."""
        c = CONSTANTS[const]
        return c["A"] ** (1 / c["beta"]) * self.tau ** (c["alpha"] / c["beta"]) * self.dt

    def subset(self, keep_paths):
        """New Snapshot containing only the selected paths (bool mask)."""
        keep_paths = np.asarray(keep_paths, bool)
        node_keep = keep_paths[self.pid]
        lens = np.diff(self.offsets)[keep_paths]
        new_pid = np.repeat(np.arange(len(lens)), lens)
        sub = Snapshot(self.name, self.phase, new_pid, self.xyz[node_keep], self.tau[node_keep],
                       self.vel[node_keep], self.dt[node_keep], self.t[node_keep],
                       self.floored[node_keep], np.r_[0, np.cumsum(lens)])
        sub.paths = self.paths[keep_paths].reset_index(drop=True)
        return sub

    def time_to_exit(self):
        """Time remaining until the end of each node's path [s]."""
        T = np.repeat(self.paths["T_s"].to_numpy(), np.diff(self.offsets))
        return T - self.t


def build_snapshot(path: Path, floor_scale: float = 1.0) -> Snapshot:
    nid, vals, pairs = _read_export(path)
    lookup = pd.Series(np.arange(len(nid)), index=nid)
    seqs = _reconstruct(pairs)
    rows = []
    for s in seqs:
        s = s[np.isin(s, nid)]
        if len(s) >= 2:
            rows.append(lookup.loc[s].to_numpy())
    lens = np.array([len(r) for r in rows])
    order = np.concatenate(rows)
    offsets = np.r_[0, np.cumsum(lens)]
    pid = np.repeat(np.arange(len(rows)), lens)
    xyz, tau, vel = vals[order, 0:3], vals[order, 3], vals[order, 4]

    first = np.zeros(len(order), bool)
    first[offsets[:-1]] = True
    seg_len = np.zeros(len(order))
    seg_len[1:] = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    seg_len[first] = 0.0
    seg_vel = np.zeros(len(order))
    seg_vel[1:] = 0.5 * (vel[1:] + vel[:-1])
    seg_vel[first] = 0.0

    pos = vel > 0
    sum_pos = np.bincount(pid, weights=np.where(pos, vel, 0.0), minlength=len(rows))
    n_pos = np.bincount(pid, weights=pos.astype(float), minlength=len(rows))
    mean_v = np.where(n_pos > 0, sum_pos / np.maximum(n_pos, 1), 0.1)
    vfloor = np.maximum(0.01 * floor_scale * mean_v, 0.001 * floor_scale)[pid]
    floored = (seg_vel < vfloor) & ~first
    seg_vel = np.maximum(seg_vel, vfloor)
    dt = np.where(first, 0.0, seg_len / seg_vel)
    csum = np.cumsum(dt)
    t = csum - np.repeat(csum[offsets[:-1]], lens)

    snap = Snapshot(path.stem, phase_from_name(path), pid, xyz, tau, vel, dt, t,
                    floored, offsets)
    snap.paths = _path_table(snap)
    return snap


def _path_table(s: Snapshot) -> pd.DataFrame:
    n = s.n_paths
    bc = lambda w: np.bincount(s.pid, weights=w, minlength=n)
    T = bc(s.dt)
    out = {"T_s": T, "transit_ms": T * 1e3}
    for k, c in CONSTANTS.items():
        lin = bc(s.inner(k))
        out[f"{k}_lin"] = lin
        out[f"{k}_HI3"] = lin ** c["beta"]
        # effective stress: HI3 = A T^b tau_eff^a exactly
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"{k}_tau_eff"] = (lin / (c["A"] ** (1 / c["beta"]) * T)) ** (c["beta"] / c["alpha"])
    c = CONSTANTS["GW"]
    with np.errstate(divide="ignore", invalid="ignore"):
        tp = np.where(s.t > 0, s.t ** (c["beta"] - 1), 0.0)
    d2 = np.nan_to_num(c["beta"] * c["A"] * tp * s.tau ** c["alpha"] * s.dt)
    out["GW_HI2"] = bc(d2)
    out["SA"] = bc(s.tau * s.dt) * PA_S_TO_DYNE
    out["tau_mean_time"] = bc(s.tau * s.dt) / np.where(T > 0, T, np.nan)
    out["tau_peak"] = np.maximum.reduceat(s.tau, s.offsets[:-1])
    out["v_inlet"] = s.vel[s.offsets[:-1]]
    out["n_points"] = np.diff(s.offsets)
    out["floored_time_frac"] = bc(np.where(s.floored, s.dt, 0.0)) / np.where(T > 0, T, np.nan)
    last = s.offsets[1:] - 1
    for j, ax in enumerate("xyz"):
        out[f"start_{ax}"] = s.xyz[s.offsets[:-1], j]
        out[f"end_{ax}"] = s.xyz[last, j]
    df = pd.DataFrame(out)
    df["end_plane_dist_mm"], df["completed"] = classify_completion(df)
    return df


def classify_completion(paths: pd.DataFrame, tol_mm: float = 1.0):
    """Flag paths whose end point lies on the (planar) outlet.

    Fits a plane to the end points iteratively (start from the 70 % closest to
    the median, refit on the 70 % closest to the plane). CFD-Post streamlines
    that reach the outlet stop exactly on it; paths that stop inside the
    domain (wall contact, stagnation, step/time limits) lie off the plane.
    """
    E = paths[["end_x", "end_y", "end_z"]].to_numpy() * 1e3
    d = np.linalg.norm(E - np.median(E, axis=0), axis=1)
    sel = d <= np.percentile(d, 70)
    for _ in range(4):
        c = E[sel].mean(axis=0)
        n = np.linalg.svd(E[sel] - c)[2][2]
        dp = np.abs((E - c) @ n)
        sel = dp <= np.percentile(dp, 70)
    return dp, dp <= tol_mm


def file_hash(path: Path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def load_snapshot(path: Path, cache_dir: Path | None = None, floor_scale=1.0) -> Snapshot:
    """Build (or load cached) snapshot. Cache key = file hash + floor scale."""
    path = Path(path)
    if cache_dir is None:
        return build_snapshot(path, floor_scale)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = f"{path.stem}_{file_hash(path)}_f{floor_scale:g}.npz"
    cp = cache_dir / key
    if cp.exists():
        z = np.load(cp, allow_pickle=False)
        phase = phase_from_name(path)          # always from the name (parser may change)
        s = Snapshot(path.stem, phase, z["pid"], z["xyz"], z["tau"], z["vel"],
                     z["dt"], z["t"], z["floored"], z["offsets"])
        s.paths = _path_table(s)
        return s
    s = build_snapshot(path, floor_scale)
    np.savez_compressed(cp, pid=s.pid, xyz=s.xyz, tau=s.tau, vel=s.vel, dt=s.dt,
                        t=s.t, floored=s.floored, offsets=s.offsets,
                        phase=np.nan if s.phase is None else s.phase)
    return s


def verify_against_haem5(path: Path, n_check=None):
    """Compare per-path HI2/HI3/SA/transit with haem5.run_pipeline (GW)."""
    import haem5
    ref = haem5.run_pipeline(str(path))["summary"]
    s = build_snapshot(Path(path))
    p = s.paths
    assert len(ref) == len(p), (len(ref), len(p))
    rel = lambda a, b: np.nanmax(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))
    return {
        "n_paths": len(p),
        "max_rel_err_HI3": rel(p["GW_HI3"].to_numpy(), ref["HI3_percent"].to_numpy()),
        "max_rel_err_HI2": rel(p["GW_HI2"].to_numpy(), ref["HI2_percent"].to_numpy()),
        "max_rel_err_SA": rel(p["SA"].to_numpy(), ref["SA_dyne_s_cm2"].to_numpy()),
        "max_rel_err_transit": rel(p["transit_ms"].to_numpy(), ref["transit_ms"].to_numpy()),
    }
