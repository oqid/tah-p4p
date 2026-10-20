#!/usr/bin/env python3
"""
cfx_out_summary.py - collate ANSYS CFX solver settings and run outcomes from .out files.

Walks one or more folders (e.g. Workbench "<project>_files" folders, or the folders that
contain them), finds every CFX solver output file (*.out), and writes one row per run to a
CSV: case, mesh size, material, boundary conditions, solver controls, iterations,
final residuals, termination reason, mass flow through each boundary and timing.

Standard library only. Usage (Windows example):
    python cfx_out_summary.py "E:\\P4P TORTURE\\50CC SIMS" "E:\\35cc Retry" -o sim_settings.csv
Options:
    --latest-only   keep only the newest .out in each <project>_files folder (by solver start time)
"""
import argparse
import csv
import os
import re
import sys
from datetime import datetime

KG_S_TO_L_MIN = None  # set from density per file


def last(pattern, text, flags=0, group=1):
    m = None
    for m in re.finditer(pattern, text, flags):
        pass
    return m.group(group).strip() if m else ""


def first(pattern, text, flags=0, group=1):
    m = re.search(pattern, text, flags)
    return m.group(group).strip() if m else ""


def ccl_block(text):
    """Return the CCL echoed for the (last) run in the file."""
    starts = [m.end() for m in re.finditer(r"CFX Command Language for Run.*?\n", text)]
    if not starts:
        return ""
    s = starts[-1]
    e = text.find("+-----", s + 200)
    return text[s:e if e > 0 else s + 60000]


def boundaries(ccl):
    """Parse BOUNDARY blocks -> list of dicts."""
    out = []
    dom = ccl.split("DOMAIN MODELS:")[0]
    for m in re.finditer(r"BOUNDARY: (.+?)\n(.*?)(?=BOUNDARY: |\Z)", dom, re.S):
        name, body = m.group(1).strip(), m.group(2)
        body = re.sub(r"\\\n\s*", "", body)  # join continued lines
        b = {"name": name,
             "type": first(r"Boundary Type = (\S+)", body),
             "location": first(r"Location = (.+)", body)}
        for key in ["Normal Speed", "Relative Pressure", "Mass Flow Rate", "Option"]:
            vals = re.findall(rf"{key} = (.+)", body)
            if vals:
                b[key] = "; ".join(v.strip() for v in vals) if key == "Option" else vals[0].strip()
        out.append(b)
    return out


def boundary_massflows(text):
    m = None
    for m in re.finditer(r"\|\s+P-Mass\s+\|\s*\n\s*\+-+\+\s*\n(.*?)Domain Imbalance :\s+(\S+)", text, re.S):
        pass
    if not m:
        return {}, ""
    flows = dict((n.strip(), float(v)) for n, v in re.findall(r"Boundary\s+:\s+(.+?)\s+(-?\d\.\d+E[+-]\d+)", m.group(1)))
    return flows, m.group(2)


def parse(path):
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    ccl = ccl_block(text)
    r = {"out_file": path}
    # identify case / workbench project
    parts = re.split(r"[\\/]", path)
    proj = next((p for p in parts if p.endswith("_files")), "")
    r["wb_project"] = proj[:-6] if proj else ""
    r["system"] = next((p for p in parts if re.match(r"CFX(-\d+)?$", p)), "")
    r["run_name"] = first(r'-fullname "([^"]+)"', text)
    r["solver_version"] = first(r"This run of the CFX (.+?) Solver", text)
    started = last(r"Solver started at (\d\d:\d\d:\d\d on \d\d \w{3}\s*\n?\s*\d{4})", text)
    r["started"] = re.sub(r"\s+", " ", started)
    r["user_host"] = last(r"by user (\S+ on \S+)", text)
    r["analysis_type"] = first(r"ANALYSIS TYPE:\s*\n\s*Option = (.+)", ccl)
    r["material"] = first(r"MATERIAL: (.+)", ccl)
    r["density"] = first(r"Density = (.+)", ccl)
    r["viscosity"] = first(r"Dynamic Viscosity = (.+)", ccl)
    r["turbulence_model"] = first(r"TURBULENCE MODEL:\s*\n\s*Option = (.+)", ccl)
    r["wall_function"] = first(r"TURBULENT WALL FUNCTIONS:\s*\n\s*Option = (.+)", ccl)
    r["advection_scheme"] = first(r"ADVECTION SCHEME:\s*\n\s*Option = (.+)", ccl)
    r["turbulence_numerics"] = first(r"Turbulence Numerics = (.+)", ccl)
    r["timescale_control"] = first(r"Timescale Control = (.+)", ccl)
    r["physical_timescale"] = first(r"Physical Timescale = (.+)", ccl)
    r["max_iterations"] = first(r"Maximum Number of Iterations = (.+)", ccl)
    r["residual_target"] = first(r"Residual Target = (.+)", ccl)
    r["residual_type"] = first(r"Residual Type = (.+)", ccl)
    r["double_precision"] = first(r"Double Precision = (.+)", ccl)
    r["partitions"] = first(r"Maximum Number of Processes = (.+)", ccl)
    r["nodes"] = last(r"Total Number of Nodes\s+=\s+(\d+)", text)
    r["elements"] = last(r"Total Number of Elements\s+=\s+(\d+)", text)
    for b in boundaries(ccl):
        tag = re.sub(r"\W+", "_", b["name"]).strip("_")
        val = b.get("Normal Speed") or b.get("Relative Pressure") or b.get("Mass Flow Rate") or ""
        r[f"BC_{tag}"] = f"{b['type']} @ {b['location']}" + (f" | {val}" if val else "")
    # outcome
    r["iterations"] = last(r"OUTER LOOP ITERATION =\s+(\d+)", text)
    for eq in ["U-Mom", "V-Mom", "W-Mom", "P-Mass"]:
        r[f"final_rms_{eq}"] = last(rf"\|\s*{eq}\s+\|\s*[\d.]+\s*\|\s*(\S+)\s*\|", text)
    term = first(r"Termination and Interrupt Condition Summary\s*=+\s*\n(.*?)\n\s*=+", text, re.S)
    r["termination"] = re.sub(r"\s+", " ", term).strip() or ("NOT FINISHED / crashed" if "has finished" not in text else "")
    flows, imb = boundary_massflows(text)
    rho = float(first(r"([\d.]+)", r["density"]) or 1060)
    for n, v in flows.items():
        tag = re.sub(r"\W+", "_", n).strip("_")
        r[f"massflow_{tag}_kg_s"] = v
        r[f"flow_{tag}_L_min"] = round(abs(v) / rho * 60000, 3)
    r["mass_imbalance"] = imb
    r["wall_clock_s"] = last(r"Total wall clock time: (\S+) seconds", text)
    r["finished_ok"] = "has finished" in text[-3000:]
    return r


def start_key(r):
    try:
        return datetime.strptime(r["started"], "%H:%M:%S on %d %b %Y")
    except ValueError:
        return datetime.min


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+")
    ap.add_argument("-o", "--output", default="cfx_sim_settings.csv")
    ap.add_argument("--latest-only", action="store_true")
    a = ap.parse_args()
    rows = []
    for root in a.roots:
        for dp, _, files in os.walk(root):
            for f in files:
                if f.lower().endswith(".out"):
                    p = os.path.join(dp, f)
                    try:
                        rows.append(parse(p))
                    except Exception as e:  # keep going
                        print(f"!! {p}: {e}", file=sys.stderr)
    if a.latest_only:
        best = {}
        for r in rows:
            k = r["wb_project"] or os.path.dirname(r["out_file"])
            if k not in best or start_key(r) > start_key(best[k]):
                best[k] = r
        rows = list(best.values())
    rows.sort(key=lambda r: (r["wb_project"], start_key(r)))
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} runs -> {a.output}")


if __name__ == "__main__":
    main()
