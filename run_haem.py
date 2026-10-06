"""One entry point for selected haem analyses, grouped in a new run folder."""

import argparse
import contextlib
from datetime import datetime
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import haem5


class Tee:
    def __init__(self, terminal, log):
        self.terminal, self.log = terminal, log

    def write(self, text):
        self.terminal.write(text)
        return self.log.write(text)

    def flush(self):
        self.terminal.flush()
        self.log.flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exports", nargs="*", type=Path)
    parser.add_argument("--manifest", type=Path, help="Same case manifest as haem_diagnostics.py")
    parser.add_argument("--analyses", nargs="+", choices=["summary", "diagnostics", "empirical", "plots"],
                        default=["summary", "diagnostics"])
    parser.add_argument("--outdir", type=Path, help="New run folder; existing folders are refused")
    parser.add_argument("--compare-constants", action="store_true", help="Compare GW/HO/TZ in the core analysis")
    parser.add_argument("--bpm", type=float, default=120.)
    parser.add_argument("--expected-seeds", type=int)
    parser.add_argument("--outlet-plane", nargs=6, type=float)
    parser.add_argument("--outlet-tolerance", type=float, default=1e-4)
    args = parser.parse_args()
    cases = [{"path": str(p.resolve())} for p in args.exports]
    if args.manifest:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8-sig"))
        for entry in manifest["cases"]:
            case = dict(entry)
            for key in ("path", "weights"):
                if case.get(key):
                    case[key] = str((args.manifest.parent / case[key]).resolve())
            cases.append(case)
    if not cases:
        parser.error("Provide exports or --manifest.")
    defaults = {"bpm": args.bpm, "expected_seeds": args.expected_seeds,
                "outlet_plane": args.outlet_plane, "outlet_tolerance": args.outlet_tolerance,
                "flow_model": "SST", "stress_model": "molecular viscosity * scalar shear strain rate"}
    cases = [{**defaults, **case} for case in cases]
    for case in cases:
        if not Path(case["path"]).is_file():
            parser.error(f"Input does not exist: {case['path']}")
        if not np.isfinite(case["bpm"]) or case["bpm"] <= 0:
            parser.error("BPM must be finite and positive.")
        if not np.isfinite(case["outlet_tolerance"]) or case["outlet_tolerance"] < 0:
            parser.error("Outlet tolerance must be finite and nonnegative.")
    analyses = set(args.analyses)
    root = args.outdir or Path(__file__).resolve().parent / "outputs" / "runs" / datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
    if root.exists():
        parser.error(f"Run folder already exists; choose a new --outdir: {root}")
    root.mkdir(parents=True)
    provenance = {"analyses": sorted(analyses), "compare_constants": args.compare_constants,
                  "cases": cases, "command": sys.argv, "status": "running"}
    record = root / "run.json"
    record.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    index = ["# Haem analysis run", "", "Configuration and completion status: [run.json](run.json).",
             "", "Core summary/empirical/plots use all reconstructed paths and inlet-speed proxy weighting.",
             "Diagnostics applies the manifest's outlet selection and supplied flux weights, when provided.",
             "These populations and weights may differ; consult the diagnostics report before comparing.", ""]
    core_rows, diagnostic_rows = [], []
    try:
        for number, case in enumerate(cases, 1):
            path = Path(case["path"])
            folder = root / f"{number:02d}_{path.stem}"
            folder.mkdir()
            case["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            index.extend([f"## {case.get('label', path.stem)}", "",
                          f"- [Console log]({folder.name}/console.txt)", ""])
            with (folder / "console.txt").open("w", encoding="utf-8") as log:
                with contextlib.redirect_stdout(Tee(sys.stdout, log)), contextlib.redirect_stderr(Tee(sys.stderr, log)):
                    print(f"Input: {path}\nAnalyses: {', '.join(sorted(analyses))}")
                    result = None
                    if analyses & {"summary", "empirical", "plots"}:
                        if args.compare_constants:
                            results, comparison = haem5.run_pipeline_all_constants(path)
                            comparison.to_csv(folder / "constants.csv")
                            result = results[haem5.ACTIVE_CONSTANTS]
                        else:
                            result = haem5.run_pipeline(path)
                        result["summary"].to_csv(folder / "streamlines.csv", index=False)
                    if "summary" in analyses:
                        haem5.print_summary(result)
                        core_rows.append({"case": case.get("label", path.stem), "input": str(path),
                                          "bpm": case["bpm"], "size_cc": case.get("size_cc"), "phase": case.get("phase"),
                                          "mean_transit_ms": result["summary"].transit_ms.mean(),
                                          "max_stress_Pa": result["summary"].max_stress_Pa.max(),
                                          **{k: v for k, v in result.items() if k.startswith("device_")}})
                    if "diagnostics" in analyses:
                        import haem_diagnostics
                        diagnostic_rows.append(haem_diagnostics.analyze(path, case, folder / "diagnostics"))
                        index.append(f"- [Diagnostics report]({folder.name}/diagnostics/report.md)")
                    if "empirical" in analyses:
                        import haem5_empiricalthresholds as empirical
                        target = folder / "empirical"
                        target.mkdir()
                        lines, device, envelope = empirical.analyze_empirical_thresholds(result)
                        lines.to_csv(target / "streamlines.csv", index=False)
                        device.to_csv(target / "device.csv", index=False)
                        envelope.to_csv(target / "exposure_envelope.csv", index=False)
                        empirical.print_device_summary(device)
                        if "plots" in analyses:
                            plt.close(empirical.plot_empirical_thresholds(envelope, target / "stress_duration.png"))
                    if "plots" in analyses:
                        target = folder / "plots"
                        target.mkdir()
                        plt.close(haem5.plot_streamlines_3d(result["streamlines"], save_path=target / "streamlines_3d.png"))
                        plt.close(haem5.plot_shear_vs_time(result["streamlines"], save_path=target / "shear_vs_time.png"))
                        plt.close(haem5.plot_hi_histogram(result["summary"], save_path=target / "hi_histogram.png"))
                        plt.close(haem5.plot_sa_pdf(result["summary"], save_path=target / "sa_pdf.png", label=path.stem))
            index.append("")
        if core_rows:
            pd.DataFrame(core_rows).to_csv(root / "summary.csv", index=False)
            index.extend(["[Core case summaries](summary.csv)", ""])
        if diagnostic_rows:
            pd.DataFrame(diagnostic_rows).to_csv(root / "comparison.csv", index=False)
            index.extend(["[Diagnostics case comparison](comparison.csv)", ""])
        provenance["status"] = "complete"
    except Exception as exc:
        provenance.update(status="failed", error=str(exc))
        raise
    finally:
        record.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
        (root / "README.md").write_text("\n".join(index), encoding="utf-8")
    print(f"\nRun saved: {root.resolve()}\nStart with README.md; console.txt contains each case's printed output.")


if __name__ == "__main__":
    main()
