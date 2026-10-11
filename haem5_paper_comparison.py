"""Paper-inspired snapshot comparisons; equal-path HI is the primary result.

python haem5_paper_comparison.py z_inputs --label "50 cc" --bpm 120
python haem5_paper_comparison.py z_inputs --comparison-folder z_inputs_35cc --comparison-label "35 cc" --bpm 120 --comparison-bpm 120
python haem5_paper_comparison.py --replot outputs_paper_comparison/RUN

Taskin et al. (2012), Figures 4, 6, 7, 9: path histories and model sensitivity.
Marom et al. (2014), Figures 4-6: spatial stress, SA distributions and their tails.
"""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import numpy as np
import pandas as pd

import haem5
from haem5_snapshots import natural_key, phase_from_name

COLORS = ["#548BB0", "#8E7AB4", "#BA76A4", "#D58F9D", "#C89161"]
MODELS = ("GW", "HO", "TZ")
BANDS = np.array([0., 1., 10., 50., 100., np.inf])
BAND_LABELS = ["<1 Pa", "1-10 Pa", "10-50 Pa", "50-100 Pa", ">=100 Pa"]

def coefficients_on_path(stress, dt, constants):
    """Same right-endpoint quadrature and exponent convention as haem5."""
    stress, dt = np.asarray(stress, float), np.asarray(dt, float)
    if (not np.isfinite(stress).all() or not np.isfinite(dt).all()
            or np.any(stress < 0) or np.any(dt < 0)):
        raise ValueError("Stress and segment duration must be finite and nonnegative.")
    A, alpha, beta = (constants[k] for k in ("A", "alpha", "beta"))
    t = np.cumsum(dt)
    power = np.zeros_like(t)
    positive = t > 0
    power[positive] = t[positive] ** (beta - 1)
    hi2 = np.sum(beta * A * power * stress**alpha * dt)
    dose = A**(1/beta) * stress**(alpha/beta) * dt
    return float(hi2), float(dose.sum()**beta), dose

def aggregate_hi(hi, speed):
    hi, speed = np.asarray(hi, float), np.asarray(speed, float)
    if not len(hi) or not np.isfinite(hi).all() or np.any(hi < 0):
        raise ValueError("HI must be finite and nonnegative.")
    if not np.isfinite(speed).all() or np.any(speed < 0):
        raise ValueError("Inlet speeds must be finite and nonnegative.")
    weighted = float(np.average(hi, weights=speed)) if speed.sum() > 0 else np.nan
    return float(hi.mean()), weighted

def reflection_kde(values, grid, bandwidth):
    """Ordinary SA density on a linear axis, with reflection at zero."""
    values = np.asarray(values)
    result = np.zeros_like(grid)
    for chunk in np.array_split(values, max(1, int(np.ceil(len(values)/500)))):
        d1 = (grid[:, None] - chunk[None, :])/bandwidth
        d2 = (grid[:, None] + chunk[None, :])/bandwidth
        result += (np.exp(-.5*d1*d1) + np.exp(-.5*d2*d2)).sum(axis=1)
    return result/(len(values)*bandwidth*np.sqrt(2*np.pi))

def top_contribution(values):
    values = np.sort(np.asarray(values))[::-1]
    return np.r_[0, np.arange(1, len(values)+1)/len(values)*100], (
        np.r_[0, np.cumsum(values)/values.sum()*100] if values.sum() > 0 else np.zeros(len(values)+1))

def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def extract_case(path, meta, directory):
    gw = haem5.POWER_LAW_CONSTANTS["GW"]
    result = haem5.run_pipeline(path, A=gw["A"], alpha=gw["alpha"], beta=gw["beta"])
    summary = result["summary"].copy()
    paths = result["streamlines"]
    extra = []
    durations, doses = np.zeros(5), np.zeros(5)
    for df in paths:
        tau, dt = df.shear.to_numpy(), df.dt.to_numpy()
        transit = dt.sum()
        if transit <= 0:
            raise ValueError("Zero-duration reconstructed path in "+str(path))
        row = {"mean_stress_Pa": float(np.sum(tau*dt)/transit),
               "path_mean_speed_m_s": float(df.arc_length.sum()/transit),
               "path_length_mm": float(df.arc_length.sum()*1000),
               "peak_time_fraction": float(df.t_cumulative.iloc[int(np.argmax(tau))]/transit)}
        for model in MODELS:
            hi2, hi3, dose = coefficients_on_path(tau, dt, haem5.POWER_LAW_CONSTANTS[model])
            row[model+"_HI2_percent"], row[model+"_HI3_percent"] = hi2, hi3
            if model == "GW":
                for b in range(5):
                    mask = (tau >= BANDS[b]) & (tau < BANDS[b+1])
                    durations[b] += dt[mask].sum()
                    doses[b] += dose[mask].sum()
        extra.append(row)
    data = pd.concat([summary.reset_index(drop=True), pd.DataFrame(extra)], axis=1)
    for model in ("HI2", "HI3"):
        np.testing.assert_allclose(data["GW_"+model+"_percent"], data[model+"_percent"], rtol=1e-10, atol=1e-15)
    if not np.isfinite(data.to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite path metrics in "+str(path))
    data.to_csv(directory/"paths.csv", index=False)
    # Select actual paths by HI3 rank; never construct a synthetic mean trajectory.
    order = np.argsort(data.GW_HI3_percent.to_numpy(), kind="stable")
    arrays = {}
    for q in (50, 95, 99):
        index = int(order[round(q/100*(len(order)-1))])
        df = paths[index]
        arrays["p"+str(q)] = df[["t_cumulative", "shear", "HI3_percent_cumulative"]].to_numpy()
        meta["path_id_p"+str(q)] = index
    rng = np.random.default_rng(20261010)
    selected = np.sort(rng.choice(len(paths), min(200, len(paths)), replace=False))
    clouds = []
    for index in selected:
        df = paths[index]
        # Every vertex from the selected paths, not a CFD field interpolation.
        clouds.append(df[["x", "y", "z", "shear"]].to_numpy())
    arrays["cloud"] = np.concatenate(clouds)
    arrays["time_bands_ms"] = durations/len(paths)*1000
    arrays["dose_bands_percent"] = doses/doses.sum()*100 if doses.sum() else np.zeros(5)
    np.savez_compressed(directory/"plot_data.npz", **arrays)
    meta.update(n_paths=len(paths), sampled_spatial_paths=len(selected),
                input=str(path.resolve()), sha256=fingerprint(path))
    del result
    return data, arrays

def load_cases(record, base):
    cases = []
    phases = sorted(set(c["phase"] for c in record["cases"]))
    for meta in record["cases"]:
        directory = base/meta["id"]
        case = dict(meta)
        case["data"] = pd.read_csv(directory/"paths.csv")
        with np.load(directory/"plot_data.npz") as arrays:
            case["arrays"] = {k: arrays[k] for k in arrays.files}
        case["color"] = COLORS[phases.index(case["phase"]) % len(COLORS)]
        case["style"] = "-" if case["group"] == 0 else "--"
        case["display"] = f'{case["label"]} | {case["phase"]:g}%'
        cases.append(case)
    return cases

def finish(fig, output, name, description, gallery):
    fig.savefig(output/(name+".png"), dpi=170, bbox_inches="tight")
    plt.close(fig)
    gallery.append((name, description))

def tidy(ax, xlabel=None, ylabel=None):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(alpha=.2)

def lines_by_group(cases, ax, key, statistic="mean"):
    for group in sorted(set(c["group"] for c in cases)):
        subset = sorted([c for c in cases if c["group"] == group], key=lambda c:c["phase"])
        values = [getattr(c["data"][key], statistic)() for c in subset]
        ax.plot([c["phase"] for c in subset], values, marker="o",
                color=COLORS[group*2], ls="-" if group == 0 else "--", label=subset[0]["label"])

def sa_figures(cases, output, gallery):
    values = np.concatenate([c["data"].SA_dyne_s_cm2 for c in cases])
    bulk_end = max(1., float(np.quantile(values, .99)))
    end = max(float(values.max())*1.05, bulk_end*1.1)
    # A common bandwidth prevents separate smoothing choices from creating apparent differences.
    bandwidth = max(1.06*np.std(values)*len(values)**(-.2), .01)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    for ax, limit in zip(axes, [bulk_end, end]):
        grid = np.linspace(0, limit, 700)
        for c in cases:
            v = c["data"].SA_dyne_s_cm2.to_numpy()
            ax.plot(grid, reflection_kde(v, grid, bandwidth), color=c["color"], ls=c["style"], label=c["display"]+f" | median {np.median(v):.2f}")
        ax.set_xlim(0, limit)
        ax.set_ylim(bottom=0)
        tidy(ax, "Stress accumulation [dyne s/cm²]", "Probability density [per dyne s/cm²]")
    for c in cases:
        axes[0].axvline(c["data"].SA_dyne_s_cm2.median(), color=c["color"], alpha=.6, ls=":", lw=.8)
    axes[0].set_title("Bulk to pooled P99 | dotted lines: medians")
    axes[1].set_title("Tail | pooled P99 to observed maximum")
    tail_grid = np.linspace(bulk_end, end, 200)
    tail_max = max(reflection_kde(c["data"].SA_dyne_s_cm2.to_numpy(), tail_grid, bandwidth).max() for c in cases)
    axes[1].set_xlim(bulk_end, end)
    axes[1].set_ylim(0, max(tail_max*1.15, 1e-4))
    axes[0].legend(fontsize=7, ncol=1 if len(cases)<7 else 2)
    fig.suptitle("Marom Fig. 5 inspired | ordinary linear axes and a shared smoothing bandwidth")
    fig.supxlabel("Equal-path population. Density is not a percentage; the area under a curve gives probability.", fontsize=9)
    finish(fig, output, "01_sa_linear_density", "Compare the bulk and tail of SA, using Marom's linear-axis distribution approach.", gallery)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    for c in cases:
        v = c["data"].SA_dyne_s_cm2.to_numpy()
        for ax, limit in zip(axes, [bulk_end, max(35., end)]):
            grid = np.linspace(0, limit, 600)
            sorted_v = np.sort(v)
            fraction = (len(v)-np.searchsorted(sorted_v, grid, side="right"))/len(v)*100
            ax.step(grid, fraction, where="post", color=c["color"], ls=c["style"], label=c["display"])
    axes[0].set(xlim=(0, bulk_end), ylim=(0, 100), title="How many paths exceed this SA?")
    axes[1].set(xlim=(0, max(35., end)), ylim=(0, 5), title="Tail zoom | 0-5% of paths")
    axes[1].axvline(10, color=".5", ls=":", lw=1)
    axes[1].axvline(35, color=".5", ls=":", lw=1)
    for ax in axes:
        tidy(ax, "SA cutoff [dyne s/cm²]", "Paths ABOVE the cutoff [%]")
    axes[0].legend(fontsize=7, ncol=1 if len(cases)<7 else 2)
    fig.suptitle("SA exceedance | direct percentages, no density or logarithmic axis")
    fig.supxlabel("10: tail cutoff discussed by Marom; 35: historical SA reference, not a universal pass/fail safety test.", fontsize=9)
    finish(fig, output, "02_sa_exceedance_percent", "Read the fraction of paths above an SA cutoff directly; inspect the highest-exposure minority.", gallery)

def model_figures(cases, output, gallery):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharey=True, layout="constrained")
    for ax, model in zip(axes, MODELS):
        for group in sorted(set(c["group"] for c in cases)):
            cs = sorted([c for c in cases if c["group"] == group], key=lambda c:c["phase"])
            for method, marker, ls in (("HI2","o","--"), ("HI3","s","-")):
                ax.plot([c["phase"] for c in cs], [c["data"][model+"_"+method+"_percent"].mean() for c in cs],
                        color=COLORS[group*2], marker=marker, ls=ls, label=cs[0]["label"]+" "+method)
        ax.set_title(model+" coefficients")
        ax.set_yscale("log")
        tidy(ax, "Systolic phase [%]", "Equal-path mean HI [%]")
        ax.legend(fontsize=8)
    fig.suptitle("Taskin Figs. 6 and 9 inspired | does the comparison depend on the damage model?")
    fig.supxlabel("Same trajectories for every model. Logarithmic HI axis; no experimental validation or Eulerian solution implied.", fontsize=9)
    finish(fig, output, "03_taskin_model_sensitivity", "Check whether model/constant choices preserve or change phase and device rankings.", gallery)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    for ax, method in zip(axes, ("HI2", "HI3")):
        for group in sorted(set(c["group"] for c in cases)):
            cs = sorted([c for c in cases if c["group"] == group], key=lambda c:c["phase"])
            ratios = []
            for c in cs:
                mean, weighted = aggregate_hi(c["data"]["GW_"+method+"_percent"], c["data"].mean_inlet_velocity)
                ratios.append(weighted/mean if mean else np.nan)
            ax.plot([c["phase"] for c in cs], ratios, "o-", color=COLORS[group*2], label=cs[0]["label"])
        ax.axhline(1, color=".4", ls="--", lw=1)
        ax.set_title(method+" | GW coefficients")
        tidy(ax, "Systolic phase [%]", "Inlet-speed-weighted mean / equal-path mean")
        ax.legend(fontsize=8)
    fig.suptitle("Weighting sensitivity only | the main figures use equal-path HI")
    fig.supxlabel("1 = no change. This compares two population summaries; individual path HI is unchanged.", fontsize=9)
    finish(fig, output, "09_weighting_sensitivity", "Quantify how much the inlet-speed proxy changes the reported mean, without making it the primary estimate.", gallery)

def relationships(cases, output, gallery):
    fig, axes = plt.subplots(len(cases), 3, figsize=(13, 2.5*len(cases)), squeeze=False, layout="constrained")
    predictors = [("mean_inlet_velocity","Inlet speed [m/s]"), ("path_mean_speed_m_s","Path distance / transit time [m/s]"), ("transit_ms","Transit time [ms]")]
    max_hi = max(c["data"].GW_HI3_percent.max() for c in cases)
    min_hi = min(c["data"].GW_HI3_percent[c["data"].GW_HI3_percent > 0].min() for c in cases)
    for c, row in zip(cases, axes):
        d = c["data"]
        for ax, (key, label) in zip(row, predictors):
            ax.scatter(d[key], d.GW_HI3_percent, s=2, color=c["color"], alpha=.13, rasterized=True)
            valid = (d[key] > 0) & (d.GW_HI3_percent > 0)
            rho = d.loc[valid, key].rank().corr(d.loc[valid, "GW_HI3_percent"].rank())
            ax.set(yscale="log", ylim=(min_hi/1.2, max_hi*1.2), title=c["display"]+f" | rank r = {rho:.2f}")
            tidy(ax, label, "Individual path HI3 [%]")
            ax.xaxis.set_major_locator(MaxNLocator(4))
    fig.suptitle("Velocity, residence time and damage | all paths; GW HI3; x-ranges vary by row; correlations are descriptive", fontsize=11)
    finish(fig, output, "04_velocity_residence_hi", "Distinguish inlet-speed weighting from the physical speed/residence-time relationship on individual paths.", gallery)
    cols = min(3, len(cases))
    fig, axes = plt.subplots(int(np.ceil(len(cases)/cols)), cols, figsize=(4.5*cols, 4*int(np.ceil(len(cases)/cols))), layout="constrained", squeeze=False)
    norm = LogNorm(min_hi, max_hi)
    for ax, c in zip(axes.flat, cases):
        d = c["data"]
        order = np.argsort(d.GW_HI3_percent)
        artist = ax.scatter(d.transit_ms.iloc[order], d.mean_stress_Pa.iloc[order],
                            c=d.GW_HI3_percent.iloc[order], cmap="magma", norm=norm, s=4, alpha=.5, rasterized=True)
        ax.set(xscale="log", yscale="log", title=c["display"])
        tidy(ax, "Transit time [ms, log scale]", "Time-averaged stress [Pa, log scale]")
    for ax in list(axes.flat)[len(cases):]:
        ax.set_visible(False)
    xlim=(min(c["data"].transit_ms.min() for c in cases)*.9, max(c["data"].transit_ms.max() for c in cases)*1.1)
    ylim=(min(c["data"].mean_stress_Pa[c["data"].mean_stress_Pa>0].min() for c in cases)*.9, max(c["data"].mean_stress_Pa.max() for c in cases)*1.1)
    for ax in list(axes.flat)[:len(cases)]:
        ax.set(xlim=xlim, ylim=ylim)
    fig.colorbar(artist, ax=list(axes.flat)[:len(cases)], label="Individual path HI3 [%] | shared logarithmic colour scale", shrink=.7)
    fig.suptitle("Stress-duration mechanism map | does damage arise from long exposure, strong stress, or both?")
    finish(fig, output, "05_stress_duration_map", "Separate residence-time effects from stress magnitude; colour shows path damage.", gallery)

def histories_and_space(cases, output, gallery):
    fig, axes = plt.subplots(len(cases), 2, figsize=(12, 2.7*len(cases)), layout="constrained", squeeze=False)
    max_time = max(c["arrays"]["p"+str(q)][-1,0]*1000 for c in cases for q in (50,95,99))
    max_stress = max(c["arrays"]["p"+str(q)][:,1].max() for c in cases for q in (50,95,99))
    max_hi = max(c["arrays"]["p"+str(q)][-1,2] for c in cases for q in (50,95,99))
    for c, row in zip(cases, axes):
        for q, color in zip((50,95,99), (COLORS[0],COLORS[1],COLORS[2])):
            a = c["arrays"]["p"+str(q)]
            label=f"HI3 rank P{q} | path {c['path_id_p'+str(q)]}"
            row[0].plot(a[:,0]*1000,a[:,1],color=color,lw=1,label=label)
            row[1].plot(a[:,0]*1000,a[:,2],color=color,lw=1)
        row[0].set(title=c["display"], xlim=(0,max_time), ylim=(0,max_stress*1.05))
        row[1].set(title=c["display"], xlim=(0,max_time), ylim=(0,max_hi*1.05))
        row[0].legend(fontsize=7)
        tidy(row[0], "Time along path [ms]", "Local shear stress [Pa]")
        tidy(row[1], "Time along path [ms]", "Cumulative path HI3 [%]")
    fig.suptitle("Taskin Figs. 4 and 7 inspired | actual median-, P95- and P99-damage paths; GW coefficients")
    fig.supxlabel("Paths selected by final HI3 rank within each snapshot. These are examples, not population percentile curves.", fontsize=9)
    finish(fig, output, "06_representative_path_histories", "Locate the stress bursts that produce typical versus high path damage, with objective rank-based selection.", gallery)
    cloud = np.concatenate([c["arrays"]["cloud"] for c in cases])
    # The largest two spatial extents give an informative, shared projection.
    dims = np.sort(np.argsort(np.ptp(cloud[:,:3],axis=0))[-2:])
    positive = cloud[:,3][cloud[:,3] > 0]
    norm = LogNorm(max(positive.min(), .001), positive.max())
    cols = min(3,len(cases))
    fig, axes=plt.subplots(int(np.ceil(len(cases)/cols)),cols,figsize=(4.6*cols,4.2*int(np.ceil(len(cases)/cols))),layout="constrained",squeeze=False)
    for ax,c in zip(axes.flat,cases):
        a=c["arrays"]["cloud"]
        order=np.argsort(a[:,3])
        artist=ax.scatter(a[order,dims[0]]*1000,a[order,dims[1]]*1000,
                          c=np.maximum(a[order,3],norm.vmin),cmap="magma",norm=norm,s=1,alpha=.5,rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(c["display"])
        ax.set_xlim(cloud[:,dims[0]].min()*1000-2,cloud[:,dims[0]].max()*1000+2)
        ax.set_ylim(cloud[:,dims[1]].min()*1000-2,cloud[:,dims[1]].max()*1000+2)
        tidy(ax, "XYZ"[dims[0]]+" [mm]", "XYZ"[dims[1]]+" [mm]")
    for ax in list(axes.flat)[len(cases):]:
        ax.set_visible(False)
    fig.colorbar(artist,ax=list(axes.flat)[:len(cases)],label="Local shear [Pa] | logarithmic colour scale",shrink=.7)
    fig.suptitle("Taskin / Marom inspired spatial stress view | 200 reproducibly sampled paths per snapshot")
    fig.supxlabel("Projected path vertices, not an Eulerian contour or a transient particle-position snapshot. No anatomical alignment inferred.",fontsize=9)
    finish(fig,output,"07_spatial_path_stress","Find sampled high-stress regions; keep equal spatial scales when comparing sizes.",gallery)

def concentration_and_bands(cases,output,gallery):
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
    for c in cases:
        x,y=top_contribution(c["data"].GW_HI3_percent)
        for ax in axes:
            ax.plot(x,y,color=c["color"],ls=c["style"],label=c["display"])
    axes[0].plot([0,100],[0,100],color=".5",ls=":",label="Equal contribution")
    axes[0].set(xlim=(0,100),ylim=(0,100),title="Complete population")
    axes[1].set(xlim=(0,10),ylim=(0,100),title="Most damaging 10% of paths")
    for ax in axes:
        tidy(ax,"Paths included, ranked from highest HI3 [%]","Share of summed path HI3 [%]")
    axes[0].legend(fontsize=7,ncol=1 if len(cases)<7 else 2)
    fig.suptitle("Damage concentration | does a small minority dominate the equal-path mean?")
    finish(fig,output,"08_damage_concentration","Quantify the contribution of high-damage paths; this is a share of modelled path HI, not measured haemoglobin mass.",gallery)
    fig,axes=plt.subplots(1,2,figsize=(12,max(4.5,len(cases)*.55)),layout="constrained")
    y=np.arange(len(cases))
    for ax,key in zip(axes,["time_bands_ms","dose_bands_percent"]):
        left=np.zeros(len(cases))
        for b in range(5):
            values=np.array([c["arrays"][key][b] for c in cases])
            ax.barh(y,values,left=left,color=COLORS[b],label=BAND_LABELS[b])
            left+=values
        ax.set_yticks(y,[c["display"] for c in cases])
        ax.invert_yaxis()
        ax.grid(axis="x",alpha=.2)
    axes[0].set(xlabel="Mean time per path in stress band [ms]",title="Where is time spent?")
    axes[1].set(xlabel="Share of total linearized HI3 dose [%]",xlim=(0,100),title="Where does the nonlinear dose arise?")
    axes[1].legend(fontsize=8,loc="center left",bbox_to_anchor=(1.02,.5))
    fig.suptitle("Exposure budget versus HI3 dose | GW coefficients")
    fig.supxlabel("Bands are descriptive, not injury thresholds. Linearized dose is additive; final HI3 is its fractional power.",fontsize=9)
    finish(fig,output,"10_exposure_and_dose_bands","Explain why short strong-stress exposures can dominate damage despite low mean stress.",gallery)

def subset_plot(cases,output,gallery):
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
    rng=np.random.default_rng(73)
    for c in cases:
        data=c["data"]
        sizes=sorted(set([min(n,len(data)) for n in [100,250,500,1000,2000,len(data)]]))
        for ax,key,kind in zip(axes,["GW_HI3_percent","SA_dyne_s_cm2"],["mean","median"]):
            v=data[key].to_numpy()
            full=np.mean(v) if kind=="mean" else np.median(v)
            results=[]
            for n in sizes:
                estimates=[(np.mean if kind=="mean" else np.median)(v[rng.choice(len(v),n,replace=False)]) for _ in range(100)]
                results.append(np.percentile(estimates,[2.5,50,97.5])/full*100)
            results=np.array(results)
            ax.plot(sizes,results[:,1],color=c["color"],ls=c["style"],label=c["display"])
            ax.fill_between(sizes,results[:,0],results[:,2],color=c["color"],alpha=.1)
    for ax,title in zip(axes,["Mean HI3","Median SA"]):
        ax.axhline(100,color=".4",ls=":",lw=1)
        ax.set(title=title,xlabel="Paths sampled without replacement",ylabel="Estimate / complete export estimate [%]")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7,ncol=1 if len(cases)<7 else 2)
    fig.suptitle("Taskin Fig. 3 inspired sampling check | 100 random subsets per sample size")
    fig.supxlabel("Bands: central 95% of subset estimates. Conditional sampling variability, not CFD/mesh/seed-convergence proof.",fontsize=9)
    finish(fig,output,"11_subset_stability","Check whether observed device differences exceed variability from subsampling the existing paths.",gallery)

def summary_table(cases):
    rows=[]
    for c in cases:
        d=c["data"]
        row={k:c.get(k) for k in ("label","phase","bpm","cardiac_output_l_min","n_paths")}
        for model in MODELS:
            for method in ("HI2","HI3"):
                key=model+"_"+method+"_percent"
                mean,weighted=aggregate_hi(d[key],d.mean_inlet_velocity)
                row[key+"_mean"]=mean
                row[key+"_weighted"]=weighted
                row[key+"_p95"]=d[key].quantile(.95)
        for col in ("SA_dyne_s_cm2","transit_ms","max_stress_Pa","mean_stress_Pa","path_mean_speed_m_s","mean_inlet_velocity"):
            for q in (.5,.95,.99):
                row[col+"_p"+str(round(q*100))]=d[col].quantile(q)
        for cutoff in (1,2,5,10,35):
            row[f"SA_above_{cutoff}_percent"]=float((d.SA_dyne_s_cm2>cutoff).mean()*100)
        row["highest_1percent_share_HI3"]=float(np.sort(d.GW_HI3_percent)[-max(1,int(np.ceil(len(d)*.01))):].sum()/d.GW_HI3_percent.sum()*100)
        rows.append(row)
    return pd.DataFrame(rows)

def phase_summary(cases,output,gallery):
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
    for ax,key,stat,label in zip(axes.flat,
            ["GW_HI3_percent","SA_dyne_s_cm2","transit_ms","max_stress_Pa"],
            ["mean","median","median","median"],
            ["Mean HI3 [%] | GW","Median SA [dyne s/cm²]","Median transit [ms]","Median per-path peak stress [Pa]"]):
        lines_by_group(cases,ax,key,stat)
        tidy(ax,"Systolic phase [%]",label)
        ax.legend(fontsize=8)
    fig.suptitle("Comparison dashboard | equal-path statistics at matched systolic phases")
    finish(fig,output,"12_phase_comparison","Separate changes in predicted haemolysis, accumulated stress, residence time and stress peaks.",gallery)

def render(cases,output,record):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    gallery=[]
    sa_figures(cases,output,gallery)
    model_figures(cases,output,gallery)
    relationships(cases,output,gallery)
    histories_and_space(cases,output,gallery)
    concentration_and_bands(cases,output,gallery)
    subset_plot(cases,output,gallery)
    phase_summary(cases,output,gallery)
    table=summary_table(cases)
    table.to_csv(output/"discussion_metrics.csv",index=False)
    groups=sorted(set(c["group"] for c in cases))
    if len(groups)==2:
        a,b=[table[table.label==next(c["label"] for c in cases if c["group"]==g)].set_index("phase") for g in groups]
        common=a.index.intersection(b.index)
        if not len(common):
            raise ValueError("No matching systolic phases between devices; supply matched exports or explicit --phases.")
        record["matched_phases_percent"] = common.tolist()
        record["unmatched_phases_percent"] = {"reference": a.index.difference(common).tolist(), "comparison": b.index.difference(common).tolist()}
        metrics=["GW_HI3_percent_mean","SA_dyne_s_cm2_p50","SA_dyne_s_cm2_p95","transit_ms_p50","max_stress_Pa_p99"]
        compare=pd.DataFrame({"phase":common})
        for key in metrics:
            compare[key+"_comparison_over_reference"]=b.loc[common,key].to_numpy()/a.loc[common,key].replace(0,np.nan).to_numpy()
        compare.to_csv(output/"matched_phase_ratios.csv",index=False)
    gallery.sort()
    import html
    conditions="; ".join(f"{c['label']}: {c['bpm']:g} BPM, cardiac output " + (f"{c['cardiac_output_l_min']:g} L/min" if c.get("cardiac_output_l_min") is not None else "not supplied") for c in [next(c for c in cases if c["group"]==g) for g in groups])
    cards="".join('<section><h2>'+html.escape(name.replace("_"," "))+'</h2><p>'+html.escape(desc)+'</p><a href="'+name+'.png"><img loading="lazy" src="'+name+'.png"></a></section>' for name,desc in gallery)
    (output/"gallery.html").write_text('<!doctype html><meta charset="utf-8"><title>Paper-inspired haem comparison</title><style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;background:#faf8fb;color:#302938}section{background:white;padding:24px;margin:24px 0;border:1px solid #e6deec;border-radius:12px}img{max-width:100%;height:auto}h1{font-size:30px}p{line-height:1.6}</style><h1>Paper-inspired haem comparison</h1><p>Primary HI summaries give each reconstructed path equal weight. '+html.escape(record["comparison_context"])+'. '+html.escape(conditions)+'. Frozen-flow snapshots; no transient path stitching or whole-cycle integration.</p>'+cards,encoding="utf-8")
    readme = """# Paper-inspired comparison

Open gallery.html to browse all twelve figures. discussion_metrics.csv contains
the values used for quantitative discussion. Per-case paths.csv and plot_data.npz
allow --replot without running CFD-export processing again.

## What the papers motivate

Sources: [Taskin et al.](../../papers/taskin-etal.pdf) and
[Marom et al.](../../papers/marom-etal.pdf).

Taskin et al. 2012, Figs. 4 and 7: actual path stress histories, shown here with
objective P50/P95/P99 HI3-rank path selection and accumulated HI3.
Figs. 6 and 9: model and coefficient sensitivity. We compare HI2 and HI3 with
GW, HO and TZ on identical trajectories; there is no Eulerian model or experimental
validation in these exports. Fig. 3 motivates the subset check, which is not
a substitute for new-seed or CFD discretization convergence.

Marom et al. 2014, Fig. 5: SA distribution shape, median and tail. The new SA
density uses linear axes and a shared bandwidth, with a separate tail panel.
The exceedance figure states the percentage of paths above each SA cutoff.
Figs. 4 and 6 motivate the spatial stress display. Our display projects vertices
along frozen-flow paths, not transient particle locations or a solved field contour.
Marom's transient chamber-seeded particles, flow models and operating conditions
differ from these inlet-seeded frozen-flow paths; published medians are not
validation targets for this pipeline.

## HI weighting

The primary estimate is arithmetic mean of individual path HI: sum(HI_i)/N.
The sensitivity figure additionally reports sum(v_inlet_i*HI_i)/sum(v_inlet_i).
Weighting does not change any path's damage calculation. Speed already enters
the transit-time calculation; changing operating flow can change both stress and
exposure. Taskin's flow-versus-HI trend does not establish a universal inverse
local-velocity law or determine the appropriate sampling weights.

An equal-path mean describes the sampled path population. It represents blood
throughput only if the seeding already represents equal transported amounts.
For equal-area inlet seeding, normal flux times area provides physical throughput
weights; scalar inlet speed is only a proxy under forward, normal-flow assumptions.
Double-weighting flux-proportional particle releases would be inappropriate.
All paths are included here; outlet arrival has not been checked.

HI3 is [sum(A^(1/beta)*tau^(alpha/beta)*dt)]^beta in this code's notation.
The right-hand dose-band plot partitions that additive inner dose, not final HI3.
The existing velocity floor and scalar-shear model are inherited from haem5.
Coefficients have limited calibration ranges; sensitivity is not validation.
In particular, TZ was fitted at 50-320 Pa; these paths include much lower stresses.

## 35 cc versus 50 cc discussion

First compare at 120 BPM in both geometries. Then compare a separate matched-
cardiac-output condition with the corresponding size-specific BPMs.
Compare like systolic phases and use identical seeding, stress definitions,
outlet criteria, mesh/segment resolution and model coefficients where possible.
Geometry-only causation cannot be inferred if flow/boundary conditions also differ.
Changing --bpm is metadata only: new CFD exports must embody the changed condition.
Cycle-mean cardiac output is metadata supplied by the user; nominal chamber
capacity alone does not establish delivered stroke volume or flow.

Questions for the figures:
- Is the bulk SA distribution shifted, or only its high-exposure tail?
- Does a longer transit or stronger stress explain the shift?
- Are a few paths dominating mean HI3?
- Do geometry rankings persist under HI2/HI3 and all three coefficient sets?
- Does seeding/weighting change the conclusion?
- Where are sampled high-stress regions, and do those regions move with phase?

Spatial axes retain physical millimetres. Different device coordinate systems
must be aligned upstream for anatomical location comparisons.
Rank correlations are descriptive, not causal tests. Subset bands describe the
fixed exported paths, not biological or CFD-model uncertainty.
"""
    (output/"README.md").write_text(readme,encoding="utf-8")

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("folder",nargs="?",type=Path,default=Path("z_inputs"))
    p.add_argument("--label",default="50 cc")
    p.add_argument("--bpm",type=float,default=120.)
    p.add_argument("--cardiac-output-l-min",type=float)
    p.add_argument("--comparison-folder",type=Path)
    p.add_argument("--comparison-label",default="35 cc")
    p.add_argument("--comparison-bpm",type=float,default=120.)
    p.add_argument("--comparison-cardiac-output-l-min",type=float)
    p.add_argument("--comparison-context",help="Describe the operating-condition comparison; default is inferred from supplied BPM metadata.")
    p.add_argument("--phases",nargs="+",type=float,help="Phases in natural filename order, used for each folder.")
    p.add_argument("--output-dir",type=Path)
    p.add_argument("--replot",type=Path,help="Render cached per-case data from a previous run into a new output folder.")
    args=p.parse_args()
    output=args.output_dir or Path("outputs_paper_comparison")/datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
    output.mkdir(parents=True,exist_ok=False)
    if args.replot:
        record=json.loads((args.replot/"run.json").read_text())
        record.update(status="rendering",replot_source=str(args.replot.resolve()))
        base=Path(record.get("data_root",args.replot.resolve()))
    else:
        record={"status":"processing","comparison_context":args.comparison_context or ("Same-BPM comparison" if not args.comparison_folder or args.bpm==args.comparison_bpm else "Different-BPM operating conditions"),"cases":[],
                "data_root":str(output.resolve()),"HI_primary_weighting":"equal path",
                "source_papers":{"Taskin":"papers/taskin-etal.pdf","Marom":"papers/marom-etal.pdf"},
                "constants":haem5.POWER_LAW_CONSTANTS}
        base=output
    try:
        if not args.replot:
            groups=[(args.folder,args.label,args.bpm,args.cardiac_output_l_min)]
            if args.comparison_folder:
                groups.append((args.comparison_folder,args.comparison_label,args.comparison_bpm,args.comparison_cardiac_output_l_min))
                if args.label==args.comparison_label:
                    raise ValueError("Use different group labels.")
            for g,(folder,label,bpm,co) in enumerate(groups):
                if not np.isfinite(bpm) or bpm<=0 or (co is not None and (not np.isfinite(co) or co<=0)):
                    raise ValueError("BPM and supplied cardiac output must be finite and positive.")
                files=sorted(folder.glob("*.csv"),key=natural_key)
                phases=args.phases if args.phases is not None else [phase_from_name(f) for f in files]
                if len(files)<2 or len(phases)!=len(files) or not np.isfinite(phases).all() or np.any(np.diff(phases)<=0) or any(x<0 or x>100 for x in phases):
                    raise ValueError("Need at least two files and increasing phases in [0,100].")
                for number,(path,phase) in enumerate(zip(files,phases)):
                    meta=dict(id=f"g{g}_{number:02d}",group=g,label=label,phase=phase,bpm=bpm,cardiac_output_l_min=co)
                    directory=output/meta["id"]
                    directory.mkdir()
                    print(f"Processing {label}, {phase:g}%: {path.name}",flush=True)
                    extract_case(path,meta,directory)
                    record["cases"].append(meta)
                    (output/"run.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
        cases=load_cases(record,base)
        print("Rendering twelve paper-inspired figures...",flush=True)
        render(cases,output,record)
        record["status"]="complete"
    except Exception as exc:
        record.update(status="failed",error=str(exc))
        raise
    finally:
        (output/"run.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
    print("Saved: "+str(output.resolve()),flush=True)

if __name__=="__main__":
    main()
