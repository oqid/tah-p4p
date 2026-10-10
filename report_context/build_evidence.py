"""Rebuild report evidence from saved results; does not rerun or change CFD."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
sources = [
    'sweeps/sweep_results.csv', 'sweeps/300k_sweep_results.csv',
    'sweeps/a1to10k_sweep_results.csv',
    'sweeps/segment_sweeps/segment_sweep_results.csv',
    'sweeps/tolerance_sweep_results.csv',
    'outputs/runs/runner_check_2026-10-06/01_1000p_100_alotpa_3djawn/diagnostics/streamlines.csv',
]
original, large, dense, seg, tol, paths = [pd.read_csv(ROOT / p) for p in sources]
population = pd.concat([original, large, dense]).drop_duplicates(['N_requested', 'constants'], keep='last')
population = population.sort_values(['constants', 'N_requested'])
for name, frame, group, metric in [
    ('path_count_evidence.csv', population, 'N_requested', 'HI3_NIH_mg_per_100L'),
    ('segment_count_evidence.csv', seg, 'segments_requested', 'HI3_NIH_mg_per_100L'),
    ('tolerance_evidence.csv', tol, 'tolerance_grid_relative', 'HI3_NIH_mg_per_100L'),
]:
    frame = frame.copy()
    for constants in frame.constants.unique():
        mask = frame.constants.eq(constants)
        subset = frame.loc[mask]
        reference = subset.loc[subset[group].eq(.01), metric].iloc[0] if name.startswith('tolerance') else subset.sort_values(group)[metric].iloc[-1]
        frame.loc[mask, 'HI3_change_from_reference_percent'] = 100 * (subset[metric] / reference - 1)
    frame.to_csv(OUT / name, index=False)

paths = paths.sort_values('GW_HI3_percent', ascending=False).copy()
paths['normalized_proxy_weight'] = paths.weight / paths.weight.sum()
paths['weighted_HI3_contribution'] = paths.normalized_proxy_weight * paths.GW_HI3_percent
paths['cumulative_proxy_weight_percent'] = paths.normalized_proxy_weight.cumsum() * 100
paths['cumulative_HI3_contribution_percent'] = paths.weighted_HI3_contribution.cumsum() / paths.weighted_HI3_contribution.sum() * 100
paths.to_csv(OUT / 'preliminary_100_path_contributions.csv', index=False)

hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}
segment_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT / 'sweeps/segment_sweeps').glob('sweep_N5000_S*.csv')}
(OUT / 'evidence_provenance.json').write_text(json.dumps({'source_sha256': hashes, 'segment_export_sha256': segment_hashes}, indent=2), encoding='utf-8')

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for c, color, marker in [('GW', '#8865a0', 'o'), ('HO', '#5485ad', 's')]:
    for ax, filename, xcol in zip(axes, ['path_count_evidence.csv', 'segment_count_evidence.csv', 'tolerance_evidence.csv'], ['N_requested', 'segments_requested', 'tolerance_grid_relative']):
        d = pd.read_csv(OUT / filename).query('constants == @c').sort_values(xcol)
        ax.plot(d[xcol], d.HI3_change_from_reference_percent, color=color, marker=marker, ms=4, label=c)
        ax.set_xscale('log')
        ax.axhline(0, color='#666666', linewidth=.7)
        ax.grid(alpha=.15)
axes[0].set(title='Release count', xlabel='Requested paths', ylabel='HI3 change from reference (%)')
axes[1].set(title='Configured segment count', xlabel='Requested segments')
axes[2].set(title='Tracking tolerance', xlabel='Grid-relative tolerance')
axes[2].set_xticks([.0025, .005, .01, .02], ['0.0025', '0.005', '0.01', '0.02'])
axes[2].xaxis.set_minor_locator(NullLocator())
axes[0].legend(frameon=False)
fig.suptitle('Preliminary export sensitivity — one development flow field', fontsize=13)
fig.text(.5, .01, 'References: 300,000 requested paths; 5,000 configured segments; tolerance 0.01. These are separate sweeps, not combined error bounds.', ha='center', fontsize=9)
fig.tight_layout(rect=(0,.055,1,.94))
fig.savefig(OUT / 'export_sensitivity.png', dpi=200)
fig.savefig(OUT / 'export_sensitivity.pdf')
plt.close(fig)
print('Evidence tables, provenance and sensitivity figure written.')
print('Identical S1000/S2000/S5000:', len({segment_hashes[f'sweep_N5000_S{s}.csv'] for s in [1000,2000,5000]}) == 1)
