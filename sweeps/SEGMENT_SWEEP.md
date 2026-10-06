# Segment-count sweep

Export the same requested number of pathlines at each segment-count setting,
using names such as `sweep_N1000_S100.csv`, `sweep_N1000_S200.csv`, and
`sweep_N1000_S500.csv`. Place them beside `run_segment_sweep.py` or pass an
input folder. The script processes existing exports; it does not resample them.

```powershell
python sweeps/run_segment_sweep.py --pathlines 1000 --dry-run
python sweeps/run_segment_sweep.py --pathlines 1000 --resume --constants GW HO
python sweeps/run_segment_sweep.py "D:\CFD exports" --pathlines 1000 --min-segments 100 --max-segments 500
```

Results go to `segment_sweep_results.csv`, with requested pathlines,
requested segments, actual processed streamlines, and the same haemolysis
metrics as the original sweep. Resume uses pathlines, segments, and constants
together, so multiple fixed pathline counts can share a results file.
`--out` selects another results file. Without `--resume`, results are replaced.

The portable folder contains the same runner and its own pipeline copy:
`.venv\Scripts\python.exe run_segment_sweep.py --pathlines 1000 --resume`.
