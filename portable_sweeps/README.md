# Portable sweep scripts

Copy this entire folder onto the lab computer (USB is fine). It contains
both scripts and the `haem5.py` pipeline they need, with no dependency on
the original repository. The calculation code is unchanged.

This is a portable source bundle, not a standalone executable: the computer
still needs **64-bit Python 3.11 or newer**. Python 3.11 is the tested version.
The setup creates a local environment inside this folder; no administrator
access or globally installed scientific packages are needed.

## Windows setup and use

1. Double-click `setup_windows.cmd` once. Installing packages requires internet
   access. If Python is missing, install Python first (or ask your lab technician).
2. Put your `sweep_N<count>.csv` exports beside `run_sweep.py`.
   Only `sweep_N10.csv` is included as a small working sample; copy the other
   exports from the original `sweeps` folder if you want the full sweep.
3. Double-click `run_sweep.cmd`. It creates `sweep_results.csv` here.
4. Double-click `plot_sweep.cmd`. It creates `hi_vs_pathlines.png` here.

The launchers keep the window open so you can read errors or completion messages.
Existing results are not bundled, so the included sample produces fresh results.
Without `--resume`, running a sweep replaces the results file.

For options, open a terminal in this folder:

```bat
run_sweep.cmd --dry-run
run_sweep.cmd --resume --constants GW HO --max-n 3000
run_sweep.cmd "D:\CFD exports" --resume
plot_sweep.cmd --colour
plot_sweep.cmd --metric nih --out nih_vs_pathlines.png
```

Default inputs and outputs stay beside the scripts even when launched from
another directory. Explicit relative paths are relative to your terminal's
current directory. Use full paths when in doubt. The sample has only one
particle-count level; add more exports to obtain convergence curves.

## Other systems / direct Python use

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_sweep.py --resume
.venv/bin/python plot_sweep.py
```

On Windows, direct Python commands use `.venv\Scripts\python.exe` instead.
The packages are pinned to the versions used to verify this copy.

## If the lab computer has no internet

On an internet-connected Windows computer with the same Python version and
architecture as the lab computer, run this from the portable folder:

```bat
python -m pip download --only-binary=:all: -r requirements.txt -d wheels
```

Copy the folder including `wheels` to the lab computer, then run:

```bat
python -m venv .venv
.venv\Scripts\python.exe -m pip install --no-index --find-links=wheels -r requirements.txt
```

You can then use the launchers normally. Python itself must already be
installed, or brought separately with an installer. Create `.venv` on the
destination computer; virtual environments should not be copied between PCs.
