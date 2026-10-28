"""Run the paper's analysis pipeline in order and write all figures.

  1. build_stellar_offsets.py  host-minus-neighbour Teff, logg, [Fe/H]
  2. robustness.py             family offsets, chemical clocks, robustness tests
  3. figures_v3.py             Figures 1-5
  4. ext1_iso_age.py           isochrone ages against [Y/Mg] ages
  5. ext1_fig.py               the isochrone figure

Each script prints its results; outputs go to output/ and figures to output/figures/.

Usage:  python run_pipeline.py
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = ['build_stellar_offsets.py', 'robustness.py', 'figures_v3.py', 'ext1_iso_age.py', 'ext1_fig.py']
env = dict(os.environ, MPLBACKEND='Agg')   # write figures without opening windows

for s in SCRIPTS:
    print(f'\n--- {s} ' + '-' * max(0, 70 - len(s)), flush=True)
    if subprocess.run([sys.executable, str(HERE / s)], env=env).returncode:
        sys.exit(f'{s} failed')
print(f'\nfigures in {HERE / "output" / "figures"}')
