"""Run the paper's analysis pipeline in order and write all figures.

  1. build_stellar_offsets.py  host-minus-neighbour Teff, logg, [Fe/H]
  2. robustness.py             family offsets, chemical clocks, robustness tests
  3. figures_v3.py             Figures 1-5
  4. ext1_iso_age.py           isochrone ages against [Y/Mg] ages
  5. ext1_fig.py               the isochrone figure

Each script prints its results. With the v2 neighbour table (the default), outputs go
to output/ and figures to output/figures/.

--neighbours FILE redoes everything with another neighbour table (columns SwP, First,
Second, Third, ..., e.g. a neighbours_<tag>.csv from match_neighbours.py with other
--k or --scales). The host-minus-neighbour table is built from epsilon_v0.csv first,
and all outputs go to output/<tag>/.

Usage:  python run_pipeline.py                                  # v2
        python run_pipeline.py --neighbours FILE [--tag NAME]   # any other neighbour table
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = ['build_stellar_offsets.py', 'robustness.py', 'figures_v3.py', 'ext1_iso_age.py', 'ext1_fig.py']
ap = argparse.ArgumentParser()
ap.add_argument('--neighbours', default=None, help='neighbour table to use instead of v2')
ap.add_argument('--tag', default=None, help='name of the output folder (default: the file name)')
a = ap.parse_args()

env = dict(os.environ, MPLBACKEND='Agg')   # write figures without opening windows
out = HERE / 'output'
if a.neighbours:
    nb = Path(a.neighbours).resolve()
    out = out / (a.tag or nb.stem.removeprefix('neighbours_'))
    table = out / 'final_test.csv'
    build = [sys.executable, str(HERE.parents[1] / 'chemclock-matching' / 'build_table.py'),
             '--neighbours', str(nb), '--out', str(table)]
    if subprocess.run(build).returncode:
        sys.exit('could not build the table for ' + str(nb))
    env.update(CHEMCLOCK_NEIGHBOURS=str(nb), CHEMCLOCK_FINAL_TEST=str(table), CHEMCLOCK_OUT=str(out))

for s in SCRIPTS:
    print(f'\n--- {s} ' + '-' * max(0, 70 - len(s)), flush=True)
    if subprocess.run([sys.executable, str(HERE / s)], env=env).returncode:
        sys.exit(f'{s} failed')
print(f'\nfigures in {out / "figures"}')
