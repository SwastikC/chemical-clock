"""Run the four neighbour-reuse checks in order (about 1.5 minutes).

  1. reuse_bootstrap.py   sharing structure, bootstraps, leave-one-neighbour-out
  2. calibrate.py         simulation on the real host-neighbour graph; reuse inflation
  3. rematch.py           rank-mean matching redone, with and without replacement
  4. matching_spread.py   400 random draws of 3 from each host's 10 best candidates

Each script prints its results; CSVs go to output/.

Usage:  python run_checks.py
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for s in ['reuse_bootstrap.py', 'calibrate.py', 'rematch.py', 'matching_spread.py']:
    print(f'\n--- {s} ' + '-' * max(0, 70 - len(s)), flush=True)
    if subprocess.run([sys.executable, str(HERE / s)]).returncode:
        sys.exit(f'{s} failed')
print(f'\nCSVs in {HERE / "output"}')
