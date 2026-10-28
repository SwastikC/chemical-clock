"""Run the whole chemical-clock analysis in order, from a fresh clone.

  0. check the data files against chemclock-data/SHA256SUMS
  1. crossmatch the HARPS-GTO stars with the NASA Exoplanet Archive snapshot
  2. choose the neighbours (scaled distance, 10 per host) and check them against
     the committed neighbour_table_v2.csv
  3. build the host-minus-neighbour table and check it against final_test_v2.csv
  4. the analysis: family offsets, chemical clocks, robustness tests, isochrone
     ages, Figures 1-6
  5. the uncertainty of the main result including the noise of the control stars

Stops at the first step that fails. Outputs go to chemclock-matching/output/,
chemclock-analysis/pipeline/output/ (figures) and chemclock-analysis/checks/output/.

Usage:  python reproduce.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'chemclock-data'
MATCH = ROOT / 'chemclock-matching'
ANALYSIS = ROOT / 'chemclock-analysis'

STEPS = [
    ('1. host crossmatch', [MATCH / 'crossmatch_hosts.py']),
    ('2. neighbours (v2)', [MATCH / 'match_neighbours.py']),
    ('3. host-minus-neighbour table', [MATCH / 'build_table.py']),
    ('4. analysis and figures', [ANALYSIS / 'pipeline' / 'run_pipeline.py']),
    ('5. uncertainty including control-star noise', [ANALYSIS / 'checks' / 'control_bootstrap.py']),
]


def banner(text):
    print(f'\n=== {text} ' + '=' * max(0, 70 - len(text)), flush=True)


banner('0. data checksums')
bad = []
for line in (DATA / 'SHA256SUMS').read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    ok = hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    print(f'{name}: {"OK" if ok else "FAILED"}')
    if not ok:
        bad.append(name)
if bad:
    sys.exit(f'{len(bad)} data file(s) differ from the committed versions; clone the repository again')

for title, args in STEPS:
    banner(title)
    if subprocess.run([sys.executable] + [str(x) for x in args]).returncode:
        sys.exit(f'step "{title}" failed')

banner('done')
print(f'figures in {ANALYSIS / "pipeline" / "output" / "figures"}')
