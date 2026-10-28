"""Run the whole chemical-clock pipeline in order, from a fresh clone.

  0. check the data files against chemclock-data/SHA256SUMS
  1. crossmatch the HARPS-GTO stars with the NASA Exoplanet Archive snapshot
  2. rebuild the paper v6 input table (final_test_v0.csv) and check it
  3. neighbour matching with the default settings and two variants
  4. the paper's analysis pipeline and its figures
  5. how the results depend on the neighbours (reuse, rematching)
  6. the original notebook's slope-sign count, then the corrected count

Stops at the first step that fails. Outputs go to chemclock-matching/output/,
chemclock-analysis/pipeline/output/ (figures) and chemclock-analysis/output/.

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
    ('2. rebuild paper input table', [MATCH / 'rebuild_paper_table.py']),
    ('3a. matching: defaults', [MATCH / 'match_neighbours.py']),
    ('3b. matching: without replacement', [MATCH / 'match_neighbours.py', '--no-replacement']),
    ('3c. matching: HARPS parameters, archive hosts removed',
     [MATCH / 'match_neighbours.py', '--params', 'harps', '--drop-archive-hosts', '--tag', 'harps_clean']),
    ('4. paper analysis and figures', [ANALYSIS / 'pipeline' / 'run_pipeline.py']),
    ('5. neighbour-reuse and matching checks', [ANALYSIS / 'checks' / 'neighbour_reuse' / 'run_checks.py']),
    ('6a. slope-sign count: original notebook', [ANALYSIS / 'tc_slopes.py', '--paper']),
    ('6b. slope-sign count: corrected', [ANALYSIS / 'tc_slopes.py']),
]


def banner(text):
    print(f'\n=== {text} ' + '=' * max(0, 70 - len(text)), flush=True)


banner('0. data checksums (tag v0)')
bad = []
for line in (DATA / 'SHA256SUMS').read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    ok = hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    print(f'{name}: {"OK" if ok else "FAILED"}')
    if not ok:
        bad.append(name)
if bad:
    sys.exit(f'{len(bad)} data file(s) differ from tag v0; re-download them before going on')

for title, args in STEPS:
    banner(title)
    if subprocess.run([sys.executable] + [str(x) for x in args]).returncode:
        sys.exit(f'step "{title}" failed')

banner('done')
print(f'outputs in {MATCH / "output"}, {ANALYSIS / "output"}')
print(f'figures in {ANALYSIS / "pipeline" / "output" / "figures"}')
