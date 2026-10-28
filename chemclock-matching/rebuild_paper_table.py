"""Rebuild the paper v6 input table (final_test_v0.csv) from its ingredients.

For every host in final_test_v0.csv, takes the three neighbours from
neighbour_table_v6.csv and recomputes, from epsilon_v0.csv,
  value  = epsilon[host] - epsilon[neighbour]
  error  = e_host and e_neighbour added in quadrature
for [Fe/H], [Fe/H]_a and the 18 [X/Fe] columns. The result has the same layout
as final_test_v0.csv and is compared with it cell by cell.

This checks that the stored table follows from the v0 abundances and the v6
neighbours. It does not regenerate the v6 neighbours themselves (see the
known issues in chemclock-data/README.md).

Output: output/final_test_rebuilt.csv. Exit status 1 if the rebuilt table does
not match the stored one.

With --neighbours, builds the same table for any other neighbour table (columns
SwP, First, Second, Third, ... with any number of neighbours, as written by
match_neighbours.py) and writes it to
--table-out, without comparing; this is the input the paper pipeline needs.

Usage:  python rebuild_paper_table.py [--data DIR] [--out DIR] [--tol 1e-9]
        python rebuild_paper_table.py --neighbours output/neighbours_harps_clean.csv --table-out FILE
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=str(HERE.parent / 'chemclock-data' / 'data'))
ap.add_argument('--out', default=str(HERE / 'output'))
ap.add_argument('--tol', type=float, default=1e-9, help='largest allowed absolute difference')
ap.add_argument('--neighbours', default=None, help='build the table for this neighbour table instead of v6')
ap.add_argument('--table-out', default=None, help='with --neighbours: where to write the table')
a = ap.parse_args()
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
norm = lambda s: str(s).replace(' ', '').strip()

# ---- inputs
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
EPS = eps.drop_duplicates('n').set_index('n')
stored = pd.read_csv(D / 'reference/final_test_v0.csv')
nb = pd.read_csv(a.neighbours or D / 'reference/neighbour_table_v6.csv')
slots = [c for c in nb.columns if c != 'SwP']          # First, Second, Third, N4, ...
nb = {norm(r['SwP']): [norm(r[c]) for c in slots] for _, r in nb.iterrows()}
qs = [c[2:] for c in stored.columns if c.startswith('1_') and not c.startswith('1_e_')]
hosts = list(nb) if a.neighbours else list(stored['id'].map(norm))

# ---- rebuild, in the column order of the stored table
rows = []
for h in hosts:
    row = {'id': h}
    for i, (slot, n) in enumerate(zip(slots, nb[h]), 1):
        row[slot] = n
        for q in qs:
            row[f'{i}_{q}'] = EPS.at[h, q] - EPS.at[n, q]
        for q in qs:
            row[f'{i}_e_{q}'] = float(np.hypot(EPS.at[h, f'e_{q}'], EPS.at[n, f'e_{q}']))
    rows.append(row)
layout = ['id'] + [c for i, slot in enumerate(slots, 1)
                   for c in [slot] + [f'{i}_{q}' for q in qs] + [f'{i}_e_{q}' for q in qs]]
rebuilt = pd.DataFrame(rows)[layout]                 # same layout as final_test_v0.csv, k blocks
if a.neighbours:
    dest = Path(a.table_out or OUT / f'final_test_{Path(a.neighbours).stem}.csv')
    dest.parent.mkdir(parents=True, exist_ok=True)
    rebuilt.to_csv(dest, index=False)
    print(f'built {len(rebuilt)} hosts x 3 neighbours from {a.neighbours}; wrote {dest}')
    raise SystemExit(0)
rebuilt.to_csv(OUT / 'final_test_rebuilt.csv', index=False)

# ---- compare
names = ['id', 'First', 'Second', 'Third']
same_names = (stored[names].map(norm) == rebuilt[names]).all().all()
num = [c for c in stored.columns if c not in names]
diff = (stored[num] - rebuilt[num]).abs().to_numpy()
worst = np.nanmax(diff)
print(f'rebuilt {len(rebuilt)} hosts x 3 neighbours x {len(qs)} quantities from epsilon_v0 + neighbour_table_v6')
print(f'  star names identical: {same_names};  largest |stored - rebuilt|: {worst:.1e}')
print(f'wrote {OUT / "final_test_rebuilt.csv"}')
if not same_names or worst > a.tol:
    raise SystemExit('MISMATCH: the rebuilt table differs from final_test_v0.csv')
print('OK: final_test_v0.csv is reproduced exactly')
