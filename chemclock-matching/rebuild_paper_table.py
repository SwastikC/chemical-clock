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

Usage:  python rebuild_paper_table.py [--data DIR] [--out DIR] [--tol 1e-9]
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
a = ap.parse_args()
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
norm = lambda s: str(s).replace(' ', '').strip()

# ---- inputs
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
EPS = eps.drop_duplicates('n').set_index('n')
stored = pd.read_csv(D / 'reference/final_test_v0.csv')
v6 = pd.read_csv(D / 'reference/neighbour_table_v6.csv')
v6 = {norm(r['SwP']): [norm(r[c]) for c in ('First', 'Second', 'Third')] for _, r in v6.iterrows()}
qs = [c[2:] for c in stored.columns if c.startswith('1_') and not c.startswith('1_e_')]

# ---- rebuild, in the column order of the stored table
rows = []
for h in stored['id'].map(norm):
    row = {'id': h}
    for i, (slot, n) in enumerate(zip(('First', 'Second', 'Third'), v6[h]), 1):
        row[slot] = n
        for q in qs:
            row[f'{i}_{q}'] = EPS.at[h, q] - EPS.at[n, q]
        for q in qs:
            row[f'{i}_e_{q}'] = float(np.hypot(EPS.at[h, f'e_{q}'], EPS.at[n, f'e_{q}']))
    rows.append(row)
rebuilt = pd.DataFrame(rows)[list(stored.columns)]
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
