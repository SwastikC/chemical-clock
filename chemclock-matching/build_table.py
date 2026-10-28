"""Build the host-minus-neighbour table that the analysis pipeline reads.

For every host in a neighbour table (columns SwP, First, Second, Third, N4, ...,
any number of neighbours, as written by match_neighbours.py) and each of its
neighbours, computes from epsilon_v0.csv
  value  = epsilon[host] - epsilon[neighbour]
  error  = e_host and e_neighbour added in quadrature
for [Fe/H], [Fe/H]_a and the 18 [X/Fe] columns. One row per host, one block of
columns per neighbour (prefix 1_, 2_, ...).

With the default input (the committed v2 neighbour table) the result is compared
cell by cell with the committed final_test_v2.csv; exit status 1 if they differ.

Usage:  python build_table.py                                   # v2, checked
        python build_table.py --neighbours output/neighbours_k3.csv --out FILE
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=str(HERE.parent / 'chemclock-data' / 'data'))
ap.add_argument('--neighbours', default=None, help='neighbour table (default: the committed v2 table)')
ap.add_argument('--out', default=None, help='where to write the table (default: output/final_test_<name>.csv)')
ap.add_argument('--tol', type=float, default=1e-9, help='largest allowed difference from the committed table')
a = ap.parse_args()
D = Path(a.data)
norm = lambda s: str(s).replace(' ', '').strip()
QS = ['[Fe/H]', '[Fe/H]_a'] + [f'[{e}/Fe]' for e in
      ['Zn', 'Y', 'Sc', 'Mg', 'Ca', 'Na', 'Si', 'V', 'Mn', 'Ti', 'Ba', 'Ni', 'Cu', 'Sr', 'Ce', 'Co', 'Cr', 'Al']]

eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
EPS = eps.drop_duplicates('n').set_index('n')
src = Path(a.neighbours) if a.neighbours else D / 'reference/neighbour_table_v2.csv'
nb = pd.read_csv(src)
slots = [c for c in nb.columns if c != 'SwP']          # First, Second, Third, N4, ...
rows = []
for _, r in nb.iterrows():
    h = norm(r['SwP']); row = {'id': h}
    for i, slot in enumerate(slots, 1):
        n = norm(r[slot]); row[slot] = n
        for q in QS:
            row[f'{i}_{q}'] = EPS.at[h, q] - EPS.at[n, q]
        for q in QS:
            row[f'{i}_e_{q}'] = float(np.hypot(EPS.at[h, f'e_{q}'], EPS.at[n, f'e_{q}']))
    rows.append(row)
layout = ['id'] + [c for i, slot in enumerate(slots, 1) for c in [slot] + [f'{i}_{q}' for q in QS] + [f'{i}_e_{q}' for q in QS]]
table = pd.DataFrame(rows)[layout]
dest = Path(a.out) if a.out else HERE / 'output' / f'final_test_{src.stem.replace("neighbour_table_", "").replace("neighbours_", "")}.csv'
dest.parent.mkdir(parents=True, exist_ok=True)
table.to_csv(dest, index=False)
print(f'built {len(table)} hosts x {len(slots)} neighbours x {len(QS)} quantities from {src.name}; wrote {dest}')

if a.neighbours is None:                               # the committed v2 table: compare
    old = pd.read_csv(D / 'reference/final_test_v2.csv'); new = pd.read_csv(dest)
    names = ['id'] + slots
    x, y = (t.drop(columns=names).to_numpy(float) for t in (new, old))
    same = (list(new.columns) == list(old.columns) and new[names].equals(old[names])
            and (np.isnan(x) == np.isnan(y)).all())    # missing values in the same cells
    worst = np.nanmax(np.abs(x - y)) if same else np.inf
    print(f'  same stars and missing values as final_test_v2.csv: {same};  largest difference: {worst:.1e}')
    if worst > a.tol:
        raise SystemExit('MISMATCH: the built table differs from final_test_v2.csv')
    print('OK: final_test_v2.csv is reproduced exactly')
