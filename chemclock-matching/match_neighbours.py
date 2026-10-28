"""Neighbour matching for the chemical-clock sample (the v2 neighbour table).

For each of the 131 hosts in chemclock-data/data/reference/hosts.csv:
  - distance to every control star, each parameter difference divided by its
    typical HARPS-GTO uncertainty (--scales, default 50 K, 0.10 dex, 0.03 dex):
        d = sqrt( (dTeff/50)^2 + (dlogg/0.10)^2 + (d[Fe/H]/0.03)^2 )
  - the --k (default 10) controls with the smallest distance are the neighbours.

Parameters are those of the General Properties table. Control stars are the
HARPS-GTO stars without a planet flag, minus every star that is a planet host in
the archive snapshot (needs output/host_crossmatch.csv from crossmatch_hosts.py).
A control can serve several hosts.

Writes
  output/neighbours_<tag>.csv  one row per host: SwP, First, Second, Third, N4, ... (closest first)
  output/pairs_<tag>.csv       one row per host-neighbour pair: d[X/Fe] and its error for
                               18 elements from epsilon_v0.csv, dTeff, dlogg, d[Fe/H], distance
and prints the match quality. With the defaults the result is identical to
chemclock-data/data/reference/neighbour_table_v2.csv, which is checked.

Usage:  python match_neighbours.py                       # v2
        python match_neighbours.py --k 3 --tag k3          # other choices, for comparison
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=str(HERE.parent / 'chemclock-data' / 'data'))
ap.add_argument('--k', type=int, default=10, help='neighbours per host')
ap.add_argument('--scales', type=float, nargs=3, default=[50.0, 0.10, 0.03], metavar=('T', 'G', 'F'),
                help='uncertainties of Teff (K), logg and [Fe/H] (dex) that divide the differences')
ap.add_argument('--out', default=str(HERE / 'output'))
ap.add_argument('--tag', default='v2')
a = ap.parse_args()
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
norm = lambda s: str(s).replace(' ', '').strip()
ELS = ['Zn', 'Cu', 'Na', 'Mn', 'V', 'Cr', 'Co', 'Ni', 'Mg', 'Si', 'Ca', 'Ti', 'Sc', 'Al', 'Sr', 'Y', 'Ba', 'Ce']

# ---- inputs
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
EPS = eps.drop_duplicates('n').set_index('n')
if not (OUT / 'host_crossmatch.csv').exists():
    raise SystemExit(f'needs {OUT / "host_crossmatch.csv"}: run crossmatch_hosts.py first')
xm = pd.read_csv(OUT / 'host_crossmatch.csv')
archive_hosts = set(xm.loc[xm['archive_host'].notna(), 'star'])
pool_names = set(eps.loc[eps['SwP?'] == 'No', 'n'])
print(f'removed {len(pool_names & archive_hosts)} archive planet hosts from the control pool')
pool_names -= archive_hosts

gp = pd.read_csv(D / 'harps_gto/general_properties.csv'); gp['n'] = gp['Star'].map(norm)
P = gp.drop_duplicates('n').set_index('n')[['Teff (K)', 'logg [cm/s2]', '[Fe/H]']]
P.columns = ['T', 'g', 'F']
hosts = list(pd.read_csv(D / 'reference/hosts.csv')['host'].map(norm))
pool = P[P.index.isin(pool_names)]

# ---- scaled distance, k nearest
scale = pd.Series(a.scales, index=list('TgF'))
R = {}
for h in hosts:
    R[h] = np.sqrt((((pool[list('TgF')] - P.loc[h, list('TgF')]) / scale) ** 2).sum(axis=1)).sort_values(kind='stable')
assign = {h: list(R[h].index[:a.k]) for h in hosts}

# ---- outputs
names = ['First', 'Second', 'Third'] + [f'N{i}' for i in range(4, a.k + 1)]
pd.DataFrame([[h] + ns for h, ns in assign.items()], columns=['SwP'] + names[:a.k]).to_csv(OUT / f'neighbours_{a.tag}.csv', index=False)
pairs = []
for h, ns in assign.items():
    for rank, n in enumerate(ns, 1):
        row = dict(matching_id=a.tag, host=h, neighbour=n, rank=rank, scaled_distance=R[h][n],
                   dTeff=P.at[h, 'T'] - P.at[n, 'T'], dlogg=P.at[h, 'g'] - P.at[n, 'g'], dFeH=P.at[h, 'F'] - P.at[n, 'F'])
        for e in ELS:
            row[f'd[{e}/Fe]'] = EPS.at[h, f'[{e}/Fe]'] - EPS.at[n, f'[{e}/Fe]']
            row[f'e_d[{e}/Fe]'] = float(np.hypot(EPS.at[h, f'e_[{e}/Fe]'], EPS.at[n, f'e_[{e}/Fe]']))
        pairs.append(row)
pairs = pd.DataFrame(pairs); pairs.to_csv(OUT / f'pairs_{a.tag}.csv', index=False)

# ---- report
use = pairs['neighbour'].value_counts()
d = pairs[['dTeff', 'dlogg', 'dFeH']].abs()
print(f'matching {a.tag}: {len(hosts)} hosts, {len(pairs)} pairs, {pairs.neighbour.nunique()} unique controls, '
      f'max use {use.max()} ({use.index[0]})')
print(f'  median |dTeff| {d.dTeff.median():.0f} K, |dlogg| {d.dlogg.median():.3f}, |d[Fe/H]| {d.dFeH.median():.3f};  '
      f'max |dTeff| {d.dTeff.max():.0f} K')
print(f'wrote {OUT / f"neighbours_{a.tag}.csv"} and {OUT / f"pairs_{a.tag}.csv"}')
ref = D / 'reference/neighbour_table_v2.csv'
if a.k == 10 and list(a.scales) == [50.0, 0.10, 0.03] and ref.exists():
    same = pd.read_csv(OUT / f'neighbours_{a.tag}.csv').equals(pd.read_csv(ref))
    print(f'  identical to {ref.name}: {same}')
    if not same:
        raise SystemExit(f'MISMATCH: the neighbours differ from {ref.name}')
