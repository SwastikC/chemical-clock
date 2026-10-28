"""Rank-mean neighbour matching for the chemical-clock sample.

For each host:
  - candidates = control stars with |Teff_host - Teff| <= --window K
  - rank the candidates separately on |dTeff|, |dlogg| and |d[Fe/H]|
  - average the three ranks; the --k lowest mean ranks are the neighbours.

With replacement (default), each host takes its own top k, so a control can serve
several hosts. With --no-replacement, every control is used at most once and the
assignment minimises the total mean rank over all hosts at once (Hungarian
algorithm), so the result does not depend on the order in which hosts are processed.

Writes
  output/neighbours_<tag>.csv  one row per host: SwP, First, Second, ... (best first)
  output/pairs_<tag>.csv       one row per host-neighbour pair, in the agreed schema:
                               d[X/Fe] and its error for 18 elements from epsilon_v0.csv,
                               plus dTeff, dlogg, d[Fe/H] from the parameters used to match
and prints the match quality and the overlap with the stored v6 neighbour table.

Options left open by the paper's Methods text are exposed as flags:
  --params   general | harps   parameter source (General Properties table, or the
                               HARPS-GTO columns TeffK, logg[cm/s2], [Fe/H]---_x with
                               General Properties as fallback for hosts missing there)
  --ties     ordinal | average how tied distances are ranked
  --include-host               rank the host together with its candidates (as the
                               original notebook did)
  --drop-archive-hosts         remove controls that are planet hosts in the archive
                               snapshot (needs output/host_crossmatch.csv)
  --hosts    v6 | final        host list: the 133 of neighbour_table_v6 or the 131 of final_test

Usage:  python match_neighbours.py [--params general] [--no-replacement] [--tag mytest]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.optimize import linear_sum_assignment

ap = argparse.ArgumentParser()
ap.add_argument('--data', default='../chemclock-data/data')
ap.add_argument('--params', choices=['general', 'harps'], default='general')
ap.add_argument('--hosts', choices=['v6', 'final'], default='final')
ap.add_argument('--window', type=float, default=160.0, help='Teff window, K')
ap.add_argument('--k', type=int, default=3)
ap.add_argument('--no-replacement', action='store_true')
ap.add_argument('--ties', choices=['ordinal', 'average'], default='ordinal')
ap.add_argument('--include-host', action='store_true')
ap.add_argument('--drop-archive-hosts', action='store_true')
ap.add_argument('--out', default='output')
ap.add_argument('--tag', default=None)
a = ap.parse_args()
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(exist_ok=True)
tag = a.tag or f'{a.params}_{"nr" if a.no_replacement else "wr"}_w{int(a.window)}_k{a.k}'
norm = lambda s: str(s).replace(' ', '').strip()
ELS = ['Zn', 'Cu', 'Na', 'Mn', 'V', 'Cr', 'Co', 'Ni', 'Mg', 'Si', 'Ca', 'Ti', 'Sc', 'Al', 'Sr', 'Y', 'Ba', 'Ce']

# ---- inputs
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
EPS = eps.drop_duplicates('n').set_index('n')
pool_names = set(eps.loc[eps['SwP?'] == 'No', 'n'])
if a.drop_archive_hosts:
    xm = pd.read_csv(OUT / 'host_crossmatch.csv')
    archive_hosts = set(xm.loc[xm['archive_host'].notna(), 'star'])
    dropped = sorted(pool_names & archive_hosts)
    pool_names -= archive_hosts
    print(f'dropped {len(dropped)} archive hosts from the control pool')

gp = pd.read_csv(D / 'harps_gto/general_properties.csv'); gp['n'] = gp['Star'].map(norm)
G = gp.drop_duplicates('n').set_index('n')[['Teff (K)', 'logg [cm/s2]', '[Fe/H]']]
G.columns = ['T', 'g', 'F']
if a.params == 'general':
    P = G
else:
    cols = ['TeffK', 'logg[cm/s2]', '[Fe/H]---_x']; rec = {}
    for fn, keys in (('starswithoutplanet.csv', ['Star---']), ('starswithplanet.csv', ['Star---', 'hostname_1'])):
        for _, r in pd.read_csv(D / 'harps_gto' / fn).dropna(subset=cols).iterrows():
            for k in keys: rec.setdefault(norm(r[k]), tuple(float(v) for v in r[cols]))
    H = pd.DataFrame.from_dict(rec, orient='index', columns=['T', 'g', 'F'])
    P = pd.concat([H, G[~G.index.isin(H.index)]])

if a.hosts == 'v6':
    hosts = list(pd.read_csv(D / 'reference/neighbour_table_v6.csv')['SwP'].map(norm))
else:
    hosts = list(pd.read_csv(D / 'reference/final_test_v0.csv')['id'].map(norm))
pool = P[P.index.isin(pool_names)]

# ---- rank-mean per host
R = {}
for h in hosts:
    hp = P.loc[h]
    c = pool[(pool['T'] - hp['T']).abs() <= a.window]
    if a.include_host:
        c = pd.concat([hp.to_frame().T.set_axis(['__host__']), c])
    mr = np.mean([rankdata((c[q] - hp[q]).abs(), method=a.ties) for q in 'TgF'], axis=0)
    s = pd.Series(mr, index=c.index).drop('__host__', errors='ignore')
    R[h] = s.sort_values(kind='stable')
    if len(R[h]) < a.k:
        raise SystemExit(f'{h}: only {len(R[h])} candidates within {a.window} K')

if not a.no_replacement:
    assign = {h: list(R[h].index[:a.k]) for h in hosts}
else:
    cand = sorted(set().union(*[set(r.index) for r in R.values()])); ci = {n: i for i, n in enumerate(cand)}
    BIG = 1e9
    C = np.full((a.k * len(hosts), len(cand)), BIG)
    for j, h in enumerate(hosts):
        for n, v in R[h].items():
            C[a.k * j:a.k * (j + 1), ci[n]] = v
    rows, cols_ = linear_sum_assignment(C)
    if (C[rows, cols_] >= BIG).any():
        raise SystemExit('no assignment without replacement fits inside the Teff window')
    assign = {h: [] for h in hosts}
    for r_, c_ in zip(rows, cols_): assign[hosts[r_ // a.k]].append(cand[c_])
    assign = {h: sorted(ns, key=lambda n: R[h][n]) for h, ns in assign.items()}

# ---- outputs
names = ['First', 'Second', 'Third'] + [f'N{i}' for i in range(4, a.k + 1)]
pd.DataFrame([[h] + ns for h, ns in assign.items()], columns=['SwP'] + names[:a.k]).to_csv(OUT / f'neighbours_{tag}.csv', index=False)
pairs = []
for h, ns in assign.items():
    for rank, n in enumerate(ns, 1):
        row = dict(matching_id=tag, host=h, neighbour=n, rank=rank, mean_rank=R[h][n],
                   dTeff=P.at[h, 'T'] - P.at[n, 'T'], dlogg=P.at[h, 'g'] - P.at[n, 'g'], dFeH=P.at[h, 'F'] - P.at[n, 'F'])
        for e in ELS:
            row[f'd[{e}/Fe]'] = EPS.at[h, f'[{e}/Fe]'] - EPS.at[n, f'[{e}/Fe]']
            row[f'e_d[{e}/Fe]'] = float(np.hypot(EPS.at[h, f'e_[{e}/Fe]'], EPS.at[n, f'e_[{e}/Fe]']))
        pairs.append(row)
pairs = pd.DataFrame(pairs); pairs.to_csv(OUT / f'pairs_{tag}.csv', index=False)

# ---- report
use = pairs['neighbour'].value_counts()
d = pairs[['dTeff', 'dlogg', 'dFeH']].abs()
print(f'matching {tag}: {len(hosts)} hosts, {len(pairs)} pairs, {pairs.neighbour.nunique()} unique controls, '
      f'max use {use.max()} ({use.index[0]})')
print(f'  median |dTeff| {d.dTeff.median():.0f} K, |dlogg| {d.dlogg.median():.3f}, |d[Fe/H]| {d.dFeH.median():.3f};  '
      f'max |dTeff| {d.dTeff.max():.0f} K')
st = pd.read_csv(D / 'reference/neighbour_table_v6.csv')
st = {norm(r['SwP']): [norm(r[c]) for c in ('First', 'Second', 'Third')] for _, r in st.iterrows()}
common = [h for h in hosts if h in st]
same_set = sum(set(assign[h][:3]) == set(st[h]) for h in common)
same_order = sum(assign[h][:3] == st[h] for h in common)
shared = np.mean([len(set(assign[h][:3]) & set(st[h])) / 3 for h in common])
print(f'  vs stored v6 table ({len(common)} hosts): same three neighbours {same_set}, same order {same_order}, '
      f'neighbour overlap {shared:.0%}')
print(f'wrote {OUT / f"neighbours_{tag}.csv"} and {OUT / f"pairs_{tag}.csv"}')
