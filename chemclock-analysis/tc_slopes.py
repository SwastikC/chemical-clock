"""Condensation-temperature slopes of host-minus-neighbour abundances, and the
small/giant planet count.

For every host and each of its neighbours:
  - fit a straight line (ordinary least squares) to d[X/Fe] against the 50%
    condensation temperature Tc of element X;
  - keep the slope (dex/K).
Each host then gets a vote: +1 if most of its neighbour slopes are positive,
-1 if most are negative, 0 otherwise (with 3 neighbours: at least 2 of 3).
The votes are counted separately for small and giant planets.

--paper reproduces the numbers of the original notebook
(original/condensation_temperature_original.ipynb) exactly:
  --table legacy --pairing reversed --count planets
"reversed" is what the notebook did: the abundances were listed in order of
increasing Tc and the temperatures in order of decreasing Tc, so each element
was fitted at the temperature of another element. Use it only to
reproduce the notebook; "correct" puts every element at its own Tc.

Count modes:
  hosts    one vote per host; class small/giant from the largest Mp sin i of
           the host (output of chemclock-matching/crossmatch_hosts.py)
  planets  as the notebook: one row per planet, joined on the archive hd_name
           only; a host with three planets counts three times, and a host
           with both small and giant planets counts in both groups

Input is either a reference table (--table) or a pair table written by
chemclock-matching/match_neighbours.py (--pairs).

Output: output/tc_slopes_<tag>.csv, one row per host with its slopes and vote.

Usage:  python tc_slopes.py                 # v0 data, correct Tc, one vote per host
        python tc_slopes.py --paper         # the notebook's numbers
        python tc_slopes.py --pairs ../chemclock-matching/output/pairs_general_wr_w160_k3.csv
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=str(HERE.parent / 'chemclock-data' / 'data'))
ap.add_argument('--table', choices=['v0', 'legacy'], default='v0',
                help='reference/final_test_v0.csv, or the pre-Cr-fix legacy/final_test_original.csv')
ap.add_argument('--pairs', default=None, help='pair table from match_neighbours.py (overrides --table)')
ap.add_argument('--pairing', choices=['correct', 'reversed'], default='correct',
                help='which Tc each element is fitted at; "reversed" reproduces the notebook')
ap.add_argument('--exclude', nargs='*', default=[], metavar='EL', help='elements to leave out of the fit, e.g. Zn')
ap.add_argument('--count', choices=['hosts', 'planets'], default='hosts')
ap.add_argument('--cut', type=float, default=0.3, help='small/giant split in Mp sin i, Jupiter masses')
ap.add_argument('--crossmatch', default=str(HERE.parent / 'chemclock-matching' / 'output' / 'host_crossmatch.csv'),
                help='host classes for --count hosts')
ap.add_argument('--paper', action='store_true', help='same as --table legacy --pairing reversed --count planets')
ap.add_argument('--out', default=str(HERE / 'output'))
ap.add_argument('--tag', default=None)
a = ap.parse_args()
if a.paper:
    a.table, a.pairing, a.count, a.pairs = 'legacy', 'reversed', 'planets', None
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
source = Path(a.pairs).stem.removeprefix('pairs_') if a.pairs else a.table
tag = a.tag or ('paper' if a.paper else f'{source}_{a.pairing}_{a.count}')
norm = lambda s: str(s).replace(' ', '').strip()

# ---- d[X/Fe] per host, one array (neighbours x elements) each
ELS = ['Zn', 'Cu', 'Na', 'Mn', 'V', 'Cr', 'Co', 'Ni', 'Mg', 'Si', 'Ca', 'Ti', 'Sc', 'Al', 'Sr', 'Y', 'Ba', 'Ce']
tc = pd.read_csv(D / 'condensation_temps.csv', encoding='utf-8-sig')
tc = dict(zip(tc['Element'].str.strip(), tc['Temperature'].astype(float)))
els = sorted((e for e in ELS if e not in a.exclude), key=lambda e: tc[e])   # increasing Tc
T = np.array([tc[e] for e in els])
if a.pairing == 'reversed':
    T = T[::-1]

if a.pairs:
    p = pd.read_csv(a.pairs).sort_values(['host', 'rank'], kind='stable')
    diff = {norm(h): g[[f'd[{e}/Fe]' for e in els]].to_numpy(float) for h, g in p.groupby('host', sort=False)}
else:
    path = 'reference/final_test_v0.csv' if a.table == 'v0' else 'legacy/final_test_original.csv'
    t = pd.read_csv(D / path)
    diff = {norm(r['id']): np.array([[r[f'{i}_[{e}/Fe]'] for e in els] for i in (1, 2, 3)], float)
            for _, r in t.iterrows()}

# ---- slopes and votes
rows = []
for h, y in diff.items():
    slopes = [np.polyfit(T, yi, 1)[0] for yi in y]
    pos, neg = sum(s > 0 for s in slopes), sum(s < 0 for s in slopes)
    vote = 1 if pos > len(slopes) / 2 else -1 if neg > len(slopes) / 2 else 0
    rows.append(dict(host=h, **{f'slope_{i}': s for i, s in enumerate(slopes, 1)}, vote=vote))
res = pd.DataFrame(rows)

# ---- count small/giant
if a.count == 'hosts':
    if not Path(a.crossmatch).exists():
        raise SystemExit(f'--count hosts needs {a.crossmatch}: run chemclock-matching/crossmatch_hosts.py first')
    xm = pd.read_csv(a.crossmatch).drop_duplicates('star').set_index('star')
    res['mp_sini_mj'] = res['host'].map(xm['mp_sini_mj'])
    res['class'] = np.where(res['mp_sini_mj'].isna(), 'not in archive',
                            np.where(res['mp_sini_mj'] < a.cut, 'small', 'giant'))
    groups = {g: res[res['class'] == g] for g in ('small', 'giant')}
    unit = 'hosts'
else:
    arc = pd.read_csv(D / 'nasa/nasa_archive_snapshot.csv', low_memory=False)
    arc['hd'] = arc['hd_name'].map(lambda s: norm(s) if pd.notna(s) else '')
    m = arc.merge(res, left_on='hd', right_on='host')
    groups = {'small': m[m['pl_bmassj'] < a.cut], 'giant': m[m['pl_bmassj'] > a.cut]}
    unit = 'planets'
res.to_csv(OUT / f'tc_slopes_{tag}.csv', index=False)

# ---- report
print(f'Tc slopes {tag}: {"pairs " + a.pairs if a.pairs else "table " + a.table}, {len(els)} elements, '
      f'Tc pairing {a.pairing}, counting {a.count}')
v = res['vote'].value_counts()
print(f'  {len(res)} hosts; votes +1: {v.get(1, 0)}, -1: {v.get(-1, 0)}, 0: {v.get(0, 0)}')
for g, d in groups.items():
    print(f'  {g:5s} ({"<" if g == "small" else ">" if a.count == "planets" else ">="} {a.cut} Mj): '
          f'positive {(d["vote"] > 0).sum():3d}, negative {(d["vote"] < 0).sum():3d}  ({unit}, {len(d)} in group)')
if a.pairing == 'reversed':
    print('  NOTE: reversed pairing fits each element at the wrong Tc. Use it only to reproduce the notebook.')
print(f'wrote {OUT / f"tc_slopes_{tag}.csv"}')
