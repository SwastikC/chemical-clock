"""Crossmatch NASA Exoplanet Archive hosts to the HARPS-GTO sample.

For every HARPS-GTO star in general_properties.csv:
  1. match by identifier: the archive hd_name or hostname equal to the HARPS name
     (spaces removed), also allowing a trailing component letter A on the HARPS side;
  2. otherwise, nearest archive host within --radius arcsec.
The archive positions are at the Gaia epoch and the HARPS ones at J2000, so
high-proper-motion stars (tau Cet, eps Eri, 61 Vir) can sit 15-50 arcsec apart.
That is why the identifier match comes first.

Each matched star gets the largest Mp sin i among its planets that have a mass,
and a class: small (< --cut Mj), giant (>= --cut), or none (no mass).

Output: output/host_crossmatch.csv, plus a summary comparing with the HARPS
host flag (epsilon SwP?) and the 131 hosts of the sample.

Usage:  python crossmatch_hosts.py [--radius 60] [--cut 0.3] [--data DIR] [--out DIR]
        (by default reads ../chemclock-data/data and writes output/ next to this script)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--data', default=str(HERE.parent / 'chemclock-data' / 'data'))
ap.add_argument('--radius', type=float, default=60.0, help='positional match radius, arcsec')
ap.add_argument('--cut', type=float, default=0.3, help='small/giant split in Mp sin i, Jupiter masses')
ap.add_argument('--out', default=str(HERE / 'output'))
a = ap.parse_args()
D = Path(a.data); OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
norm = lambda s: str(s).replace(' ', '').strip()

# ---- archive hosts: one row per host
arc = pd.read_csv(D / 'nasa/nasa_archive_snapshot.csv', low_memory=False)
arc['hd'] = arc['hd_name'].map(lambda s: norm(s) if pd.notna(s) else '')
arc['hn'] = arc['hostname'].map(norm)
hosts = (arc.groupby('hostname')
            .agg(hd=('hd', 'first'), hn=('hn', 'first'), ra=('ra', 'first'), dec=('dec', 'first'),
                 n_planets=('pl_name', 'size'), mp_sini_mj=('pl_bmassj', 'max'))
            .reset_index())

# ---- HARPS-GTO stars
gp = pd.read_csv(D / 'harps_gto/general_properties.csv'); gp['n'] = gp['Star'].map(norm)
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
flag = dict(zip(eps['n'], eps['SwP?']))

by_name = {}
for i, r in hosts.iterrows():
    for key in (r['hd'], r['hn']):
        if key: by_name.setdefault(key, i)

def sep_arcsec(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600

rows = []
for _, s in gp.iterrows():
    d = sep_arcsec(s['RA (degree)'], s['Dec (degree)'], hosts['ra'].values, hosts['dec'].values)
    key = s['n'] if s['n'] in by_name else (s['n'][:-1] if s['n'].endswith('A') and s['n'][:-1] in by_name else None)
    if key:   # HARPS names a component (HD26965A) where the archive names the system (HD 26965)
        i, method = by_name[key], 'identifier'
    else:
        i = int(np.argmin(d)); method = 'position' if d[i] <= a.radius else None
    if method is None:
        rows.append(dict(star=s['n'], harps_flag=flag.get(s['n']), archive_host=None, method=None,
                         sep_arcsec=np.nan, n_planets=0, mp_sini_mj=np.nan))
        continue
    h = hosts.loc[i]
    rows.append(dict(star=s['n'], harps_flag=flag.get(s['n']), archive_host=h['hostname'], method=method,
                     sep_arcsec=round(float(d[i]), 2), n_planets=int(h['n_planets']), mp_sini_mj=h['mp_sini_mj']))
x = pd.DataFrame(rows)
x['class'] = np.where(x['mp_sini_mj'].isna(), np.where(x['archive_host'].notna(), 'none', ''),
                      np.where(x['mp_sini_mj'] < a.cut, 'small', 'giant'))
x.to_csv(OUT / 'host_crossmatch.csv', index=False)

# ---- summary
m = x[x['archive_host'].notna()]
print(f'HARPS-GTO stars: {len(x)};  archive hosts matched: {len(m)} '
      f'({(m.method == "identifier").sum()} by identifier, {(m.method == "position").sum()} by position)')
print(f'  classes at {a.cut} Mj: small {(m["class"] == "small").sum()}, giant {(m["class"] == "giant").sum()}, '
      f'no mass {(m["class"] == "none").sum()}')
flagged = set(eps.loc[eps['SwP?'] == 'Yes', 'n'])
print(f'HARPS host flag (epsilon SwP? = Yes): {len(flagged)} stars')
miss = sorted(flagged - set(m['star']))
print(f'  flagged as hosts but no planet in this archive snapshot ({len(miss)}): {", ".join(miss)}')
extra = m[~m['star'].isin(flagged)].sort_values('mp_sini_mj')
print(f'  archive hosts flagged as NON-hosts in HARPS ({len(extra)}), i.e. in the control pool:')
for _, r in extra.iterrows():
    print(f'    {r.star:12s} {r.archive_host:14s} Mp sin i = {r.mp_sini_mj:7.3f} Mj  ({r.method}, {r.sep_arcsec:.1f}")')
ref = set(pd.read_csv(D / 'reference/hosts.csv')['host'].map(norm))
mm = m[m['star'].isin(ref)]
print(f'the {len(ref)} hosts of the sample (hosts.csv): {len(mm)} matched to the archive; small {(mm["class"] == "small").sum()}, '
      f'giant {(mm["class"] == "giant").sum()}, no mass {(mm["class"] == "none").sum()}; '
      f'not in archive: {", ".join(sorted(ref - set(m["star"])))}')
print(f'\nwrote {OUT / "host_crossmatch.csv"}')
