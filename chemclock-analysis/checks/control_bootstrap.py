"""Uncertainty of the small-vs-giant result including the noise of the control stars.

The pipeline's sigma treats each host's neighbour average as exact. With the
scaled matching and 10 neighbours, controls are shared between hosts and their own
scatter adds to the uncertainty. This script redoes the scaled matching
(match_neighbours.py --method scaled --drop-archive-hosts, General Properties
parameters) inside a bootstrap:

  hosts     resample the small and the giant hosts (separately, with replacement)
  controls  resample the control pool (with replacement) and rematch every host
  both      both at once

and reports the standard error of the small-minus-giant difference for the
s-process family and the [Y/Mg] clock under each scheme. The "hosts" scheme is the
usual standard error; "both" adds the control noise. The checks in
neighbour_reuse/ show that resampling both tends to overstate the error, so the
two bracket the true value.

The unresampled estimate is the same quantity the pipeline computes (neighbours
combined by inverse variance, element by element), and is checked against it.

Usage:  python control_bootstrap.py [--k 10] [--scales 50 0.10 0.03] [--n-boot 500]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / 'chemclock-data' / 'data'
ap = argparse.ArgumentParser()
ap.add_argument('--k', type=int, default=10)
ap.add_argument('--scales', type=float, nargs=3, default=[50.0, 0.10, 0.03], metavar=('T', 'G', 'F'))
ap.add_argument('--n-boot', type=int, default=500)
ap.add_argument('--seed', type=int, default=2026)
ap.add_argument('--crossmatch', default=str(ROOT / 'chemclock-matching' / 'output' / 'host_crossmatch.csv'))
ap.add_argument('--out', default=str(Path(__file__).resolve().parent / 'output'))
a = ap.parse_args()
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)
norm = lambda s: str(s).replace(' ', '').strip()
SP = ['Sr', 'Y', 'Ba', 'Ce']
ELS = SP + ['Mg']                     # s-process family and the [Y/Mg] clock
YMG_SLOPE = -0.041                    # dex/Gyr, Tucci Maia et al. 2016, as in the pipeline

# ---- inputs, as in match_neighbours.py --method scaled --drop-archive-hosts --params general
eps = pd.read_csv(D / 'harps_gto/epsilon_v0.csv'); eps['n'] = eps['Star'].map(norm)
E = eps.drop_duplicates('n').set_index('n')
gp = pd.read_csv(D / 'harps_gto/general_properties.csv'); gp['n'] = gp['Star'].map(norm)
P = gp.drop_duplicates('n').set_index('n')[['Teff (K)', 'logg [cm/s2]', '[Fe/H]']].set_axis(list('TgF'), axis=1)
if not Path(a.crossmatch).exists():
    raise SystemExit(f'needs {a.crossmatch}: run chemclock-matching/crossmatch_hosts.py first')
xm = pd.read_csv(a.crossmatch)
pool = set(eps.loc[eps['SwP?'] == 'No', 'n']) - set(xm.loc[xm['archive_host'].notna(), 'star'])
C = P[P.index.isin(pool)]
hosts = list(pd.read_csv(D / 'reference/final_test_v0.csv')['id'].map(norm))

# planet class as in the pipeline: largest Mp sin i, joined on the archive hd_name
arc = pd.read_csv(D / 'nasa/nasa_archive_snapshot.csv', low_memory=False)
arc['hd'] = arc['hd_name'].astype(str).map(norm)
mass = arc.dropna(subset=['pl_bmassj']).groupby('hd')['pl_bmassj'].max()
cls = np.array([('small' if mass[h] < 0.3 else 'giant') if h in mass.index else '' for h in hosts])

# ---- pair-level arrays for every host x control
Hp = P.loc[hosts].to_numpy(float); Cp = C.to_numpy(float)
dist = np.sqrt((((Hp[:, None, :] - Cp[None, :, :]) / np.array(a.scales)) ** 2).sum(-1))   # hosts x controls
hv = E.loc[hosts, [f'[{e}/Fe]' for e in ELS]].to_numpy(float); he = E.loc[hosts, [f'e_[{e}/Fe]' for e in ELS]].to_numpy(float)
cv = E.loc[C.index, [f'[{e}/Fe]' for e in ELS]].to_numpy(float); ce = E.loc[C.index, [f'e_[{e}/Fe]' for e in ELS]].to_numpy(float)
V = hv[:, None, :] - cv[None, :, :]                                   # d[X/Fe] per pair
with np.errstate(divide="ignore"):                                    # zero errors are dropped below
    W = 1 / (he[:, None, :] ** 2 + ce[None, :, :] ** 2)               # inverse variance per pair
W = np.where(np.isfinite(V) & np.isfinite(W) & (W > 0), W, 0.0); V = np.nan_to_num(V)
isp = [ELS.index(e) for e in SP]; iy, img = ELS.index('Y'), ELS.index('Mg')


def per_host(cidx):
    """s-process and [Y/Mg] age per host, matching each host to its k nearest controls in cidx."""
    order = np.argsort(dist[:, cidx], axis=1, kind='stable')[:, :a.k]
    sel = np.asarray(cidx)[order]                                      # hosts x k control indices
    r = np.arange(len(hosts))[:, None]
    v, w = V[r, sel], W[r, sel]                                        # hosts x k x elements
    s1, s0 = (v * w).sum(1), w.sum(1)                                  # per element
    sproc = s1[:, isp].sum(1) / s0[:, isp].sum(1)
    y = s1 / np.where(s0 > 0, s0, np.nan)
    return sproc, (y[:, iy] - y[:, img]) / YMG_SLOPE


def diff(x, hidx):
    x, c = x[hidx], cls[hidx]
    return x[c == 'small'].mean() - x[c == 'giant'].mean()


allc = np.arange(len(C)); allh = np.where(cls != '')[0]
base = per_host(allc)
d0 = {q: diff(x, allh) for q, x in zip(('s-process', '[Y/Mg]'), base)}


def welch(x):
    s, g = x[cls == 'small'], x[cls == 'giant']
    return (s.mean() - g.mean()) / np.sqrt(s.var(ddof=1) / len(s) + g.var(ddof=1) / len(g))


def host_draw():
    return np.concatenate([rng.choice(np.where(cls == c)[0], (cls == c).sum()) for c in ('small', 'giant')])


boot = {s: {q: [] for q in d0} for s in ('hosts', 'controls', 'both')}
for _ in range(a.n_boot):
    hb, cb = host_draw(), rng.integers(0, len(C), len(C))
    for scheme, (hi, ci) in (('hosts', (hb, None)), ('controls', (allh, cb)), ('both', (hb, cb))):
        vals = base if ci is None else per_host(ci)
        for q, x in zip(d0, vals):
            boot[scheme][q].append(diff(x, hi))

print(f'scaled matching, k = {a.k}, scales {a.scales}; {len(allh)} classified hosts '
      f'({(cls == "small").sum()} small, {(cls == "giant").sum()} giant), {len(C)} controls, {a.n_boot} bootstrap draws')
rows = []
for q in d0:
    se = {s: np.std(boot[s][q], ddof=1) for s in boot}
    rows.append(dict(quantity=q, diff=d0[q], welch_sigma=welch(base[list(d0).index(q)]),
                     se_hosts=se['hosts'], se_controls=se['controls'], se_both=se['both'],
                     sigma_hosts=d0[q] / se['hosts'], sigma_both=d0[q] / se['both']))
    unit = 'dex' if q == 's-process' else 'Gyr'
    print(f'  {q:9s} small - giant {d0[q]:+.4f} {unit}: Welch {rows[-1]["welch_sigma"]:+.2f} sigma | '
          f'bootstrap SE hosts {se["hosts"]:.4f} ({d0[q] / se["hosts"]:+.2f} sigma), controls {se["controls"]:.4f}, '
          f'both {se["both"]:.4f} ({d0[q] / se["both"]:+.2f} sigma)')
pd.DataFrame(rows).to_csv(OUT / f'control_bootstrap_k{a.k}.csv', index=False)
print(f'wrote {OUT / f"control_bootstrap_k{a.k}.csv"}')
