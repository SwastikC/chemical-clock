"""Neighbour-reuse tests on the stored v6 neighbour table (131 hosts, 393 pairs)."""
import sys; sys.path.insert(0, '.')
from core import *
from scipy.stats import norm as N01
rng = np.random.default_rng(20261025)
B = 10000
P = from_final_test()
nP = len(P.ph)

# ---- sharing structure
use = np.bincount(P.pn, minlength=P.U)
pcls = P.cls[P.ph]
by = {c: np.bincount(P.pn[pcls == c], minlength=P.U) for c in (1, 2)}
print('== Sharing structure (classified hosts only)')
print('neighbour use counts (all 131 hosts):', dict(zip(*np.unique(use, return_counts=True))))
cross = (by[1] > 0) & (by[2] > 0)
print(f'neighbours used by >1 small host: {(by[1]>1).sum()},  >1 giant host: {(by[2]>1).sum()},  by both classes: {cross.sum()}')
print('most used:', [(P.nb_names[i], int(use[i])) for i in np.argsort(-use)[:6]])
# connected components of the host-neighbour graph
parent = list(range(P.H + P.U))
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
for h, u in zip(P.ph, P.pn): parent[find(h)] = find(P.H + u)
comp = np.array([find(h) for h in range(P.H)])
sizes = np.sort(np.unique(comp, return_counts=True)[1])[::-1]
print(f'host clusters (hosts linked through shared neighbours): {len(sizes)}; largest sizes {sizes[:6].tolist()}; singletons {(sizes==1).sum()}')

# ---- bootstraps
obs = {k: welch(v[0], P.cls) for k, v in P.host_values(np.ones((1, nP))).items()}
def host_draw(B):
    w = np.zeros((B, P.H))
    for c in (0, 1, 2):
        idx = np.where(P.cls == c)[0]
        draws = rng.integers(0, len(idx), size=(B, len(idx)))
        for b in range(B): np.add.at(w[b], idx[draws[b]], 1)
    return w
def nb_draw(B):
    w = np.zeros((B, P.U))
    draws = rng.integers(0, P.U, size=(B, P.U))
    for b in range(B): np.add.at(w[b], draws[b], 1)
    return w
results = {}
CH = 1000
for scheme in ('hosts only', 'neighbours only', 'two-way'):
    diffs = {k: [] for k in obs}
    for start in range(0, B, CH):
        hw = host_draw(CH) if scheme != 'neighbours only' else np.ones((CH, P.H))
        nw = nb_draw(CH) if scheme != 'hosts only' else np.ones((CH, P.U))
        hv = P.host_values(nw[:, P.pn])
        for k, x in hv.items(): diffs[k].append(group_diff(x, hw, P.cls))
    results[scheme] = {k: np.concatenate(v) for k, v in diffs.items()}

print(f'\n== Bootstrap standard errors of the small-minus-giant difference (B = {B})')
print(f"{'quantity':10s} {'diff':>8s} {'Welch SE':>9s} {'sig':>6s} | {'SE host':>8s} {'SE nb':>8s} {'SE 2way':>8s} {'sig 2way':>9s} {'infl':>5s} {'frac same sign':>15s}")
rows = []
for k, (d, se, s, ns, ng) in obs.items():
    se_h, se_n, se_2 = (np.nanstd(results[sc][k], ddof=1) for sc in ('hosts only', 'neighbours only', 'two-way'))
    fs = np.mean(np.sign(results['two-way'][k]) == np.sign(d))
    lo, hi = np.nanpercentile(results['two-way'][k], [2.5, 97.5])
    print(f'{k:10s} {d:+8.4f} {se:9.4f} {s:+6.2f} | {se_h:8.4f} {se_n:8.4f} {se_2:8.4f} {d/se_2:+9.2f} {se_2/se:5.2f} {fs:15.4f}')
    rows.append(dict(quantity=k, diff=d, welch_se=se, welch_sigma=s, se_hosts=se_h, se_neighbours=se_n,
                     se_twoway=se_2, sigma_twoway=d/se_2, ci95_lo=lo, ci95_hi=hi, frac_same_sign=fs))
pd.DataFrame(rows).to_csv('reuse_bootstrap_summary.csv', index=False)

# ---- leave-one-neighbour-out
print('\n== Leave-one-neighbour-out (each of the', P.U, 'neighbours removed from every host that uses it)')
W = np.ones((P.U, nP)); W[P.pn, np.arange(nP)] = 0
hv = P.host_values(W)
lono = {k: np.array([welch(x[i], P.cls)[2] for i in range(P.U)]) for k, x in hv.items()}
for k, s in lono.items():
    i_min, i_max = np.argmin(np.abs(s)), np.argmax(np.abs(s))
    print(f'{k:10s} sigma range {s.min():+.2f} to {s.max():+.2f}; weakest after dropping {P.nb_names[i_min]} ({s[i_min]:+.2f})')
i = list(P.nb_names).index('HD139590')
print('dropping HD139590 (used', int(use[i]), 'times):', {k: round(float(v[i]), 2) for k, v in lono.items()})
top = np.argsort(-use)[:10]
Wt = np.ones((1, nP)); Wt[0, np.isin(P.pn, top)] = 0
print('dropping the 10 most-reused neighbours at once:', {k: round(welch(x[0], P.cls)[2], 2) for k, x in P.host_values(Wt).items()})
