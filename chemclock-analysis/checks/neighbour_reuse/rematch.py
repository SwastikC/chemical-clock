"""Redo the rank-mean matching with and without replacement and compare the headline numbers.
The stored neighbour table cannot be regenerated exactly (see README), so both arms use the same
re-implementation; the with-replacement arm is the like-for-like control for the no-replacement arm."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import *
from scipy.stats import rankdata
from scipy.optimize import linear_sum_assignment
eps = pd.read_csv(FILES['epsilon']); eps['n'] = eps['Star'].map(norm)
POOL = set(eps.loc[eps['SwP?'] == 'No', 'n'])
HOSTS = [norm(x) for x in pd.read_csv(FILES['final_test'])['id']]

def params(src):
    g = pd.read_csv(FILES['properties']); g['n'] = g['Star'].map(norm)
    G = g.drop_duplicates('n').set_index('n')[['Teff (K)', 'logg [cm/s2]', '[Fe/H]']].set_axis(['T', 'g', 'F'], axis=1)
    if src == 'general': return G
    P = {}
    cols = ['TeffK', 'logg[cm/s2]', '[Fe/H]---_x']
    for fn, keys in ((FILES['swop'], ['Star---']), (FILES['swp'], ['Star---', 'hostname_1'])):
        for _, r in pd.read_csv(fn).dropna(subset=cols).iterrows():
            for k in keys: P.setdefault(norm(r[k]), tuple(r[cols]))
    H = pd.DataFrame.from_dict(P, orient='index', columns=['T', 'g', 'F'])
    return pd.concat([H, G[~G.index.isin(H.index)]])      # HARPS columns, General Properties fallback

def rank_table(P):
    pool = P[P.index.isin(POOL)]
    out = {}
    for h in HOSTS:
        hp = P.loc[h]; c = pool[(pool['T'] - hp['T']).abs() <= 160]
        R = np.mean([rankdata((c[k] - hp[k]).abs(), method='ordinal') for k in 'TgF'], axis=0)
        out[h] = pd.Series(R, index=c.index).sort_values(kind='stable')
    return out, pool

def with_replacement(R):
    return {h: list(r.index[:3]) for h, r in R.items()}

def without_replacement(R, pool):
    cand = list(pool.index); ci = {n: i for i, n in enumerate(cand)}
    BIG = 1e6
    C = np.full((3*len(HOSTS), len(cand)), BIG)
    for j, h in enumerate(HOSTS):
        for n, r in R[h].items():
            C[3*j:3*j+3, ci[n]] = r
    rows, cols = linear_sum_assignment(C)
    assert (C[rows, cols] < BIG).all(), 'some host could not be matched within 160 K'
    out = {h: [] for h in HOSTS}
    for r, c in zip(rows, cols): out[HOSTS[r//3]].append(cand[c])
    return out

def quality(assign, P):
    d = np.array([(P.loc[h] - P.loc[n]).abs().values for h, ns in assign.items() for n in ns])
    return np.median(d, 0), d.max(0)

def summarise(assign):
    Pr = from_epsilon(assign)
    hv = Pr.host_values(np.ones((1, len(Pr.ph))))
    return Pr, {k: welch(v[0], Pr.cls) for k, v in hv.items()}

if __name__ == '__main__':   # matching_spread.py imports the functions above without running this
    rows = []
    stored = from_final_test()
    st = {k: welch(v[0], stored.cls) for k, v in stored.host_values(np.ones((1, len(stored.ph)))).items()}
    rows.append(('stored v6 table', 'n/a', stored, st, None))
    for src in ('general', 'harps'):
        P = params(src)
        R, pool = rank_table(P)
        for name, fn in (('with replacement', with_replacement), ('without replacement', without_replacement)):
            a = fn(R, pool) if fn is without_replacement else fn(R)
            Pr, res = summarise(a)
            rows.append((name, src, Pr, res, quality(a, P)))
            pd.DataFrame([(h, *ns) for h, ns in a.items()], columns=['SwP', 'First', 'Second', 'Third']).to_csv(
                f'{OUT}/neighbours_{src}_{"wr" if fn is with_replacement else "nr"}.csv', index=False)

    Q = list(st)
    print(f"{'matching':22s} {'params':8s} {'uniq nb':>7s} {'max use':>7s} {'med dT':>6s} {'med dg':>6s} {'med dF':>6s} {'max dT':>6s} | " + ' '.join(f'{q:>9s}' for q in Q))
    out = []
    for name, src, Pr, res, qual in rows:
        use = np.bincount(Pr.pn)
        qs = f'{qual[0][0]:6.0f} {qual[0][1]:6.3f} {qual[0][2]:6.3f} {qual[1][0]:6.0f}' if qual else ' '*27
        print(f'{name:22s} {src:8s} {Pr.U:7d} {use.max():7d} {qs} | ' + ' '.join(f'{res[q][2]:+9.2f}' for q in Q))
        for q in Q:
            out.append(dict(matching=name, params=src, quantity=q, diff=res[q][0], se=res[q][1], sigma=res[q][2], n_small=res[q][3], n_giant=res[q][4]))
    print('\ndifferences (s-process in dex, clocks in Gyr):')
    print(f"{'matching':22s} {'params':8s} | " + ' '.join(f'{q:>9s}' for q in Q))
    for name, src, Pr, res, qual in rows:
        print(f'{name:22s} {src:8s} | ' + ' '.join(f'{res[q][0]:+9.4f}' if q == 's-process' else f'{res[q][0]:+9.2f}' for q in Q))
    pd.DataFrame(out).to_csv(f'{OUT}/rematch_summary.csv', index=False)
