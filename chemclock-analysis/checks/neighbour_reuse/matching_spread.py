"""How much does the result move across reasonable neighbour choices?
For each host, draw 3 of its 10 best rank-mean candidates at random; repeat. Same epsilon differences as v6."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rematch import *
rng = np.random.default_rng(11)
NDRAW = 400
eps_ = eps.drop_duplicates('n').set_index('n')
idx = [ELS.index(e) for e in SP]
for src in ('general', 'harps'):
    P = params(src); R, pool = rank_table(P)
    top10 = {h: list(r.index[:10]) for h, r in R.items()}
    # precompute pair values for host x candidate
    base = from_epsilon({h: top10[h] for h in HOSTS})     # 10 pairs per host, in order
    cls = base.cls
    sig = {k: [] for k in CLOCKS}; sig['s-process'] = []
    for _ in range(NDRAW):
        w = np.zeros((1, len(base.ph)))
        for j in range(len(HOSTS)):
            w[0, 10*j + rng.choice(10, 3, replace=False)] = 1
        for k, x in base.host_values(w).items(): sig[k].append(welch(x[0], cls)[2])
    st = welch(from_final_test().host_values(np.ones((1, 393)))['s-process'][0], cls)[2]
    s = np.array(sig['s-process'])
    print(f'{src}: s-process sigma over {NDRAW} random top-10 draws: median {np.median(s):+.2f}, 16-84% [{np.percentile(s,16):+.2f}, {np.percentile(s,84):+.2f}], min {s.min():+.2f}, max {s.max():+.2f}; '
          f'fraction at or beyond stored {st:+.2f}: {np.mean(s <= st):.4f}')
    for k in CLOCKS:
        a = np.array(sig[k]); print(f'    {k:8s} median {np.median(a):+.2f} [{np.percentile(a,16):+.2f}, {np.percentile(a,84):+.2f}]')
    pd.DataFrame(sig).to_csv(f'{OUT}/matching_spread_{src}.csv', index=False)
