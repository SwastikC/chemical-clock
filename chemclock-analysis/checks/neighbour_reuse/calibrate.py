"""Is the two-way bootstrap SE right for this design? Simulate on the real host-neighbour graph.
Model: y_h = mu_class + a_h - mean_{n in N(h)} b_n.  sigma_a, sigma_b estimated from the data."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import *
rng = np.random.default_rng(7)
P = from_final_test()
idx = [ELS.index(e) for e in SP]
# pair-level s-process family values
num = (P.IV[:, idx]*P.V[:, idx]).sum(1); den = P.IV[:, idx].sum(1)
fp = np.where(den > 0, num/np.where(den > 0, den, 1), np.nan)
within = pd.Series(fp).groupby(P.ph).var(ddof=1)
sb2 = np.nanmean(within)
hv = P.host_values(np.ones((1, len(P.ph))))['s-process'][0]
resid = np.concatenate([hv[P.cls == c] - np.nanmean(hv[P.cls == c]) for c in (1, 2)])
tot = np.nanvar(resid, ddof=2)
sa2 = max(tot - sb2/3, 1e-8)
print(f'sigma_b (neighbour-to-neighbour, from within-host pair scatter) = {np.sqrt(sb2):.4f} dex')
print(f'sigma_total per host = {np.sqrt(tot):.4f};  implied sigma_a = {np.sqrt(sa2):.4f} dex')

def design_from(ph, pn, H):
    M = np.zeros((H, pn.max()+1)); np.add.at(M, (ph, pn), 1)
    return M / M.sum(1, keepdims=True)          # host x neighbour averaging matrix
cls = P.cls; keep = cls > 0
Mreal = design_from(P.ph, P.pn, P.H)
# nested design with the same counts: every pair gets its own neighbour
Mnest = design_from(P.ph, np.arange(len(P.ph)), P.H)
c = np.zeros(P.H); c[cls == 1] = 1/(cls == 1).sum(); c[cls == 2] = -1/(cls == 2).sum()
def exact_sd(M):
    # Var(c.y) = sa2*|c|^2 + sb2*|M^T c|^2
    return np.sqrt(sa2*(c**2).sum() + sb2*((M.T @ c)**2).sum())

def sim(M, nsim=600, nboot=300):
    H, U = M.shape
    pairs = np.argwhere(M > 0)
    welch_se, boot_h, boot_2 = [], [], []
    d_all = []
    for _ in range(nsim):
        av = rng.normal(0, np.sqrt(sa2), H); bv = rng.normal(0, np.sqrt(sb2), U)
        y = av - M @ bv
        s, g = y[cls == 1], y[cls == 2]
        d_all.append(s.mean() - g.mean())
        welch_se.append(np.sqrt(s.var(ddof=1)/len(s) + g.var(ddof=1)/len(g)))
        # bootstraps
        bd_h, bd_2 = [], []
        cnt = (M > 0).astype(float)
        for _b in range(nboot):
            hw = np.zeros(H)
            for cc in (1, 2):
                ii = np.where(cls == cc)[0]; np.add.at(hw, rng.choice(ii, len(ii)), 1)
            nw = np.bincount(rng.integers(0, U, U), minlength=U).astype(float)
            # host-only
            bd_h.append((hw*y*(cls == 1)).sum()/(hw*(cls == 1)).sum() - (hw*y*(cls == 2)).sum()/(hw*(cls == 2)).sum())
            # two-way: reweight neighbours inside each host's average
            Wn = cnt * nw[None]; den = Wn.sum(1)
            ok = den > 0
            yb = av - np.where(ok, (Wn @ bv)/np.where(ok, den, 1), np.nan)
            w = hw*ok
            bd_2.append((w*np.nan_to_num(yb)*(cls == 1)).sum()/(w*(cls == 1)).sum() - (w*np.nan_to_num(yb)*(cls == 2)).sum()/(w*(cls == 2)).sum())
        boot_h.append(np.std(bd_h, ddof=1)); boot_2.append(np.std(bd_2, ddof=1))
    return np.std(d_all, ddof=1), np.mean(welch_se), np.mean(boot_h), np.mean(boot_2)

for name, M in (('real graph (with reuse)', Mreal), ('nested, no reuse', Mnest)):
    t, w, bh, b2 = sim(M)
    print(f'{name:26s} exact SD {exact_sd(M):.5f} | simulated SD {t:.5f} | mean Welch SE {w:.5f} | host bootstrap {bh:.5f} | two-way bootstrap {b2:.5f}')
print(f'\nVariance inflation from reuse, exact: SD_real/SD_nested = {exact_sd(Mreal)/exact_sd(Mnest):.3f}')
# apply to observed difference
d_obs, se_obs = welch(hv, cls)[:2]
infl = exact_sd(Mreal)/exact_sd(Mnest)
print(f'Observed s-process difference {d_obs:+.4f}: Welch {d_obs/se_obs:+.2f} sigma -> reuse-corrected {d_obs/(se_obs*infl):+.2f} sigma')

print('\nExact reuse inflation for every quantity (sigma_b from within-host pair scatter of that quantity):')
single = P.host_values(np.eye(len(P.ph)))           # each pair alone -> pair-level values on its host
allv = P.host_values(np.ones((1, len(P.ph))))
rows = []
for k in allv:
    pv = single[k][np.arange(len(P.ph)), P.ph]
    sb2k = np.nanmean(pd.Series(pv).groupby(P.ph).var(ddof=1))
    yk = allv[k][0]
    r = np.concatenate([yk[(cls == cc) & np.isfinite(yk)] - np.nanmean(yk[cls == cc]) for cc in (1, 2)])
    sa2k = max(np.var(r, ddof=2) - sb2k/3, 1e-12)
    ok = np.isfinite(yk) & keep
    cc_ = np.zeros(P.H); cc_[(cls == 1) & ok] = 1/((cls == 1) & ok).sum(); cc_[(cls == 2) & ok] = -1/((cls == 2) & ok).sum()
    v = lambda M: sa2k*(cc_**2).sum() + sb2k*((M.T @ cc_)**2).sum()
    f = np.sqrt(v(Mreal)/v(Mnest))
    d, se, sg = welch(yk, cls)[:3]
    rows.append(dict(quantity=k, diff=d, welch_sigma=sg, inflation=f, corrected_sigma=sg/f))
    print(f'  {k:10s} inflation {f:.3f}   {sg:+.2f} sigma -> {sg/f:+.2f} sigma')
pd.DataFrame(rows).to_csv(f'{OUT}/reuse_inflation.csv', index=False)
