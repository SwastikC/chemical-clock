"""Pair-level data and fast statistics for the neighbour-reuse tests."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'pipeline'))
from paths import FILES
import numpy as np, pandas as pd
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')   # CSVs written by the checks
os.makedirs(OUT, exist_ok=True)
norm = lambda s: str(s).replace(' ', '').strip()
ELS = ['Zn','Cu','Na','Mn','V','Cr','Co','Ni','Mg','Si','Ca','Ti','Sc','Al','Sr','Y','Ba','Ce']
SP = ['Sr','Y','Ba','Ce']
CLOCKS = {'[Y/Mg]':('Y','Mg',-0.0410), '[Y/Al]':('Y','Al',-0.0459), '[Ba/Mg]':('Ba','Mg',-0.0450),
          '[Sr/Mg]':('Sr','Mg',-0.0420), '[Ce/Mg]':('Ce','Mg',-0.0400)}

def planet_class(hosts):
    comp = pd.read_csv(FILES['archive']); comp['hd'] = comp['hd_name'].astype(str).map(norm)
    m = (comp.dropna(subset=['pl_bmassj']).sort_values('pl_bmassj', ascending=False)
             .drop_duplicates('hd').set_index('hd')['pl_bmassj'])
    mass = np.array([m.get(h, np.nan) for h in hosts])
    cls = np.where(np.isnan(mass), 0, np.where(mass < 0.3, 1, 2))   # 1 small, 2 giant, 0 unclassified
    return cls, mass

class Pairs:
    """hosts[h], nb_names[u]; pair arrays: ph (host idx), pn (neighbour idx), V, IV (inverse variance, 0 if unusable)."""
    def __init__(self, hosts, host_of_pair, nb_of_pair, V, E):
        self.hosts = list(hosts); self.H = len(hosts)
        self.nb_names, self.pn = np.unique(np.asarray(nb_of_pair), return_inverse=True)
        self.U = len(self.nb_names)
        self.ph = np.asarray(host_of_pair)
        ok = np.isfinite(V) & np.isfinite(E) & (E > 0)
        self.V = np.where(ok, V, 0.0); self.IV = np.where(ok, 1.0/np.where(ok, E, 1)**2, 0.0)
        self.cls, self.mass = planet_class(self.hosts)

    def host_sums(self, pair_w):
        """pair_w: (B, P) pair multiplicities -> per-host sums S1=sum w*iv*v, S0=sum w*iv, shape (B, H, nel)."""
        B = pair_w.shape[0]
        S1 = np.zeros((B, self.H, len(ELS))); S0 = np.zeros_like(S1)
        wv = pair_w[:, :, None] * (self.IV * self.V)[None]; w0 = pair_w[:, :, None] * self.IV[None]
        for b in range(B):
            np.add.at(S1[b], self.ph, wv[b]); np.add.at(S0[b], self.ph, w0[b])
        return S1, S0

    def host_values(self, pair_w):
        """Return dict of per-host quantities (B, H): s-process family and clock ages (NaN where undefined)."""
        S1, S0 = self.host_sums(pair_w)
        out = {}
        idx = [ELS.index(e) for e in SP]
        num = S1[:, :, idx].sum(-1); den = S0[:, :, idx].sum(-1)
        out['s-process'] = np.where(den > 0, num/np.where(den > 0, den, 1), np.nan)
        Y = np.where(S0 > 0, S1/np.where(S0 > 0, S0, 1), np.nan)
        for name, (a, b, slope) in CLOCKS.items():
            out[name] = (Y[:, :, ELS.index(a)] - Y[:, :, ELS.index(b)]) / slope
        return out

def group_diff(x, host_w, cls):
    """x (B,H) per-host values, host_w (B,H) host multiplicities. Returns small-minus-giant weighted mean difference."""
    ok = np.isfinite(x)
    res = []
    for c in (1, 2):
        w = host_w * ok * (cls == c)[None]
        res.append((np.where(ok, x, 0)*w).sum(1)/w.sum(1))
    return res[0] - res[1]

def welch(x, cls):
    s = x[(cls == 1) & np.isfinite(x)]; g = x[(cls == 2) & np.isfinite(x)]
    d = s.mean() - g.mean(); se = np.sqrt(s.var(ddof=1)/len(s) + g.var(ddof=1)/len(g))
    return d, se, d/se, len(s), len(g)

def from_final_test():
    f = pd.read_csv(FILES['final_test'])
    hosts = [norm(x) for x in f['id']]
    ph, pn, V, E = [], [], [], []
    for h, (_, r) in enumerate(f.iterrows()):
        for k, c in zip((1,2,3), ('First','Second','Third')):
            ph.append(h); pn.append(norm(r[c]))
            V.append([r[f'{k}_[{e}/Fe]'] for e in ELS]); E.append([r[f'{k}_e_[{e}/Fe]'] for e in ELS])
    return Pairs(hosts, ph, pn, np.array(V, float), np.array(E, float))

def from_epsilon(assign):
    """assign: dict host -> list of neighbour names. Differences from epsilon.csv (exact recipe of final_test)."""
    eps = pd.read_csv(FILES['epsilon']); eps['n'] = eps['Star'].map(norm); Ep = eps.drop_duplicates('n').set_index('n')
    hosts = list(assign); ph, pn, V, E = [], [], [], []
    for h, name in enumerate(hosts):
        for n in assign[name]:
            ph.append(h); pn.append(n)
            V.append([Ep.at[name, f'[{e}/Fe]'] - Ep.at[n, f'[{e}/Fe]'] for e in ELS])
            E.append([np.hypot(Ep.at[name, f'e_[{e}/Fe]'], Ep.at[n, f'e_[{e}/Fe]']) for e in ELS])
    return Pairs(hosts, ph, pn, np.array(V, float), np.array(E, float))
