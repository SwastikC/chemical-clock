"""
Extension 1: Isochrone-age vs chemical-clock-age cross-validation.

Tests whether the SwP-minus-neighbour differential ages we derive from
the [Y/Mg] chemical clock are consistent with independently-measured
isochrone ages from the same HARPS-GTO source catalogue.

This is the strongest possible internal cross-check: same stars,
two independent age methods (photospheric chemistry via Tucci Maia 2019,
vs HR-diagram position via Adibekyan 2012 isochrone fits).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import FILES, OUT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, pearsonr, spearmanr
import pickle

rng = np.random.default_rng(42)
plt.rcParams.update({
    'font.family':'serif','font.size':10,'axes.labelsize':11,'axes.titlesize':11,
    'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'axes.linewidth':0.8,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'figure.dpi':130,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT = OUT
import os
os.makedirs(OUT, exist_ok=True)
os.makedirs(f'{OUT}/figures', exist_ok=True)

# =====================================================================
# Load everything
# =====================================================================
nb_table  = pd.read_csv(FILES['neighbours'])
swp_full  = pd.read_csv(FILES['swp'])
swop_full = pd.read_csv(FILES['swop'])
df_f      = pd.read_csv(FILES['final_test'])
prop      = pd.read_csv(FILES['properties'])
comp      = pd.read_csv(FILES['archive'])

def norm(s):
    return str(s).replace(' ','').strip()

def parse_age(s):
    if pd.isna(s): return np.nan, np.nan
    try:
        parts = str(s).split()
        return float(parts[0]), float(parts[1])
    except: return np.nan, np.nan

# Parse ages from both tables
swop_full['Star_n'] = swop_full['Star---'].apply(norm)
swp_full[['age_iso','age_iso_err']]  = swp_full['Age (e)Gyr'].apply(lambda s: pd.Series(parse_age(s)))
swop_full[['age_iso','age_iso_err']] = swop_full['Age (e)Gyr'].apply(lambda s: pd.Series(parse_age(s)))

# Hosts are found by common name (hostname_1) or HD identifier (Star---), as in
# build_stellar_offsets.py. hostname_1 alone misses tau Cet, eps Eri, 61 Vir, ...
swp_keys = pd.concat([swp_full.assign(Star_n=swp_full[c].apply(norm)) for c in ('hostname_1', 'Star---')])
swp_lookup  = swp_keys.drop_duplicates('Star_n').set_index('Star_n')[['age_iso','age_iso_err']]
swop_lookup = swop_full.drop_duplicates('Star_n').set_index('Star_n')[['age_iso','age_iso_err']]

# =====================================================================
# Compute SwP-minus-neighbour isochrone differential ages
# =====================================================================
records = []
for _, r in nb_table.iterrows():
    swp = norm(r['SwP'])
    if swp not in swp_lookup.index: continue
    a_swp, e_swp = swp_lookup.loc[swp, ['age_iso','age_iso_err']]
    if not np.isfinite(a_swp): continue

    nb_ages, nb_errs = [], []
    for nb in [norm(r['First']), norm(r['Second']), norm(r['Third'])]:
        if nb in swop_lookup.index:
            a_nb, e_nb = swop_lookup.loc[nb, ['age_iso','age_iso_err']]
            if np.isfinite(a_nb):
                nb_ages.append(a_nb); nb_errs.append(e_nb)
    if not nb_ages: continue
    
    # Inverse-variance weighted neighbour age
    nb_ages = np.array(nb_ages); nb_errs = np.array(nb_errs)
    valid = nb_errs > 0
    if not valid.any(): continue
    w = 1/nb_errs[valid]**2
    a_nb_mean = (nb_ages[valid]*w).sum()/w.sum()
    a_nb_err  = 1/np.sqrt(w.sum())
    
    d_age_iso = a_swp - a_nb_mean
    d_age_iso_err = np.sqrt(e_swp**2 + a_nb_err**2)
    
    records.append({
        'Star': swp,
        'age_swp': a_swp, 'age_swp_err': e_swp,
        'age_nb_mean': a_nb_mean, 'age_nb_err': a_nb_err,
        'd_age_iso': d_age_iso, 'd_age_iso_err': d_age_iso_err,
        'n_nb_with_age': len(nb_ages),
    })

iso_df = pd.DataFrame(records)
print(f"Stars with isochrone-age differential: {len(iso_df)}")
print(f"  with 3 neighbours having ages: {(iso_df['n_nb_with_age']==3).sum()}")
print(f"  with 2 neighbours having ages: {(iso_df['n_nb_with_age']==2).sum()}")
print(f"  with 1 neighbour having age:   {(iso_df['n_nb_with_age']==1).sum()}")

# =====================================================================
# Now compute [Y/Mg] differential ages for the SAME stars
# =====================================================================
elem_cols = [c for c in df_f.columns
             if c.startswith('1_[') and c.endswith('/Fe]')
             and 'e_' not in c and c != '1_[Fe/H]']
elements = sorted(set(c.replace('1_[','').replace('/Fe]','') for c in elem_cols))
elements = ['Zn','Cu','Na','Mn','Cr','Si','Mg','Co','Ni','V','Ba','Ce','Ca','Sc','Sr','Y','Ti','Al']
# Use proper order
tc = pd.read_csv(FILES['tc'])
tc_map = dict(zip(tc['Element'], tc['Temperature']))
elements = sorted([e for e in elements if e in tc_map], key=lambda e: tc_map[e])

def get_yE(row, n):
    y = np.array([row[f'{n}_[{e}/Fe]'] for e in elements], dtype=float)
    e = np.array([row[f'{n}_e_[{e}/Fe]'] for e in elements], dtype=float)
    return y, e

def combine_neighbours(row):
    Y = np.zeros((3, len(elements))); E = np.zeros((3, len(elements)))
    for k, n in enumerate((1,2,3)):
        Y[k], E[k] = get_yE(row, n)
    invvar = np.where(np.isfinite(E) & (E>0), 1.0/np.where(E>0, E, np.nan)**2, 0.0)
    w_sum = invvar.sum(0)
    y_comb = np.where(w_sum>0, np.nansum(Y*invvar,0)/np.where(w_sum>0,w_sum,1), np.nan)
    e_comb = np.where(w_sum>0, 1.0/np.sqrt(np.where(w_sum>0,w_sum,1)), np.nan)
    return y_comb, e_comb

TM_SLOPE = -24.39  # Gyr/dex = 1/(-0.041 dex/Gyr), Tucci Maia 2019

ymg_records = []
for _, row in df_f.iterrows():
    star = row['id']
    Y, E = combine_neighbours(row)
    iy = elements.index('Y'); img = elements.index('Mg')
    if not (np.isfinite(Y[iy]) and np.isfinite(Y[img])): continue
    dymg = Y[iy] - Y[img]
    dymg_err = np.sqrt(E[iy]**2 + E[img]**2)
    ymg_records.append({
        'Star': star,
        'd_age_ymg': TM_SLOPE * dymg,
        'd_age_ymg_err': abs(TM_SLOPE) * dymg_err,
    })
ymg_df = pd.DataFrame(ymg_records)

# =====================================================================
# Merge & compare
# =====================================================================
merged = iso_df.merge(ymg_df, on='Star', how='inner')
# Planet class
comp['hd_name'] = comp['hd_name'].astype(str).str.replace(' ','')
host = (comp.dropna(subset=['pl_bmassj'])
            .sort_values('pl_bmassj', ascending=False)
            .drop_duplicates('hd_name', keep='first')
            [['hd_name','pl_bmassj']])
merged = merged.merge(host, left_on='Star', right_on='hd_name', how='left')
merged['planet_class'] = pd.NA
merged.loc[merged['pl_bmassj'] < 0.3, 'planet_class'] = 'small'
merged.loc[merged['pl_bmassj'] >=0.3, 'planet_class'] = 'giant'

print(f"\nStars with BOTH isochrone Δage and [Y/Mg] Δage: {len(merged)}")
print(f"  small: {(merged['planet_class']=='small').sum()}")
print(f"  giant: {(merged['planet_class']=='giant').sum()}")

# =====================================================================
# Test 1: Direct correlation between the two age methods
# =====================================================================
print("\n" + "="*78)
print("TEST 1: Per-star correlation, isochrone Δage vs [Y/Mg] Δage")
print("="*78)
x = merged['d_age_iso'].values
y = merged['d_age_ymg'].values
xe = merged['d_age_iso_err'].values
ye = merged['d_age_ymg_err'].values

valid = np.isfinite(x) & np.isfinite(y)
x, y, xe, ye = x[valid], y[valid], xe[valid], ye[valid]

r_pear, p_pear = pearsonr(x, y)
r_spear, p_spear = spearmanr(x, y)
slope, ic, r_, p_, se = linregress(x, y)

print(f"  N = {len(x)}")
print(f"  Pearson r  = {r_pear:+.3f}, p = {p_pear:.4f}")
print(f"  Spearman r = {r_spear:+.3f}, p = {p_spear:.4f}")
print(f"  Linear fit: [Y/Mg] Δage = ({slope:+.3f} ± {se:.3f}) · iso Δage + ({ic:+.3f})")

# =====================================================================
# Test 2: Group-level small-vs-giant test for ISOCHRONE Δage
# =====================================================================
print("\n" + "="*78)
print("TEST 2: Isochrone-age small-vs-giant differential")
print("="*78)
small = merged[merged['planet_class']=='small']
giant = merged[merged['planet_class']=='giant']

for label, group, col in [
    ('iso Δage', 'small', 'd_age_iso'),
    ('iso Δage', 'giant', 'd_age_iso'),
    ('Y/Mg Δage', 'small', 'd_age_ymg'),
    ('Y/Mg Δage', 'giant', 'd_age_ymg'),
]:
    g = small if group=='small' else giant
    vals = g[col].dropna().values
    m, se_ = vals.mean(), vals.std(ddof=1)/np.sqrt(len(vals))
    print(f"  {label:<12} {group}: mean = {m:+.3f} ± {se_:.3f} Gyr   ({m/se_:+.1f}σ)")

# Two-sample permutation
def perm(a, b, n=20000):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b]); na = len(a)
    diffs = np.empty(n)
    for i in range(n):
        rng.shuffle(pool); diffs[i] = pool[:na].mean() - pool[na:].mean()
    p = (np.abs(diffs) >= np.abs(obs)).mean()
    se_pool = np.sqrt(a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))
    return obs, obs/se_pool, p, len(a), len(b)

iso_diff, iso_sig, iso_p, na, nb = perm(small['d_age_iso'].values, giant['d_age_iso'].values)
print(f"\nIsochrone:  Δ(S-G) = {iso_diff:+.3f} Gyr  σ = {iso_sig:+.2f}  p_perm = {iso_p:.4f}")
ymg_diff, ymg_sig, ymg_p, _, _ = perm(small['d_age_ymg'].values, giant['d_age_ymg'].values)
print(f"[Y/Mg]:     Δ(S-G) = {ymg_diff:+.3f} Gyr  σ = {ymg_sig:+.2f}  p_perm = {ymg_p:.4f}")

# =====================================================================
# Test 3: Combined age estimator (variance-weighted average of the two methods)
# =====================================================================
print("\n" + "="*78)
print("TEST 3: Variance-weighted combined-age estimator")
print("="*78)
both = merged.dropna(subset=['d_age_iso','d_age_ymg','d_age_iso_err','d_age_ymg_err']).copy()
# replace zero/negative errors with median
for c in ['d_age_iso_err','d_age_ymg_err']:
    med = both[c].median()
    both.loc[both[c] <= 0, c] = med

w_iso = 1/both['d_age_iso_err']**2
w_ymg = 1/both['d_age_ymg_err']**2
both['d_age_comb'] = (w_iso*both['d_age_iso'] + w_ymg*both['d_age_ymg']) / (w_iso + w_ymg)
both['d_age_comb_err'] = 1/np.sqrt(w_iso + w_ymg)

small_c = both[both['planet_class']=='small']['d_age_comb'].values
giant_c = both[both['planet_class']=='giant']['d_age_comb'].values
comb_diff, comb_sig, comb_p, ncs, ncg = perm(small_c, giant_c)
print(f"Combined: Δ(S-G) = {comb_diff:+.3f} Gyr  σ = {comb_sig:+.2f}  p_perm = {comb_p:.4f}")
print(f"  Small N = {ncs}, mean = {small_c.mean():+.3f} ± {small_c.std(ddof=1)/np.sqrt(ncs):.3f}")
print(f"  Giant N = {ncg}, mean = {giant_c.mean():+.3f} ± {giant_c.std(ddof=1)/np.sqrt(ncg):.3f}")

# Save
merged.to_csv(f'{OUT}/iso_ymg_comparison.csv', index=False)
with open(f'{OUT}/_state_ext1.pkl','wb') as f:
    pickle.dump({
        'merged': merged, 'both': both,
        'r_pearson': r_pear, 'p_pearson': p_pear, 'slope': slope, 'slope_se': se,
        'iso_diff': iso_diff, 'iso_sig': iso_sig, 'iso_p': iso_p,
        'ymg_diff': ymg_diff, 'ymg_sig': ymg_sig, 'ymg_p': ymg_p,
        'comb_diff': comb_diff, 'comb_sig': comb_sig, 'comb_p': comb_p,
    }, f)
print("\n✓ Done")
