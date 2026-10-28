"""
Comprehensive robustness analysis for the chemical-clock paper.
Addresses 10 anticipated referee objections with explicit tests.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import FILES, OUT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, wilcoxon, ttest_1samp, ks_2samp
from scipy.special import gammaln
import pickle

rng = np.random.default_rng(42)
N_PERM = 20000
N_BOOT = 5000

plt.rcParams.update({
    'font.family':'serif','font.size':10,'axes.labelsize':11,'axes.titlesize':11,
    'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'axes.linewidth':0.8,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'figure.dpi':130,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT  = OUT
import os
os.makedirs(f'{OUT}/figures', exist_ok=True)

# Load
df_f = pd.read_csv(FILES['final_test'])
tc   = pd.read_csv(FILES['tc'])
prop = pd.read_csv(FILES['properties'])
comp = pd.read_csv(FILES['archive'])
stellar_offsets = pd.read_csv(OUT+'/paper_stellar_offsets.csv')
# neighbour blocks in the host-minus-neighbour table: 1_, 2_, 3_, ... (any number)
NB = sorted({int(c.split('_')[0]) for c in df_f.columns if c.split('_')[0].isdigit()})

elem_cols = [c for c in df_f.columns
             if c.startswith('1_[') and c.endswith('/Fe]')
             and 'e_' not in c and c != '1_[Fe/H]']
elements_present = sorted(set(c.replace('1_[','').replace('/Fe]','') for c in elem_cols))
tc_map = dict(zip(tc['Element'], tc['Temperature']))
elements = sorted([e for e in elements_present if e in tc_map],
                  key=lambda e: tc_map[e])
T_arr = np.array([tc_map[e] for e in elements])

FAMILIES = {
    'Volatile'      : ['Zn','Cu','Na','Mn'],
    'Iron-peak ref' : ['V','Cr','Co','Ni'],
    'Alpha + Al ref': ['Mg','Si','Ca','Ti','Sc','Al'],
    's-process'     : ['Sr','Y','Ba','Ce'],
}

def get_yE(row, n):
    y = np.array([row[f'{n}_[{e}/Fe]'] for e in elements], dtype=float)
    e = np.array([row[f'{n}_e_[{e}/Fe]'] for e in elements], dtype=float)
    return y, e

def combine_neighbours(row):
    Y = np.zeros((len(NB), len(elements))); E = np.zeros((len(NB), len(elements)))
    for k, n in enumerate(NB):
        Y[k], E[k] = get_yE(row, n)
    invvar = np.where(np.isfinite(E) & (E>0), 1.0/np.where(E>0, E, np.nan)**2, 0.0)
    w_sum = invvar.sum(0)
    y_comb = np.where(w_sum>0, np.nansum(Y*invvar,0)/np.where(w_sum>0,w_sum,1), np.nan)
    e_comb = np.where(w_sum>0, 1.0/np.sqrt(np.where(w_sum>0,w_sum,1)), np.nan)
    return y_comb, e_comb

# Build per-star table
records = []
for _, row in df_f.iterrows():
    star = row['id']
    y, e = combine_neighbours(row)
    rec = {'Star': star}
    for j, el in enumerate(elements):
        rec[f'd_{el}'] = y[j]; rec[f'e_{el}'] = e[j]
    for fam, es in FAMILIES.items():
        idx = [elements.index(el) for el in es]
        yf, ef = y[idx], e[idx]
        m = np.isfinite(yf) & np.isfinite(ef) & (ef>0)
        if m.sum():
            w = 1/ef[m]**2
            rec[f'fam_{fam}'] = (yf[m]*w).sum()/w.sum()
        else:
            rec[f'fam_{fam}'] = np.nan
    records.append(rec)
per_star = pd.DataFrame(records)

# Merge stellar params + planet info
per_star = per_star.merge(prop[['Star','Teff (K)','logg [cm/s2]','[Fe/H]']],
                          on='Star', how='left')
comp['hd_name'] = comp['hd_name'].astype(str).str.replace(' ','')
host = (comp.dropna(subset=['pl_bmassj'])
            .sort_values('pl_bmassj', ascending=False)
            .drop_duplicates('hd_name', keep='first')
            [['hd_name','pl_bmassj','pl_radj']])
per_star = per_star.merge(host, left_on='Star', right_on='hd_name', how='left')
per_star['planet_class'] = pd.NA
per_star.loc[per_star['pl_bmassj'] < 0.3, 'planet_class'] = 'small'
per_star.loc[per_star['pl_bmassj'] >=0.3, 'planet_class'] = 'giant'
# Merge stellar offsets
per_star = per_star.merge(stellar_offsets, left_on='Star', right_on='SwP', how='left')

print(f"Per-star table: {len(per_star)} rows")
print(f"With planet class: {per_star['planet_class'].notna().sum()}")
print(f"With stellar offsets: {per_star['dTeff'].notna().sum()}")

small_df = per_star[per_star['planet_class']=='small']
giant_df = per_star[per_star['planet_class']=='giant']
print(f"Small: {len(small_df)}, Giant: {len(giant_df)}")

# =====================================================================
# OBJECTION 1: Could this be Δ[Fe/H] driving s-process, not age?
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 1: Is the s-process signal driven by residual Δ[Fe/H]?")
print("="*78)

both = per_star.dropna(subset=['fam_s-process','dFeH','planet_class']).copy()
print(f"Stars with both s-process and Δ[Fe/H]: {len(both)}")

# Correlation between Δ[Fe/H] and Δ[s-process]
sl, ic, r, p, se = linregress(both['dFeH'], both['fam_s-process'])
print(f"\nLinear regression: Δ[s-process] = β · Δ[Fe/H] + α")
print(f"  β (slope)  = {sl:+.3f} ± {se:.3f} dex/dex")
print(f"  Pearson r  = {r:+.3f}, p = {p:.4f}")

# Partial: regress out Δ[Fe/H] from Δ[s-process], then redo small-vs-giant
both['fam_sp_resid'] = both['fam_s-process'] - (sl * both['dFeH'] + ic)

s_resid = both[both['planet_class']=='small']['fam_sp_resid'].values
g_resid = both[both['planet_class']=='giant']['fam_sp_resid'].values
diff_resid = s_resid.mean() - g_resid.mean()
se_resid = np.sqrt(s_resid.var(ddof=1)/len(s_resid) + g_resid.var(ddof=1)/len(g_resid))
sigma_resid = diff_resid / se_resid

# Permutation
obs = diff_resid
pool = np.concatenate([s_resid, g_resid]); na = len(s_resid)
diffs = np.empty(N_PERM)
for i in range(N_PERM):
    rng.shuffle(pool); diffs[i] = pool[:na].mean() - pool[na:].mean()
p_resid = (np.abs(diffs) >= np.abs(obs)).mean()

print(f"\nAfter regressing out Δ[Fe/H]:")
print(f"  Small: {s_resid.mean():+.4f} dex, Giant: {g_resid.mean():+.4f} dex")
print(f"  Difference (S-G): {diff_resid:+.4f} dex   σ = {sigma_resid:+.2f}   p = {p_resid:.5f}")
obj1 = dict(beta=sl, beta_err=se, r=r, p_corr=p, sigma_resid=sigma_resid, p_resid=p_resid)

# =====================================================================
# OBJECTION 2: Selection bias - are RV-quiet stars older?
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 2: Δlogg as proxy for age - does it absorb signal?")
print("="*78)

both = per_star.dropna(subset=['fam_s-process','dlogg','planet_class']).copy()

# logg increases as stars evolve UP MS slightly; more relevant: Teff drops
# Multivariate regression: Δ[s-process] = β_T·ΔT + β_g·Δlogg + β_M·ΔFeH
from numpy.linalg import lstsq
X = np.column_stack([np.ones(len(both)),
                     both['dTeff'].values,
                     both['dlogg'].values,
                     both['dFeH'].values])
y = both['fam_s-process'].values
coefs, *_ = lstsq(X, y, rcond=None)
y_pred = X @ coefs
resid = y - y_pred
print(f"Multivariate regression: Δ[s-process] ~ ΔTeff + Δlogg + Δ[Fe/H]")
print(f"  β_T   = {coefs[1]:+.5e} dex/K")
print(f"  β_g   = {coefs[2]:+.4f} dex/dex")
print(f"  β_M   = {coefs[3]:+.4f} dex/dex")
print(f"  R²    = {1 - resid.var()/y.var():.3f}")

# Now test small vs giant on the residuals
both['sp_clean'] = resid
s_clean = both[both['planet_class']=='small']['sp_clean'].values
g_clean = both[both['planet_class']=='giant']['sp_clean'].values
diff_clean = s_clean.mean() - g_clean.mean()
se_clean = np.sqrt(s_clean.var(ddof=1)/len(s_clean) + g_clean.var(ddof=1)/len(g_clean))
sigma_clean = diff_clean/se_clean
obs = diff_clean
pool = np.concatenate([s_clean, g_clean]); na = len(s_clean)
diffs = np.empty(N_PERM)
for i in range(N_PERM):
    rng.shuffle(pool); diffs[i] = pool[:na].mean() - pool[na:].mean()
p_clean = (np.abs(diffs) >= np.abs(obs)).mean()

print(f"\nAfter regressing out ΔTeff, Δlogg, Δ[Fe/H]:")
print(f"  Small: {s_clean.mean():+.4f} dex, Giant: {g_clean.mean():+.4f} dex")
print(f"  Difference: {diff_clean:+.4f} dex   σ = {sigma_clean:+.2f}   p = {p_clean:.5f}")
obj2 = dict(coefs=coefs, sigma_clean=sigma_clean, p_clean=p_clean)

# =====================================================================
# OBJECTION 3: Y/Mg clock could be wrong — do other clocks agree?
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 3: Cross-check chemical clocks - Y/Mg, Y/Al, Ba/Mg, Sr/Mg, Ce/Mg")
print("="*78)

# Tucci Maia 2019 + Spina 2018 calibrations (in dex per Gyr)
# All have similar slopes; let's compute for each
CLOCKS = {
    # Slopes in dex/Gyr. [Y/Mg]: Tucci Maia et al. 2016 (= Spina et al. 2016).
    # [Y/Al]: Spina et al. 2016 direct fit.
    # [Ba/Mg], [Ce/Mg]: Spina et al. 2016 [X/Fe]-age slopes minus their [Mg/Fe] slope (+0.0100).
    # [Sr/Mg]: Spina et al. 2018 [Sr/Fe] slope (-0.032) minus Spina et al. 2016 [Mg/Fe] slope.
    '[Y/Mg]'  : ('Y', 'Mg', -0.0410),
    '[Y/Al]'  : ('Y', 'Al', -0.0459),
    '[Ba/Mg]' : ('Ba','Mg', -0.0450),
    '[Sr/Mg]' : ('Sr','Mg', -0.0420),
    '[Ce/Mg]' : ('Ce','Mg', -0.0400),
}

clock_results = {}
for clock_name, (e1, e2, slope_per_gyr) in CLOCKS.items():
    s_dage, g_dage = [], []
    for _, row in df_f.iterrows():
        star = row['id']
        # Inverse-variance combine each element across the 3 neighbours
        # (same weighting as the family decomposition), then difference.
        def ivw(elem):
            v = np.array([row.get(f'{n}_[{elem}/Fe]', np.nan) for n in NB], float)
            e = np.array([row.get(f'{n}_e_[{elem}/Fe]', np.nan) for n in NB], float)
            m = np.isfinite(v) & np.isfinite(e) & (e > 0)
            if not m.any(): return np.nan
            w = 1/e[m]**2
            return (v[m]*w).sum()/w.sum()
        c1, c2 = ivw(e1), ivw(e2)
        if not (np.isfinite(c1) and np.isfinite(c2)):
            continue
        d_mean = c1 - c2
        # Translate to age: Δ[X/Y] = slope * Δage  →  Δage = Δ[X/Y]/slope
        d_age = d_mean / slope_per_gyr
        cls = per_star.loc[per_star['Star']==star, 'planet_class']
        if len(cls) == 0: continue
        cls = cls.values[0]
        if pd.isna(cls): continue
        if cls == 'small': s_dage.append(d_age)
        elif cls == 'giant': g_dage.append(d_age)
    s_dage = np.array(s_dage); g_dage = np.array(g_dage)
    clock_results[clock_name] = {
        'small_mean': s_dage.mean(), 'small_se': s_dage.std(ddof=1)/np.sqrt(len(s_dage)),
        'giant_mean': g_dage.mean(), 'giant_se': g_dage.std(ddof=1)/np.sqrt(len(g_dage)),
        'n_small': len(s_dage), 'n_giant': len(g_dage),
    }
    diff = s_dage.mean() - g_dage.mean()
    se = np.sqrt(s_dage.var(ddof=1)/len(s_dage) + g_dage.var(ddof=1)/len(g_dage))
    sigma = diff/se
    clock_results[clock_name]['diff'] = diff
    clock_results[clock_name]['sigma'] = sigma
    print(f"  {clock_name:<10}: Δage_small = {s_dage.mean():+.2f}±{s_dage.std(ddof=1)/np.sqrt(len(s_dage)):.2f}, "
          f"Δage_giant = {g_dage.mean():+.2f}±{g_dage.std(ddof=1)/np.sqrt(len(g_dage)):.2f}, "
          f"Δ(S-G) = {diff:+.2f} Gyr ({sigma:+.1f}σ)")

# =====================================================================
# OBJECTION 4: Cherry-picked family boundaries? Element-by-element
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 4: Per-element 1-sample sigma for small/giant")
print("="*78)

per_elem_data = []
for el in elements:
    s = small_df[f'd_{el}'].dropna().values
    g = giant_df[f'd_{el}'].dropna().values
    s_mean = s.mean(); s_se = s.std(ddof=1)/np.sqrt(len(s))
    g_mean = g.mean(); g_se = g.std(ddof=1)/np.sqrt(len(g))
    diff = s_mean - g_mean
    sediff = np.sqrt(s.var(ddof=1)/len(s) + g.var(ddof=1)/len(g))
    sigma_diff = diff/sediff
    per_elem_data.append({
        'element': el, 'TC': tc_map[el],
        's_mean': s_mean, 's_se': s_se, 's_sigma': s_mean/s_se,
        'g_mean': g_mean, 'g_se': g_se, 'g_sigma': g_mean/g_se,
        'diff': diff, 'diff_sigma': sigma_diff,
    })
per_elem_df = pd.DataFrame(per_elem_data).sort_values('TC')
print(per_elem_df[['element','TC','s_sigma','g_sigma','diff_sigma']].to_string(index=False))

# Leave-one-out: drop each s-process element, redo family test
print("\nLeave-one-element-out (s-process family stability):")
sproc = FAMILIES['s-process']
for skip in [None] + sproc:
    elems_use = [e for e in sproc if e != skip]
    # recompute small s-process mean dropping this element
    s_vals, g_vals = [], []
    for _, row in df_f.iterrows():
        star = row['id']
        cls = per_star.loc[per_star['Star']==star, 'planet_class'].values
        if len(cls)==0 or pd.isna(cls[0]): continue
        Y, E = combine_neighbours(row)
        idx = [elements.index(el) for el in elems_use]
        yf, ef = Y[idx], E[idx]
        m = np.isfinite(yf) & np.isfinite(ef) & (ef>0)
        if m.sum()==0: continue
        w = 1/ef[m]**2
        val = (yf[m]*w).sum()/w.sum()
        if cls[0]=='small': s_vals.append(val)
        elif cls[0]=='giant': g_vals.append(val)
    s_vals, g_vals = np.array(s_vals), np.array(g_vals)
    diff = s_vals.mean() - g_vals.mean()
    se = np.sqrt(s_vals.var(ddof=1)/len(s_vals) + g_vals.var(ddof=1)/len(g_vals))
    label = f"all 4 ({','.join(sproc)})" if skip is None else f"drop {skip}"
    print(f"  {label:<28}  Δ(S-G) = {diff:+.4f}   σ = {diff/se:+.2f}")

# =====================================================================
# OBJECTION 5: Robustness to T_C cut
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 5: T_C split robustness — grid from 900 to 1500 K")
print("="*78)

tc_grid = np.arange(900, 1501, 50)
robustness = []
for tc_split in tc_grid:
    vol_idx = [i for i, t in enumerate(T_arr) if t < tc_split]
    ref_idx = [i for i, t in enumerate(T_arr) if t >= tc_split]
    s_rv, g_rv = [], []
    for _, row in df_f.iterrows():
        star = row['id']
        cls = per_star.loc[per_star['Star']==star, 'planet_class'].values
        if len(cls)==0 or pd.isna(cls[0]): continue
        Y, E = combine_neighbours(row)
        def wmean(idx):
            yi, ei = Y[idx], E[idx]
            m = np.isfinite(yi) & np.isfinite(ei) & (ei>0)
            if m.sum()==0: return np.nan
            w = 1/ei[m]**2
            return (yi[m]*w).sum()/w.sum()
        if not vol_idx or not ref_idx: continue
        rv = wmean(ref_idx) - wmean(vol_idx)
        if cls[0]=='small': s_rv.append(rv)
        elif cls[0]=='giant': g_rv.append(rv)
    s_rv, g_rv = np.array(s_rv), np.array(g_rv)
    diff = s_rv.mean() - g_rv.mean()
    se = np.sqrt(s_rv.var(ddof=1)/len(s_rv) + g_rv.var(ddof=1)/len(g_rv))
    robustness.append({'tc_split': tc_split, 'diff': diff, 'sigma': diff/se,
                       'n_vol': len(vol_idx), 'n_ref': len(ref_idx)})
robust_df = pd.DataFrame(robustness)
print(robust_df.to_string(index=False))

# =====================================================================
# OBJECTION 6: Comparison with a Bedell-style Tc-slope fit
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 6: Bedell-style Tc-slope fit (per-star linear regression)")
print("="*78)

slopes_full, slopes_pure = [], []  # full = include s-process; pure = drop
sproc_indices = [elements.index(e) for e in FAMILIES['s-process']]
for _, row in df_f.iterrows():
    star = row['id']
    cls = per_star.loc[per_star['Star']==star, 'planet_class'].values
    if len(cls)==0 or pd.isna(cls[0]): continue
    Y, E = combine_neighbours(row)
    m = np.isfinite(Y) & np.isfinite(E) & (E>0)
    if m.sum() < 4: continue
    w = 1/E[m]**2
    # Full
    sl, ic, r, p, _ = linregress(T_arr[m], Y[m])
    slopes_full.append({'star': star, 'class': cls[0], 'slope': sl})
    # Pure (drop s-process)
    m_pure = m.copy()
    for i in sproc_indices: m_pure[i] = False
    if m_pure.sum() >= 4:
        sl, *_ = linregress(T_arr[m_pure], Y[m_pure])
        slopes_pure.append({'star': star, 'class': cls[0], 'slope': sl})

slopes_full = pd.DataFrame(slopes_full)
slopes_pure = pd.DataFrame(slopes_pure)

print("Tc-slope results (units: 10⁻⁴ dex/K):")
for label, df_sl in [("Full (with s-process)", slopes_full),
                     ("Pure (no s-process)",   slopes_pure)]:
    s = df_sl[df_sl['class']=='small']['slope'].values * 1e4
    g = df_sl[df_sl['class']=='giant']['slope'].values * 1e4
    diff = s.mean() - g.mean()
    se = np.sqrt(s.var(ddof=1)/len(s) + g.var(ddof=1)/len(g))
    print(f"  {label:<24}: small={s.mean():+.3f}±{s.std(ddof=1)/np.sqrt(len(s)):.2f}, "
          f"giant={g.mean():+.3f}±{g.std(ddof=1)/np.sqrt(len(g)):.2f}, "
          f"Δ={diff:+.3f} ({diff/se:+.2f}σ)")

# =====================================================================
# OBJECTION 7: Power for the s-process detection
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 7: Power analysis specific to s-process detection")
print("="*78)
from scipy.stats import norm
s_sp = small_df['fam_s-process'].dropna().values
g_sp = giant_df['fam_s-process'].dropna().values
delta = s_sp.mean() - g_sp.mean()
pooled = np.sqrt(((len(s_sp)-1)*s_sp.var(ddof=1) + (len(g_sp)-1)*g_sp.var(ddof=1))
                 / (len(s_sp) + len(g_sp) - 2))
d = delta / pooled
print(f"Effect size (Cohen's d) = {abs(d):.3f}  [observed delta = {delta:+.4f}, pooled SD = {pooled:.4f}]")
for alpha, label in [(0.05,'p<0.05'), (0.0027,'3σ'), (5.7e-7,'5σ')]:
    z_a = norm.ppf(1-alpha/2); z_b = norm.ppf(0.80)
    N = int(np.ceil(2*(z_a+z_b)**2 / d**2))
    print(f"  N per group for 80% power, {label}: {N}")
print(f"  Current sample: small N={len(s_sp)}, giant N={len(g_sp)}")
# Achieved power at current N
def achieved_power(n, d, alpha):
    z_a = norm.ppf(1-alpha/2)
    z = abs(d) * np.sqrt(n/2) - z_a
    return norm.cdf(z)
n_eff = (2*len(s_sp)*len(g_sp))/(len(s_sp)+len(g_sp))
print(f"  Achieved power at p<0.05 (current N): {achieved_power(n_eff, d, 0.05)*100:.0f}%")
print(f"  Achieved power at 3σ (current N):     {achieved_power(n_eff, d, 0.0027)*100:.0f}%")

# =====================================================================
# OBJECTION 8: Bayes factor for two-population vs one-population model
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 8: Bayes factor for two-population vs one-population s-process")
print("="*78)
all_sp = np.concatenate([s_sp, g_sp])
def gauss_loglik(x, mu, sigma):
    return -0.5*np.sum(((x-mu)/sigma)**2 + np.log(2*np.pi*sigma**2))
# 1-population: Gaussian fit
mu1 = all_sp.mean(); s1 = all_sp.std(ddof=0)
ll1 = gauss_loglik(all_sp, mu1, s1)
# 2-population: separate Gaussians
mu_s = s_sp.mean(); ss_s = s_sp.std(ddof=0)
mu_g = g_sp.mean(); ss_g = g_sp.std(ddof=0)
ll2 = gauss_loglik(s_sp, mu_s, ss_s) + gauss_loglik(g_sp, mu_g, ss_g)
# BIC: penalty for params (1-pop has 2, 2-pop has 4)
n_total = len(all_sp)
bic1 = 2*np.log(n_total) - 2*ll1
bic2 = 4*np.log(n_total) - 2*ll2
delta_bic = bic1 - bic2
print(f"  1-population log-lik: {ll1:.2f}, BIC = {bic1:.2f}")
print(f"  2-population log-lik: {ll2:.2f}, BIC = {bic2:.2f}")
print(f"  ΔBIC (1-pop − 2-pop) = {delta_bic:+.2f}  →  Bayes factor ≈ exp(ΔBIC/2) = {np.exp(delta_bic/2):.2e}")
print(f"  Interpretation: 2-population model is favoured by {np.exp(delta_bic/2):.1e}:1")

# =====================================================================
# OBJECTION 9: Δlogg as additional regressor in family decomposition
# (already done above in Obj 2, summarized here)
# =====================================================================

# =====================================================================
# OBJECTION 10: Continuous mass trend - regression with slope CI
# =====================================================================
print("\n" + "="*78)
print("OBJECTION 10: Continuous mass trend - linear regression with CI")
print("="*78)

# Use [Y/Mg] differential age via TM2019 calibration
TM_SLOPE = -24.39  # Gyr/dex = 1/(-0.041 dex/Gyr), Tucci Maia 2019
mass_age_df = per_star.dropna(subset=['pl_bmassj','d_Y','d_Mg']).copy()
mass_age_df['dYMg'] = mass_age_df['d_Y'] - mass_age_df['d_Mg']
mass_age_df['delta_age'] = TM_SLOPE * mass_age_df['dYMg']
mass_age_df['log_mass'] = np.log10(mass_age_df['pl_bmassj'])

x = mass_age_df['log_mass'].values
y = mass_age_df['delta_age'].values
sl, ic, r, p, se = linregress(x, y)
# Bootstrap slope
slope_boots = []
n = len(x)
for _ in range(N_BOOT):
    idx = rng.integers(0, n, n)
    s_b, *_ = linregress(x[idx], y[idx])
    slope_boots.append(s_b)
sb = np.array(slope_boots)
lo, hi = np.percentile(sb, [2.5, 97.5])
print(f"Linear regression: Δage = β · log10(M_p sin i / M_J) + α")
print(f"  β  = {sl:+.3f} ± {se:.3f} Gyr/dex")
print(f"  95% bootstrap CI on β: [{lo:+.3f}, {hi:+.3f}] Gyr/dex")
print(f"  Pearson r = {r:+.3f}, p = {p:.4f}")
print(f"  N = {n}")
obj10 = dict(slope=sl, slope_se=se, slope_lo=lo, slope_hi=hi, r=r, p=p, n=n)

# Save state
with open(f'{OUT}/_state_v3.pkl','wb') as f:
    pickle.dump({
        'per_star': per_star, 'small_df': small_df, 'giant_df': giant_df,
        'elements': elements, 'T_arr': T_arr, 'tc_map': tc_map,
        'FAMILIES': FAMILIES,
        'per_elem_df': per_elem_df,
        'robust_df': robust_df,
        'clock_results': clock_results,
        'slopes_full': slopes_full, 'slopes_pure': slopes_pure,
        'obj1': obj1, 'obj2': obj2, 'obj10': obj10,
        'mass_age_df': mass_age_df,
        's_sp': s_sp, 'g_sp': g_sp,
        'd': d, 'delta_bic': delta_bic,
    }, f)

print("\n✓ All robustness analyses done. State pickled.")
