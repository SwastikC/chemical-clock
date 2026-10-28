"""All figures for the strengthened paper."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.stats import linregress

plt.rcParams.update({
    'font.family':'serif','font.size':10,'axes.labelsize':11,'axes.titlesize':11,
    'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'axes.linewidth':0.8,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'figure.dpi':130,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT = OUT
with open(f'{OUT}/_state_v3.pkl','rb') as f:
    S = pickle.load(f)

per_star = S['per_star']
small_df = S['small_df']
giant_df = S['giant_df']
elements = S['elements']
T_arr    = S['T_arr']
tc_map   = S['tc_map']
FAMILIES = S['FAMILIES']
per_elem_df = S['per_elem_df']
robust_df   = S['robust_df']
clock_results = S['clock_results']
slopes_full = S['slopes_full']
slopes_pure = S['slopes_pure']
obj1     = S['obj1']
obj2     = S['obj2']
obj10    = S['obj10']
mass_age_df = S['mass_age_df']
s_sp = S['s_sp']; g_sp = S['g_sp']
d        = S['d']
delta_bic = S['delta_bic']

FAMILY_COLOR = {'Volatile':'C0','Iron-peak ref':'C2','Alpha + Al ref':'C1','s-process':'C3'}
FAMILY_FOR_ELEM = {e:f for f,es in FAMILIES.items() for e in es}

# =====================================================================
# FIGURE 1 — Per-element offsets (smoking gun)
# =====================================================================
fig, axes = plt.subplots(2, 1, figsize=(11, 6.8), sharex=True)
for ax, grp_label, df_grp, color, marker in [
    (axes[0], 'Small-planet hosts', small_df, 'steelblue', 'o'),
    (axes[1], 'Giant-planet hosts', giant_df, 'firebrick', 's')]:
    means = np.array([df_grp[f'd_{e}'].mean() for e in elements])
    ses   = np.array([df_grp[f'd_{e}'].std(ddof=1)/np.sqrt(df_grp[f'd_{e}'].notna().sum())
                       for e in elements])
    for i, (el, t) in enumerate(sorted(zip(elements, T_arr), key=lambda p: p[1])):
        fam = FAMILY_FOR_ELEM[el]
        ax.axvspan(t-22, t+22, color=FAMILY_COLOR[fam], alpha=0.07, lw=0)

    ax.errorbar(T_arr, means, yerr=ses, fmt=marker, ms=7, color=color,
                ecolor='grey', capsize=2.5, lw=1.2, mec='k', mew=0.4, zorder=3)
    for i, el in enumerate(elements):
        weight = 'bold' if FAMILY_FOR_ELEM[el]=='s-process' else 'normal'
        ax.annotate(el, (T_arr[i], means[i]),
                    xytext=(0, 11 if means[i]>=0 else -14),
                    textcoords='offset points',
                    fontsize=8, ha='center', color='dimgray', weight=weight)
    ax.axhline(0, color='grey', lw=0.6)
    ax.axvline(1200, color='goldenrod', ls=':', lw=1)
    ax.set_ylabel(r'$\overline{\Delta[X/\mathrm{Fe}]}$  (SwP $-$ neighbours, dex)')
    ax.text(0.02, 0.94, f'{grp_label}  (N = {len(df_grp)})',
            transform=ax.transAxes, fontsize=11, weight='bold',
            bbox=dict(facecolor='white', edgecolor='grey', alpha=0.85))
    ax.set_ylim(-0.07, 0.06)

legend_elems = [Patch(facecolor=FAMILY_COLOR[f], alpha=0.4, label=f) for f in FAMILIES]
axes[1].legend(handles=legend_elems, loc='lower right', frameon=True, ncol=4)   # empty part of the plot
axes[1].set_xlabel(r'Condensation temperature  $T_{\rm C}$  (K)')
plt.tight_layout()
plt.savefig(f'{OUT}/figures/fig1_per_element.pdf')
plt.savefig(f'{OUT}/figures/fig1_per_element.png')
plt.close()

# =====================================================================
# FIGURE 2 — Family decomposition (bars + sigma)
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
ax = axes[0]
fam_list = list(FAMILIES.keys())
x_pos = np.arange(len(fam_list))
w = 0.35
rng = np.random.default_rng(7)
for i, fam in enumerate(fam_list):
    s_vals = small_df[f'fam_{fam}'].dropna().values
    g_vals = giant_df[f'fam_{fam}'].dropna().values
    s_means = np.array([rng.choice(s_vals, len(s_vals), replace=True).mean() for _ in range(5000)])
    g_means = np.array([rng.choice(g_vals, len(g_vals), replace=True).mean() for _ in range(5000)])
    s_lo, s_hi = np.percentile(s_means, [2.5,97.5])
    g_lo, g_hi = np.percentile(g_means, [2.5,97.5])
    ax.bar(i-w/2-0.02, s_vals.mean(), w, color='steelblue', alpha=0.85,
           yerr=[[s_vals.mean()-s_lo],[s_hi-s_vals.mean()]], capsize=4,
           edgecolor='k', lw=0.6,
           label='Small-planet hosts' if i==0 else None)
    ax.bar(i+w/2+0.02, g_vals.mean(), w, color='firebrick', alpha=0.85,
           yerr=[[g_vals.mean()-g_lo],[g_hi-g_vals.mean()]], capsize=4,
           edgecolor='k', lw=0.6,
           label='Giant-planet hosts' if i==0 else None)
ax.axhline(0, color='grey', lw=0.6)
ax.set_xticks(x_pos); ax.set_xticklabels(fam_list, rotation=12, ha='right')
ax.set_ylabel(r'$\langle\overline{\Delta[X/\mathrm{Fe}]}\rangle_{\rm family}$  (dex)')
ax.set_title('(a) Family means (bootstrap mean $\\pm$ 95% CI)')
ax.legend(loc='lower left', frameon=True)

ax = axes[1]
sigmas, ps = [], []
for fam in fam_list:
    s = small_df[f'fam_{fam}'].dropna().values
    g = giant_df[f'fam_{fam}'].dropna().values
    diff = s.mean() - g.mean()
    se = np.sqrt(s.var(ddof=1)/len(s) + g.var(ddof=1)/len(g))
    sigmas.append(diff/se)
    # quick perm
    pool = np.concatenate([s,g]); na = len(s)
    ds = np.empty(2000)
    for j in range(2000):
        rng.shuffle(pool); ds[j] = pool[:na].mean()-pool[na:].mean()
    ps.append((np.abs(ds) >= np.abs(diff)).mean())
colors = [FAMILY_COLOR[f] for f in fam_list]
bars = ax.barh(x_pos, sigmas, color=colors, edgecolor='k', alpha=0.85, lw=0.6)
for i, (s, p) in enumerate(zip(sigmas, ps)):
    label = f'$\\sigma$={s:+.1f}\np={p:.4f}'.replace('p=0.0000', f'p<{1/len(ds):.4f}')
    ax.text(s + (0.4 if s>=0 else -0.4), i, label,
            va='center', ha='left' if s>=0 else 'right', fontsize=9)
ax.axvline(0, color='grey', lw=0.6)
ax.axvline(-3, color='C3', ls=':', lw=0.8, alpha=0.5)
ax.axvline(+3, color='C3', ls=':', lw=0.8, alpha=0.5)
ax.set_yticks(x_pos); ax.set_yticklabels(fam_list)
ax.set_xlabel('Small $-$ Giant offset  ($\\sigma$, two-sample permutation)')
ax.set_title('(b) Where the small/giant difference lives')
ax.set_xlim(-7, 4)

plt.tight_layout()
plt.savefig(f'{OUT}/figures/fig2_family_decomposition.pdf')
plt.savefig(f'{OUT}/figures/fig2_family_decomposition.png')
plt.close()

# =====================================================================
# FIGURE 3 — Five chemical clocks (UPDATED)
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

# Panel A: 5 clocks bar chart
ax = axes[0]
clocks = list(clock_results.keys())
y_pos = np.arange(len(clocks))
small_means = [clock_results[c]['small_mean'] for c in clocks]
small_ses   = [clock_results[c]['small_se']   for c in clocks]
giant_means = [clock_results[c]['giant_mean'] for c in clocks]
giant_ses   = [clock_results[c]['giant_se']   for c in clocks]

w = 0.4
ax.barh(y_pos-w/2, small_means, w, xerr=small_ses, capsize=3,
        color='steelblue', alpha=0.85, edgecolor='k', lw=0.6,
        label='Small-planet hosts', error_kw=dict(lw=1.2, ecolor='black'))
ax.barh(y_pos+w/2, giant_means, w, xerr=giant_ses, capsize=3,
        color='firebrick', alpha=0.85, edgecolor='k', lw=0.6,
        label='Giant-planet hosts', error_kw=dict(lw=1.2, ecolor='black'))
ax.axvline(0, color='grey', lw=0.6, ls='--')
ax.set_yticks(y_pos); ax.set_yticklabels(clocks)
ax.set_xlabel(r'$\Delta_{\rm age}$  (Gyr)')
ax.set_title('(a) Five chemical clocks all agree')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
# Annotate sigma to the right of the longest bar; the axis follows the data
right = max(max(m + e for m, e in zip(small_means, small_ses)), max(m + e for m, e in zip(giant_means, giant_ses)))
left = min(min(m - e for m, e in zip(small_means, small_ses)), min(m - e for m, e in zip(giant_means, giant_ses)))
for i, c in enumerate(clocks):
    sig = clock_results[c]['sigma']
    ax.text(right + 0.1, i, f' $\\sigma$ = {sig:+.1f}', va='center', fontsize=8, color='dimgray')
ax.set_xlim(min(-0.7, left - 0.2), right + 0.8)

# Panel B: continuous mass trend with regression line
ax = axes[1]
both = mass_age_df.dropna(subset=['log_mass','delta_age']).copy()
ax.errorbar(both['pl_bmassj'], both['delta_age'],
            yerr=None, fmt='o', ms=4, color='dimgray',
            alpha=0.55)
# Plot regression
x_grid = np.linspace(both['log_mass'].min(), both['log_mass'].max(), 50)
y_fit = obj10['slope']*x_grid + (both['delta_age'].mean() - obj10['slope']*both['log_mass'].mean())
ax.plot(10**x_grid, y_fit, 'r-', lw=2, label=f'Linear fit')
# Bootstrap CI band
rng2 = np.random.default_rng(7)
slopes_boot, intercepts_boot = [], []
xa = both['log_mass'].values; ya = both['delta_age'].values
n = len(xa)
for _ in range(2000):
    idx = rng2.integers(0, n, n)
    sl, ic, *_ = linregress(xa[idx], ya[idx])
    slopes_boot.append(sl); intercepts_boot.append(ic)
slopes_boot = np.array(slopes_boot); intercepts_boot = np.array(intercepts_boot)
y_predictions = np.array([s*x_grid+i for s,i in zip(slopes_boot, intercepts_boot)])
y_lo = np.percentile(y_predictions, 16, axis=0)
y_hi = np.percentile(y_predictions, 84, axis=0)
ax.fill_between(10**x_grid, y_lo, y_hi, color='red', alpha=0.18,
                label=r'68% bootstrap CI')

ax.axhline(0, color='grey', lw=0.6)
ax.axvline(0.3, color='goldenrod', ls='--', lw=1, label='small/giant cut')
ax.set_xscale('log')
ax.set_xlabel(r'Planet mass  $M_p \sin i$  ($M_{\rm J}$)')
ax.set_ylabel(r'$\Delta_{\rm age}$  (Gyr)')
slope_str = (f"$\\beta = {obj10['slope']:+.2f}^{{{obj10['slope_hi']-obj10['slope']:+.2f}}}"
             f"_{{{obj10['slope_lo']-obj10['slope']:+.2f}}}$ Gyr/dex\n"
             f"$r = {obj10['r']:+.2f}$, $p = {obj10['p']:.3f}$")
ax.text(0.04, 0.95, slope_str, transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.85))
ax.set_title('(b) Continuous trend: $\\Delta_{\\rm age}$ vs planet mass')
ax.legend(loc='lower right', frameon=True, fontsize=8.5)
ax.set_ylim(-3, 5)

plt.tight_layout()
plt.savefig(f'{OUT}/figures/fig3_clocks_and_trend.pdf')
plt.savefig(f'{OUT}/figures/fig3_clocks_and_trend.png')
plt.close()

# =====================================================================
# FIGURE 4 — Naive vs Purified (kept the same as v2)
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for ax, col, title in [
    (axes[0], 'fam_s-process',
     '(a) Refractory$-$volatile slope (per-star Tc fit)\nINCLUDING s-process refractories'),
    (axes[1], None,
     '(b) Same fit, EXCLUDING s-process refractories')]:
    pass

# Use the slopes from per-star Tc fits (Bedell-style)
ax = axes[0]
s_vals = slopes_full[slopes_full['class']=='small']['slope'].values * 1e4
g_vals = slopes_full[slopes_full['class']=='giant']['slope'].values * 1e4
parts = ax.violinplot([s_vals, g_vals], positions=[1,2],
                      widths=0.7, showmeans=False, showmedians=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(['steelblue','firebrick'][i])
    pc.set_alpha(0.5); pc.set_edgecolor('black'); pc.set_linewidth(0.6)
for partname in ('cbars','cmins','cmaxes'):
    parts[partname].set_color('grey'); parts[partname].set_linewidth(0.6)
rng3 = np.random.default_rng(7)
for i, vals in enumerate([s_vals, g_vals]):
    means = np.array([rng3.choice(vals, len(vals), replace=True).mean() for _ in range(5000)])
    ax.errorbar(i+1, vals.mean(),
                yerr=[[vals.mean()-np.percentile(means,2.5)],
                      [np.percentile(means,97.5)-vals.mean()]],
                fmt='o', color='black', ms=8, capsize=4, lw=1.4, mfc='white', zorder=5)
diff = s_vals.mean() - g_vals.mean()
se_pool = np.sqrt(s_vals.var(ddof=1)/len(s_vals) + g_vals.var(ddof=1)/len(g_vals))
sigma = diff/se_pool
ax.axhline(0, color='grey', lw=0.6)
ax.set_xticks([1,2])
ax.set_xticklabels([f'Small\n(n={len(s_vals)})', f'Giant\n(n={len(g_vals)})'])
ax.set_ylabel(r'Tc-slope  ($10^{-4}$ dex/K)')
ax.set_title('(a) Tc-slope INCLUDING s-process elements\n(standard literature method)', fontsize=10.5)
ax.text(0.5, 0.04,
        f'Small$-$Giant difference: {diff:+.3f}\n$\\sigma$ = {sigma:+.2f}',
        transform=ax.transAxes, ha='center', fontsize=10,
        bbox=dict(facecolor='lightyellow', edgecolor='C3', alpha=0.9))

ax = axes[1]
s_vals = slopes_pure[slopes_pure['class']=='small']['slope'].values * 1e4
g_vals = slopes_pure[slopes_pure['class']=='giant']['slope'].values * 1e4
parts = ax.violinplot([s_vals, g_vals], positions=[1,2],
                      widths=0.7, showmeans=False, showmedians=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(['steelblue','firebrick'][i])
    pc.set_alpha(0.5); pc.set_edgecolor('black'); pc.set_linewidth(0.6)
for partname in ('cbars','cmins','cmaxes'):
    parts[partname].set_color('grey'); parts[partname].set_linewidth(0.6)
for i, vals in enumerate([s_vals, g_vals]):
    means = np.array([rng3.choice(vals, len(vals), replace=True).mean() for _ in range(5000)])
    ax.errorbar(i+1, vals.mean(),
                yerr=[[vals.mean()-np.percentile(means,2.5)],
                      [np.percentile(means,97.5)-vals.mean()]],
                fmt='o', color='black', ms=8, capsize=4, lw=1.4, mfc='white', zorder=5)
diff = s_vals.mean() - g_vals.mean()
se_pool = np.sqrt(s_vals.var(ddof=1)/len(s_vals) + g_vals.var(ddof=1)/len(g_vals))
sigma = diff/se_pool
ax.axhline(0, color='grey', lw=0.6)
ax.set_xticks([1,2])
ax.set_xticklabels([f'Small\n(n={len(s_vals)})', f'Giant\n(n={len(g_vals)})'])
ax.set_ylabel(r'Tc-slope  ($10^{-4}$ dex/K)')
ax.set_title('(b) Tc-slope EXCLUDING s-process elements\n(planet-formation-relevant signal)', fontsize=10.5)
ax.text(0.5, 0.04,
        f'Small$-$Giant difference: {diff:+.3f}\n$\\sigma$ = {sigma:+.2f}',
        transform=ax.transAxes, ha='center', fontsize=10,
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9))

fig.text(0.5, 0.02,
         "Removing the s-process elements collapses the apparent Tc-slope difference between small and giant planet hosts.",
         ha='center', fontsize=10, style='italic')
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(f'{OUT}/figures/fig4_tcslope_full_vs_pure.pdf')
plt.savefig(f'{OUT}/figures/fig4_tcslope_full_vs_pure.png')
plt.close()

# =====================================================================
# FIGURE 5 — NEW Robustness panel
# =====================================================================
fig = plt.figure(figsize=(13, 8))
gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.32)

# Small-giant sigma of the s-process family, before any regression
sigma_raw = (s_sp.mean() - g_sp.mean()) / np.sqrt(s_sp.var(ddof=1)/len(s_sp) + g_sp.var(ddof=1)/len(g_sp))

# (a) Stellar-parameter regression: Δ[s-process] vs Δ[Fe/H]
ax = fig.add_subplot(gs[0, 0])
both = per_star.dropna(subset=['fam_s-process','dFeH','planet_class']).copy()
colors = both['planet_class'].map({'small':'steelblue','giant':'firebrick'})
ax.scatter(both['dFeH'], both['fam_s-process'], c=colors, alpha=0.6, s=22,
           edgecolor='k', lw=0.4)
sl, ic, r, p, _ = linregress(both['dFeH'], both['fam_s-process'])
xfit = np.array([both['dFeH'].min(), both['dFeH'].max()])
ax.plot(xfit, sl*xfit + ic, 'k--', lw=1.2,
        label=f'$r={r:+.2f}$, $p={p:.2f}$')
ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)
ax.set_xlabel(r'$\Delta$[Fe/H]  (dex)')
ax.set_ylabel(r'$\Delta$[$s$-process]  (dex)')
ax.set_title('(a) No correlation with $\\Delta$[Fe/H]')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.legend(loc='upper left', fontsize=8)

# (b) Leave-one-out s-process stability
ax = fig.add_subplot(gs[0, 1])
sproc = FAMILIES['s-process']
labels, vals = ['all 4'], [sigma_raw]
# Recompute leave-one-out on the fly
for skip in sproc:
    use = [e for e in sproc if e != skip]
    s_vals, g_vals = [], []
    for _, row in per_star.iterrows():
        if pd.isna(row['planet_class']): continue
        ys = [row[f'd_{el}'] for el in use if not pd.isna(row[f'd_{el}'])]
        es = [row[f'e_{el}'] for el in use if not pd.isna(row[f'e_{el}']) and row[f'e_{el}']>0]
        if len(ys) == 0 or len(es) != len(ys): continue
        ys = np.array(ys); es = np.array(es)
        w = 1/es**2
        v = (ys*w).sum()/w.sum()
        if row['planet_class']=='small': s_vals.append(v)
        elif row['planet_class']=='giant': g_vals.append(v)
    s_vals = np.array(s_vals); g_vals = np.array(g_vals)
    diff = s_vals.mean() - g_vals.mean()
    se = np.sqrt(s_vals.var(ddof=1)/len(s_vals) + g_vals.var(ddof=1)/len(g_vals))
    labels.append(f'drop {skip}')
    vals.append(diff/se)
y_pos = np.arange(len(labels))
ax.barh(y_pos, vals, color='C3', alpha=0.7, edgecolor='k', lw=0.5)
for i, v in enumerate(vals):
    ax.text(v + 0.15, i, f'{v:+.1f}σ', va='center', ha='left',
            color='white' if abs(v) > 1 else 'black', fontsize=9)
ax.axvline(-3, color='gray', ls=':', lw=0.8)
ax.axvline(0, color='grey', lw=0.6)
ax.set_yticks(y_pos); ax.set_yticklabels(labels)
ax.set_xlabel('Small$-$Giant $\\sigma$')
ax.set_title('(b) Leave-one-out: $s$-process stability')
ax.set_xlim(-6, 0.5)

# (c) T_C cut robustness
ax = fig.add_subplot(gs[0, 2])
ax.plot(robust_df['tc_split'], robust_df['sigma'], 'o-',
        color='C3', markersize=7, lw=1.4, mec='k', mew=0.4)
ax.axhline(0, color='grey', lw=0.6)
ax.axvline(1200, color='goldenrod', ls='--', lw=1, label='conventional cut')
ax.fill_between([1100, 1300], -2, 0, alpha=0.15, color='C3', label='robust regime')
ax.set_xlabel(r'$T_{\rm C}$ split (K)')
ax.set_ylabel('Small$-$Giant $\\sigma$ in $\\Delta_{\\rm ref-vol}$')
ax.set_title('(c) Robustness to $T_{\\rm C}$ split')
ax.legend(fontsize=8, loc='upper right')

# (d) Power curves for the s-process detection
ax = fig.add_subplot(gs[1, 0])
from scipy.stats import norm
Ns = np.logspace(np.log10(10), np.log10(500), 80)
def power_at_N(N, d, alpha):
    z_a = norm.ppf(1-alpha/2)
    z = abs(d) * np.sqrt(N/2) - z_a
    return norm.cdf(z)
for alpha, c, lbl in [(0.05, 'C0', '$p<0.05$'),
                      (0.0027, 'C3', '$3\\sigma$'),
                      (5.7e-7, 'darkred', '$5\\sigma$')]:
    ax.plot(Ns, [power_at_N(n, abs(d), alpha) for n in Ns],
            color=c, lw=1.7, label=lbl)
ax.axhline(0.80, color='grey', ls=':', lw=1, label='80% power')
N_eff = (2*len(s_sp)*len(g_sp))/(len(s_sp)+len(g_sp))
ax.axvline(N_eff, color='black', ls='--', lw=1)
ax.text(N_eff*1.05, 0.05, f'current\nN$_{{eff}}$={N_eff:.0f}', fontsize=8)
ax.set_xscale('log')
ax.set_xlabel('Effective sample size per group')
ax.set_ylabel('Statistical power')
ax.set_title(f'(d) Power for $s$-process effect\n(Cohen\'s $d$ = {abs(d):.2f})')
ax.set_ylim(0, 1.05)
ax.legend(loc='lower right', fontsize=8)

# (e) Bayes factor: small / giant s-process histograms with two-Gaussian fit
ax = fig.add_subplot(gs[1, 1])
bins = np.linspace(-0.20, 0.18, 28)
ax.hist(s_sp, bins=bins, color='steelblue', alpha=0.6, edgecolor='k', lw=0.4,
        density=True, label=f'Small (N={len(s_sp)})')
ax.hist(g_sp, bins=bins, color='firebrick', alpha=0.5, edgecolor='k', lw=0.4,
        density=True, label=f'Giant (N={len(g_sp)})')
xx = np.linspace(-0.20, 0.18, 200)
def gauss(x, mu, s): return np.exp(-0.5*((x-mu)/s)**2)/(s*np.sqrt(2*np.pi))
ax.plot(xx, gauss(xx, s_sp.mean(), s_sp.std(ddof=0)), 'b-', lw=1.5)
ax.plot(xx, gauss(xx, g_sp.mean(), g_sp.std(ddof=0)), 'r-', lw=1.5)
ax.axvline(0, color='grey', lw=0.6, ls='--')
ax.axvline(s_sp.mean(), color='steelblue', lw=1.2, ls=':')
ax.axvline(g_sp.mean(), color='firebrick', lw=1.2, ls=':')
bf = np.exp(delta_bic/2)
ax.text(0.02, 0.95, f'$\\Delta$BIC = {delta_bic:+.1f}\nBayes factor = {bf:.0f}:1\n(2-pop favoured)',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.9))
ax.set_xlabel(r'$\Delta$[$s$-process]  (dex)')
ax.set_ylabel('Density')
ax.set_title('(e) Bayes-factor support for two populations')
ax.legend(loc='upper right', fontsize=8)

# (f) After ΔT/Δlogg/Δ[Fe/H] regression: signal stable
ax = fig.add_subplot(gs[1, 2])
panels = ['Raw', '$-\\Delta$[Fe/H]', '$-\\Delta T,\\Delta g, \\Delta$[Fe/H]']
sigmas_panel = [sigma_raw, obj1['sigma_resid'], obj2['sigma_clean']]
y_pos = np.arange(len(panels))
ax.barh(y_pos, sigmas_panel, color='C3', alpha=0.7, edgecolor='k', lw=0.5)
for i, v in enumerate(sigmas_panel):
    ax.text(v + 0.15, i, f'{v:+.1f}σ', va='center', ha='left',
            color='white', fontsize=9)
ax.axvline(0, color='grey', lw=0.6)
ax.axvline(-3, color='gray', ls=':', lw=0.8)
ax.set_yticks(y_pos); ax.set_yticklabels(panels)
ax.set_xlabel('Small$-$Giant $\\sigma$')
ax.set_title('(f) Stable under stellar-param regression')
ax.set_xlim(-6, 0.5)

plt.savefig(f'{OUT}/figures/fig5_robustness.pdf')
plt.savefig(f'{OUT}/figures/fig5_robustness.png')
plt.close()

print("All 5 figures generated.")
