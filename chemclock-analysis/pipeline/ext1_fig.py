"""Figure: isochrone vs [Y/Mg] differential age cross-validation."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from matplotlib.gridspec import GridSpec
from scipy.stats import linregress

plt.rcParams.update({
    'font.family':'serif','font.size':10,'axes.labelsize':11,'axes.titlesize':11,
    'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'axes.linewidth':0.8,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'figure.dpi':130,'savefig.dpi':300,'savefig.bbox':'tight',
})

OUT = OUT
with open(f'{OUT}/_state_ext1.pkl','rb') as f:
    S = pickle.load(f)

merged = S['merged']
both   = S['both']

small = merged[merged['planet_class']=='small']
giant = merged[merged['planet_class']=='giant']

fig = plt.figure(figsize=(13, 8.5))
gs = GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.35,
              top=0.94, bottom=0.08, left=0.07, right=0.97)

# ============ (a) Per-star scatter: iso vs Y/Mg ============
ax = fig.add_subplot(gs[0, 0:2])
colors = merged['planet_class'].map({'small':'steelblue','giant':'firebrick'})

x = merged['d_age_iso'].values
y = merged['d_age_ymg'].values
xe = merged['d_age_iso_err'].values
ye = merged['d_age_ymg_err'].values

# Color by planet class
for cls, color, marker in [('small','steelblue','o'), ('giant','firebrick','s')]:
    sub = merged[merged['planet_class']==cls]
    ax.errorbar(sub['d_age_iso'], sub['d_age_ymg'],
                xerr=sub['d_age_iso_err'], yerr=sub['d_age_ymg_err'],
                fmt=marker, color=color, alpha=0.6, ms=5,
                ecolor='lightgrey', elinewidth=0.7, capsize=0,
                label=f'{cls.capitalize()}-planet hosts (N={len(sub)})')

# 1:1 line
lim = max(np.nanmax(np.abs(x)), np.nanmax(np.abs(y))) * 0.9
ax.plot([-lim, lim], [-lim, lim], 'k:', lw=1, alpha=0.5, label='1:1 line')

# Linear fit
valid = np.isfinite(x) & np.isfinite(y)
slope, ic, r_, p_, se = linregress(x[valid], y[valid])
xf = np.linspace(-lim, lim, 50)
ax.plot(xf, slope*xf + ic, 'C2-', lw=2, label=f'Linear fit: slope = {slope:+.2f} ± {se:.2f}')

# Bootstrap CI
rng = np.random.default_rng(7)
n = valid.sum(); slopes = []
xx, yy = x[valid], y[valid]
for _ in range(2000):
    idx = rng.integers(0, n, n)
    s_, *_ = linregress(xx[idx], yy[idx])
    slopes.append(s_)
slopes = np.array(slopes)
preds = np.array([s*xf + (yy.mean() - s*xx.mean()) for s in slopes])
ax.fill_between(xf, np.percentile(preds, 16, axis=0), np.percentile(preds, 84, axis=0),
                color='C2', alpha=0.15, label='68% bootstrap CI')

ax.axhline(0, color='grey', lw=0.5)
ax.axvline(0, color='grey', lw=0.5)
ax.set_xlabel(r'Isochrone $\Delta_{\rm age}$  (SwP $-$ neighbours, Gyr)')
ax.set_ylabel(r'[Y/Mg] chemical-clock $\Delta_{\rm age}$  (Gyr)')
ax.set_title(f'(a) Per-star agreement of two independent age methods   ('
             f'Pearson $r={S["r_pearson"]:+.2f}$, $p={S["p_pearson"]:.3f}$)')
ax.legend(loc='upper left', frameon=True, fontsize=8.5)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

# ============ (b) Group means: iso vs Y/Mg vs combined ============
ax = fig.add_subplot(gs[0, 2])
methods = ['Isochrone\n($\\Delta_{\\rm age}^{\\rm iso}$)', 
           '[Y/Mg]\n($\\Delta_{\\rm age}^{\\rm YMg}$)',
           'Combined\n(var-weighted)']
small_means = [
    small['d_age_iso'].mean(),
    small['d_age_ymg'].mean(),
    both[both['planet_class']=='small']['d_age_comb'].mean(),
]
small_ses = [
    small['d_age_iso'].std(ddof=1)/np.sqrt(small['d_age_iso'].notna().sum()),
    small['d_age_ymg'].std(ddof=1)/np.sqrt(small['d_age_ymg'].notna().sum()),
    both[both['planet_class']=='small']['d_age_comb'].std(ddof=1)/
        np.sqrt(both[both['planet_class']=='small']['d_age_comb'].notna().sum()),
]
giant_means = [
    giant['d_age_iso'].mean(),
    giant['d_age_ymg'].mean(),
    both[both['planet_class']=='giant']['d_age_comb'].mean(),
]
giant_ses = [
    giant['d_age_iso'].std(ddof=1)/np.sqrt(giant['d_age_iso'].notna().sum()),
    giant['d_age_ymg'].std(ddof=1)/np.sqrt(giant['d_age_ymg'].notna().sum()),
    both[both['planet_class']=='giant']['d_age_comb'].std(ddof=1)/
        np.sqrt(both[both['planet_class']=='giant']['d_age_comb'].notna().sum()),
]

x_pos = np.arange(len(methods))
w = 0.35
ax.bar(x_pos-w/2-0.02, small_means, w, yerr=small_ses, capsize=4,
       color='steelblue', alpha=0.85, edgecolor='k', lw=0.6,
       label='Small-host', error_kw=dict(lw=1.2, ecolor='black'))
ax.bar(x_pos+w/2+0.02, giant_means, w, yerr=giant_ses, capsize=4,
       color='firebrick', alpha=0.85, edgecolor='k', lw=0.6,
       label='Giant-host', error_kw=dict(lw=1.2, ecolor='black'))

# Annotate sigmas
sigmas = [S['iso_sig'], S['ymg_sig'], S['comb_sig']]
ps     = [S['iso_p'],   S['ymg_p'],   S['comb_p']]
for i, (sg, p) in enumerate(zip(sigmas, ps)):
    y_max = max(small_means[i]+small_ses[i], giant_means[i]+giant_ses[i])
    ax.annotate(f'S-G:\n{sg:+.1f}$\\sigma$\np={p:.3f}',
                xy=(i, y_max + 0.15),
                ha='center', va='bottom', fontsize=8.5,
                bbox=dict(facecolor='lightyellow' if abs(sg)>2 else 'white',
                         edgecolor='grey', alpha=0.8))

ax.axhline(0, color='grey', lw=0.6)
ax.set_xticks(x_pos); ax.set_xticklabels(methods, fontsize=9)
ax.set_ylabel(r'$\Delta_{\rm age}$  (Gyr)')
ax.set_title('(b) Small-vs-giant in each method')
ax.legend(loc='lower left', frameon=True, fontsize=8.5)
ax.set_ylim(-0.8, 1.5)

# ============ (c) Δage_iso histograms ============
ax = fig.add_subplot(gs[1, 0])
bins = np.linspace(-5, 5, 22)
s_iso = small['d_age_iso'].dropna().values
g_iso = giant['d_age_iso'].dropna().values
ax.hist(s_iso, bins=bins, color='steelblue', alpha=0.6, edgecolor='k', lw=0.4,
        label=f'Small (N={len(s_iso)})')
ax.hist(g_iso, bins=bins, color='firebrick', alpha=0.5, edgecolor='k', lw=0.4,
        label=f'Giant (N={len(g_iso)})')
ax.axvline(0, color='grey', lw=0.5, ls='--')
ax.axvline(s_iso.mean(), color='steelblue', lw=1.4, ls=':')
ax.axvline(g_iso.mean(), color='firebrick', lw=1.4, ls=':')
ax.set_xlabel(r'Isochrone $\Delta_{\rm age}$  (Gyr)')
ax.set_ylabel('Number of stars')
ax.set_title('(c) Isochrone $\\Delta_{\\rm age}$ distribution')
ax.legend(loc='upper right', frameon=True, fontsize=8.5)

# ============ (d) Δage_YMg histograms (same axes) ============
ax = fig.add_subplot(gs[1, 1])
s_y = small['d_age_ymg'].dropna().values
g_y = giant['d_age_ymg'].dropna().values
ax.hist(s_y, bins=bins, color='steelblue', alpha=0.6, edgecolor='k', lw=0.4,
        label=f'Small (N={len(s_y)})')
ax.hist(g_y, bins=bins, color='firebrick', alpha=0.5, edgecolor='k', lw=0.4,
        label=f'Giant (N={len(g_y)})')
ax.axvline(0, color='grey', lw=0.5, ls='--')
ax.axvline(s_y.mean(), color='steelblue', lw=1.4, ls=':')
ax.axvline(g_y.mean(), color='firebrick', lw=1.4, ls=':')
ax.set_xlabel(r'[Y/Mg] $\Delta_{\rm age}$  (Gyr)')
ax.set_ylabel('Number of stars')
ax.set_title('(d) [Y/Mg] $\\Delta_{\\rm age}$ distribution')
ax.legend(loc='upper right', frameon=True, fontsize=8.5)

# ============ (e) Combined Δage histograms ============
ax = fig.add_subplot(gs[1, 2])
s_c = both[both['planet_class']=='small']['d_age_comb'].dropna().values
g_c = both[both['planet_class']=='giant']['d_age_comb'].dropna().values
bins_c = np.linspace(-3, 3, 22)
ax.hist(s_c, bins=bins_c, color='steelblue', alpha=0.6, edgecolor='k', lw=0.4,
        label=f'Small (N={len(s_c)})')
ax.hist(g_c, bins=bins_c, color='firebrick', alpha=0.5, edgecolor='k', lw=0.4,
        label=f'Giant (N={len(g_c)})')
ax.axvline(0, color='grey', lw=0.5, ls='--')
ax.axvline(s_c.mean(), color='steelblue', lw=1.4, ls=':')
ax.axvline(g_c.mean(), color='firebrick', lw=1.4, ls=':')
ax.set_xlabel(r'Combined $\Delta_{\rm age}$  (Gyr)')
ax.set_ylabel('Number of stars')
ax.set_title(f'(e) Variance-weighted combined $\\Delta_{{\\rm age}}$  '
             f'(S$-$G: {S["comb_sig"]:+.1f}$\\sigma$)')
ax.legend(loc='upper right', frameon=True, fontsize=8.5)

plt.savefig(f'{OUT}/figures/fig_iso_vs_ymg.pdf')
plt.savefig(f'{OUT}/figures/fig_iso_vs_ymg.png')
plt.close()

print("Figure saved.")
