# chemclock-analysis

| | |
|---|---|
| [`pipeline/`](pipeline/) | The analysis: element-family offsets, chemical-clock ages, robustness tests, isochrone cross-check, Figures 1–6. |
| [`checks/control_bootstrap.py`](checks/control_bootstrap.py) | Uncertainty of the main result including the noise of the control stars. |

```
python pipeline/run_pipeline.py      # results and figures (v2 neighbours)
python checks/control_bootstrap.py   # uncertainty including control-star noise (seconds)
```

Outputs go to `pipeline/output/` and `checks/output/` (not committed). The expected numbers are in
the [top-level README](../README.md#step-4-the-analysis-and-figures).

## pipeline/

Reads the host list, the v2 neighbour table and its host-minus-neighbour table
(`final_test_v2.csv`), the HARPS-GTO parameters and isochrone ages, and the archive snapshot; all
paths are in `pipeline/paths.py`. `run_pipeline.py` runs the five scripts in order:

| Script | What it does | Writes |
|---|---|---|
| `build_stellar_offsets.py` | host minus mean-of-neighbours ΔTeff, Δlog g, Δ[Fe/H] | `paper_stellar_offsets.csv` |
| `robustness.py` | combines the neighbours of each host (inverse-variance, element by element), family means, small − giant tests, five chemical clocks, robustness checks | `_state_v3.pkl` |
| `figures_v3.py` | Figures 1–5 from `_state_v3.pkl` | `figures/fig1`–`fig5` |
| `ext1_iso_age.py` | isochrone Δage against [Y/Mg] Δage, same hosts | `iso_ymg_comparison.csv`, `_state_ext1.pkl` |
| `ext1_fig.py` | the isochrone figure | `figures/fig_iso_vs_ymg` |

Hosts are classified by the largest Mp sin i of their planets (small < 0.3 Mj ≤ giant), joined
on the archive `hd_name`: 121 of the 131 hosts, 51 small and 70 giant. Every element is fitted at
its own condensation temperature. `run_pipeline.py --neighbours FILE` runs everything on another
neighbour table (for example from `match_neighbours.py --k 3 --tag k3`), with outputs in
`pipeline/output/<tag>/`.

## How the neighbours are chosen

Each host is compared with the 10 stars without planets closest in
d = √[(ΔTeff/50 K)² + (Δlog g/0.10)² + (Δ[Fe/H]/0.03)²], each difference divided by its typical
uncertainty (`chemclock-matching/match_neighbours.py`). This replaced a rank-mean matching
(candidates within ±160 K, ranked separately on the three differences, ranks averaged, 3 taken),
which has five weaknesses:

1. Ranks discard how large a difference is, and the measurement uncertainties are never used.
2. log g is quoted to 0.01 dex, so many candidates tie, and the tie rule changes which win.
3. The hard Teff window is arbitrary.
4. Three neighbours are noisy: picking 3 of many nearly equivalent candidates moved the result by
   about ±0.5σ on its own.
5. It weights the three parameters equally, although s-process [X/Fe] changes strongly with
   [Fe/H] and only weakly with Teff (about 0.01 dex per 100 K among the controls).

In tests on the same data, rank-mean variants gave an s-process difference of −0.019 to
−0.027 dex (−2.0σ to −3.1σ). Scaled-distance matchings gave about −0.04 dex (−0.042 dex with 10
neighbours, −0.044 dex with 3) and stayed between −3.8σ and −5.2σ (before control noise) for Teff
scales of 27–160 K, either parameter set, and 3 or 10 neighbours. The scaled distance matches
[Fe/H] more closely, which is where the s-process abundances vary most. Those comparisons, and an
older neighbour table that could not be regenerated, are in the git history (commit `d2b41f2`,
`chemclock-analysis/compare_matchings.py`).

## checks/control_bootstrap.py

With 10 neighbours, control stars are shared between hosts, so the pipeline's σ, which treats the
neighbour averages as exact, is too optimistic. The script redoes the matching inside a bootstrap
over the hosts, the controls, or both:

| Small − giant | difference | pipeline σ | bootstrap: hosts | bootstrap: hosts + controls |
|---|---|---|---|---|
| s-process | −0.0425 dex | −4.72 | −5.12 | **−4.45** |
| [Y/Mg] clock | +1.05 Gyr | +2.66 | +2.82 | **+2.49** |

Resampling both tends to overstate the uncertainty a little, so the true value lies between the
last two columns. `--k 3` gives −4.40 / −4.69 / −3.72 for the s-process.

## Fixes made to the code as received

1. **Isochrone ages: hosts looked up by name only.** `ext1_iso_age.py` found host ages by common
   name (`hostname_1`) only, the problem `build_stellar_offsets.py` already describes and fixes.
   Six hosts with ages were dropped (τ Cet, 61 Vir, ε Eri, GJ 3021, ι Hor, HD 26965).
2. **Figure 5 used a typed-in −4.4** for the "all 4" bar of panel (b) and the "Raw" bar of panel
   (f). Both are now computed.
3. **Figure 5 labels were invisible**: white text drawn outside the bars. They are now inside.
4. **Figure 2 claimed p < 0.0001** from 2000 permutations, whose smallest possible p is 1/2000. The
   label now follows the number of permutations.
5. **Figure 3**: the axis of panel (a) was fixed at 2 Gyr and cut off bars; it now follows the data,
   and the legend no longer covers a bar. Panel (b) printed β with a doubled `+` sign.
6. **Figure 1**: the legend covered the Ba label; it is now in the empty lower panel.
7. The pipeline assumed exactly 3 neighbours; it now takes any number.
8. A divide-by-zero warning was printed several hundred times; the computation was already
   correct and is unchanged.

## Open questions (not changed)

- **Selection.** Small planets are found by radial velocity mostly around quiet, inactive stars,
  which tend to be older, and magnetic activity affects Ba abundances. No matching removes this.
- **Planet classes use `hd_name` only** (data README, issue 3): 121 hosts classified, against 124
  with the identifier-and-position crossmatch of `crossmatch_hosts.py`.
- **Figure 5(d)** computes power from the observed effect size, which only restates the observed
  significance; the "robust regime" band in Figure 5(c) is drawn by hand at 1100–1300 K.
- The [Y/Mg] calibration (−0.041 dex/Gyr) is cited as Tucci Maia 2019 in four comments and as
  Tucci Maia et al. 2016 in one (`robustness.py`, `ext1_iso_age.py`); the paper cites 2016.
