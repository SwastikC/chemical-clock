# chemclock-analysis

| | |
|---|---|
| [`pipeline/`](pipeline/) | The paper's analysis: element-family offsets, chemical-clock ages, robustness tests, isochrone cross-check, Figures 1–5. |
| [`checks/neighbour_reuse/`](checks/neighbour_reuse/) | How the results depend on shared neighbours and on the choice of neighbours. |
| [`compare_matchings.py`](compare_matchings.py) | The whole pipeline redone for v6, v2 (recommended), v1 and six other matchings. |
| [`checks/control_bootstrap.py`](checks/control_bootstrap.py) | Uncertainty of the v2 result including the noise of the control stars. |
| [`tc_slopes.py`](tc_slopes.py) | The earlier slope-sign count from the original notebook, with a mode that reproduces it exactly. |
| [`original/`](original/) | The original notebook, for reference. |

```
python pipeline/run_pipeline.py    # the paper's results and figures
python checks/neighbour_reuse/run_checks.py   # sensitivity to the neighbours (about 1.5 minutes)
python compare_matchings.py        # everything redone with v2 and other neighbour tables (about 5 minutes)
python checks/control_bootstrap.py # v2 uncertainty including control-star noise (seconds)
python tc_slopes.py --paper        # the numbers of the original notebook, exactly
python tc_slopes.py                # corrected slope-sign count: each element at its own Tc, one vote per host
```

Outputs go to `pipeline/output/` and `output/` (not committed).

## pipeline/: the paper's analysis

Reads the paper's input table (`final_test_v0.csv`), the v6 neighbour table, the HARPS-GTO
parameters and ages, and the archive snapshot; all paths are in `pipeline/paths.py`.
`run_pipeline.py` runs the five scripts in order:

| Script | What it does | Writes |
|---|---|---|
| `build_stellar_offsets.py` | host minus mean-of-neighbours ΔTeff, Δlog g, Δ[Fe/H] | `paper_stellar_offsets.csv` |
| `robustness.py` | combines the neighbours of each host (inverse-variance), family means, small − giant tests, five chemical clocks, ten robustness checks | `_state_v3.pkl` |
| `figures_v3.py` | Figures 1–5 from `_state_v3.pkl` | `figures/fig1`–`fig5` |
| `ext1_iso_age.py` | isochrone Δage against [Y/Mg] Δage, same hosts | `iso_ymg_comparison.csv`, `_state_ext1.pkl` |
| `ext1_fig.py` | the isochrone figure | `figures/fig_iso_vs_ymg` |

Hosts are classified by the largest Mp sin i of their planets (small < 0.3 Mj ≤ giant), joined
on the archive `hd_name`: 121 of the 131 hosts, 51 small and 70 giant. The expected main numbers
are in the [top-level README](../README.md#step-4-the-papers-analysis-and-figures). Elements are
paired with their own Tc throughout; the problem of the original notebook does not occur here.
Running on the legacy table (before the Cr correction) changes no s-process number.

### Changes made to the pipeline as received

The first commit that adds `pipeline/` holds the scripts as received, with only the input paths
changed; its output is identical to the original scripts. The changes after that:

1. **Isochrone ages: hosts looked up by name only.** `ext1_iso_age.py` found host ages by common
   name (`hostname_1`) only, the same problem `build_stellar_offsets.py` describes and fixes. Six
   hosts with ages were dropped (τ Cet, 61 Vir, ε Eri, GJ 3021, ι Hor, HD 26965). With them the
   comparison has 99 hosts instead of 94: Pearson r +0.22 → +0.27, isochrone small − giant
   +2.03σ → +1.97σ, combined +3.30σ → +3.38σ.
2. **Figure 5 used a typed-in −4.4.** The "all 4" bar of panel (b) and the "Raw" bar of panel (f)
   were the number −4.4 written in the code, not computed, so they would not follow a change of
   input. They are now computed (−4.40 with the paper's input, so the figure is unchanged).
3. **Figure 5 labels were invisible.** The σ labels of panels (b) and (f) were white text drawn
   outside the bars; they are now inside.
4. **Figure 2 claimed p < 0.0001** when no permutation beat the observed value, but it runs 2000
   permutations, so the smallest p it can show is 1/2000. The label now says p < 0.0005.
5. A divide-by-zero warning was printed several hundred times; the computation was already
   correct and is unchanged.

### Neighbour reuse and the choice of neighbours

Tested in [`checks/neighbour_reuse/`](checks/neighbour_reuse/) (`python checks/neighbour_reuse/run_checks.py`):

- **Reuse is a small effect.** 83 of the 270 neighbours serve 2–4 of the 131 hosts, and 38 serve
  both a small- and a giant-planet host. Accounting for it lowers the s-process result from
  −4.40σ to −4.30σ (variance inflation 1.02–1.03 for every quantity). Leaving out any one
  neighbour keeps it between −4.05σ and −4.52σ.
- **The choice of neighbours is a large effect.** The rank-mean matching the paper describes,
  redone here, gives −2.97σ with the General Properties parameters and −2.27σ with the HARPS
  columns, with or without replacement. Random choices of 3 among each host's 10 best candidates
  give a median of −3.80σ (−3.47σ with HARPS columns); the stored table's −4.40σ is at the 13th
  (3rd) percentile of those draws. The clock ages move the same way.
- **The stored neighbours are looser matches than the stated method gives:** median |ΔTeff| 43–45 K
  against 15–17 K. Together with data issue 6 (the stored table cannot be regenerated), this is
  the main open question for the paper: the headline significance depends on a neighbour table
  whose construction is not known.

### Choosing the neighbours

The paper's rank-mean method has five weaknesses, each seen in the tests above:

1. **Ranks discard how large a difference is.** The best-ranked candidate is first whether it is
   5 K or 60 K away, depending on how crowded the Teff window is, and the measurement
   uncertainties of the three parameters are never used.
2. **Ties decide neighbours.** log g is quoted to 0.01 dex, so many candidates tie, and the tie
   rule changes which neighbours win.
3. **The hard Teff window and the parameter set are arbitrary,** and the result follows them:
   rank-mean variants give −2.0σ to −3.1σ, the looser stored v6 table −4.4σ (table below).
4. **Three neighbours are noisy.** Picking 3 of many nearly equivalent candidates moves the result
   by about ±0.5σ on its own (random-draw test above).
5. **It weights the parameters equally,** although s-process [X/Fe] changes strongly with [Fe/H]
   and only weakly with Teff (about 0.01 dex per 100 K among the controls).

**v2 uses a scaled distance instead**
(`match_neighbours.py --method scaled --k 10 --drop-archive-hosts`): every control star gets
d = √[(ΔTeff/50 K)² + (Δlog g/0.10)² + (Δ[Fe/H]/0.03)²], each difference divided by the typical
HARPS-GTO uncertainty quoted in the paper, and the 10 smallest distances are the neighbours. No
window, no ranks, almost no ties, and 10 neighbours average out much of the noise of a single
control star. General Properties parameters; archive planet hosts removed from the controls. Other
scales can be tried with `--scales`; the method is described in the
[data README](../chemclock-data/README.md#v2-a-better-matching-recommended).

With 10 neighbours, controls are shared between hosts, so the pipeline's σ, which treats the
neighbour averages as exact, is too optimistic. `checks/control_bootstrap.py` redoes the matching
inside a bootstrap over the hosts, the controls, or both:

| s-process, small − giant | difference | pipeline σ | bootstrap: hosts | bootstrap: hosts + controls |
|---|---|---|---|---|
| v2 (10 neighbours) | −0.0425 dex | −4.72 | −5.12 | **−4.45** |
| scaled, 3 neighbours | −0.0443 dex | −4.40 | −4.69 | −3.72 |

Resampling both tends to overstate the uncertainty (see `checks/neighbour_reuse/`), so the true
value lies between the last two columns. For the [Y/Mg] clock with v2 the same range is
+2.5σ to +2.8σ.

### Results for every matching

`python compare_matchings.py` runs the whole pipeline on each neighbour table and writes
`output/matching_comparison.csv`. v2 and v1 are committed in `chemclock-data/data/reference/`,
and the script checks that both are regenerated exactly.

| | v6 (paper) | **v2** | scaled, k=3 | v1 (rank-mean) | rank-mean HARPS | … no repl. | rank-mean General | … no repl. | … hosts removed |
|---|---|---|---|---|---|---|---|---|---|
| s-process diff (dex) | −0.039 | **−0.043** | −0.044 | −0.021 | −0.021 | −0.019 | −0.026 | −0.027 | −0.027 |
| s-process σ | −4.40 | **−4.72** | −4.40 | −2.31 | −2.27 | −1.99 | −2.97 | −3.04 | −3.10 |
| … after Δ[Fe/H] | −4.54 | **−4.72** | −4.55 | −2.43 | −2.40 | −2.19 | −3.23 | −3.39 | −3.49 |
| … after ΔTeff, Δlog g, Δ[Fe/H] | −4.72 | **−4.42** | −4.09 | −2.49 | −2.50 | −2.37 | −3.22 | −3.38 | −3.48 |
| [Y/Mg] Δage (Gyr) | +0.84 | **+1.05** | +0.84 | +0.30 | +0.25 | +0.11 | +0.49 | +0.50 | +0.55 |
| [Y/Mg] σ | +2.19 | **+2.66** | +1.89 | +0.76 | +0.63 | +0.26 | +1.31 | +1.24 | +1.48 |
| Tc slope σ, all elements | −2.21 | **−2.23** | −1.99 | −0.73 | −0.68 | −0.40 | −1.37 | −1.36 | −1.41 |
| Tc slope σ, without s-process | −0.32 | **+0.45** | −0.02 | +0.39 | +0.37 | +0.56 | +0.06 | +0.12 | +0.11 |
| isochrone vs [Y/Mg] age, r | +0.27 | **+0.48** | +0.51 | +0.39 | +0.42 | +0.30 | +0.42 | +0.36 | +0.45 |
| combined age σ | +3.38 | **+3.55** | +3.09 | +3.02 | +2.98 | +1.89 | +2.89 | +2.01 | +3.47 |

With the scaled matching, the paper's results come back at about the paper's size: small-planet
hosts are lower in s-process elements than giant-planet hosts by 0.043 dex (−4.4σ to −4.7σ with
control noise and stellar-parameter regression included); the [Y/Mg] clock gives +1.05 Gyr
(+2.5σ to +2.8σ); removing the s-process elements removes the Tc-slope difference
(−2.2σ → +0.4σ); and the isochrone and [Y/Mg] ages agree better than with v6 (r = 0.48). The
rank-mean matchings give about half the s-process difference; with the scaled distance the
neighbours match [Fe/H] more closely, which is where the s-process abundances vary most.

What no matching can remove: small planets are found by radial velocity mostly around quiet,
inactive stars, which tend to be older, and magnetic activity is known to affect Ba abundances.
A detection bias of this kind would make small-planet hosts look older for reasons unrelated
to planet formation, and should be discussed alongside the result.

The figure texts in `figures_v3.py` (for example "Five chemical clocks all agree", "No correlation
with Δ[Fe/H]") were written for the v6 results and are not updated for other matchings.

### Other open questions (not changed)

- **Planet classes use `hd_name` only** (data README, issue 3): 121 hosts classified, against 124
  with the identifier-and-position crossmatch of `crossmatch_hosts.py`.
- **Archive hosts in the control pool** (data README, issue 2): 10 of the v6 neighbours are
  themselves planet hosts.
- **Figure 5(d)** computes power from the observed effect size, which only restates the observed
  significance; and the "robust regime" band in Figure 5(c) is drawn by hand at 1100–1300 K.
- The [Y/Mg] calibration (−0.041 dex/Gyr) is cited as Tucci Maia 2019 in four comments and as
  Tucci Maia et al. 2016 in one (`robustness.py`, `ext1_iso_age.py`); check which the paper means.

## tc_slopes.py

`--count hosts` (the default) needs `../chemclock-matching/output/host_crossmatch.csv`, so run
`chemclock-matching/crossmatch_hosts.py` first.

For every host and each of its three neighbours, a straight line is fitted (ordinary least
squares) to Δ[X/Fe] against Tc for the 18 elements (Zn, Cu, Na, Mn, Cr, Si, Mg, Co, Ni, V, Ba,
Ce, Ca, Sc, Sr, Y, Ti, Al; Tc from `condensation_temps.csv`). Each host gets a vote: +1 if at least
two of its three slopes are positive, −1 if at least two are negative. The votes are then counted
for small (< 0.3 Mj) and giant planets.

| Flag | Choices | Default | `--paper` |
|---|---|---|---|
| `--table` | `v0`, `legacy` (before the Cr correction) | `v0` | `legacy` |
| `--pairs FILE` | a `pairs_<tag>.csv` from `match_neighbours.py`, instead of `--table` | | |
| `--pairing` | `correct`, `reversed` (see below) | `correct` | `reversed` |
| `--count` | `hosts` (one vote per host), `planets` (one row per planet) | `hosts` | `planets` |
| `--exclude` | elements to leave out, e.g. `--exclude Zn` | none | none |
| `--cut` | small/giant split, Mj | 0.3 | 0.3 |

Expected output:

| Command | Small + | Small − | Giant + | Giant − | Counted |
|---|---|---|---|---|---|
| `--paper` | 63 | 52 | 44 | 45 | planets |
| `--table legacy --pairing correct --count planets` | 51 | 64 | 43 | 46 | planets |
| `--count planets` | 51 | 64 | 44 | 45 | planets |
| (defaults) | 26 | 26 | 38 | 34 | hosts |

## Problems found in the original notebook

These are why `--paper` and the defaults differ. Each one is a flag, so its effect can be seen on
its own.

1. **Elements fitted at the wrong Tc (`--pairing`).** The notebook listed the abundances in order of
   increasing Tc (Zn, Cu, Na, …, Al) and the temperatures as they appear in the CSV, in order of
   decreasing Tc (Al, Ti, Y, …, Zn). Each abundance was fitted at the temperature of a different
   element: Zn at 1652 K (Al), Al at 704 K (Zn); no element was at its own Tc. The x axis is therefore
   close to reversed, and the slopes come out with nearly the opposite sign (the correlation
   between the two sets of slopes is −0.88). The majority sign changes for 112 of the 131 hosts,
   and the small-planet result reverses (63+ / 52− becomes 51+ / 64−).
2. **Planets counted instead of hosts (`--count`).** The vote is per host, but the notebook joined it
   to the archive one row per planet. A host with three planets is counted three times, and a host
   with small and giant planets is counted in both groups. The 204 planets counted belong to only
   121 hosts.
3. **Hosts joined on `hd_name` only.** With `--count planets`, hosts without an HD number in the
   archive (HIP 5158, BD-08 2823) are not counted, as in the notebook. `--count hosts` uses the
   identifier-and-position crossmatch instead.
4. **Cr.** The notebook used the abundances before the Cr correction (`--table legacy`); see the
   [data README](../chemclock-data/README.md#v0-changes-from-the-original-files).

The notebook also holds code that was not used for its stored results: 3σ clipping (commented
"2σ"), setting slopes with p > 0.05 to zero, and a later version with 17 elements (no Zn). The
stored results come from an earlier version without these, with Zn included; `--paper`
reproduces them exactly.

## original/

`condensation_temperature_original.ipynb` is the notebook the paper's counts were made with, kept
for reference with its figures removed (they made it 20 MB). It is not runnable as it stands:
it reads files under their old names (`epsilon.csv`, `Neighbour Dataset.csv`,
`Composite_data.csv`), its cells were run out of order, and cells 58 onward are unfinished
scratch work.
