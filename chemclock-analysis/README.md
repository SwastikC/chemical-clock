# chemclock-analysis

| | |
|---|---|
| [`pipeline/`](pipeline/) | The paper's analysis: element-family offsets, chemical-clock ages, robustness tests, isochrone cross-check, Figures 1–5. |
| [`tc_slopes.py`](tc_slopes.py) | The earlier slope-sign count from the original notebook, with a mode that reproduces it exactly. |
| [`original/`](original/) | The original notebook, for reference. |

```
python pipeline/run_pipeline.py    # the paper's results and figures
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
| `robustness.py` | combines the 3 neighbours per host (inverse-variance), family means, small − giant tests, five chemical clocks, ten robustness checks | `_state_v3.pkl` |
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

### Open questions (not changed)

- **Hosts are not independent.** 85 of the 273 neighbour stars serve more than one host (up to
  five), and 38 serve both a small- and a giant-planet host. The σ values treat the hosts as
  independent, so they are probably somewhat too large. Comparing with a matching without
  replacement (`match_neighbours.py --no-replacement`) would show by how much.
- **Planet classes use `hd_name` only** (data README, issue 3): 121 hosts classified, against 124
  with the identifier-and-position crossmatch of `crossmatch_hosts.py`.
- **Archive hosts in the control pool** (data README, issue 2): 10 of the v6 neighbours are
  themselves planet hosts.
- **Figure 5(d)** computes power from the observed effect size, which only restates the observed
  significance; and the "robust regime" band in Figure 5(c) is drawn by hand at 1100–1300 K.
- The [Y/Mg] calibration (−0.041 dex/Gyr) is cited as Tucci Maia 2019 in four comments and as
  Tucci Maia et al. 2016 in one (`robustness.py`, `ext1_iso_age.py`); check which the paper means.
- The pipeline reads the v6 table only. To run it on a new matching, `final_test` and
  `neighbours` in `paths.py` would have to point at tables in the same layout.

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
