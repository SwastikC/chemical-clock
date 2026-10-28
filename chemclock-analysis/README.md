# chemclock-analysis

Condensation-temperature (Tc) slopes of host-minus-neighbour abundances, and the count of
positive and negative slopes for small and giant planets.

```
python tc_slopes.py --paper        # the numbers of the original notebook, exactly
python tc_slopes.py                # corrected analysis: v0 data, each element at its own Tc, one vote per host
```

`--count hosts` (the default) needs `../chemclock-matching/output/host_crossmatch.csv`, so run
`chemclock-matching/crossmatch_hosts.py` first. Outputs go to `output/` (not committed).

## tc_slopes.py

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
