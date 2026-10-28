# Neighbour-reuse tests (open issue 1)

Do the paper's results depend on neighbours being shared between hosts, or on which neighbours
were chosen? `python run_checks.py` runs the four scripts in order (about 1.5 minutes); each
prints its results and writes CSVs to `output/`. They can be run from any folder and read the
tagged data through `../../pipeline/paths.py`.

```
python reuse_bootstrap.py    # sharing structure, host/neighbour/two-way bootstraps, leave-one-neighbour-out
python calibrate.py          # simulation on the real host-neighbour graph; exact reuse inflation per quantity
python rematch.py            # rank-mean matching redone with and without replacement (Hungarian assignment)
python matching_spread.py    # 400 random draws of 3 from each host's 10 best candidates
```

`core.py` holds the pair-level engine. On the stored table it reproduces the v6 values exactly
(s-process -0.0394 dex, 4.40 sigma; [Y/Mg] +0.84 Gyr, 2.19 sigma; etc.).

## Checked in this repository

- Every number under **Results** below is reproduced exactly with tag v0. The checks were written
  against the files before the Cr correction (`epsilon.csv`, `final_test.csv`, now in
  `chemclock-data/data/legacy/`); Cr enters neither the s-process family nor the clocks, so v0
  gives the same numbers.
- `rematch.py` and `chemclock-matching/match_neighbours.py` are independent implementations of the
  same rank-mean matching. With replacement they choose the same three neighbours for all 131
  hosts, with both parameter sets. Without replacement they agree for 126 of 131; the two
  assignments have the same total mean rank (9211.33), so the other 5 are ties between equally
  good solutions.
- The stored pairs' median |dTeff| is 43 K (General Properties) or 45 K (HARPS columns), against
  15-17 K for the rank-mean method.

## Data provenance established here

The file names in this section are those of the files before the Cr correction (now
`chemclock-data/data/legacy/`).

- `final_test.csv` pair values are epsilon[host] - epsilon[neighbour] and errors are the quadrature sum,
  exactly, for all 18 elements and all 393 pairs (`epsilon.csv`, 157 hosts + 868 non-hosts).
- The `[Cr/Fe]` column in epsilon.csv, and therefore in final_test, holds [Cr/H]. Every one of the 357
  pairs checkable against the HARPS catalogues matches delta[CrI/H], not delta[CrI/Fe]. Cr enters the
  iron-peak family, the Tc slopes and the refractory-volatile index. It does not enter the s-process
  family or any clock.
- The stored neighbour table cannot be regenerated with the notebook algorithm
  (`Sequential Approach- T then G, M.ipynb`) from either parameter set: at most 14 of 133 hosts get the
  same three neighbours. Stored pairs have median |dTeff| 44 K; a strict rank-mean top-3 gives 15-17 K.

## Results

Sharing: 270 unique neighbours over 393 slots in the 131-host set; 83 neighbours serve 2-4 hosts;
38 serve both a small and a giant host. Hosts form 47 linked clusters, the largest with 42 hosts.

Reuse itself is a small effect. Exact variance inflation relative to an identical design with no
reuse is 1.02-1.03 for every quantity (s-process 4.40 -> 4.30 sigma). Leave-one-neighbour-out keeps the
s-process within 4.05-4.52 sigma; dropping HD139590 gives 4.25; dropping the 10 most-reused at once, 4.37.

The two-way (pigeonhole) bootstrap gives 3.14 sigma, but the simulation shows it overstates the SE by
about 35% even with zero reuse, because each host value already contains its neighbours' scatter. It
should not be quoted.

Matching choice is a large effect. Like-for-like, with vs without replacement changes nothing
(General Properties params: -2.97 vs -2.94 sigma; HARPS columns: -2.27 vs -1.99). But both are far
below the stored table's -4.40. Random 3-of-10 draws give a median of -3.80 sigma
(16-84%: -4.28 to -3.32) with General Properties params and -3.47 (-3.99 to -2.99) with HARPS columns;
the stored table sits at the 13th and 3rd percentile. See `rematch_summary.csv` and
`matching_spread_*.csv` for the clocks.
