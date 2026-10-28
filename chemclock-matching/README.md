# chemclock-matching

Host crossmatch and neighbour matching for the chemical-clock project. The step-by-step guide
with expected numbers is in the [top-level README](../README.md); `python ../reproduce.py` runs
everything.

The scripts read `../chemclock-data/data` and write to `output/` next to this file, whatever
folder you run them from (pass `--data` or `--out` to change that). `output/` is not committed.

```
pip install -r ../requirements.txt
python crossmatch_hosts.py                  # archive hosts -> HARPS-GTO stars, classes at 0.3 Mj
python rebuild_paper_table.py               # rebuild final_test_v0.csv and check it
python match_neighbours.py                  # rank-mean, General Properties parameters, with replacement
python match_neighbours.py --no-replacement
python match_neighbours.py --params harps --drop-archive-hosts --tag harps_clean
```

## crossmatch_hosts.py

Matches every HARPS-GTO star to the archive snapshot: first by identifier (HD, HIP, BD, with a
trailing component letter allowed on the HARPS side), then by position within 60 arcsec. The
positions differ in epoch (J2000 against Gaia), which is why identifiers come first. Each matched
star gets the largest Mp sin i among its planets and a small/giant class. The summary lists stars
whose HARPS host flag disagrees with the archive.

## rebuild_paper_table.py

Rebuilds the paper v6 input table `final_test_v0.csv` from `neighbour_table_v6.csv` and
`epsilon_v0.csv` (host minus neighbour, errors in quadrature), writes
`output/final_test_rebuilt.csv` and compares it cell by cell with the stored table. It exits with
an error if they differ. With tag v0 the largest difference is about 1e-16.

## match_neighbours.py

Rank-mean matching: controls within ±160 K, ranked separately on |dTeff|, |dlogg| and |d[Fe/H]|,
the three ranks averaged, the k lowest taken. Without replacement, the assignment is solved
globally (Hungarian algorithm), so the result does not depend on host order.

Every choice the paper's Methods leave open is a flag: parameter set, Teff window, k, replacement,
tie ranking, ranking the host with its candidates, removing archive hosts from the control pool,
and the host list. Run `python match_neighbours.py -h` for the list.

Each run writes `neighbours_<tag>.csv` and `pairs_<tag>.csv`. The pair table is the agreed schema:
one row per host-neighbour pair with d[X/Fe] and its error for 18 elements (from
`epsilon_v0.csv`), dTeff, dlogg and d[Fe/H]. Combining the neighbours per host and all statistics
are left to the analysis code.

The run also reports how many hosts get the same three neighbours as the stored v6 table. With the
defaults this is 11 of 131, with about half of all neighbours shared.

## original/

`rank_mean_original.ipynb` is the notebook the v6 neighbours were first built with. It is kept
for reference only: it expects the CSV files in its own folder and asks for one host at a time,
and running its method over all hosts does not reproduce the v6 table either. It ranks with
`argsort` (ties broken by sort order), includes the host in its own ranking, uses the HARPS-GTO
parameter columns and the full non-host table as the pool, and keys hosts on `hostname_1`.
