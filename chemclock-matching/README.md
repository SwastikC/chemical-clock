# chemclock-matching

Host crossmatch, neighbour matching and the host-minus-neighbour table for the chemical-clock
project. The step-by-step guide with expected numbers is in the [top-level README](../README.md);
`python ../reproduce.py` runs everything.

The scripts read `../chemclock-data/data` and write to `output/` next to this file, whatever
folder you run them from (pass `--data` or `--out` to change that). `output/` is not committed.

```
pip install -r ../requirements.txt
python crossmatch_hosts.py      # archive hosts -> HARPS-GTO stars, classes at 0.3 Mj
python match_neighbours.py      # the v2 neighbours, checked against the committed table
python build_table.py           # the host-minus-neighbour table, checked against final_test_v2.csv
```

## crossmatch_hosts.py

Matches every HARPS-GTO star to the archive snapshot: first by identifier (HD, HIP, BD, with a
trailing component letter allowed on the HARPS side), then by position within 60 arcsec. The
positions differ in epoch (J2000 against Gaia), which is why identifiers come first. Each matched
star gets the largest Mp sin i among its planets and a small/giant class. The summary lists stars
whose HARPS host flag disagrees with the archive. `match_neighbours.py` uses the output to remove
archive planet hosts from the controls.

## match_neighbours.py

For each host in `hosts.csv`, the distance to every control star is

d = √[ (ΔTeff / 50 K)² + (Δlog g / 0.10)² + (Δ[Fe/H] / 0.03)² ]

with General Properties parameters, and the 10 closest controls are the neighbours. Controls are
the non-host stars of `epsilon_v0.csv` minus the archive planet hosts; a control can serve
several hosts. With the defaults the result is the committed `neighbour_table_v2.csv`, and the
script stops if it is not. `--k` and `--scales` change the number of neighbours and the scales,
for comparisons (`--tag` names the output).

Each run writes `neighbours_<tag>.csv` (one row per host, closest first) and `pairs_<tag>.csv`
(one row per host-neighbour pair with Δ[X/Fe] and its error for 18 elements, ΔTeff, Δlog g,
Δ[Fe/H] and the distance).

## build_table.py

Builds the table the analysis reads: one row per host, one block of Δ[X/Fe] and errors per
neighbour (prefixes `1_` … `10_`), from `epsilon_v0.csv`. With the default input it checks the
result cell by cell against the committed `final_test_v2.csv`. `--neighbours FILE --out FILE`
builds the table for another neighbour table.
