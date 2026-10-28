# chemical-clock

Code and input data for the chemical-clock project.

| Folder | Contents |
|---|---|
| [`chemclock-data/`](chemclock-data/) | Input data (tag **v0**): HARPS-GTO parameters and abundances, NASA Exoplanet Archive snapshot, reference tables. Read its README for column notes and known issues. |
| [`chemclock-matching/`](chemclock-matching/) | Host crossmatch and rank-mean neighbour matching scripts. |

## Getting started

```
git clone https://github.com/SwastikC/chemical-clock.git
cd chemical-clock/chemclock-data
shasum -a 256 -c SHA256SUMS        # check your copy of the data (Linux: sha256sum -c SHA256SUMS)
cd ../chemclock-matching
pip install numpy scipy pandas
python crossmatch_hosts.py
python match_neighbours.py
```

The scripts read the data from `../chemclock-data/data` by default, so run them from inside
`chemclock-matching/`. Outputs go to `chemclock-matching/output/`.
