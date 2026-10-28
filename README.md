# chemical-clock

Code and input data for the chemical-clock project: compare the abundances of planet-host stars
in the HARPS-GTO sample with those of their closest "twin" stars without planets.

| Folder | Contents |
|---|---|
| [`chemclock-data/`](chemclock-data/) | Input data (tag **v0**): HARPS-GTO parameters and abundances, NASA Exoplanet Archive snapshot, the paper's reference tables. **Read its README**: it has the column notes and the list of known issues. |
| [`chemclock-matching/`](chemclock-matching/) | Scripts: host crossmatch, rebuild of the paper table, neighbour matching. |
| [`reproduce.py`](reproduce.py) | Runs every step below in order. |

---

## Quick start (about 1 minute)

You need Python 3.9 or newer and `git`.

```bash
git clone https://github.com/SwastikC/chemical-clock.git
```

```bash
cd chemical-clock
```

Create a virtual environment so the packages do not touch your system Python:

```bash
python3 -m venv .venv
```

Activate it. On Linux or macOS:

```bash
source .venv/bin/activate
```

On Windows (PowerShell), run `.venv\Scripts\Activate.ps1` instead.

Install the packages and run everything:

```bash
pip install -r requirements.txt
```

```bash
python reproduce.py
```

The last line should read `outputs in .../chemclock-matching/output`. If any step fails, the
script stops and says which one.

---

## What each step does, and what you should see

You can also run the steps one by one. The scripts can be run from any folder; they find the data
in `chemclock-data/data/` and write to `chemclock-matching/output/`.

### Step 0. Check the data

`reproduce.py` checks every data file against `chemclock-data/SHA256SUMS`. You can also do it by
hand from inside `chemclock-data/` with `sha256sum -c SHA256SUMS` (on macOS:
`shasum -a 256 -c SHA256SUMS`). Every line must say `OK`. If one does not, your copy of the
data is not tag v0 and your numbers will not match the ones below.

### Step 1. Which HARPS-GTO stars host planets?

```bash
python chemclock-matching/crossmatch_hosts.py
```

Matches the 1028 HARPS-GTO stars to the NASA Exoplanet Archive snapshot and splits hosts into
small (< 0.3 Mj) and giant planets. Writes `output/host_crossmatch.csv`. Expected summary:

| Quantity | Expected |
|---|---|
| Archive hosts matched | 178 (177 by identifier, 1 by position) |
| Small / giant | 60 / 118 |
| Flagged as hosts in HARPS but no planet in the archive | 7 |
| Planet hosts sitting in the HARPS non-host (control) pool | 28 |
| Paper's 131 hosts (`final_test_v0`) | 124 in the archive: 52 small, 72 giant |

### Step 2. Rebuild the paper's input table

```bash
python chemclock-matching/rebuild_paper_table.py
```

Takes the three neighbours per host that the paper used (`neighbour_table_v6.csv`) and recomputes
every host-minus-neighbour Δ[X/Fe] and its error from `epsilon_v0.csv`. It compares the result
with the stored `final_test_v0.csv` and must end with:

```
OK: final_test_v0.csv is reproduced exactly
```

This is the table the paper's statistics were computed from. It is written to
`output/final_test_rebuilt.csv`.

### Step 3. Neighbour matching

```bash
python chemclock-matching/match_neighbours.py
```

Picks, for each of the 131 hosts, the 3 non-host stars closest in Teff, log g and [Fe/H] (the
rank-mean method from the paper's Methods section). Writes `neighbours_<tag>.csv` (three
neighbours per host) and `pairs_<tag>.csv` (one row per host-neighbour pair with Δ[X/Fe] for 18
elements). Expected:

| Run | Tag | Unique controls | Same 3 neighbours as v6 | Neighbour overlap with v6 |
|---|---|---|---|---|
| `match_neighbours.py` | `general_wr_w160_k3` | 268 | 11 / 131 | 50% |
| `... --no-replacement` | `general_nr_w160_k3` | 393 | 9 / 131 | 45% |
| `... --params harps --drop-archive-hosts --tag harps_clean` | `harps_clean` | 264 | 11 / 131 | 41% |

Run `python chemclock-matching/match_neighbours.py -h` for all options (Teff window, k,
parameter set, tie handling, and so on). Every option the paper's Methods leave open is a flag.

---

## What can and cannot be reproduced at the moment

- **Reproduced exactly:** the paper's input table `final_test_v0.csv`, from the v6 neighbours and
  the v0 abundances (Step 2).
- **Not reproduced:** the v6 neighbours themselves. Neither these scripts nor the original
  notebook (`chemclock-matching/original/`) gives the same three neighbours for more than about
  10 of the 131 hosts. The procedure that produced `neighbour_table_v6.csv` was not recorded well
  enough to regenerate it. This is open issue 6 in the
  [data README](chemclock-data/README.md#known-issues-in-the-data-open-to-be-resolved-by-the-students).
- **Not in this repository yet:** the analysis code that turns the pair table into the paper's
  figures and statistics (combining the neighbours per host, trends against condensation
  temperature). `pairs_<tag>.csv` and `final_test_v0.csv` are its inputs.

## Rules for reporting a result

1. Quote the data tag (**v0**) and the matching tag (for example `general_wr_w160_k3`).
2. Say which parameter set you used (`general` or `harps`); they are not the same numbers.
3. Do not edit files in `chemclock-data/`. If you think the data is wrong, open an issue.
4. Do not commit anything in `output/`; it is regenerated by the scripts.

## Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'scipy'` (or pandas) | The virtual environment is not active, or you skipped `pip install -r requirements.txt`. |
| `AttributeError: 'DataFrame' object has no attribute 'map'` | pandas is older than 2.1. Run `pip install -U pandas`. |
| `--drop-archive-hosts needs .../host_crossmatch.csv` | Run Step 1 first. |
| A checksum says `FAILED` | Your data differs from tag v0. Clone the repository again. |
| `python3` not found on Windows | Use `python` or `py` instead. |
