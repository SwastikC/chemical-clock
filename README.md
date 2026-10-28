# chemical-clock

Code and input data for the chemical-clock project: compare the abundances of planet-host stars
in the HARPS-GTO sample with those of their closest "twin" stars without planets.

| Folder | Contents |
|---|---|
| [`chemclock-data/`](chemclock-data/) | Input data (tag **v0**): HARPS-GTO parameters and abundances, NASA Exoplanet Archive snapshot, the paper's reference tables. **Read its README**: it has the column notes and the list of known issues. |
| [`chemclock-matching/`](chemclock-matching/) | Scripts: host crossmatch, rebuild of the paper table, neighbour matching. |
| [`chemclock-analysis/`](chemclock-analysis/) | The paper's analysis pipeline (Figures 1–5) and the earlier slope-sign count. **Read its README**: it lists the problems found and fixed. |
| [`reproduce.py`](reproduce.py) | Runs every step below in order. |

---

## Quick start (about 10 minutes)

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

The last lines should read `=== done` and say where the outputs are. If any step fails, the
script stops and says which one.

---

## What each step does, and what you should see

You can also run the steps one by one. The scripts can be run from any folder; they find the data
in `chemclock-data/data/` and write to the `output/` folder next to the script.

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

### Step 4. The paper's analysis and figures

```bash
python chemclock-analysis/pipeline/run_pipeline.py
```

Runs the paper's pipeline on the paper's input table (`final_test_v0.csv`). For every host it
combines the three neighbours (inverse-variance weighted), averages Δ[X/Fe] within four element
families (volatile, iron-peak, α + Al, s-process), and compares small-planet hosts with
giant-planet hosts. It writes Figures 1–5 and the isochrone figure to
`chemclock-analysis/pipeline/output/figures/`. Expected main results (121 classified hosts: 51
small, 70 giant):

| Result | Expected |
|---|---|
| s-process, small − giant | −0.039 dex, **−4.40σ** |
| … after regressing out Δ[Fe/H] | −4.54σ |
| … after regressing out ΔTeff, Δlog g, Δ[Fe/H] | −4.72σ |
| [Y/Mg] clock, Δage small − giant | +0.84 Gyr (+2.2σ) |
| Tc slope small − giant, with / without s-process | −2.21σ / −0.32σ |
| Isochrone ages vs [Y/Mg] ages (99 hosts) | Pearson r = +0.27 |

The full printout of every test is in the terminal output.
[chemclock-analysis/README.md](chemclock-analysis/README.md#pipeline-the-papers-analysis) lists
the fixes made to the pipeline and the questions it leaves open.

### Step 5. How much do the results depend on the neighbours?

```bash
python chemclock-analysis/checks/neighbour_reuse/run_checks.py
```

Takes about 1.5 minutes. Expected, for the s-process small − giant result:

| Test | s-process σ |
|---|---|
| Stored neighbour table (Step 4) | −4.40 |
| … corrected for neighbours shared between hosts | −4.30 |
| … leaving out any one neighbour | −4.05 to −4.52 |
| Rank-mean matching redone, General Properties parameters | −2.97 (−2.94 without replacement) |
| Rank-mean matching redone, HARPS columns | −2.27 (−1.99 without replacement) |
| Random 3 of each host's 10 best candidates (median, 16–84%) | −3.80 (−4.28 to −3.32) |

Sharing neighbours barely matters; which neighbours are chosen matters a lot. The stored table
gives a stronger result than the matching method it is said to come from, and it cannot be
regenerated (see "What can and cannot be reproduced" below). Details:
[chemclock-analysis/README.md](chemclock-analysis/README.md#neighbour-reuse-and-the-choice-of-neighbours).

### Step 6. Redo everything with a reproducible neighbour table (v2)

```bash
python chemclock-analysis/compare_matchings.py
```

```bash
python chemclock-analysis/checks/control_bootstrap.py
```

The first takes about 5 minutes: it builds the neighbour table again for each matching choice,
runs the whole paper pipeline on each, and prints the main numbers side by side
(`chemclock-analysis/output/matching_comparison.csv`; figures for each in
`chemclock-analysis/pipeline/output/<matching>/figures/`). The second (a few seconds) adds the
noise of the control stars to the uncertainty of the main result.

**v2 is the recommended neighbour table.** Each host gets the 10 control stars closest in
(Teff, log g, [Fe/H]), with every difference measured in units of its typical uncertainty
(50 K, 0.10 dex, 0.03 dex, the values quoted in the paper):
d = √[(ΔTeff/50)² + (Δlog g/0.10)² + (Δ[Fe/H]/0.03)²]. General Properties parameters; archive planet
hosts removed from the controls. Unlike the paper's rank-mean method it needs no Teff window,
does not throw away how large a difference is, has almost no ties, and averages 10 neighbours
instead of 3. It is committed as `chemclock-data/data/reference/neighbour_table_v2.csv` (with
`final_test_v2.csv`), and the script checks that both are regenerated exactly. Why this method:
[chemclock-analysis/README.md](chemclock-analysis/README.md#choosing-the-neighbours).

| | v6 (paper) | **v2 (recommended)** | scaled, 3 neighbours | v1 (rank-mean) | rank-mean, General Prop. |
|---|---|---|---|---|---|
| s-process, small − giant | −0.039 dex, −4.40σ | **−0.043 dex, −4.72σ** | −0.044 dex, −4.40σ | −0.021 dex, −2.31σ | −0.026 dex, −2.97σ |
| … including control-star noise | | **−4.45σ** | −3.72σ | | |
| … after ΔTeff, Δlog g, Δ[Fe/H] | −4.72σ | **−4.42σ** | −4.09σ | −2.49σ | −3.22σ |
| [Y/Mg] Δage, small − giant | +0.84 Gyr (+2.2σ) | **+1.05 Gyr (+2.7σ)** | +0.84 Gyr (+1.9σ) | +0.30 Gyr (+0.8σ) | +0.49 Gyr (+1.3σ) |
| Tc slope, with / without s-process | −2.2σ / −0.3σ | **−2.2σ / +0.4σ** | −2.0σ / −0.0σ | −0.7σ / +0.4σ | −1.4σ / +0.1σ |
| Isochrone vs [Y/Mg] age, r | +0.27 | **+0.48** | +0.51 | +0.39 | +0.42 |
| Isochrone + [Y/Mg] combined age | +3.4σ | **+3.5σ** | +3.1σ | +3.0σ | +2.9σ |

To run the pipeline on a single neighbour table:
`python chemclock-analysis/pipeline/run_pipeline.py --neighbours chemclock-data/data/reference/neighbour_table_v2.csv --tag v2`.

### Step 7. The early slope-sign count (original notebook)

```bash
python chemclock-analysis/tc_slopes.py --paper
```

Reproduces the original notebook, an earlier and simpler analysis: one Tc slope per host-neighbour
pair, and a count of hosts with mostly positive or mostly negative slopes.

| | Small planets (< 0.3 Mj) | Giant planets (> 0.3 Mj) |
|---|---|---|
| Positive slope | 63 | 44 |
| Negative slope | 52 | 45 |

> [!WARNING]
> The notebook fitted each element at the condensation temperature of a different element, which
> nearly reverses the sign of the slopes, and it counted planets rather than hosts. The paper
> pipeline (Step 4) does not have these problems. See
> [chemclock-analysis/README.md](chemclock-analysis/README.md#problems-found-in-the-original-notebook).

The corrected count (`python chemclock-analysis/tc_slopes.py`): small 26 positive / 26 negative,
giant 38 / 34 hosts. To use your own matching, pass a pair table from Step 3 with `--pairs`.

---

## What can and cannot be reproduced at the moment

- **Reproduced exactly:** the paper's input table `final_test_v0.csv`, from the v6 neighbours and
  the v0 abundances (Step 2); the paper's results and figures from that table (Step 4); and the
  notebook's slope-sign counts (Step 7, `--paper`).
- **Not reproduced:** the v6 neighbours themselves. Neither these scripts nor the original
  notebook (`chemclock-matching/original/`) gives the same three neighbours for more than about
  10 of the 131 hosts. The procedure that produced `neighbour_table_v6.csv` was not recorded well
  enough to regenerate it. This is open issue 6 in the
  [data README](chemclock-data/README.md#known-issues-in-the-data-open-to-be-resolved-by-the-students).
  It matters: the stored neighbours are looser matches (median |ΔTeff| 43–45 K, against 15–17 K
  for the rank-mean method) and give a stronger result (−4.40σ) than the rank-mean method does
  when redone (−2.0σ to −3.0σ; Step 5). **v2** (Step 6) is a neighbour table made entirely by
  code with a better matching method, and the whole analysis redone on it: the s-process result
  is −0.043 dex at −4.4σ to −4.7σ.
- **Errors found:** in the original notebook (Step 7) and small ones in the pipeline (Step 4); see
  [chemclock-analysis/README.md](chemclock-analysis/README.md).

## Rules for reporting a result

1. Quote the data tag (**v0**), the matching tag (for example `general_wr_w160_k3`) and the
   analysis options (for example `--pairing correct --count hosts`).
2. Say which parameter set you used (`general` or `harps`); they are not the same numbers.
3. Do not edit files in `chemclock-data/`. If you think the data is wrong, open an issue.
4. Do not commit anything in `output/`; it is regenerated by the scripts.

## Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'scipy'` (or pandas) | The virtual environment is not active, or you skipped `pip install -r requirements.txt`. |
| `AttributeError: 'DataFrame' object has no attribute 'map'` | pandas is older than 2.1. Run `pip install -U pandas`. |
| `--drop-archive-hosts needs .../host_crossmatch.csv` or `--count hosts needs ...` | Run Step 1 first. |
| A checksum says `FAILED` | Your data differs from tag v0. Clone the repository again. |
| `python3` not found on Windows | Use `python` or `py` instead. |
