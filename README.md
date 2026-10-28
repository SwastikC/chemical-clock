# chemical-clock

Code and input data for the chemical-clock project: compare the abundances of planet-host stars
in the HARPS-GTO sample with those of their closest "twin" stars without planets, and ask
whether hosts of small planets and hosts of giant planets differ.

| Folder | Contents |
|---|---|
| [`chemclock-data/`](chemclock-data/) | Input data: HARPS-GTO parameters and abundances, NASA Exoplanet Archive snapshot, the host list and the v2 neighbour table. **Read its README**: it has the column notes and the list of known issues. |
| [`chemclock-matching/`](chemclock-matching/) | Host crossmatch, neighbour matching, and the host-minus-neighbour table. |
| [`chemclock-analysis/`](chemclock-analysis/) | The analysis (Figures 1–6) and the uncertainty check. |
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

The last lines should read `=== done` and say where the figures are. If any step fails, the
script stops and says which one.

---

## What each step does, and what you should see

You can also run the steps one by one. The scripts can be run from any folder; they find the data
in `chemclock-data/data/` and write to the `output/` folder next to the script.

### Step 0. Check the data

`reproduce.py` checks every data file against `chemclock-data/SHA256SUMS`. You can also do it by
hand from inside `chemclock-data/` with `sha256sum -c SHA256SUMS` (on macOS:
`shasum -a 256 -c SHA256SUMS`). Every line must say `OK`. If one does not, your copy of the data
differs from the repository and your numbers will not match the ones below.

### Step 1. Which HARPS-GTO stars host planets?

```bash
python chemclock-matching/crossmatch_hosts.py
```

Matches the 1028 HARPS-GTO stars to the NASA Exoplanet Archive snapshot, by identifier and then by
position, and splits hosts into small (< 0.3 Mj) and giant planets. Writes
`output/host_crossmatch.csv`. Expected:

| Quantity | Expected |
|---|---|
| Archive hosts matched | 178 (177 by identifier, 1 by position) |
| Small / giant | 60 / 118 |
| Flagged as hosts in HARPS but no planet in the archive | 7 |
| Planet hosts among the HARPS non-hosts | 28 (removed from the controls in Step 2) |
| The 131 hosts of the sample | 124 in the archive: 52 small, 72 giant |

### Step 2. Choose the neighbours

```bash
python chemclock-matching/match_neighbours.py
```

For each of the 131 hosts, picks the 10 stars without planets that are closest in
(Teff, log g, [Fe/H]). Each difference is measured in units of its typical uncertainty:

d = √[ (ΔTeff / 50 K)² + (Δlog g / 0.10)² + (Δ[Fe/H] / 0.03)² ]

so a difference counts in proportion to how well it is measured. The 28 non-host stars that do
host planets in the archive are not used. Expected: 1310 pairs, 538 different control stars,
median |ΔTeff| 45 K, |Δlog g| 0.05, |Δ[Fe/H]| 0.03, and the line
`identical to neighbour_table_v2.csv: True`. Why this method:
[chemclock-analysis/README.md](chemclock-analysis/README.md#how-the-neighbours-are-chosen).

### Step 3. Build the host-minus-neighbour table

```bash
python chemclock-matching/build_table.py
```

Computes, for every host and each of its 10 neighbours, Δ[X/Fe] = host − neighbour for 18
elements, with the two errors added in quadrature. Must end with
`OK: final_test_v2.csv is reproduced exactly`.

### Step 4. The analysis and figures

```bash
python chemclock-analysis/pipeline/run_pipeline.py
```

For every host, combines its 10 neighbours (inverse-variance weighted), averages Δ[X/Fe] within
four element families, and compares small-planet hosts with giant-planet hosts (121 hosts with a
planet mass: 51 small, 70 giant). Writes Figures 1–5 and the isochrone figure to
`chemclock-analysis/pipeline/output/figures/`. Expected:

| Family | Small hosts − neighbours | Giant hosts − neighbours | Small − giant |
|---|---|---|---|
| Volatile (Zn, Cu, Na, Mn) | +0.019 ± 0.006 dex | +0.004 ± 0.005 | +0.014 (+2.0σ) |
| Iron-peak (V, Cr, Co, Ni) | +0.014 ± 0.003 | +0.003 ± 0.004 | +0.011 (+2.1σ) |
| α + Al (Mg, Si, Ca, Ti, Sc, Al) | +0.022 ± 0.006 | +0.007 ± 0.004 | +0.015 (+2.0σ) |
| **s-process (Sr, Y, Ba, Ce)** | **−0.040 ± 0.007** | +0.003 ± 0.006 | **−0.043 (−4.7σ)** |

| Test | Expected |
|---|---|
| s-process, after regressing out ΔTeff, Δlog g, Δ[Fe/H] | −4.4σ |
| s-process, leaving out one element | −4.1σ (no Ba) to −4.9σ (no Y) |
| Chemical clocks, small − giant age | +1.05 to +1.37 Gyr, +2.6σ to +3.4σ ([Y/Mg] +1.05 Gyr, +2.7σ) |
| Tc slope small − giant, with / without the s-process elements | −2.2σ / +0.4σ |
| Age difference vs planet mass | −0.67 ± 0.22 Gyr per dex, r = −0.27, p = 0.003 |
| Isochrone ages vs [Y/Mg] ages (101 hosts) | r = +0.48; small − giant: isochrone +0.46 Gyr (+1.6σ), combined +1.29 Gyr (+3.5σ) |
| Two-population model | ΔBIC = +11.1 |

### Step 5. Uncertainty including the control stars

```bash
python chemclock-analysis/checks/control_bootstrap.py
```

The σ values in Step 4 treat each host's 10-neighbour average as exact, but control stars are
shared between hosts and have their own scatter. This step redoes the matching inside a bootstrap
over the hosts and the control stars (a few seconds). Expected for the s-process difference
(−0.0425 dex): −5.1σ resampling the hosts only, **−4.45σ** resampling hosts and controls. For the
[Y/Mg] clock: +2.8σ and +2.5σ. Quote the s-process result as −0.043 dex at about −4.5σ.

---

## Things to keep in mind

- **Selection.** Small planets are found by radial velocity mostly around quiet, inactive stars,
  which tend to be older, and magnetic activity is known to affect Ba abundances. No choice of
  neighbours removes this; it has to be discussed with the result.
- **The other families.** Small-planet hosts are also slightly richer than their neighbours in the
  volatile, iron-peak and α elements (+0.014 to +0.022 dex). The s-process deficit is the only
  negative offset and the only small − giant difference beyond 2.1σ.
- **History.** An earlier draft of the paper used a different neighbour table, which could not be
  regenerated from any code, and an earlier slope-count notebook with an element–temperature
  mix-up. Both have been removed. They remain in the git history (commit `d2b41f2`) for anyone
  who needs to trace the old numbers.

## Rules for reporting a result

1. Quote the neighbour table (**v2**) and any option you changed (for example `--k 3`).
2. Quote the s-process result with the uncertainty from Step 5, not only the σ from Step 4.
3. Do not edit files in `chemclock-data/`. If you think the data is wrong, open an issue.
4. Do not commit anything in `output/`; it is regenerated by the scripts.

## Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: No module named 'scipy'` (or pandas) | The virtual environment is not active, or you skipped `pip install -r requirements.txt`. |
| `AttributeError: 'DataFrame' object has no attribute 'map'` | pandas is older than 2.1. Run `pip install -U pandas`. |
| `needs .../host_crossmatch.csv: run crossmatch_hosts.py first` | Run Step 1 first. |
| A checksum says `FAILED` | Your data differs from the repository. Clone the repository again. |
| `MISMATCH` in Step 2 or 3 | Something changed in the data or the code; clone the repository again. |
| `python3` not found on Windows | Use `python` or `py` instead. |
