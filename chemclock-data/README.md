# chemclock-data

Input data for the chemical-clock project. Read-only for students: new versions (v1, the
matching ensemble) arrive by pull request and are tagged. Check your copy against
`SHA256SUMS` (`sha256sum -c SHA256SUMS`) and quote the tag you used in every result.

Current tag: **v0**. Proposed next: adds two reproducible neighbour tables, **v2** (recommended) and v1; see below.

## Files

| File | Rows | Contents |
|---|---|---|
| `data/harps_gto/general_properties.csv` | 1028 | Teff, log g, [Fe/H] with errors, RA, Dec for every HARPS-GTO star |
| `data/harps_gto/starswithplanet.csv` | 220 | HARPS-GTO hosts, one row per planet (147 unique hosts), joined to the archive |
| `data/harps_gto/starswithoutplanet.csv` | 947 | HARPS-GTO stars without a known planet at the time of the original build |
| `data/harps_gto/epsilon_v0.csv` | 1025 | Per-star [X/Fe] and errors for the 18 elements, plus the host flag `SwP?` (157 Yes, 868 No). Source of every abundance difference |
| `data/nasa/nasa_archive_snapshot.csv` | 5595 | NASA Exoplanet Archive composite table (`pscomppars`), one row per planet |
| `data/reference/neighbour_table_v6.csv` | 133 | The three neighbours per host used in paper v6 |
| `data/reference/final_test_v0.csv` | 131 | Paper v6 input table: host-minus-neighbour d[X/Fe] and errors for each of the three neighbours (`1_`, `2_`, `3_` prefixes) |
| `data/reference/neighbour_table_v2.csv` | 131 | **Recommended.** Ten neighbours per host from the scaled-distance matching (v2, see below). Regenerated exactly by the code |
| `data/reference/final_test_v2.csv` | 131 | Same layout as `final_test_v0.csv`, with ten neighbour blocks (`1_` … `10_`) |
| `data/reference/neighbour_table_v1.csv` | 131 | Three neighbours per host from the documented rank-mean matching (v1, superseded by v2) |
| `data/reference/final_test_v1.csv` | 131 | Same layout as `final_test_v0.csv`, for the v1 neighbours |
| `data/condensation_temps.csv` | 25 | 50% condensation temperatures |
| `data/legacy/` | | The two files before the Cr correction. Do not use; kept for traceability |

Every pair value in `final_test_v0.csv` equals epsilon_v0[host] minus epsilon_v0[neighbour], and
every error is the two star errors added in quadrature.

## v2: a better matching (recommended)

```
python chemclock-matching/crossmatch_hosts.py
python chemclock-matching/match_neighbours.py --method scaled --k 10 --drop-archive-hosts --tag v2
python chemclock-matching/rebuild_paper_table.py --neighbours chemclock-matching/output/neighbours_v2.csv --table-out final_test_v2.csv
```

Each host gets the 10 control stars with the smallest
d = √[(ΔTeff/50 K)² + (Δlog g/0.10)² + (Δ[Fe/H]/0.03)²], every difference divided by the typical
HARPS-GTO uncertainty quoted in the paper. General Properties parameters; archive planet hosts
removed from the control pool (issue 2); a control can serve several hosts; the 131 hosts of
`final_test_v0.csv`. No Teff window and no ranks. Median |dTeff| 45 K, |dlogg| 0.05,
|d[Fe/H]| 0.03 (all within the measurement uncertainties); 20 of the 1310 pairs differ by more
than 150 K. Why this replaces the rank-mean method, and the results, are in
`chemclock-analysis/README.md` ("Choosing the neighbours"). `chemclock-analysis/compare_matchings.py`
regenerates both files and stops if they differ from the committed ones.

## v1: a reproducible neighbour table (rank-mean, superseded by v2)

`neighbour_table_v6.csv` cannot be regenerated (issue 6), and its neighbours are looser matches
than the method described in the paper gives (median |dTeff| 43-45 K against 15-17 K). v1
replaces it with a table made by code:

```
python chemclock-matching/crossmatch_hosts.py
python chemclock-matching/match_neighbours.py --params harps --drop-archive-hosts --tag v1
python chemclock-matching/rebuild_paper_table.py --neighbours chemclock-matching/output/neighbours_v1.csv --table-out final_test_v1.csv
```

Choices, fixed before any result was looked at: the rank-mean method (Teff window 160 K, 3
neighbours), the HARPS-GTO parameter columns as in the original notebook, archive planet hosts
removed from the control pool (issue 2), neighbours may serve more than one host. The 131 hosts of
`final_test_v0.csv`. Median |dTeff| 16 K; 11 hosts have the same three neighbours as v6.
`chemclock-analysis/compare_matchings.py` regenerates both files and stops if they differ from
the committed ones.

## Column notes

- Star names: remove spaces before joining (`HD 10700` and `HD10700` are the same star). Never join
  on common names (tau Cet, eps Eri, 61 Vir): join on HD/HIP/BD identifiers or on position.
- In the two `starswith*` files, columns marked `(e) (o)` hold a string `"value error n_lines"`;
  `9.99` means not measured.
- Parameters: `general_properties.csv` and the `TeffK`, `logg[cm/s2]`, `[Fe/H]---_x` columns of the
  `starswith*` files are two different parameter sets. Say which one you use.
- `epsilon_v0.csv` has two iron columns. `[Fe/H]_a` equals `[Fe/H]---_x` of the HARPS-GTO tables and
  is the iron used to convert [X/H] to [X/Fe].

## v0 changes from the original files

- **Cr.** In the original `epsilon.csv`, the `[Cr/Fe]` column held [CrI/H] − 0.03 instead of
  [CrI/H] − [Fe/H]. v0 holds [CrI/H] − [Fe/H], checked against all 1078 catalogue rows. The
  per-pair Cr values in `final_test_v0.csv` are corrected the same way. Nothing else changed.

## Known issues in the data (open, to be resolved by the students)

1. 7 stars flagged as hosts in HARPS-GTO have no planet in this archive snapshot:
   HD111031, HD131664, HD150433, HD157172, HD215456, HD38858, HD9578.
2. 28 stars flagged as non-hosts in HARPS-GTO are hosts in this archive snapshot (masses from
   0.02 to 27 Mj). Ten of them are neighbours in the v6 table.
3. The paper classified hosts by joining on the archive `hd_name` only. Hosts without an HD number
   (BD-08 2823, HIP 5158) and HARPS names with a component letter (HD26965A) were left unclassified.
4. `neighbour_table_v6.csv` has 133 hosts; `final_test_v0.csv` has 131. HD10647 and HD108147 are
   missing from the second; the reason is not documented.
5. The download date of the archive snapshot is unknown (before July 2025).
6. The procedure that produced `neighbour_table_v6.csv` is not documented well enough to regenerate it.
   v2 (above) is a reproducible replacement made with a better matching method; the results with
   it are in `chemclock-analysis/README.md`.

## Sources

- Adibekyan et al. 2012, A&A 545, A32 (parameters, alpha and iron-peak abundances)
- Delgado Mena et al. 2017, A&A 606, A94 (neutron-capture abundances)
- NASA Exoplanet Archive, Planetary Systems Composite Parameters
