# chemclock-data

Input data for the chemical-clock project. Read-only for students: changes arrive by pull request.
Check your copy against `SHA256SUMS` (`sha256sum -c SHA256SUMS`) before quoting a result.

## Files

| File | Rows | Contents |
|---|---|---|
| `data/harps_gto/general_properties.csv` | 1028 | Teff, log g, [Fe/H] with errors, RA, Dec for every HARPS-GTO star. The parameters used for matching |
| `data/harps_gto/starswithplanet.csv` | 220 | HARPS-GTO hosts, one row per planet (147 unique hosts), joined to the archive; isochrone ages |
| `data/harps_gto/starswithoutplanet.csv` | 947 | HARPS-GTO stars without a known planet at the time of the original build; isochrone ages |
| `data/harps_gto/epsilon_v0.csv` | 1025 | Per-star [X/Fe] and errors for the 18 elements, plus the host flag `SwP?` (157 Yes, 868 No). Source of every abundance difference |
| `data/nasa/nasa_archive_snapshot.csv` | 5595 | NASA Exoplanet Archive composite table (`pscomppars`), one row per planet |
| `data/condensation_temps.csv` | 25 | 50% condensation temperatures |
| `data/reference/hosts.csv` | 131 | The planet hosts of the sample |
| `data/reference/neighbour_table_v2.csv` | 131 | The 10 neighbours of each host (columns `SwP`, `First`, `Second`, `Third`, `N4` … `N10`, closest first) |
| `data/reference/final_test_v2.csv` | 131 | Host-minus-neighbour Δ[X/Fe] and errors for each of the 10 neighbours (prefixes `1_` … `10_`) |

The two `reference/*_v2.csv` tables are made by code, and `reproduce.py` checks that they are
regenerated exactly:

```
python chemclock-matching/crossmatch_hosts.py
python chemclock-matching/match_neighbours.py
python chemclock-matching/build_table.py
```

## How the v2 neighbours are chosen

Each host gets the 10 control stars with the smallest
d = √[(ΔTeff/50 K)² + (Δlog g/0.10)² + (Δ[Fe/H]/0.03)²], every difference divided by the typical
HARPS-GTO uncertainty. Parameters from `general_properties.csv`. The controls are the stars with
`SwP? = No` in `epsilon_v0.csv`, minus the 28 that are planet hosts in the archive snapshot. A
control can serve several hosts. Median |ΔTeff| 45 K, |Δlog g| 0.05, |Δ[Fe/H]| 0.03, all within the
measurement uncertainties; 20 of the 1310 pairs differ by more than 150 K. Why this method:
`chemclock-analysis/README.md`.

Every value in `final_test_v2.csv` equals epsilon_v0[host] minus epsilon_v0[neighbour], and every
error is the two star errors added in quadrature.

## Column notes

- Star names: remove spaces before joining (`HD 10700` and `HD10700` are the same star). Never join
  on common names (tau Cet, eps Eri, 61 Vir): join on HD/HIP/BD identifiers or on position.
- In the two `starswith*` files, columns marked `(e) (o)` hold a string `"value error n_lines"`;
  `9.99` means not measured. `Age (e)Gyr` holds `"age error"`.
- Parameters: `general_properties.csv` and the `TeffK`, `logg[cm/s2]`, `[Fe/H]---_x` columns of the
  `starswith*` files are two different parameter sets. The matching uses the first.
- `epsilon_v0.csv` has two iron columns. `[Fe/H]_a` equals `[Fe/H]---_x` of the HARPS-GTO tables and
  is the iron used to convert [X/H] to [X/Fe].
- **Cr.** In the original abundance file, the `[Cr/Fe]` column held [CrI/H] − 0.03 instead of
  [CrI/H] − [Fe/H]. `epsilon_v0.csv` holds [CrI/H] − [Fe/H], checked against all 1078 catalogue
  rows. Nothing else changed.

## Known issues in the data (open, to be resolved by the students)

1. 7 stars flagged as hosts in HARPS-GTO have no planet in this archive snapshot:
   HD111031, HD131664, HD150433, HD157172, HD215456, HD38858, HD9578. They are in `hosts.csv` but
   have no planet mass, so the small/giant comparison leaves them out.
2. 28 stars flagged as non-hosts in HARPS-GTO are hosts in this archive snapshot (masses from
   0.02 to 27 Mj). They are removed from the controls.
3. The analysis classifies hosts by joining on the archive `hd_name` only. Hosts without an HD
   number (BD-08 2823, HIP 5158) and HARPS names with a component letter (HD26965A) are left
   unclassified: 121 classified hosts, against 124 with the crossmatch of `crossmatch_hosts.py`.
4. `hosts.csv` is the host list of the earlier version of the analysis. Two hosts of that version's
   source list, HD10647 and HD108147, are not in it; the reason is not documented.
5. The download date of the archive snapshot is unknown (before July 2025).

## Sources

- Adibekyan et al. 2012, A&A 545, A32 (parameters, alpha and iron-peak abundances)
- Delgado Mena et al. 2017, A&A 606, A94 (neutron-capture abundances)
- NASA Exoplanet Archive, Planetary Systems Composite Parameters
