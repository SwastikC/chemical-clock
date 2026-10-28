"""Input files and output folder for the pipeline. Edit here if you move things.

The inputs are the tagged files in chemclock-data/data. The names in brackets are
the ones the scripts used before they were moved into this repository.
"""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
D = ROOT / 'chemclock-data' / 'data'
FILES = {
    'neighbours': D / 'reference/neighbour_table_v6.csv',    # (Neighbour Dataset.csv)
    'final_test': D / 'reference/final_test_v0.csv',         # (final_test.csv)
    'properties': D / 'harps_gto/general_properties.csv',    # (Properties_f.csv)
    'swp':        D / 'harps_gto/starswithplanet.csv',
    'swop':       D / 'harps_gto/starswithoutplanet.csv',
    'tc':         D / 'condensation_temps.csv',              # (Condensation Temps.csv)
    'archive':    D / 'nasa/nasa_archive_snapshot.csv',      # (Comp_data.csv)
}
OUT = str(Path(__file__).resolve().parent / 'output')
(Path(OUT) / 'figures').mkdir(parents=True, exist_ok=True)
