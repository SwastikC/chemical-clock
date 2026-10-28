"""Redo the paper's analysis for several neighbour tables and compare the results.

For each matching below: build the neighbour table with
chemclock-matching/match_neighbours.py, run the paper pipeline on it
(pipeline/run_pipeline.py --neighbours), and collect the main numbers. The stored
v6 table (the paper) is the first row. v1 is the matching proposed as the
reproducible replacement for v6; its regenerated neighbours and host-minus-neighbour
table are checked against the committed neighbour_table_v1.csv and final_test_v1.csv
in chemclock-data/data/reference/.

Output: output/matching_comparison.csv, and a table on screen. Takes about 3 minutes.

Usage:  python compare_matchings.py
"""
import pickle
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MATCH = ROOT / 'chemclock-matching'
PIPE = HERE / 'pipeline'
OUT = HERE / 'output'; OUT.mkdir(parents=True, exist_ok=True)
V1 = ROOT / 'chemclock-data' / 'data' / 'reference' / 'neighbour_table_v1.csv'

# tag: (description, match_neighbours.py options)
MATCHINGS = {
    'v1':                 ('HARPS params, archive hosts removed, with replacement (proposed v1)',
                           ['--params', 'harps', '--drop-archive-hosts']),
    'harps_wr':           ('HARPS params, with replacement', ['--params', 'harps']),
    'harps_nr':           ('HARPS params, without replacement', ['--params', 'harps', '--no-replacement']),
    'general_wr':         ('General Properties params, with replacement', ['--params', 'general']),
    'general_nr':         ('General Properties params, without replacement', ['--params', 'general', '--no-replacement']),
    'general_clean':      ('General Properties params, archive hosts removed, with replacement',
                           ['--params', 'general', '--drop-archive-hosts']),
}


def run(cmd):
    r = subprocess.run([sys.executable] + [str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f'failed: {" ".join(map(str, cmd))}\n{r.stdout}\n{r.stderr}')
    return r.stdout


def sigma(s, g):
    return (s.mean() - g.mean()) / np.sqrt(s.var(ddof=1) / len(s) + g.var(ddof=1) / len(g))


def summary(out_dir):
    S = pickle.load(open(out_dir / '_state_v3.pkl', 'rb'))
    E = pickle.load(open(out_dir / '_state_ext1.pkl', 'rb'))
    slope = {k: S[k].groupby('class')['slope'] for k in ('slopes_full', 'slopes_pure')}
    ymg = S['clock_results']['[Y/Mg]']
    return {
        'hosts small/giant': f"{len(S['small_df'])}/{len(S['giant_df'])}",
        's-process diff (dex)': S['s_sp'].mean() - S['g_sp'].mean(),
        's-process sigma': sigma(S['s_sp'], S['g_sp']),
        'sigma after d[Fe/H]': S['obj1']['sigma_resid'],
        'sigma after dTeff,dlogg,d[Fe/H]': S['obj2']['sigma_clean'],
        '[Y/Mg] dage (Gyr)': ymg['diff'],
        '[Y/Mg] sigma': ymg['sigma'],
        'Tc slope sigma (all)': sigma(*(np.array(list(slope['slopes_full'].get_group(c))) for c in ('small', 'giant'))),
        'Tc slope sigma (no s-proc)': sigma(*(np.array(list(slope['slopes_pure'].get_group(c))) for c in ('small', 'giant'))),
        'iso vs [Y/Mg] r': E['r_pearson'],
        'combined age sigma': E['comb_sig'],
    }


rows = []
print('paper (stored v6 table) ...', flush=True)
run([PIPE / 'run_pipeline.py'])
rows.append(dict(matching='v6 (paper)', description='stored table, procedure unknown', **summary(PIPE / 'output')))

run([MATCH / 'crossmatch_hosts.py'])          # needed by --drop-archive-hosts
for tag, (desc, opts) in MATCHINGS.items():
    print(f'{tag} ...', flush=True)
    run([MATCH / 'match_neighbours.py', *opts, '--tag', tag])
    nb = MATCH / 'output' / f'neighbours_{tag}.csv'
    if tag == 'v1' and V1.exists():
        same = pd.read_csv(nb).equals(pd.read_csv(V1))
        print(f'  regenerated v1 table identical to {V1.name}: {same}')
        if not same:
            sys.exit('the regenerated v1 neighbours differ from the committed table')
    run([PIPE / 'run_pipeline.py', '--neighbours', nb, '--tag', tag])
    if tag == 'v1' and V1.exists():
        table, stored = PIPE / 'output' / tag / 'final_test.csv', V1.with_name('final_test_v1.csv')
        diff = (pd.read_csv(table).drop(columns=['id', 'First', 'Second', 'Third'])
                - pd.read_csv(stored).drop(columns=['id', 'First', 'Second', 'Third'])).abs().to_numpy().max()
        print(f'  rebuilt v1 table vs {stored.name}: largest difference {diff:.1e}')
        if diff > 1e-9:
            sys.exit('the rebuilt v1 table differs from the committed final_test_v1.csv')
    rows.append(dict(matching=tag, description=desc, **summary(PIPE / 'output' / tag)))

res = pd.DataFrame(rows)
res.to_csv(OUT / 'matching_comparison.csv', index=False)
show = res.drop(columns='description').set_index('matching').T
with pd.option_context('display.float_format', '{:+.3f}'.format, 'display.width', 200):
    print('\n' + show.to_string())
print(f'\nwrote {OUT / "matching_comparison.csv"}; figures for each matching in {PIPE / "output"}/<matching>/figures')
