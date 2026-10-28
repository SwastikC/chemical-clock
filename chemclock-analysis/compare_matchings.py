"""Redo the paper's analysis for several neighbour tables and compare the results.

For each matching below: build the neighbour table with
chemclock-matching/match_neighbours.py, run the paper pipeline on it
(pipeline/run_pipeline.py --neighbours), and collect the main numbers. The stored
v6 table (the paper) is the first row. v2 (scaled distance, 10 neighbours) is the
recommended reproducible replacement for v6, v1 (rank-mean) the earlier proposal.
For both, the regenerated neighbours and host-minus-neighbour table are checked
against the committed neighbour_table_v*.csv and final_test_v*.csv in
chemclock-data/data/reference/.

Output: output/matching_comparison.csv, and a table on screen. Takes about 4 minutes.

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
REF = ROOT / 'chemclock-data' / 'data' / 'reference'
COMMITTED = ('v1', 'v2')               # tables committed in REF, checked to be regenerated exactly

# tag: (description, match_neighbours.py options)
MATCHINGS = {
    'v2':                 ('scaled distance (50 K, 0.10, 0.03), General Properties params, archive hosts removed, '
                           '10 neighbours (recommended v2)',
                           ['--method', 'scaled', '--k', '10', '--drop-archive-hosts']),
    'scaled_k3':          ('scaled distance, as v2 but 3 neighbours',
                           ['--method', 'scaled', '--k', '3', '--drop-archive-hosts']),
    'v1':                 ('rank-mean, HARPS params, archive hosts removed, with replacement (v1, superseded)',
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
    check = tag in COMMITTED and (REF / f'neighbour_table_{tag}.csv').exists()
    if check:
        same = pd.read_csv(nb).equals(pd.read_csv(REF / f'neighbour_table_{tag}.csv'))
        print(f'  regenerated {tag} neighbours identical to neighbour_table_{tag}.csv: {same}')
        if not same:
            sys.exit(f'the regenerated {tag} neighbours differ from the committed table')
    run([PIPE / 'run_pipeline.py', '--neighbours', nb, '--tag', tag])
    if check:
        new, old = (pd.read_csv(f) for f in (PIPE / 'output' / tag / 'final_test.csv', REF / f'final_test_{tag}.csv'))
        names = [c for c in old.columns if not c[0].isdigit()]          # id, First, Second, ...
        x, y = (t.drop(columns=names).to_numpy(float) for t in (new, old))
        same = (new.shape == old.shape and list(new.columns) == list(old.columns) and new[names].equals(old[names])
                and (np.isnan(x) == np.isnan(y)).all())          # missing values in the same cells
        diff = np.nanmax(np.abs(x - y)) if same else np.inf
        print(f'  rebuilt {tag} table vs final_test_{tag}.csv: same stars and missing values {same}, largest difference {diff:.1e}')
        if diff > 1e-9:
            sys.exit(f'the rebuilt {tag} table differs from the committed final_test_{tag}.csv')
    rows.append(dict(matching=tag, description=desc, **summary(PIPE / 'output' / tag)))

res = pd.DataFrame(rows)
res.to_csv(OUT / 'matching_comparison.csv', index=False)
show = res.drop(columns='description').set_index('matching').T
with pd.option_context('display.float_format', '{:+.3f}'.format, 'display.width', 200):
    print('\n' + show.to_string())
print(f'\nwrote {OUT / "matching_comparison.csv"}; figures for each matching in {PIPE / "output"}/<matching>/figures')
