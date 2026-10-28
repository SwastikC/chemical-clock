"""Per-SwP stellar-parameter offsets (host minus mean of its neighbours).

Parameters are taken from the same HARPS-GTO columns used by the rank-mean
neighbour matching (TeffK, logg[cm/s2], [Fe/H]---_x in starswithplanet.csv /
starswithoutplanet.csv). Hosts are matched on either their HD identifier
(Star---) or their common name (hostname_1). The 12 hosts absent from
starswithplanet.csv use Properties_f.csv for host and neighbours alike, so each
pair is drawn from a single parameter source.

The v5 version matched hosts only on hostname_1, which silently dropped 19 of
133 hosts (common-name mismatches such as tau Cet, eps Eri, 61 Vir, plus hosts
absent from starswithplanet.csv); 12 of them were in the 121-star classified
sample. That subset alone raised the small-vs-giant s-process significance
from 4.40 to 4.75 sigma before any regression.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import FILES, OUT
import numpy as np, pandas as pd

norm=lambda s: str(s).replace(' ','').strip()
cols=['TeffK','logg[cm/s2]','[Fe/H]---_x']
nb=pd.read_csv(FILES['neighbours']); prop=pd.read_csv(FILES['properties'])
swp=pd.read_csv(FILES['swp']); swop=pd.read_csv(FILES['swop'])
H={}
for c in ['Star---','hostname_1']:
    for _,r in swp.dropna(subset=cols).iterrows(): H.setdefault(norm(r[c]),tuple(r[cols]))
swop['n']=swop['Star---'].apply(norm)
N={r['n']:tuple(r[cols]) for _,r in swop.dropna(subset=cols).drop_duplicates('n').iterrows()}
P={norm(r['Star']):(r['Teff (K)'],r['logg [cm/s2]'],r['[Fe/H]']) for _,r in prop.iterrows()}
rec=[]
for _,r in nb.iterrows():
    s=norm(r['SwP']); nbs=[norm(r[c]) for c in nb.columns if c!='SwP']   # any number of neighbours
    if s in H and all(n in N for n in nbs): h,a=H[s],[N[n] for n in nbs]; src='harps'
    elif s in P and all(n in P for n in nbs): h,a=P[s],[P[n] for n in nbs]; src='properties_f'
    else: continue
    h=np.array(h,float); a=np.array(a,float)
    rec.append(dict(SwP=s,dTeff=h[0]-a[:,0].mean(),dlogg=h[1]-a[:,1].mean(),dFeH=h[2]-a[:,2].mean(),source=src))
pd.DataFrame(rec).to_csv(OUT+'/paper_stellar_offsets.csv',index=False)
print(len(rec),'rows written')
