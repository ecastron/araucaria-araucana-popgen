#!/usr/bin/env python
"""Counts plotted in main Figure 3B: variants private to one range versus shared.

Deliberately a DIFFERENT quantity from the `private_ratio` summary statistic used for model choice
(which uses >0.05 / <0.01 thresholds and is a ratio). Here a variant is 'present' in a range at
>1% mean frequency and 'absent' at <=1%, and raw counts are reported, because a count is
interpretable to a reader in a way a ratio is not. Both are reported so the difference is explicit.
"""
import sys, numpy as np
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from load_observed import load
from sumstats import freqs_by_deme
G,DEME,COORDS,COASTAL,stands=load()
P,N=freqs_by_deme(G,DEME,12)
co=[i for i in range(12) if COASTAL[i]]; an=[i for i in range(12) if not COASTAL[i]]
pc=np.nanmean(P[co],axis=0); pa=np.nanmean(P[an],axis=0)
nc=np.nansum(N[co],axis=0);  na=np.nansum(N[an],axis=0)
ok=np.isfinite(pc)&np.isfinite(pa)&(nc>=20)&(na>=20)
polyC=(pc>0.01)&(pc<0.99); polyA=(pa>0.01)&(pa<0.99)
onlyC=int((ok&polyC&(~polyA)&(pa<=0.01)).sum())
onlyA=int((ok&polyA&(~polyC)&(pc<=0.01)).sum())
both =int((ok&polyC&polyA).sum())
print(f"usable SNPs (>=20 chromosomes scored in both ranges): {int(ok.sum()):,} of {len(pc):,}")
print(f"  only in Coastal : {onlyC:,}")
print(f"  in both ranges  : {both:,}")
print(f"  only in Andean  : {onlyA:,}")
print(f"  total plotted   : {onlyC+both+onlyA:,}")
print(f"\n  count ratio Coastal-only : Andean-only = {onlyC/onlyA:.3f}")
print(f"  (for contrast, the private_ratio summary statistic, >0.05/<0.01 thresholds = 0.7738)")
