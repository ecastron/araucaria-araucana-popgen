#!/usr/bin/env python
"""One Slurm array task of the clock ABC: draw from the prior, simulate, summarise, save."""
import sys, time, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim9 as S
from sumstats2 import STAT_ORDER

task, nper, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
MODELS = tuple(sys.argv[4].split(","))
rows, labs, pars, fails = [], [], [], 0
t0 = time.time()
for mi, mod in enumerate(MODELS):
    for k in range(nper):
        seed = 11_000_000 + task*100_000 + mi*10_000 + k     # new seed base, disjoint from sim5-8
        r = S.one(mod, seed)
        if r is None: fails += 1; continue
        _, v, p = r
        rows.append(v); labs.append(mi); pars.append(p)
        if (k+1) % 5 == 0:
            print(f"task {task} {mod} {k+1}/{nper}  {time.time()-t0:.0f}s  fails {fails}", flush=True)
X = np.array(rows); y = np.array(labs)
keys = sorted({k for p in pars for k in p})
Pm = np.array([[p.get(k, np.nan) for k in keys] for p in pars])
np.savez_compressed(out, X=X, y=y, P=Pm, param_names=np.array(keys),
                    stat_order=np.array(list(STAT_ORDER) + S.CLOCK_NAMES), models=np.array(MODELS))
print(f"task {task}: saved {X.shape[0]} sims ({fails} failed) in {time.time()-t0:.0f}s -> {out}", flush=True)
