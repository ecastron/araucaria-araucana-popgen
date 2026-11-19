#!/usr/bin/env python
"""One Slurm array task: draw from the prior, simulate, summarise, save."""
import sys, os, time, numpy as np
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim7 as S
from sumstats2 import STAT_ORDER

task  = int(sys.argv[1])          # array index
nper  = int(sys.argv[2])          # simulations per model per task
out   = sys.argv[3]
MODELS = tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else ("H1","H2","H3")
rows, labs, pars = [], [], []
t0 = time.time()
for mi, mod in enumerate(MODELS):
    for k in range(nper):
        seed = 9_000_000 + task*100_000 + mi*10_000 + k
        r = S.one(mod, seed)
        if r is None: continue
        _, v, p = r
        if not np.all(np.isfinite(v)): continue
        rows.append(v); labs.append(mi); pars.append(p)
        if (k+1) % 5 == 0:
            print(f"task {task} {mod} {k+1}/{nper}  {time.time()-t0:.0f}s", flush=True)
X = np.array(rows); y = np.array(labs)
keys = sorted({k for p in pars for k in p})
Pm = np.array([[p.get(k, np.nan) for k in keys] for p in pars])
np.savez_compressed(out, X=X, y=y, P=Pm, param_names=np.array(keys),
                    stat_order=np.array(STAT_ORDER), models=np.array(MODELS))
print(f"task {task}: saved {X.shape[0]} sims in {time.time()-t0:.0f}s -> {out}", flush=True)
