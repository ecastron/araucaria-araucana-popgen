"""Timing and noise of clock statistics: branch-mode pi/dxy on the sim5 linked sequence vs unlinked loci."""
import sys, time, numpy as np, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim5 as B
IND_COAST = np.array([B.COASTAL[d] for d in B.SIM_DEME])
samples = {f"d{B.STAND_DEME[j]}": int(B.NPER[j]) for j in range(B.ND)}
CASES = {"typical": dict(T1=1000, Tsep=1000, Ne=1000, m=0.01, Nanc=1e4),
         "hard":    dict(T1=1e4,  Tsep=1e4,  Ne=1e4,  m=0.01, Nanc=1e5)}
def bstats(ts):
    nodes = np.array([ts.individual(i).nodes for i in range(ts.num_individuals)])
    C = nodes[IND_COAST].ravel(); A = nodes[~IND_COAST].ravel()
    pi = ts.diversity([C, A], mode="branch"); dx = ts.divergence([C, A], mode="branch")
    return pi[0], pi[1], dx
for case, p in CASES.items():
    dem = B.build("H2", p)
    for label, L, r in (("linked 1Mb r1e-7", 1_000_000, 1e-7), ("linked 1Mb r1e-6", 1_000_000, 1e-6)):
        out = []
        for s in range(3):
            t = time.time()
            ts = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L,
                                      recombination_rate=r, random_seed=100+s, ploidy=2)
            out.append(bstats(ts) + (time.time()-t, ts.num_trees))
        o = np.array(out)
        print(f"{case:8s} {label:18s} time {o[:,3].mean():6.0f}s trees {o[:,4].mean():8.0f}  "
              f"branch dxy mean {o[:,2].mean():.0f} CV {o[:,2].std()/o[:,2].mean():.3f}  "
              f"dxy/piA {np.mean(o[:,2]/o[:,1]):.3f}", flush=True)
    for nl in (200, 1000):
        t = time.time(); vals = []
        for ts in msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=1,
                                       random_seed=7, ploidy=2, num_replicates=nl):
            vals.append(bstats(ts))
        v = np.array(vals)
        print(f"{case:8s} unlinked {nl:5d} loci    time {time.time()-t:6.0f}s  branch dxy mean "
              f"{v[:,2].mean():.0f} SE/mean {v[:,2].std()/np.sqrt(nl)/v[:,2].mean():.3f}", flush=True)
