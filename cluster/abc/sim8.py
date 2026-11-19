#!/usr/bin/env python
"""ABC simulator WITH A MOLECULAR CLOCK (17 Nov 2026). Replaces the SNP-count-tuned mutation rate.

Why: with mu tuned to a fixed SNP count, absolute time was not identifiable; T tracked its prior floor
(1,000 -> 100 generations) and only T/(2Ne) was estimated. Here mu per year and generation time are
drawn from literature priors, and three clock statistics (pi Coastal, pi Andean, d_xy per site) are
added to the 12 scale-invariant shape statistics. Per-site diversity anchors Ne and therefore T.

Clock statistics are computed in branch mode on the simulated genealogy (expected pi = mu_gen x mean
pairwise branch length), with no mutation noise. Benchmark (job 45127): on the sim5 1 Mb sequence,
seed-to-seed CV of d_xy is 0.6 to 0.8%, comparable to the observed precision, at 1 to 96 s per
simulation. A first design with 30,000 unlinked loci did not finish one simulation in 90 min.

Reused unchanged from sim5: chain geometry, the four histories (build()), missingness mask, the 12
shape statistics, and the SNP-count-tuned mutation step that feeds them (shape statistics are
invariant to the mutation scale). The previous unlinked-loci version is sim8_unlinked.py.
"""
import sys, numpy as np, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim5 as B
from sumstats2 import summarise, vec

IND_COAST = np.array([B.COASTAL[d] for d in B.SIM_DEME])   # individual k -> Coastal?
CLOCK_NAMES = ["log_piC", "log_piA", "log_dxy"]

def draw(rng, model):
    mu_yr = 10**rng.uniform(-10.0, np.log10(2e-9))      # conifer range, per site per year
    # generation time (mean parent age), years. Lower bound: the 61-yr lower credible limit for the congener
    # A. angustifolia (Higuchi et al. 2026, For Ecol Manage, 95 yr [61-184]); A. araucana is slower-growing.
    # Upper bound 300 allows its longer reproductive life (lifespan >1,000 yr, Aguilera-Betti et al. 2017).
    g     = rng.uniform(60.0, 300.0)
    Tyr   = 10**rng.uniform(3.0, 7.0)                    # divergence, years
    p = dict(mu_yr=mu_yr, g=g, T_years=Tyr, mu=mu_yr*g, T1=Tyr/g,
             # Ne, m and Nanc priors identical to sim5/sim6. A wider version (Ne 10^1.5-10^4.3, Nanc up to
             # 10^6) gave pi up to 0.46 against 0.0028 observed and 20-45 min per simulation (job 45128).
             Ne=10**rng.uniform(2.3, 3.7), m=10**rng.uniform(-3.0, -1.2),
             Nanc=10**rng.uniform(3.3, 5.2))
    if model in ("H1", "H1b"):
        p["founder"] = 10**rng.uniform(-2.0, -0.3)
        p["Tb"]      = 10**rng.uniform(0.7, 2.3)
        p["Tcol"]    = p["T1"] * 10**rng.uniform(-2.0, -0.5)
        # with T in years, T1 can be a few generations; the bottleneck must end before the range split,
        # otherwise msprime reactivates a split population (LibraryError, 4 of 8 H1/H1b draws, job 45128)
        p["Tb"]      = min(p["Tb"], 0.5*(p["T1"] - p["Tcol"]))
    else:
        p["Tsep"] = p["T1"] * 10**rng.uniform(-0.7, 0.7)
    if model == "H3":
        p["mcross"]   = 10**rng.uniform(-5.0, -1.5)
        p["Tcontact"] = p["T1"] * 10**rng.uniform(-2.0, -0.1)
    return p

def clockstats(anc, mu):
    """Expected per-site pi within each range and d_xy between them, from branch lengths."""
    nodes = np.array([anc.individual(i).nodes for i in range(anc.num_individuals)])
    C = nodes[IND_COAST].ravel(); A = nodes[~IND_COAST].ravel()
    pi = anc.diversity([C, A], mode="branch"); dx = anc.divergence([C, A], mode="branch")
    return mu*np.array([pi[0], pi[1], dx])

def one(model, seed, L_bp=1_000_000, r=1e-7, target=B.NSNP):
    rng = np.random.default_rng(seed)
    p = draw(rng, model)
    try:
        dem = B.build(model, p)
        samples = {f"d{B.STAND_DEME[j]}": int(B.NPER[j]) for j in range(B.ND)}
        anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L_bp,
                                   recombination_rate=r, random_seed=seed, ploidy=2)
        cs = clockstats(anc, p["mu"])
        # shape statistics: identical to sim5.simulate (mutation rate tuned to the SNP target)
        ts = None; mu = 1e-8
        for k in range(8):
            ts = msprime.sim_mutations(anc, rate=mu, random_seed=seed+1+k,
                                       model=msprime.BinaryMutationModel())
            if ts.num_sites >= target or mu >= B.MU_CAP: break
            mu = min(mu*max(2.0, 1.3*target/max(ts.num_sites, 1)), B.MU_CAP)
        if ts.num_sites < 1000: return None
        Gm = ts.genotype_matrix()
        Gd = (Gm[:, 0::2] + Gm[:, 1::2]).T.astype(float)
        if Gd.shape[1] > target:
            Gd = Gd[:, np.linspace(0, Gd.shape[1]-1, target).astype(int)]
        r2 = np.random.default_rng(seed+2)
        cols = r2.choice(B.MASK.shape[1], size=Gd.shape[1], replace=True)
        Gd[B.MASK[:, cols]] = np.nan
        s = summarise(Gd, B.SIM_DEME, B.COORDS, B.COASTAL, nd=B.ND)
        if s is None: return None
        v = np.concatenate([vec(s), np.log10(cs)])
        if not np.all(np.isfinite(v)): return None
        p["tau"] = p["T1"]/(2*p["Ne"])
        return model, v, p
    except Exception as e:
        one.last_error = repr(e); return None
