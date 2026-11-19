#!/usr/bin/env python
"""ABC simulator WITH A MOLECULAR CLOCK (17 Nov 2026). Replaces the SNP-count-tuned mutation rate.

Why: with mu tuned to a fixed SNP count, absolute time was not identifiable; T tracked its prior floor
(1,000 -> 100 generations) and only T/(2Ne) was estimated. Here mu per year is drawn from the conifer
range and the simulated sequence is the real sequenced length, so per-site diversity (pi) and
divergence (d_xy) carry the information that anchors Ne and therefore T.

Reused unchanged from sim5: chain geometry, the four histories (build()), missingness mask, the 12
shape statistics. New: priors in years, generation time, 30,000 unlinked 86-bp loci, pi/d_xy stats.
"""
import sys, numpy as np, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim5 as B                       # build(), chains, observed missingness mask, deme layout
from sumstats2 import summarise, vec

NLOC, LLOC = 30_000, 86                 # observed mean locus length = 14.27 Mb / 166,105 loci
NSNP_CAP = 20_000                       # shape statistics computed on at most this many SNPs
IND_COAST = np.array([B.COASTAL[d] for d in B.SIM_DEME])   # individual k -> Coastal?

def draw(rng, model):
    mu_yr = 10**rng.uniform(-10.0, np.log10(2e-9))      # conifer range, per site per year
    # generation time (mean parent age), years. Lower bound: the 61-yr lower credible limit for the congener
    # A. angustifolia (Higuchi et al. 2026, For Ecol Manage, 95 yr [61-184]); A. araucana is slower-growing.
    # Upper bound 300 allows its longer reproductive life (lifespan >1,000 yr, Aguilera-Betti et al. 2017).
    g     = rng.uniform(60.0, 300.0)
    Tyr   = 10**rng.uniform(3.0, 7.0)                    # divergence, years
    p = dict(mu_yr=mu_yr, g=g, T_years=Tyr, mu=mu_yr*g, T1=Tyr/g,
             Ne=10**rng.uniform(1.5, 4.3), m=10**rng.uniform(-3.0, -1.2),
             Nanc=10**rng.uniform(2.5, 6.0))
    if model in ("H1", "H1b"):
        p["founder"] = 10**rng.uniform(-2.0, -0.3)
        p["Tb"]      = 10**rng.uniform(0.7, 2.3)
        p["Tcol"]    = p["T1"] * 10**rng.uniform(-2.0, -0.5)
    else:
        p["Tsep"] = p["T1"] * 10**rng.uniform(-0.7, 0.7)
    if model == "H3":
        p["mcross"]   = 10**rng.uniform(-5.0, -1.5)
        p["Tcontact"] = p["T1"] * 10**rng.uniform(-2.0, -0.1)
    return p

# loci separated by 1-bp spacers that recombine freely and carry no mutations
_pos = np.arange(NLOC+1)*(LLOC+1)
_bnd = np.sort(np.concatenate([_pos[:-1], _pos[:-1]+LLOC, [_pos[-1]]]))
_REC = msprime.RateMap(position=_bnd, rate=np.tile([0.0, 0.5], NLOC))
def _mutmap(mu): return msprime.RateMap(position=_bnd, rate=np.tile([mu, 0.0], NLOC))
_WIN = np.concatenate([_pos[:-1], [_pos[-1]]]).astype(float)

def clockstats(ts):
    nodes = np.array([ts.individual(i).nodes for i in range(ts.num_individuals)])
    C = nodes[IND_COAST].ravel(); A = nodes[~IND_COAST].ravel()
    piw = ts.diversity(sample_sets=[C, A], windows=_WIN, mode="site", span_normalise=False)  # per locus
    dxw = ts.divergence(sample_sets=[C, A], windows=_WIN, mode="site", span_normalise=False)
    seg = ts.segregating_sites(windows=_WIN, mode="site", span_normalise=False)
    def agg(sel):
        n = sel.sum()*LLOC
        return piw[sel,0].sum()/n, piw[sel,1].sum()/n, dxw[sel].sum()/n
    a = agg(np.ones(NLOC, bool)); b = agg(seg <= 5)
    return np.array(a + b)          # piC, piA, dxy (all loci), then the same on loci with <= 5 SNPs

CLOCK_NAMES = ["piC","piA","dxy","piC_le5","piA_le5","dxy_le5"]

def one(model, seed):
    rng = np.random.default_rng(seed)
    p = draw(rng, model)
    try:
        dem = B.build(model, p)
        samples = {f"d{B.STAND_DEME[j]}": int(B.NPER[j]) for j in range(B.ND)}
        anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=_bnd[-1],
                                   recombination_rate=_REC, random_seed=seed, ploidy=2)
        ts = msprime.sim_mutations(anc, rate=_mutmap(p["mu"]), random_seed=seed+1,
                                   model=msprime.BinaryMutationModel())
        if ts.num_sites < 1000: return None
        cs = clockstats(ts)
        Gm = ts.genotype_matrix()
        Gd = (Gm[:, 0::2] + Gm[:, 1::2]).T.astype(float)
        if Gd.shape[1] > NSNP_CAP:
            Gd = Gd[:, np.linspace(0, Gd.shape[1]-1, NSNP_CAP).astype(int)]
        r2 = np.random.default_rng(seed+2)
        cols = r2.choice(B.MASK.shape[1], size=Gd.shape[1], replace=True)
        Gd[B.MASK[:, cols]] = np.nan
        s = summarise(Gd, B.SIM_DEME, B.COORDS, B.COASTAL, nd=B.ND)
        if s is None: return None
        v = np.concatenate([vec(s), np.log10(cs)])
        if not np.all(np.isfinite(v)): return None
        p["nsnp"] = ts.num_sites
        return model, v, p
    except Exception as e:
        one.last_error = repr(e); return None
