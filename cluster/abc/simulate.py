#!/usr/bin/env python
"""msprime simulator for the three Araucaria hypotheses.

H1' single refugium, postglacial colonisation : one ancestral source; the Coastal range is founded
     from the Andean source at T_col with a founder bottleneck. No cross-gap migration afterwards.
H2' two refugia, persistent isolation          : Coastal and Andean ancestors split at T_split and
     never exchange migrants. Distance-decay migration WITHIN each range only.
H3' two refugia with secondary contact         : as H2', plus symmetric cross-gap migration
     beginning at T_contact (< T_split).

All models share: 12 demes at the real coordinates, observed per-deme sample sizes, distance-decay
migration within range. Simulated genotypes receive the OBSERVED missingness pattern before the
summary statistics are computed, so the observed point is comparable to the simulated cloud.
"""
import sys, pickle, numpy as np, msprime
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
from sumstats import summarise, vec

D = pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))
COORDS, COASTAL, DEME = D["coords"], D["coastal"], D["deme"]
ND = len(COORDS); NIND = D["n_ind"]
NPER = np.array([(DEME==i).sum() for i in range(ND)])
MISS_IND, MISS_LOC = D["miss_ind"], D["miss_loc"]
NSNP = 30000

def geo_km(a,b):
    R=6371.0
    la1,lo1=np.radians(COORDS[a,1]),np.radians(COORDS[a,0])
    la2,lo2=np.radians(COORDS[b,1]),np.radians(COORDS[b,0])
    h=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*R*np.arcsin(np.sqrt(h))
KM = np.array([[geo_km(a,b) for b in range(ND)] for a in range(ND)])

def draw(rng, model):
    p = dict(
        Ne     = 10**rng.uniform(2.3, 3.9),        # 200 - 8000 per deme
        Nanc   = 10**rng.uniform(2.6, 4.2),        # ancestral
        m0     = 10**rng.uniform(-4.0, -1.3),      # migration at zero distance
        dscale = 10**rng.uniform(1.0, 2.3),        # decay scale, km
        T1     = 10**rng.uniform(2.3, 4.3),        # split / colonisation, generations
    )
    if model == "H1":
        p["founder"] = 10**rng.uniform(-1.7, -0.1)   # founder size as a fraction of Ne
        p["Tb"]      = 10**rng.uniform(0.7, 2.3)     # founder bottleneck duration, 5-200 gens
    if model == "H3":
        p["mcross"] = 10**rng.uniform(-4.5, -1.7)
        p["Tcontact"] = p["T1"] * rng.uniform(0.05, 0.7)
    return p

def build(model, p):
    dem = msprime.Demography()
    for i in range(ND): dem.add_population(name=f"d{i}", initial_size=p["Ne"])
    within = 0
    for a in range(ND):
        for b in range(ND):
            if a == b: continue
            same = (COASTAL[a] == COASTAL[b])
            if not same: continue                          # within-range only, all models
            m = p["m0"]*np.exp(-KM[a,b]/p["dscale"])
            if m > 1e-9: dem.set_migration_rate(source=f"d{a}", dest=f"d{b}", rate=m); within += 1
    co = [i for i in range(ND) if COASTAL[i]]
    an = [i for i in range(ND) if not COASTAL[i]]
    if model == "H1":
        dem.add_population(name="ANC", initial_size=p["Nanc"])
        dem.add_population(name="SRC", initial_size=p["Ne"])
        # Backwards in time: at T1 the Coastal demes contract to the founder size and stay there
        # for Tb generations (the founder event), then coalesce into the Andean source SRC.
        for i in co:
            dem.add_population_parameters_change(time=p["T1"], population=f"d{i}",
                                                 initial_size=p["Ne"]*p["founder"])
        # cut Coastal off from each other during the bottleneck so it is a true founder event
        dem.add_migration_rate_change(time=p["T1"], rate=0)
        dem.add_population_split(time=p["T1"]+p["Tb"], derived=[f"d{i}" for i in co], ancestral="SRC")
        dem.add_population_split(time=p["T1"]+p["Tb"]+1, derived=[f"d{i}" for i in an]+["SRC"], ancestral="ANC")
    else:
        dem.add_population(name="ANC_C", initial_size=p["Nanc"])
        dem.add_population(name="ANC_A", initial_size=p["Nanc"])
        dem.add_population(name="ANC",   initial_size=p["Nanc"])
        dem.add_population_split(time=p["T1"]*0.5, derived=[f"d{i}" for i in co], ancestral="ANC_C")
        dem.add_population_split(time=p["T1"]*0.5, derived=[f"d{i}" for i in an], ancestral="ANC_A")
        dem.add_population_split(time=p["T1"], derived=["ANC_C","ANC_A"], ancestral="ANC")
        if model == "H3":
            for a in co:
                for b in an:
                    dem.set_migration_rate(source=f"d{a}", dest=f"d{b}", rate=p["mcross"])
                    dem.set_migration_rate(source=f"d{b}", dest=f"d{a}", rate=p["mcross"])
            # contact begins at Tcontact: switch cross-gap migration off before then
            for a in co:
                for b in an:
                    dem.add_migration_rate_change(time=p["Tcontact"], source=f"d{a}", dest=f"d{b}", rate=0)
                    dem.add_migration_rate_change(time=p["Tcontact"], source=f"d{b}", dest=f"d{a}", rate=0)
    dem.sort_events()
    return dem

def simulate(model, p, seed, L=1_000_000, r=1e-7, mu=1e-7, target=NSNP):
    """Ancestry is simulated once; the mutation rate is then tuned so every simulation delivers the
    SAME number of SNPs (target). The SNP count is an ascertainment property of the RAD protocol,
    not a prediction of the demographic model, and holding it fixed means FST estimation noise is
    identical across simulations and comparable to the observed matrix. Without this, low-diversity
    models yield few SNPs, noisier FST, and an attenuated IBD correlation -- an artefact that would
    be read as evidence about demography."""
    dem = build(model, p)
    samples = {f"d{i}": int(NPER[i]) for i in range(ND)}
    anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L,
                               recombination_rate=r, random_seed=seed, ploidy=2)
    ts = None; mu_try = mu
    for k in range(6):
        ts = msprime.sim_mutations(anc, rate=mu_try, random_seed=seed+1+k,
                                   model=msprime.BinaryMutationModel())
        if ts.num_sites >= target: break
        mu_try *= max(2.0, 1.3*target/max(ts.num_sites, 1))
    simulate.last_nsites = ts.num_sites
    simulate.last_mu = mu_try
    if ts.num_sites < 200: return None
    G = ts.genotype_matrix()                       # sites x haplotypes (0/1)
    Gd = (G[:, 0::2] + G[:, 1::2]).T.astype(float) # individuals x sites, dosage
    if Gd.shape[1] > target:
        idx = np.linspace(0, Gd.shape[1]-1, target).astype(int)
        Gd = Gd[:, idx]
    # apply the observed missingness structure
    rng = np.random.default_rng(seed+2)
    mi = rng.permutation(MISS_IND)[:Gd.shape[0]][:,None]
    ml = rng.choice(MISS_LOC, size=Gd.shape[1])[None,:]
    pmiss = np.clip(mi + ml - mi*ml, 0, 0.995)
    Gd[rng.random(Gd.shape) < pmiss] = np.nan
    return Gd

def one(model, seed):
    rng = np.random.default_rng(seed)
    p = draw(rng, model)
    try:
        Gd = simulate(model, p, seed)
        if Gd is None: return None
        s = summarise(Gd, DEME, COORDS, COASTAL, nd=ND)
        if s is None: return None
        return model, vec(s), p
    except Exception:
        return None
