#!/usr/bin/env python
"""The three Araucaria hypotheses on 1D stepping-stone chains.

H1' Single refugium, postglacial colonisation. Andean chain is ancestral; the Coastal chain is
    founded at T_col through a bottleneck (Ne*founder, lasting Tb) and coalesces into the Andean
    deme nearest the gap. Coastal ancestry is NESTED within Andean. No cross-gap migration.
H2' Two refugia, persistent isolation. Both chains persist to T_div, then collapse into refugia
    ANC_C / ANC_A which merge only at T_div+Tsep. No cross-gap migration, ever.
H3' Two refugia with secondary contact. As H2', plus symmetric migration through the gap corridor
    (nearest Andean<->Coastal deme pair), active from T_contact to the present.

Mutation rate is tuned to deliver a fixed SNP count but CAPPED: an uncapped search drove mu to
~3.7e-5/bp in low-diversity models, saturating sites with recurrent mutation.
Missingness resamples columns of the real observed mask, preserving its joint structure.
"""
import sys, pickle, numpy as np, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
from sumstats2 import summarise, vec
from chain import build_chains

P = "/home/ecastron/araucaria"
D = pickle.load(open(f"{P}/results/Q1/abc/observed2.pkl", "rb"))
DEME, COORDS, COASTAL = D["deme"], D["coords"], D["coastal"]
ND = len(COORDS); NPER = np.array([(DEME == i).sum() for i in range(ND)])

# CRITICAL: the observed DEME vector is in VCF order (interleaved by stand), but msprime returns
# simulated individuals GROUPED by population in dict order. Passing the observed DEME to
# summarise() on simulated data assigns simulated individuals to essentially random stands
# (only 9.2% of positions coincide), which destroys isolation by distance and inflates FST
# scatter regardless of the demography. Simulations must use the grouped vector, and the
# missingness mask rows must be reordered to match so each simulated individual inherits the
# missingness of a REAL individual from the SAME stand.
SIM_DEME = np.repeat(np.arange(ND), NPER)
_ORD = np.argsort(DEME, kind="stable")
MASK = D["miss_mask"][_ORD]
assert np.array_equal(DEME[_ORD], SIM_DEME), "mask row order does not match grouped deme order"
NSNP = 20000
MU_CAP = 2e-6

L = build_chains()
K, NA, NC = L["K"], L["nA"], L["nC"]
EDGES, STAND_DEME, NODE_C = L["edges"], L["stand_deme"], L["node_coastal"]
AN = [i for i in range(K) if not NODE_C[i]]
CO = [i for i in range(K) if NODE_C[i]]
GAP_A, GAP_C = L["gap_andean_deme"], L["gap_coastal_deme"]

def draw(rng, model):
    p = dict(Ne   = 10**rng.uniform(2.3, 3.7),
             m    = 10**rng.uniform(-3.0, -1.2),
             T1   = 10**rng.uniform(1.0, 4.8),        # FLOOR 10 gen (17 Nov 2026): control for the floor-100 lower tail
             Nanc = 10**rng.uniform(3.3, 5.2))
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

def split_many(dem, time, derived, ancestral, size, tag):
    if len(derived) <= 90:
        dem.add_population_split(time=time, derived=list(derived), ancestral=ancestral); return
    chunks = [list(derived[i:i+90]) for i in range(0, len(derived), 90)]
    hubs = []
    for ci, ch in enumerate(chunks):
        h = f"{tag}_h{ci}"; dem.add_population(name=h, initial_size=size)
        dem.add_population_split(time=time, derived=ch, ancestral=h); hubs.append(h)
    dem.add_population_split(time=time+1e-3, derived=hubs, ancestral=ancestral)

def build(model, p):
    dem = msprime.Demography()
    for i in range(K): dem.add_population(name=f"d{i}", initial_size=p["Ne"])
    for a, b in EDGES:
        dem.set_migration_rate(source=f"d{a}", dest=f"d{b}", rate=p["m"])
        dem.set_migration_rate(source=f"d{b}", dest=f"d{a}", rate=p["m"])
    co = [f"d{i}" for i in CO]; an = [f"d{i}" for i in AN]
    if model == "H1":
        dem.add_population(name="FOUNDER", initial_size=max(p["Ne"]*p["founder"], 1.0))
        dem.add_population(name="ANC", initial_size=p["Nanc"])
        # The Coastal chain splits away at Tcol while the Andean chain stays active until T1.
        # Its migration entries must be zeroed at Tcol, otherwise msprime tries to move lineages
        # into populations that are no longer active (LibraryError, 79.5% of H1 draws failed).
        for a, b in EDGES:
            if NODE_C[a] and NODE_C[b]:
                dem.add_migration_rate_change(time=p["Tcol"], source=f"d{a}", dest=f"d{b}", rate=0)
                dem.add_migration_rate_change(time=p["Tcol"], source=f"d{b}", dest=f"d{a}", rate=0)
        split_many(dem, p["Tcol"], co, "FOUNDER", max(p["Ne"]*p["founder"],1.0), "co")
        dem.add_population_split(time=p["Tcol"]+p["Tb"], derived=["FOUNDER"], ancestral=f"d{GAP_A}")
        # msprime marks any population used as a split ANCESTOR as initially_active=False.
        # d{GAP_A} is a sampled Andean deme that must be active from the present, so the flag has
        # to be set back explicitly -- otherwise every lineage entering it raises
        # "Attempt to move a lineage into an inactive population" (79.5% of H1 draws failed).
        dem[f"d{GAP_A}"].initially_active = True
        split_many(dem, p["T1"], an, "ANC", p["Nanc"], "an")
    elif model == "H1b":
        # mirror of H1: the COASTAL range is the single refugium and the Andean range is founded
        # from it. Tested because H1' as stated allows nesting in either direction.
        dem.add_population(name="FOUNDER", initial_size=max(p["Ne"]*p["founder"], 1.0))
        dem.add_population(name="ANC", initial_size=p["Nanc"])
        for a, b in EDGES:
            if (not NODE_C[a]) and (not NODE_C[b]):
                dem.add_migration_rate_change(time=p["Tcol"], source=f"d{a}", dest=f"d{b}", rate=0)
                dem.add_migration_rate_change(time=p["Tcol"], source=f"d{b}", dest=f"d{a}", rate=0)
        split_many(dem, p["Tcol"], an, "FOUNDER", max(p["Ne"]*p["founder"],1.0), "an")
        dem.add_population_split(time=p["Tcol"]+p["Tb"], derived=["FOUNDER"], ancestral=f"d{GAP_C}")
        dem[f"d{GAP_C}"].initially_active = True
        split_many(dem, p["T1"], co, "ANC", p["Nanc"], "co")
    else:
        dem.add_population(name="ANC_C", initial_size=p["Nanc"])
        dem.add_population(name="ANC_A", initial_size=p["Nanc"])
        dem.add_population(name="ANC",   initial_size=p["Nanc"])
        if model == "H3":
            dem.set_migration_rate(source=f"d{GAP_A}", dest=f"d{GAP_C}", rate=p["mcross"])
            dem.set_migration_rate(source=f"d{GAP_C}", dest=f"d{GAP_A}", rate=p["mcross"])
            dem.add_migration_rate_change(time=p["Tcontact"], source=f"d{GAP_A}", dest=f"d{GAP_C}", rate=0)
            dem.add_migration_rate_change(time=p["Tcontact"], source=f"d{GAP_C}", dest=f"d{GAP_A}", rate=0)
        split_many(dem, p["T1"], co, "ANC_C", p["Nanc"], "co")
        split_many(dem, p["T1"], an, "ANC_A", p["Nanc"], "an")
        dem.add_population_split(time=p["T1"]+p["Tsep"]+0.01, derived=["ANC_C","ANC_A"], ancestral="ANC")
    dem.sort_events()
    return dem

def simulate(model, p, seed, L_bp=1_000_000, r=1e-7, mu0=1e-8, target=NSNP, apply_missing=True):
    dem = build(model, p)
    samples = {f"d{STAND_DEME[j]}": int(NPER[j]) for j in range(ND)}
    anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L_bp,
                               recombination_rate=r, random_seed=seed, ploidy=2)
    ts = None; mu = mu0
    for k in range(8):
        ts = msprime.sim_mutations(anc, rate=mu, random_seed=seed+1+k,
                                   model=msprime.BinaryMutationModel())
        if ts.num_sites >= target or mu >= MU_CAP: break
        mu = min(mu*max(2.0, 1.3*target/max(ts.num_sites, 1)), MU_CAP)
    simulate.last_nsites = ts.num_sites; simulate.last_mu = mu
    if ts.num_sites < 1000: return None
    Gm = ts.genotype_matrix()
    Gd = (Gm[:, 0::2] + Gm[:, 1::2]).T.astype(float)
    if Gd.shape[1] > target:
        Gd = Gd[:, np.linspace(0, Gd.shape[1]-1, target).astype(int)]
    if apply_missing:
        rng = np.random.default_rng(seed+2)
        cols = rng.choice(MASK.shape[1], size=Gd.shape[1], replace=True)
        Gd[MASK[:, cols]] = np.nan
    return Gd

def one(model, seed):
    rng = np.random.default_rng(seed)
    p = draw(rng, model)
    try:
        Gd = simulate(model, p, seed)
        if Gd is None: return None
        s = summarise(Gd, SIM_DEME, COORDS, COASTAL, nd=ND)
        if s is None: return None
        v = vec(s)
        if not np.all(np.isfinite(v)): return None
        return model, v, p
    except Exception:
        return None
