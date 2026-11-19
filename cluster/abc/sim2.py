#!/usr/bin/env python
"""msprime simulator for the three Araucaria hypotheses on a stepping-stone lattice.

H1' Single refugium, postglacial colonisation. The Andean lattice is ancestral; the Coastal
    lattice is founded at T_col through a bottleneck of size Ne*founder lasting Tb generations,
    then coalesces into the Andean node nearest the gap. Coastal ancestry is NESTED in Andean.
    No cross-gap migration.
H2' Two refugia, persistent isolation. Both lattices persist to T_div, then collapse into their
    refugia ANC_C / ANC_A, which merge only at T_div+Tsep. No cross-gap migration, ever.
H3' Two refugia with secondary contact. As H2', plus symmetric migration across the closest
    Coastal-Andean node pairs, active from T_contact to the present.

Missingness is imposed by resampling COLUMNS of the real observed mask, which reproduces its joint
individual x locus structure. A parametric mask built from the two marginals double-counts them and
delivered 85.7% missing against an observed 62.1%.
"""
import sys, pickle, numpy as np, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
from sumstats import summarise, vec
from lattice import build_lattice

P = "/home/ecastron/araucaria"
D = pickle.load(open(f"{P}/results/Q1/abc/observed.pkl", "rb"))
DEME, COORDS, COASTAL = D["deme"], D["coords"], D["coastal"]
ND = len(COORDS); NPER = np.array([(DEME == i).sum() for i in range(ND)])
MASK = D["miss_mask"]                      # 163 x 78026 bool, the real thing
NSNP = 30000

L = build_lattice("tri8")
GK, EDGES, NODE_C, STAND_NODE = L["coords"], L["edges"], L["node_coastal"], L["stand_node"]
K = GK.shape[0]
CO_N = np.where(NODE_C)[0]; AN_N = np.where(~NODE_C)[0]

def _km(a, b):
    R = 6371.0
    la1, lo1 = np.radians(a[1]), np.radians(a[0]); la2, lo2 = np.radians(b[1]), np.radians(b[0])
    h = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*R*np.arcsin(np.sqrt(h))

# the closest cross-gap node pairs, used only by H3' as the secondary-contact corridor
_gap = sorted(((_km(GK[a], GK[b]), a, b) for a in CO_N for b in AN_N), key=lambda t: t[0])
GAP_PAIRS = [(a, b) for _, a, b in _gap[:6]]
GAP_KM = _gap[0][0]
# the Andean node nearest the gap is the H1' colonisation source
SRC_NODE = GAP_PAIRS[0][1]

def draw(rng, model):
    p = dict(Ne     = 10**rng.uniform(1.3, 3.3),      # 20 - 2000 per lattice node
             Nanc   = 10**rng.uniform(2.3, 4.3),
             m      = 10**rng.uniform(-3.5, -0.7),    # per-edge migration
             T1     = 10**rng.uniform(2.3, 4.3))      # T_div (H2/H3) or T_anc (H1)
    if model == "H1":
        p["founder"] = 10**rng.uniform(-2.0, -0.3)
        p["Tb"]      = 10**rng.uniform(0.7, 2.3)
        p["Tcol"]    = p["T1"] * 10**rng.uniform(-1.7, -0.3)   # colonisation is recent
    else:
        p["Tsep"] = p["T1"] * 10**rng.uniform(-0.7, 0.7)
    if model == "H3":
        p["mcross"]   = 10**rng.uniform(-5.0, -1.5)
        p["Tcontact"] = p["T1"] * 10**rng.uniform(-1.7, -0.1)
    return p

def split_many(dem, time, derived, ancestral, size, tag):
    """msprime allows at most 100 populations per event; the Andean lattice has 248 nodes.
    Chunk the split through short-lived hub populations (they exist for 1e-3 generations, so
    they contribute no coalescence of their own)."""
    if len(derived) <= 90:
        dem.add_population_split(time=time, derived=list(derived), ancestral=ancestral); return
    chunks = [list(derived[i:i+90]) for i in range(0, len(derived), 90)]
    hubs = []
    for ci, ch in enumerate(chunks):
        h = f"{tag}_h{ci}"
        dem.add_population(name=h, initial_size=size)
        dem.add_population_split(time=time, derived=ch, ancestral=h)
        hubs.append(h)
    dem.add_population_split(time=time+1e-3, derived=hubs, ancestral=ancestral)

def build(model, p):
    dem = msprime.Demography()
    for i in range(K): dem.add_population(name=f"n{i}", initial_size=p["Ne"])
    for a, b in EDGES:
        dem.set_migration_rate(source=f"n{a}", dest=f"n{b}", rate=p["m"])
        dem.set_migration_rate(source=f"n{b}", dest=f"n{a}", rate=p["m"])
    co = [f"n{i}" for i in CO_N]; an = [f"n{i}" for i in AN_N]
    if model == "H1":
        dem.add_population(name="FOUNDER", initial_size=max(p["Ne"]*p["founder"], 1.0))
        dem.add_population(name="ANC", initial_size=p["Nanc"])
        split_many(dem, p["Tcol"], co, "FOUNDER", max(p["Ne"]*p["founder"],1.0), "co")
        # the founder lineage rejoins the Andean lattice at the node closest to the gap
        dem.add_population_split(time=p["Tcol"]+p["Tb"], derived=["FOUNDER"], ancestral=f"n{SRC_NODE}")
        split_many(dem, p["T1"], an, "ANC", p["Nanc"], "an")
    else:
        dem.add_population(name="ANC_C", initial_size=p["Nanc"])
        dem.add_population(name="ANC_A", initial_size=p["Nanc"])
        dem.add_population(name="ANC",   initial_size=p["Nanc"])
        if model == "H3":
            for a, b in GAP_PAIRS:
                dem.set_migration_rate(source=f"n{a}", dest=f"n{b}", rate=p["mcross"])
                dem.set_migration_rate(source=f"n{b}", dest=f"n{a}", rate=p["mcross"])
            for a, b in GAP_PAIRS:      # contact is recent: shut the corridor off before T_contact
                dem.add_migration_rate_change(time=p["Tcontact"], source=f"n{a}", dest=f"n{b}", rate=0)
                dem.add_migration_rate_change(time=p["Tcontact"], source=f"n{b}", dest=f"n{a}", rate=0)
        split_many(dem, p["T1"], co, "ANC_C", p["Nanc"], "co")
        split_many(dem, p["T1"], an, "ANC_A", p["Nanc"], "an")
        dem.add_population_split(time=p["T1"]+p["Tsep"]+0.01, derived=["ANC_C","ANC_A"], ancestral="ANC")
    dem.sort_events()
    return dem

def simulate(model, p, seed, L_bp=1_000_000, r=1e-7, mu=1e-7, target=NSNP, apply_missing=True):
    dem = build(model, p)
    samples = {f"n{STAND_NODE[j]}": int(NPER[j]) for j in range(ND)}
    anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L_bp,
                               recombination_rate=r, random_seed=seed, ploidy=2)
    ts = None; mu_try = mu
    for k in range(6):
        ts = msprime.sim_mutations(anc, rate=mu_try, random_seed=seed+1+k,
                                   model=msprime.BinaryMutationModel())
        if ts.num_sites >= target: break
        mu_try *= max(2.0, 1.3*target/max(ts.num_sites, 1))
    simulate.last_nsites = ts.num_sites
    if ts.num_sites < 500: return None
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
        s = summarise(Gd, DEME, COORDS, COASTAL, nd=ND)
        if s is None: return None
        return model, vec(s), p
    except Exception:
        return None
