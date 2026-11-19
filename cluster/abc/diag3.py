import sys, time, numpy as np, itertools as it, msprime, pickle, collections
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import freqs_by_deme, hudson_fst
stands=pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))["stands"]

# graph (lattice) distance between the stand nodes: gene flow follows the graph, not the crow
adj=collections.defaultdict(list)
for a,b in S.EDGES: adj[a].append(b); adj[b].append(a)
def bfs(src):
    d={src:0}; q=[src]
    while q:
        u=q.pop(0)
        for v in adj[u]:
            if v not in d: d[v]=d[u]+1; q.append(v)
    return d
GD=np.full((12,12),np.nan)
for j in range(12):
    d=bfs(int(S.STAND_NODE[j]))
    for k in range(12): GD[j,k]=d.get(int(S.STAND_NODE[k]),np.nan)
EU=np.array([[S._km(S.COORDS[a],S.COORDS[b]) for b in range(12)] for a in range(12)])
iu=np.triu_indices(12,1)
ok=np.isfinite(GD[iu])
print(f"graph-distance vs euclidean correlation (within-component pairs, n={ok.sum()}): "
      f"r={np.corrcoef(GD[iu][ok],EU[iu][ok])[0,1]:+.3f}", flush=True)
print("TARGET  meanFST=+0.130 sd=0.031 CV=0.238 ibd_r(euclid)=+0.647 ibd_r(within)=+0.610\n", flush=True)

samples={f"n{S.STAND_NODE[j]}":int(S.NPER[j]) for j in range(12)}
def run(Ne,m,TA,tag):
    dem=msprime.Demography()
    for i in range(S.K): dem.add_population(name=f"n{i}", initial_size=Ne)
    for a,b in S.EDGES:
        dem.set_migration_rate(source=f"n{a}",dest=f"n{b}",rate=m)
        dem.set_migration_rate(source=f"n{b}",dest=f"n{a}",rate=m)
    dem.add_population(name="ANC", initial_size=Ne*10)
    S.split_many(dem, TA, [f"n{i}" for i in range(S.K)], "ANC", Ne*10, "all")
    dem.sort_events()
    t=time.time()
    anc=msprime.sim_ancestry(samples=samples,demography=dem,sequence_length=1_000_000,
                             recombination_rate=1e-7,random_seed=5,ploidy=2)
    ts=None; mu=1e-7
    for k in range(6):
        ts=msprime.sim_mutations(anc,rate=mu,random_seed=6+k,model=msprime.BinaryMutationModel())
        if ts.num_sites>=30000: break
        mu*=max(2.0,1.3*30000/max(ts.num_sites,1))
    Gm=ts.genotype_matrix(); Gd=(Gm[:,0::2]+Gm[:,1::2]).T.astype(float)
    if Gd.shape[1]>30000: Gd=Gd[:,np.linspace(0,Gd.shape[1]-1,30000).astype(int)]
    rng=np.random.default_rng(3); cols=rng.choice(S.MASK.shape[1],size=Gd.shape[1],replace=True)
    Gmiss=Gd.copy(); Gmiss[S.MASK[:,cols]]=np.nan
    print(f"\n########## {tag}  Ne={Ne} m={m} Nm={Ne*m:.2f} sampled/deme={32/(2*Ne)*100:.1f}% "
          f"({ts.num_sites} sites, {time.time()-t:.0f}s) ##########", flush=True)
    for G_,lab in ((Gd,"NO missing"),(Gmiss,"WITH missing")):
        Pf,N=freqs_by_deme(G_,S.DEME,12); rows=[]
        for a,b in it.combinations(range(12),2):
            msk=(N[a]>=6)&(N[b]>=6)
            if msk.sum()<50: continue
            rows.append((EU[a,b],GD[a,b],hudson_fst(Pf[a][msk],Pf[b][msk],N[a][msk],N[b][msk]),
                         S.COASTAL[a]!=S.COASTAL[b],a,b))
        eu=np.array([r[0] for r in rows]); gd=np.array([r[1] for r in rows])
        f=np.array([r[2] for r in rows]); bt=np.array([r[3] for r in rows])
        r_eu=np.corrcoef(np.log(eu),f)[0,1]
        w=~bt; r_win=np.corrcoef(np.log(eu[w]),f[w])[0,1]
        gok=np.isfinite(gd)
        r_gd=np.corrcoef(np.log(np.maximum(gd[gok],1)),f[gok])[0,1]
        print(f"  [{lab}] meanFST={f.mean():+.4f} sd={f.std():.4f} CV={f.std()/f.mean():.3f} "
              f"ibd_r(euclid)={r_eu:+.3f} ibd_r(within)={r_win:+.3f} ibd_r(GRAPH)={r_gd:+.3f}", flush=True)
        if lab=="WITH missing":
            print(f"    within-range FST {f[w].min():.4f}-{f[w].max():.4f} | between {f[bt].min():.4f}-{f[bt].max():.4f}"
                  f" | overlap={'YES' if f[w].max()>f[bt].min() else 'NO'}", flush=True)
            o=np.argsort(eu)
            print("    pairwise (km, graphsteps, FST):", flush=True)
            for i in o:
                print(f"      {stands[rows[i][4]]}-{stands[rows[i][5]]:<3} {eu[i]:6.1f}km "
                      f"{'--' if not np.isfinite(gd[i]) else int(gd[i]):>3} steps  FST={f[i]:.4f}"
                      f"  {'BETWEEN' if bt[i] else ''}", flush=True)

run(500, 1e-3, 20000, "larger demes, low migration")
run(2000, 3e-4, 40000, "large demes, very low migration")
