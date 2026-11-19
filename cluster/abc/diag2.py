import sys, time, numpy as np, itertools as it, msprime, pickle
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import freqs_by_deme, hudson_fst
stands=pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))["stands"]

def table(Gd, tag):
    Pf,N=freqs_by_deme(Gd,S.DEME,12); rows=[]
    for a,b in it.combinations(range(12),2):
        m=(N[a]>=6)&(N[b]>=6)
        if m.sum()<50: continue
        f=hudson_fst(Pf[a][m],Pf[b][m],N[a][m],N[b][m])
        rows.append((S._km(S.COORDS[a],S.COORDS[b]),f,S.COASTAL[a]!=S.COASTAL[b],a,b))
    d=np.array([r[0] for r in rows]); f=np.array([r[1] for r in rows]); bt=np.array([r[2] for r in rows])
    r_all=np.corrcoef(np.log(d),f)[0,1]
    r_win=np.corrcoef(np.log(d[~bt]),f[~bt])[0,1] if (~bt).sum()>3 else float('nan')
    print(f"  {tag}", flush=True)
    print(f"    meanFST={f.mean():+.4f} sd={f.std():.4f} CV={f.std()/f.mean():.3f} "
          f"ibd_r={r_all:+.3f} ibd_r(within)={r_win:+.3f}", flush=True)
    per=np.zeros(12); cnt=np.zeros(12)
    for (dd,ff,bb,a,b) in rows: per[a]+=ff; per[b]+=ff; cnt[a]+=1; cnt[b]+=1
    per/=np.maximum(cnt,1); o=np.argsort(per)[::-1]
    print(f"    per-stand mean FST: " + " ".join(f"{stands[i]}={per[i]:.3f}" for i in o), flush=True)

deg=np.zeros(S.K,int)
for a,b in S.EDGES: deg[a]+=1; deg[b]+=1
print("stand -> node degree: " + "  ".join(f"{stands[j]}(deg{deg[S.STAND_NODE[j]]})" for j in range(12)), flush=True)
print(f"lattice degrees: min={deg.min()} median={int(np.median(deg))} max={deg.max()}\n", flush=True)
print("TARGET  meanFST=+0.130 sd=0.031 CV=0.238 ibd_r=+0.647 ibd_r(within)=+0.610\n", flush=True)
print("CONTROL: ONE homogeneous lattice, no range split, deep ancestral collapse to bound runtime.", flush=True)

samples={f"n{S.STAND_NODE[j]}":int(S.NPER[j]) for j in range(12)}
for Ne,m,TA in ((50,1e-2,20000),(50,5e-2,20000),(200,1e-2,40000)):
    dem=msprime.Demography()
    for i in range(S.K): dem.add_population(name=f"n{i}", initial_size=Ne)
    for a,b in S.EDGES:
        dem.set_migration_rate(source=f"n{a}",dest=f"n{b}",rate=m)
        dem.set_migration_rate(source=f"n{b}",dest=f"n{a}",rate=m)
    dem.add_population(name="ANC", initial_size=Ne*20)
    S.split_many(dem, TA, [f"n{i}" for i in range(S.K)], "ANC", Ne*20, "all")
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
    print(f"\n[Ne={Ne} m={m} Nm={Ne*m:.1f} TA={TA}] {ts.num_sites} sites, {time.time()-t:.0f}s", flush=True)
    table(Gd, "NO missingness")
    rng=np.random.default_rng(3); cols=rng.choice(S.MASK.shape[1],size=Gd.shape[1],replace=True)
    G2=Gd.copy(); G2[S.MASK[:,cols]]=np.nan
    table(G2, "WITH observed missingness")
