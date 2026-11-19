import sys, time, numpy as np, itertools as it, msprime
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import summarise, freqs_by_deme, hudson_fst

def table(Gd, tag):
    Pf,N=freqs_by_deme(Gd,S.DEME,12); rows=[]
    for a,b in it.combinations(range(12),2):
        m=(N[a]>=6)&(N[b]>=6)
        if m.sum()<50: continue
        f=hudson_fst(Pf[a][m],Pf[b][m],N[a][m],N[b][m])
        d=S._km(S.COORDS[a],S.COORDS[b])
        rows.append((d,f,S.COASTAL[a]!=S.COASTAL[b],a,b))
    d=np.array([r[0] for r in rows]); f=np.array([r[1] for r in rows]); bt=np.array([r[2] for r in rows])
    r_all=np.corrcoef(np.log(d),f)[0,1]
    r_win=np.corrcoef(np.log(d[~bt]),f[~bt])[0,1] if (~bt).sum()>3 else np.nan
    print(f"\n--- {tag} ---")
    print(f"  meanFST={f.mean():+.4f} sd={f.std():.4f} CV={f.std()/f.mean():.3f} ibd_r={r_all:+.3f} ibd_r(within)={r_win:+.3f}")
    print(f"  within: {f[~bt].mean():+.4f} [{f[~bt].min():.4f}-{f[~bt].max():.4f}]   between: {f[bt].mean():+.4f} [{f[bt].min():.4f}-{f[bt].max():.4f}]" if bt.sum() else "")
    # per-stand mean FST to all others -> are some stands systematic outliers?
    st=S.__dict__.get("STANDS")
    per=np.zeros(12); cnt=np.zeros(12)
    for (dd,ff,bb,a,b) in rows:
        per[a]+=ff; per[b]+=ff; cnt[a]+=1; cnt[b]+=1
    per=per/np.maximum(cnt,1)
    order=np.argsort(per)[::-1]
    import pickle
    stands=pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))["stands"]
    print("  per-stand mean FST to all others: " + "  ".join(f"{stands[i]}={per[i]:.3f}" for i in order))
    return f

# degrees of the lattice nodes the stands sit on
deg=np.zeros(S.K,int)
for a,b in S.EDGES: deg[a]+=1; deg[b]+=1
import pickle
stands=pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))["stands"]
print("stand -> lattice node degree (low degree = peripheral = extra drift):")
print("  " + "  ".join(f"{stands[j]}:n{S.STAND_NODE[j]}(deg{deg[S.STAND_NODE[j]]})" for j in range(12)))
print(f"  lattice degree distribution: min={deg.min()} median={int(np.median(deg))} max={deg.max()}")

samples={f"n{S.STAND_NODE[j]}":int(S.NPER[j]) for j in range(12)}

print("\n================ CONTROL 1: ONE homogeneous lattice, NO split, NO missingness ================")
print("If a plain stepping stone at equilibrium gives clean IBD, the problem is in the model build.")
for Ne,m in ((50,1e-2),(50,5e-2),(200,1e-2)):
    dem=msprime.Demography()
    for i in range(S.K): dem.add_population(name=f"n{i}", initial_size=Ne)
    for a,b in S.EDGES:
        dem.set_migration_rate(source=f"n{a}",dest=f"n{b}",rate=m)
        dem.set_migration_rate(source=f"n{b}",dest=f"n{a}",rate=m)
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
    table(Gd, f"plain lattice Ne={Ne} m={m} NO missing  ({time.time()-t:.0f}s, {ts.num_sites} sites)")
    # and with missingness, to isolate its contribution
    rng=np.random.default_rng(3)
    cols=rng.choice(S.MASK.shape[1],size=Gd.shape[1],replace=True)
    G2=Gd.copy(); G2[S.MASK[:,cols]]=np.nan
    table(G2, f"plain lattice Ne={Ne} m={m} WITH missing")
