import sys, time, numpy as np, itertools as it, msprime, pickle
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import freqs_by_deme, hudson_fst
stands=pickle.load(open("/home/ecastron/araucaria/results/Q1/abc/observed.pkl","rb"))["stands"]
EU=np.array([[S._km(S.COORDS[a],S.COORDS[b]) for b in range(12)] for a in range(12)])
print("DIAGNOSIS: no IBD because the ancestral collapse TA was far shallower than the lattice", flush=True)
print("coalescence time ~2*Ne*K. All between-deme pairs then coalesce in the panmictic ancestor", flush=True)
print("at the same time TA, making FST flat regardless of distance.\n", flush=True)
print("TARGET  meanFST=+0.130 sd=0.031 CV=0.238 ibd_r=+0.647 ibd_r(within)=+0.610\n", flush=True)
samples={f"n{S.STAND_NODE[j]}":int(S.NPER[j]) for j in range(12)}
def run(Ne,m,TA):
    Nt=Ne*S.K
    print(f"\n[Ne={Ne} m={m} Nm={Ne*m:.1f} TA={TA:,} | N_total={Nt:,} lattice-coal~{2*Nt:,}] TA/coal={TA/(2*Nt):.2f}", flush=True)
    dem=msprime.Demography()
    for i in range(S.K): dem.add_population(name=f"n{i}", initial_size=Ne)
    for a,b in S.EDGES:
        dem.set_migration_rate(source=f"n{a}",dest=f"n{b}",rate=m)
        dem.set_migration_rate(source=f"n{b}",dest=f"n{a}",rate=m)
    dem.add_population(name="ANC", initial_size=Ne*10)
    S.split_many(dem, TA, [f"n{i}" for i in range(S.K)], "ANC", Ne*10, "all")
    dem.sort_events()
    t=time.time()
    try:
        anc=msprime.sim_ancestry(samples=samples,demography=dem,sequence_length=1_000_000,
                                 recombination_rate=1e-7,random_seed=5,ploidy=2)
    except Exception as e:
        print(f"  ERROR {str(e)[:80]}", flush=True); return
    ts=None; mu=1e-7
    for k in range(6):
        ts=msprime.sim_mutations(anc,rate=mu,random_seed=6+k,model=msprime.BinaryMutationModel())
        if ts.num_sites>=30000: break
        mu*=max(2.0,1.3*30000/max(ts.num_sites,1))
    Gm=ts.genotype_matrix(); Gd=(Gm[:,0::2]+Gm[:,1::2]).T.astype(float)
    if Gd.shape[1]>30000: Gd=Gd[:,np.linspace(0,Gd.shape[1]-1,30000).astype(int)]
    rng=np.random.default_rng(3); cols=rng.choice(S.MASK.shape[1],size=Gd.shape[1],replace=True)
    Gd[S.MASK[:,cols]]=np.nan
    el=time.time()-t
    Pf,N=freqs_by_deme(Gd,S.DEME,12); rows=[]
    for a,b in it.combinations(range(12),2):
        msk=(N[a]>=6)&(N[b]>=6)
        if msk.sum()<50: continue
        rows.append((EU[a,b],hudson_fst(Pf[a][msk],Pf[b][msk],N[a][msk],N[b][msk]),S.COASTAL[a]!=S.COASTAL[b]))
    eu=np.array([r[0] for r in rows]); f=np.array([r[1] for r in rows]); bt=np.array([r[2] for r in rows])
    w=~bt
    print(f"  {ts.num_sites} sites, {el:.0f}s | meanFST={f.mean():+.4f} sd={f.std():.4f} CV={f.std()/f.mean():.3f} "
          f"ibd_r={np.corrcoef(np.log(eu),f)[0,1]:+.3f} ibd_r(within)={np.corrcoef(np.log(eu[w]),f[w])[0,1]:+.3f}", flush=True)
    print(f"  within {f[w].mean():.4f} [{f[w].min():.4f}-{f[w].max():.4f}] | between {f[bt].mean():.4f} "
          f"[{f[bt].min():.4f}-{f[bt].max():.4f}] | overlap={'YES' if f[w].max()>f[bt].min() else 'NO'}", flush=True)
    near=[(e,x) for e,x in zip(eu[w],f[w]) if e<10]; far=[(e,x) for e,x in zip(eu[w],f[w]) if e>150]
    if near and far:
        print(f"  within-range FST at <10km: {np.mean([x for _,x in near]):.4f}   at >150km: {np.mean([x for _,x in far]):.4f}"
              f"   ratio={np.mean([x for _,x in far])/max(np.mean([x for _,x in near]),1e-9):.2f}  (observed ratio 0.084->0.118 = 1.41)", flush=True)
for Ne,m,TA in ((150,0.02,300_000),(150,0.02,2_000_000),(500,0.006,2_000_000)):
    run(Ne,m,TA)
