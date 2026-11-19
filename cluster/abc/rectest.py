import sys, time, numpy as np, itertools as it, msprime
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import summarise, freqs_by_deme, hudson_fst
def meanfst(G):
    Pf,N=freqs_by_deme(G,S.DEME,12); f=[]
    for a,b in it.combinations(range(12),2):
        m=(N[a]>=6)&(N[b]>=6)
        if m.sum()<50: continue
        f.append(hudson_fst(Pf[a][m],Pf[b][m],N[a][m],N[b][m]))
    return float(np.mean(f)), float(np.std(f))

p = dict(Ne=50, Nanc=2000, m=1e-3, T1=5000, Tsep=5000)
dem = S.build("H2", p)
samples = {f"n{S.STAND_NODE[j]}": int(S.NPER[j]) for j in range(12)}
print("HYPOTHESIS: excess FST scatter comes from LINKAGE, not demography.")
print("rho = 4*Ne*r*L is the number of independent genealogies. At Ne=50, L=1e6, r=1e-7 -> rho=20.")
print("The RAD data carries ~thousands of independent loci.\n")
print("TARGET  meanFST=+0.130  sdFST=0.031  CV=+0.238  ibd_r=+0.647\n")
print(f"{'r':>9} {'L':>9} {'rho':>8} | {'nsit':>6} {'meanFST':>8} {'sdFST':>7} {'CV':>6} {'ibd_r':>7} {'sec':>7}")
for r_, L_ in ((1e-7,1_000_000),(1e-6,1_000_000),(1e-5,1_000_000),(1e-4,1_000_000)):
    rho = 4*p["Ne"]*r_*L_
    t=time.time()
    try:
        anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=L_,
                                   recombination_rate=r_, random_seed=7, ploidy=2)
        ts=None; mu=1e-7
        for k in range(6):
            ts = msprime.sim_mutations(anc, rate=mu, random_seed=8+k, model=msprime.BinaryMutationModel())
            if ts.num_sites >= 30000: break
            mu *= max(2.0, 1.3*30000/max(ts.num_sites,1))
        Gm=ts.genotype_matrix(); Gd=(Gm[:,0::2]+Gm[:,1::2]).T.astype(float)
        if Gd.shape[1]>30000: Gd=Gd[:,np.linspace(0,Gd.shape[1]-1,30000).astype(int)]
        rng=np.random.default_rng(11)
        cols=rng.choice(S.MASK.shape[1],size=Gd.shape[1],replace=True)
        Gd[S.MASK[:,cols]]=np.nan
        st=summarise(Gd,S.DEME,S.COORDS,S.COASTAL,nd=12)
        mf,sf=meanfst(Gd)
        print(f"{r_:9.1e} {L_:9d} {rho:8.0f} | {ts.num_sites:6d} {mf:+8.4f} {sf:7.4f} {st['fst_cv']:6.3f} "
              f"{st['ibd_r']:+7.3f} {time.time()-t:6.1f}s", flush=True)
    except Exception as e:
        print(f"{r_:9.1e} {L_:9d} {rho:8.0f} | ERROR {str(e)[:60]}  {time.time()-t:6.1f}s", flush=True)
