import sys, time, numpy as np, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim2 as S
from sumstats import summarise, freqs_by_deme, hudson_fst
def meanfst(G):
    Pf,N=freqs_by_deme(G,S.DEME,12); f=[]
    for a,b in it.combinations(range(12),2):
        m=(N[a]>=6)&(N[b]>=6)
        if m.sum()<50: continue
        f.append(hudson_fst(Pf[a][m],Pf[b][m],N[a][m],N[b][m]))
    return float(np.mean(f))
print(f"lattice {S.K} nodes {S.EDGES.shape[0]} edges | gap {S.GAP_KM:.1f} km", flush=True)
print("TARGET  meanFST=+0.130  CV=+0.238  ibd_r=+0.647  ratio=+1.470  coast=+1.190\n", flush=True)
print(f"{'Ne':>5} {'m':>8} {'T1':>7} | {'nsit':>6} {'meanFST':>8} {'CV':>6} {'ibd_r':>7} {'ratio':>6} {'coast':>6} {'priv':>6} {'sec':>7}", flush=True)
for Ne in (50, 200, 800):
    for m in (1e-3, 1e-2, 1e-1):
        for T1 in (5000, 40000):
            p=dict(Ne=Ne,Nanc=2000,m=m,T1=T1,Tsep=T1)
            t=time.time()
            try: Gd=S.simulate("H2",p,2024)
            except Exception as e:
                print(f"{Ne:5d} {m:8.1e} {T1:7d} | ERROR {str(e)[:60]}", flush=True); continue
            if Gd is None:
                print(f"{Ne:5d} {m:8.1e} {T1:7d} | too few sites  {time.time()-t:6.1f}s", flush=True); continue
            st=summarise(Gd,S.DEME,S.COORDS,S.COASTAL,nd=12)
            print(f"{Ne:5d} {m:8.1e} {T1:7d} | {S.simulate.last_nsites:6d} {meanfst(Gd):+8.4f} {st['fst_cv']:6.3f} "
                  f"{st['ibd_r']:+7.3f} {st['fst_between_within_ratio']:6.3f} {st['coastal_excess_sd']:+6.3f} "
                  f"{st['private_ratio']:6.3f} {time.time()-t:6.1f}s", flush=True)
