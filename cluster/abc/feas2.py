import sys, time, numpy as np, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim4 as S
from sumstats import summarise, freqs_by_deme, hudson_fst
def mf(G):
    Pf,N=freqs_by_deme(G,S.DEME,12); f=[]
    for a,b in it.combinations(range(12),2):
        m=(N[a]>=6)&(N[b]>=6)
        if m.sum()<50: continue
        f.append(hudson_fst(Pf[a][m],Pf[b][m],N[a][m],N[b][m]))
    return float(np.mean(f)),float(np.std(f))
print(f"chains: K={S.K} (Andean {S.NA}, Coastal {S.NC}) | corridor d{S.GAP_A}<->d{S.GAP_C}", flush=True)
print("TARGET  meanFST=+0.130 sd=0.031 CV=0.238 ibd_r=+0.647 ratio=+1.470 coast=+1.190 priv=0.774\n", flush=True)
print(f"{'Ne':>5} {'m':>7} {'Nm':>5} {'T1':>9} | {'nsit':>6} {'mu':>8} {'meanFST':>8} {'sd':>6} {'CV':>6} {'ibd_r':>7} {'ratio':>6} {'coast':>6} {'priv':>6} {'sec':>6}", flush=True)
for Ne in (300, 1000):
    for m in (0.005, 0.02, 0.08):
        for T1 in (100_000, 500_000):
            p=dict(Ne=Ne,m=m,T1=T1,Tsep=T1*0.5)
            t=time.time()
            try: Gd=S.simulate("H2",p,4242)
            except Exception as e:
                print(f"{Ne:5d} {m:7.3f} {Ne*m:5.1f} {T1:9d} | ERROR {str(e)[:50]}", flush=True); continue
            if Gd is None:
                print(f"{Ne:5d} {m:7.3f} {Ne*m:5.1f} {T1:9d} | too few sites", flush=True); continue
            st=summarise(Gd,S.DEME,S.COORDS,S.COASTAL,nd=12); a,b=mf(Gd)
            print(f"{Ne:5d} {m:7.3f} {Ne*m:5.1f} {T1:9d} | {S.simulate.last_nsites:6d} {S.simulate.last_mu:8.1e} "
                  f"{a:+8.4f} {b:6.4f} {st['fst_cv']:6.3f} {st['ibd_r']:+7.3f} {st['fst_between_within_ratio']:6.3f} "
                  f"{st['coastal_excess_sd']:+6.3f} {st['private_ratio']:6.3f} {time.time()-t:6.1f}", flush=True)
