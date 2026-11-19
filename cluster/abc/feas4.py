import sys, time, numpy as np, itertools as it
sys.path.insert(0,"/home/ecastron/araucaria/scripts/abc")
import sim4 as S
from sumstats import summarise, freqs_by_deme, hudson_fst
from chain import _km
EU=np.array([[_km(S.COORDS[a],S.COORDS[b]) for b in range(12)] for a in range(12)])
AN=[i for i in range(12) if not S.COASTAL[i]]
print("TARGET meanFST=+0.130 sd=0.031 CV=0.238 ibd=+0.647 ibdAND=+0.610 ratio=+1.470 within=0.109 between=0.160 coast=+1.190 priv=0.774\n",flush=True)
print(f"{'Ne':>5} {'m':>6} {'Nm':>5} {'T1':>7} {'Nanc':>7} | {'meanFST':>8} {'CV':>6} {'ibd':>6} {'ibdAND':>6} {'ratio':>6} {'within':>7} {'betw':>7} {'ovl':>4} {'coast':>6} {'priv':>6} {'s':>5}",flush=True)
for Ne,m in ((1000,0.003),(1000,0.008),(3000,0.003)):
    for T1 in (10000, 40000):
        for Nanc in (30000, 150000):
            p=dict(Ne=Ne,m=m,T1=T1,Tsep=T1,Nanc=Nanc)
            t=time.time()
            try: Gd=S.simulate("H2",p,4242)
            except Exception as e:
                print(f"{Ne:5d} {m:6.3f} {Ne*m:5.1f} {T1:7d} {Nanc:7d} | ERROR {str(e)[:40]}",flush=True); continue
            if Gd is None:
                print(f"{Ne:5d} {m:6.3f} {Ne*m:5.1f} {T1:7d} {Nanc:7d} | too few sites ({S.simulate.last_nsites})",flush=True); continue
            st=summarise(Gd,S.SIM_DEME,S.COORDS,S.COASTAL,nd=12)
            Pf,N=freqs_by_deme(Gd,S.SIM_DEME,12); rec=[]
            for a,b in it.combinations(range(12),2):
                msk=(N[a]>=6)&(N[b]>=6)
                if msk.sum()<50: continue
                rec.append((a,b,EU[a,b],hudson_fst(Pf[a][msk],Pf[b][msk],N[a][msk],N[b][msk])))
            d=np.array([r[2] for r in rec]); f=np.array([r[3] for r in rec])
            isA=np.array([(r[0] in AN and r[1] in AN) for r in rec])
            bt=np.array([(r[0] in AN)!=(r[1] in AN) for r in rec])
            print(f"{Ne:5d} {m:6.3f} {Ne*m:5.1f} {T1:7d} {Nanc:7d} | {f.mean():+8.4f} {st['fst_cv']:6.3f} "
                  f"{st['ibd_r']:+6.3f} {np.corrcoef(np.log(d[isA]),f[isA])[0,1]:+6.3f} {st['fst_between_within_ratio']:6.3f} "
                  f"{f[~bt].mean():7.4f} {f[bt].mean():7.4f} {'YES' if f[~bt].max()>f[bt].min() else 'NO':>4} "
                  f"{st['coastal_excess_sd']:+6.3f} {st['private_ratio']:6.3f} {time.time()-t:5.0f}",flush=True)
