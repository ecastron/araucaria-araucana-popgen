"""Is T identified on its own, or only through T/(2Ne)? Floor-100 reference table, two-refugia models."""
import glob, pickle, numpy as np
from sklearn.ensemble import RandomForestRegressor
P="/home/ecastron/araucaria"
Xs,ys,Ps=[],[],[]
for f in sorted(glob.glob(f"{P}/results/Q1/abc_floor100/s6_*.npz")):
    d=np.load(f,allow_pickle=True)
    if d["X"].size==0: continue
    Xs.append(d["X"]); ys.append(d["y"]); Ps.append(d["P"]); PN=list(d["param_names"])
X=np.vstack(Xs); y=np.concatenate(ys); PP=np.vstack(Ps)
obs=pickle.load(open(f"{P}/results/Q1/abc/observed2.pkl","rb"))["obs"]
m=np.isin(y,[2,3]); X2,P2=X[m],PP[m]
logT=np.log10(P2[:,PN.index("T1")]); logNe=np.log10(P2[:,PN.index("Ne")]); logtau=logT-np.log10(2)-logNe
def post(t,name):
    rf=RandomForestRegressor(n_estimators=500,oob_score=True,min_samples_leaf=5,n_jobs=-1,random_state=1).fit(X2,t)
    la=rf.apply(X2); lo=rf.apply(obs[None,:])[0]; w=np.zeros(len(t))
    for k in range(la.shape[1]):
        s=la[:,k]==lo[k]; w[s]+=1/s.sum()
    w/=w.sum(); o=np.argsort(t); c=np.cumsum(w[o]); q=lambda p:10**t[o][np.searchsorted(c,p)]
    print(f"{name:8s} OOB R2 {rf.oob_score_:.3f} | median {q(.5):9.3f} 95% [{q(.025):9.3f}, {q(.975):9.3f}] | prior [{10**t.min():.3f}, {10**t.max():.1f}]")
    return w
wT=post(logT,"T"); wNe=post(logNe,"Ne"); wtau=post(logtau,"T/2Ne")
# among simulations weighted as matching the data, do T and Ne move together?
def wcorr(a,b,w):
    ma,mb=np.average(a,weights=w),np.average(b,weights=w)
    return np.average((a-ma)*(b-mb),weights=w)/np.sqrt(np.average((a-ma)**2,weights=w)*np.average((b-mb)**2,weights=w))
print(f"weighted corr(logT, logNe) among data-matching sims: {wcorr(logT,logNe,wT):+.3f}  (prior corr {np.corrcoef(logT,logNe)[0,1]:+.3f})")
# R2 inflation check: does OOB R2 for T rise just because the prior got wider?
narrow=logT>=3.0
rf=RandomForestRegressor(n_estimators=500,oob_score=True,min_samples_leaf=5,n_jobs=-1,random_state=1).fit(X2[narrow],logT[narrow])
print(f"T OOB R2 restricted to the OLD prior range (T>=1000): {rf.oob_score_:.3f}  (n={narrow.sum()})")
