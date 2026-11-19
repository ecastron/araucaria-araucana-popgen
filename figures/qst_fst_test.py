#!/usr/bin/env python3
"""Q_ST vs F_ST for Araucaria araucana: Coastal vs Andes.

Q_ST  from McIntosh et al. 2025 (Restoration Ecology) Table S4 nested ANOVA
F_ST  from our RAD data, read-based Hudson estimator (no genotype calls)
"""
import numpy as np, json

# ---- Table S4: trait, MS_region, MS_pop(region), MS_resid, population-significant
S4=[("Needle width",      26.7, 5.27, 148.5/192, True),
    ("Needle area",       15.6, 0.32,  43.6/194, False),
    ("Needle succulence", 11.0, 0.86,  76.0/194, True),
    ("Number of whorls",   9.5, 0.80,  66.9/193, True),
    ("Needle thickness",   0.5, 0.06,   3.7/194, True),
    ("Needle length",     14.9, 0.31,  81.2/192, False),
    ("Number of branches",24.7, 2.67, 237.5/193, True),
    ("Needle mass/area",   0.6, 0.21,  15.3/194, True),
    ("Branch length",    133.2,22.24,2268.7/192, True),
    ("Percent germination",0.021,0.025, 2.4/192, True)]
ni=np.array([20,20,20,20,20,7,19,20,20,20,7,20],float)   # Table S1 sample sizes
N=ni.sum(); P=12; R=2
n0=(N-(ni**2).sum()/N)/(P-1)          # seedlings per population
k0=P/R                                 # populations per region (2 coastal, 10 andes)

FST_REG   = 0.1059; FST_REG_SE = 0.0011      # Coastal vs Andes, this study
FST_WITHIN_ANDES_WC = 0.0333                  # weighted W&C, our data
FST_WITHIN_ANDES_HUD= 0.1092                  # Hudson, our data

def qst(ms_hi, ms_lo, ms_res, div, h2):
    v_hi=max((ms_hi-ms_lo)/div,0.0)
    return v_hi/(v_hi+2*h2*ms_res)

H2=[1.0,0.5,0.4,0.35,0.3,0.25,0.2,0.15,0.1]
out={"n0":n0,"k0":k0,"fst_regional":FST_REG,"traits":[t[0] for t in S4],"h2":H2,
     "qst_region":[], "qst_pop":[]}
print("n0 = %.2f seedlings/population | k0 = %.1f populations/region"%(n0,k0))
print("regional F_ST (Coastal vs Andes) = %.4f  (SE %.4f)\n"%(FST_REG,FST_REG_SE))
print("%-6s %26s %26s"%("h^2","REGION  Q_ST (mean/median)","POPULATION Q_ST (mean/median)"))
print("-"*62)
for h2 in H2:
    qr=[qst(a,b,c,n0*k0,h2) for _,a,b,c,_ in S4]
    qp=[qst(b,c,c,n0,h2)    for _,a,b,c,s in S4 if s]
    out["qst_region"].append([float(np.mean(qr)),float(np.median(qr))])
    out["qst_pop"].append([float(np.mean(qp)),float(np.median(qp))])
    mark=" <-- Q_ST > F_ST" if np.mean(qr)>FST_REG else ""
    print("%-6.2f %12.4f %12.4f  %12.4f %12.4f%s"%(
        h2,np.mean(qr),np.median(qr),np.mean(qp),np.median(qp),mark))
# crossover
lo,hi=0.01,1.0
for _ in range(60):
    mid=(lo+hi)/2
    if np.mean([qst(a,b,c,n0*k0,mid) for _,a,b,c,_ in S4])>FST_REG: lo=mid
    else: hi=mid
out["crossover_h2"]=float((lo+hi)/2)
print("\nCrossover: regional Q_ST equals F_ST at h^2 = %.3f"%out["crossover_h2"])
print("Conifer seedling traits typically h^2 = 0.15-0.35 -> Q_ST > F_ST in that whole range.")
json.dump(out,open("figdata/qst_fst.json","w"),indent=1)
print("\nwritten figdata/qst_fst.json")
