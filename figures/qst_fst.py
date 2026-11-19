#!/usr/bin/env python3
"""Q_ST from McIntosh et al. 2025 Table S4 vs our F_ST.

Nested random-effects model  trait ~ Region + Population(Region) + Residual.
  V_resid = MSE_resid
  V_pop   = (MSE_pop - MSE_resid) / n0
  Q_ST    = V_pop / (V_pop + 2 V_within)

NOTE: their residual is PHENOTYPIC within-population variance, not additive genetic
(no family term in the published model). V_P >= V_A, so these Q_ST are LOWER BOUNDS.
"""
import csv, numpy as np

# trait, MSE_region, MSE_pop, MSE_resid, df_resid, significant-for-population
S4 = [
 ("Needle width",      26.7,  5.27, 148.5/192, 192, True),
 ("Needle area",       15.6,  0.32,  43.6/194, 194, False),
 ("Needle succulence", 11.0,  0.86,  76.0/194, 194, True),
 ("Number of whorls",   9.5,  0.80,  66.9/193, 193, True),
 ("Needle thickness",   0.5,  0.06,   3.7/194, 194, True),
 ("Needle length",     14.9,  0.31,  81.2/192, 192, False),
 ("Number of branches",24.7,  2.67, 237.5/193, 193, True),
 ("Needle mass/area",   0.6,  0.21,  15.3/194, 194, True),
 ("Branch length",    133.2, 22.24,2268.7/192, 192, True),
 ("Percent survival",   0.2, 0.015,   2.8/201, 201, False),
 ("Percent germin.",  0.021,0.0250,   2.4/192, 192, True),
]
# sample sizes per population (Table S1)
ni = np.array([20,20,20,20,20,7,19,20,20,20,7,20], float)
N  = ni.sum(); p = len(ni)
n0 = (N - (ni**2).sum()/N) / (p-1)          # unbalanced-design correction
print("populations %d | total seedlings %d | n0 = %.2f\n" % (p, N, n0))

print("%-20s %9s %9s %9s %8s"%("trait","V_pop","V_resid","Q_ST","sig(pop)"))
print("-"*60)
q_sig=[]; q_all=[]
for t,msr,msp,mse,dfr,sig in S4:
    vpop = max((msp-mse)/n0, 0.0)
    q = vpop/(vpop+2*mse) if (vpop+2*mse)>0 else np.nan
    q_all.append(q)
    if sig: q_sig.append(q)
    print("%-20s %9.4f %9.4f %9.4f %8s"%(t,vpop,mse,q,"yes" if sig else "no"))
print("-"*60)
print("mean Q_ST, population-significant traits (n=%d): %.4f"%(len(q_sig),np.mean(q_sig)))
print("median                                        : %.4f"%np.median(q_sig))
print("range                                         : %.4f - %.4f"%(min(q_sig),max(q_sig)))

# ---- our F_ST, partitioned the same way ----
COASTAL={'PN','TR','VA'}
rows=list(csv.DictReader(open('figdata/host_fst_pairs.tsv'),delimiter='\t'))
within_a=[float(r['fst']) for r in rows
          if r['stand_a'] not in COASTAL and r['stand_b'] not in COASTAL]
within_c=[float(r['fst']) for r in rows
          if r['stand_a'] in COASTAL and r['stand_b'] in COASTAL]
between =[float(r['fst']) for r in rows
          if (r['stand_a'] in COASTAL) != (r['stand_b'] in COASTAL)]
print("\nour F_ST (host_fst_pairs.tsv), partitioned like their design:")
print("  within Andes  (%2d pairs) mean %.4f  median %.4f"%(len(within_a),np.mean(within_a),np.median(within_a)))
print("  within Coastal(%2d pairs) mean %.4f  median %.4f"%(len(within_c),np.mean(within_c),np.median(within_c)))
print("  between regions(%2d pairs) mean %.4f  median %.4f"%(len(between),np.mean(between),np.median(between)))
print("\nCOMPARISON (population level, within region):")
print("  Q_ST  %.4f    vs   F_ST(within Andes)  %.4f"%(np.mean(q_sig),np.mean(within_a)))
