#!/usr/bin/env python
"""Scale-free summary statistics for ABC-RF, computed identically on observed and simulated data.

Deliberately EXCLUDES the site-frequency spectrum and every heterozygosity-based quantity:
collapsed paralogs inject false intermediate-frequency variants, which is H3'/secondary-contact's
signature, so an SFS-based ABC would select that model for technical reasons.

FST is Hudson's estimator (allele frequencies only, no genotype configurations).
Statistics are ratios / correlations / dimensionless counts, so a multiplicative distortion of the
FST scale — which the dp20 comparison showed is real (~3x) — cancels.
"""
import numpy as np, itertools, warnings
warnings.filterwarnings('ignore', r'Mean of empty slice')
warnings.filterwarnings('ignore', r'invalid value encountered')

def hudson_fst(p1, p2, n1, n2):
    """Hudson FST per locus, ratio-of-averages across loci (Bhatia et al. 2013)."""
    num = (p1-p2)**2 - p1*(1-p1)/np.maximum(n1-1,1) - p2*(1-p2)/np.maximum(n2-1,1)
    den = p1*(1-p2) + p2*(1-p1)
    ok = np.isfinite(num) & np.isfinite(den) & (den > 0)
    return float(num[ok].sum()/den[ok].sum()) if ok.sum() else np.nan

def freqs_by_deme(G, deme, nd):
    """G: individuals x loci alt-allele dosage (0,1,2) with np.nan missing.
       returns (freq nd x L, n_chrom nd x L)"""
    L = G.shape[1]
    P = np.full((nd, L), np.nan); N = np.zeros((nd, L))
    for d in range(nd):
        sub = G[deme == d]
        cnt = np.nansum(sub, axis=0)
        obs = np.sum(~np.isnan(sub), axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            P[d] = cnt/(2*obs)
        N[d] = 2*obs
    return P, N

def summarise(G, deme, coords, coastal, nd=12):
    """G individuals x loci (0/1/2, nan=missing); deme int array; coords nd x 2 lon/lat;
       coastal boolean array length nd."""
    P, N = freqs_by_deme(G, deme, nd)
    pairs = list(itertools.combinations(range(nd), 2))
    fst, km, between, pia, pib = [], [], [], [], []
    R = 6371.0
    for a, b in pairs:
        m = (N[a] >= 6) & (N[b] >= 6)
        if m.sum() < 50: continue
        f = hudson_fst(P[a][m], P[b][m], N[a][m], N[b][m])
        if not np.isfinite(f): continue
        pia.append(a); pib.append(b)
        la1, lo1 = np.radians(coords[a,1]), np.radians(coords[a,0])
        la2, lo2 = np.radians(coords[b,1]), np.radians(coords[b,0])
        h = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
        fst.append(f); km.append(2*R*np.arcsin(np.sqrt(h)))
        between.append(coastal[a] != coastal[b])
    fst = np.array(fst); km = np.array(km); between = np.array(between)
    pia = np.array(pia); pib = np.array(pib)
    out = {}
    if len(fst) < 10: return None
    # 1-2 IBD: correlation and fit quality on log distance (scale-free)
    x = np.log(np.maximum(km, 1e-6))
    r = np.corrcoef(x, fst)[0,1]
    out["ibd_r"] = r
    out["ibd_r2"] = r*r
    # 3-4 use FST floored at 0 for the ratio/CV: negative Hudson estimates are sampling noise
    # around zero (Bhatia 2013), and letting them through makes the ratio flip sign and the CV
    # diverge whenever a simulation is near-panmictic. Observed pairwise FST are all >> 0, so this
    # is a no-op on the observed vector and only stabilises the simulated cloud.
    fstf = np.maximum(fst, 0.0)
    RATIO_CAP, CV_CAP = 20.0, 10.0
    bw = fstf[between].mean() if between.sum() else np.nan
    wi = fstf[~between].mean() if (~between).sum() else np.nan
    if np.isfinite(bw) and np.isfinite(wi):
        out["fst_between_within_ratio"] = float(np.clip(bw/wi, 0, RATIO_CAP)) if wi > 1e-6 else RATIO_CAP
    else:
        out["fst_between_within_ratio"] = np.nan
    mfs = fstf.mean()
    out["fst_cv"] = float(np.clip(fstf.std()/mfs, 0, CV_CAP)) if mfs > 1e-6 else CV_CAP
    # 5 residual excess of between-cordillera pairs after IBD (scale-free: standardised)
    A = np.c_[np.ones(len(x)), x]
    beta, *_ = np.linalg.lstsq(A, fst, rcond=None)
    res = fst - A@beta
    out["coastal_excess_sd"] = (res[between].mean()-res[~between].mean())/res.std() if res.std()>0 else np.nan
    # 6-7 PCA eigenvalue structure (ratios)
    Pc = P.copy()
    col_ok = np.all(np.isfinite(Pc), axis=0)
    Z = Pc[:, col_ok]
    if Z.shape[1] > 20:
        Z = (Z - Z.mean(0))
        sv = np.linalg.svd(Z, compute_uv=False)
        ev = sv**2
        out["pc1_pc2_ratio"] = ev[0]/ev[1] if len(ev) > 1 and ev[1] > 0 else np.nan
        out["pc1_share"] = ev[0]/ev.sum() if ev.sum() > 0 else np.nan
    else:
        out["pc1_pc2_ratio"] = np.nan; out["pc1_share"] = np.nan
    # 8 f3 admixture: fraction of triples negative (dimensionless)
    neg = tot = 0
    idx = [d for d in range(nd) if np.isfinite(P[d]).sum() > 100]
    for c in idx:
        for a, b in itertools.combinations([i for i in idx if i != c], 2):
            m = np.isfinite(P[a]) & np.isfinite(P[b]) & np.isfinite(P[c]) & (N[c] > 2)
            if m.sum() < 100: continue
            f3 = np.mean((P[c][m]-P[a][m])*(P[c][m]-P[b][m]) - P[c][m]*(1-P[c][m])/(N[c][m]-1))
            tot += 1; neg += (f3 < 0)
    out["f3_frac_neg"] = neg/tot if tot else np.nan
    # 9 nestedness: are coastal alleles a subset of andean? (H1' signature, dimensionless)
    co = [d for d in range(nd) if coastal[d]]; an = [d for d in range(nd) if not coastal[d]]
    pc_ = np.nanmean(P[co], axis=0); pa_ = np.nanmean(P[an], axis=0)
    m = np.isfinite(pc_) & np.isfinite(pa_)
    priv_c = np.sum((pc_[m] > 0.05) & (pa_[m] < 0.01))
    priv_a = np.sum((pa_[m] > 0.05) & (pc_[m] < 0.01))
    out["private_ratio"] = priv_c/(priv_a+1)

    # ---- statistics targeting recent contact through the gap ----
    # Secondary contact is corridor-mediated: it should depress differentiation specifically for
    # the stand pairs flanking the narrowest crossing, leaving a SPATIAL signature that none of
    # the nine statistics above measures. The corridor is defined from geography alone (the
    # closest cross-gap stand pair), so it is identical for observed and simulated data.
    co_i = [d for d in range(nd) if coastal[d]]; an_i = [d for d in range(nd) if not coastal[d]]
    def _hav(u, v):
        la1, lo1 = np.radians(u[1]), np.radians(u[0]); la2, lo2 = np.radians(v[1]), np.radians(v[0])
        h = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
        return 2*6371.0*np.arcsin(np.sqrt(h))
    _b = min(((_hav(coords[x], coords[y]), x, y) for x in an_i for y in co_i), key=lambda t: t[0])
    cA, cC = _b[1], _b[2]
    xg = between
    if xg.sum() >= 4:
        fx = fst[xg]
        out["xgap_min_ratio"] = float(fx.min()/fx.mean()) if fx.mean() > 0 else np.nan
        cd = []
        for i in np.where(xg)[0]:
            a, b = pia[i], pib[i]
            aa, bb = (a, b) if not coastal[a] else (b, a)
            cd.append(_hav(coords[aa], coords[cA]) + _hav(coords[bb], coords[cC]))
        cd = np.array(cd)
        out["xgap_corridor_r"] = (float(np.corrcoef(cd, fx)[0,1])
                                  if cd.std() > 0 and fx.std() > 0 else np.nan)
    else:
        out["xgap_min_ratio"] = np.nan; out["xgap_corridor_r"] = np.nan
    # Recent migrants carry LOW-FREQUENCY variants across the gap. Under long isolation rare
    # alleles are private to one range, so their frequencies are uncorrelated between ranges.
    pcf = np.nanmean(P[co_i], axis=0); paf = np.nanmean(P[an_i], axis=0)
    okf = np.isfinite(pcf) & np.isfinite(paf)
    pooled = (pcf + paf)/2.0
    lm = okf & (pooled > 0) & (pooled < 0.10)
    out["lowfreq_cor"] = (float(np.corrcoef(pcf[lm], paf[lm])[0,1])
                          if lm.sum() > 50 and pcf[lm].std() > 0 and paf[lm].std() > 0 else np.nan)
    return out


STAT_ORDER = ["ibd_r","ibd_r2","fst_between_within_ratio","fst_cv","coastal_excess_sd",
              "pc1_pc2_ratio","pc1_share","f3_frac_neg","private_ratio",
              "xgap_min_ratio","xgap_corridor_r","lowfreq_cor"]
def vec(d):
    return np.array([d[k] if d and k in d and np.isfinite(d[k]) else np.nan for k in STAT_ORDER])
