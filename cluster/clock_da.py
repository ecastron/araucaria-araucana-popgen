#!/usr/bin/env python
"""Model-free split date from net divergence between ranges (17 Nov 2026).

d_a = d_xy - (pi_C + pi_A)/2 per sequenced site, invariant sites included, from POPULATION ALLELE
FREQUENCIES (each call contributes its two alleles; IUPAC codes split one-each). Individual
heterozygosity is never used as a statistic. T_years ~ d_a / (2 mu_year): generation time cancels.
Three locus sets: all; locus mean H_obs <= 0.20 (the S2E ceiling); that set with <= 5 SNPs per locus.
CIs: 200 bootstrap resamples over loci. Uses the 163 relatedness-pruned trees.
"""
import numpy as np, csv, json, sys
P="/home/ecastron/araucaria"
LOCI=f"{P}/rad/ipyrad_full/full_ct091_outfiles/full_ct091.loci"
site={r[0]:r[1] for r in csv.reader(open(f"{P}/metadata/popmap_177.tsv"),delimiter="\t") if len(r)>=2}
keep={l.split()[0] for l in open(f"{P}/results/Q1/feems/keep.txt")}
COAST={"VA","TR","PN"}
MINC,MINA=6,15                     # minimum called individuals per range at a site
# allele-count lookup: A C G T, two copies per diploid call
LUT=np.zeros((256,4),np.int8)
for ch,v in {"A":(2,0,0,0),"C":(0,2,0,0),"G":(0,0,2,0),"T":(0,0,0,2),"R":(1,0,1,0),"Y":(0,1,0,1),
             "S":(0,1,1,0),"W":(1,0,0,1),"K":(0,0,1,1),"M":(1,1,0,0)}.items(): LUT[ord(ch)]=v
HET=np.zeros(256,bool)
for ch in "RYSWKM": HET[ord(ch)]=True
rows=[]   # per locus: n_sites, sum piC, sum piA, sum dxy, locus Hobs, n_snps
def flush(names,seqs):
    idx=[i for i,n in enumerate(names) if n in keep and n in site]
    if not idx: return
    pop=np.array([site[names[i]] in COAST for i in idx])
    if pop.sum()<MINC or (~pop).sum()<MINA: return
    S=np.frombuffer("".join(seqs[i] for i in idx).encode(),np.uint8).reshape(len(idx),-1)
    C=LUT[S].astype(np.int32)                       # samples x sites x 4
    called=C.sum(2)>0
    nC=called[pop].sum(0); nA=called[~pop].sum(0)
    ok=(nC>=MINC)&(nA>=MINA)
    if ok.sum()==0: return
    cC=C[pop][:,ok].sum(0); cA=C[~pop][:,ok].sum(0)   # sites x 4 allele counts
    tC=cC.sum(1,keepdims=True); tA=cA.sum(1,keepdims=True)
    pC=cC/tC; pA=cA/tA
    piC=(tC[:,0]/(tC[:,0]-1))*(1-(pC**2).sum(1)); piA=(tA[:,0]/(tA[:,0]-1))*(1-(pA**2).sum(1))
    dxy=1-(pC*pA).sum(1)
    var=((cC+cA)>0).sum(1)>1
    h=HET[S][:,ok][:,var]; cl=called[:,ok][:,var]
    hobs=h.sum()/cl.sum() if cl.sum() else 0.0
    rows.append((ok.sum(),piC.sum(),piA.sum(),dxy.sum(),hobs,var.sum()))
names,seqs=[],[]
with open(LOCI) as fh:
    for line in fh:
        if line.startswith("//"):
            flush(names,seqs); names,seqs=[],[]
        else:
            n,s=line.split(); names.append(n.rsplit(".",1)[0]); seqs.append(s)
R=np.array(rows,float)
print(f"loci used {len(R):,}, sites {R[:,0].sum():,.0f}",flush=True)
MU=[1e-10,3e-10,7e-10,1.4e-9,2e-9]
rng=np.random.default_rng(1)
def summ(M,label):
    def est(m):
        s=m[:,0].sum(); piC,piA,dxy=m[:,1].sum()/s,m[:,2].sum()/s,m[:,3].sum()/s
        return piC,piA,dxy,dxy-(piC+piA)/2
    e=est(M); boot=np.array([est(M[rng.integers(0,len(M),len(M))]) for _ in range(200)])
    lo,hi=np.percentile(boot,[2.5,97.5],axis=0)
    out={"set":label,"loci":int(len(M)),"sites":int(M[:,0].sum()),
         "piC":e[0],"piA":e[1],"dxy":e[2],"da":e[3],"da_ci":[lo[3],hi[3]],
         "T_years":{f"{m:.1e}":[e[3]/(2*m),lo[3]/(2*m),hi[3]/(2*m)] for m in MU}}
    print(f"\n[{label}] loci {len(M):,} sites {M[:,0].sum():,.0f}")
    print(f"  piC {e[0]:.5f}  piA {e[1]:.5f}  dxy {e[2]:.5f}  da {e[3]:.6f}  (95% {lo[3]:.6f} to {hi[3]:.6f})")
    for m in MU: print(f"  mu {m:.1e}/site/yr -> T {e[3]/(2*m):>10,.0f} yr  ({lo[3]/(2*m):,.0f} to {hi[3]/(2*m):,.0f})")
    return out
res=[summ(R,"all loci"), summ(R[R[:,4]<=0.20],"locus Hobs <= 0.20"),
     summ(R[(R[:,4]<=0.20)&(R[:,5]<=5)],"Hobs <= 0.20 and <= 5 SNPs")]
json.dump(res,open(f"{P}/results/clock/da_dates.json","w"),indent=1)
