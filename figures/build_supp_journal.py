#!/usr/bin/env python3
"""Supplementary figures S1-S4, journal style (white background, lettered panels, no titles).

Restructured 17 Nov 2026 from eight figures to four. Old -> new:
  S1 ABC validation  + S3 chain architecture           -> S1  (A-E, F-G)
  S2 allele balance  + S8 deficit diagnosis + S6B       -> S2  (A-E)
  S6A de novo vs reference F_ST                         -> S3
  S5 no core                                            -> S4
  S4 taxonomic resolution: CUT (it defended the kinship R2, which the main text no longer uses;
     F5D tests host genetics with SNP PCs conditioned on stand)
  S7 fungal description: CUT (its three panels are older versions of F5 A-C)
Old scripts are kept as build_S6_journal.retired-20261117.py, build_S7_journal.retired-20261117.py
and build_supp_journal.bak-8fig-20261117.py. Every plotted value is read from figdata/ or is a
constant carried over unchanged from those scripts.
"""
import csv, json, re, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

BG,PANEL="#ffffff","#ffffff"   # journal: plain white
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
GRAY="#b9b5a8"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,"axes.labelcolor":INK,
    "xtick.color":MUTED,"ytick.color":MUTED,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.labelsize":10.5,"pdf.fonttype":42})
OUT="figures/supplementary/journal"

def style(ax):
    ax.set_facecolor(PANEL)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3,width=0.9)
def tag(ax,t,extra=""):
    ax.set_title(f"{t}{('   '+extra) if extra else ''}",loc="left",fontsize=12,fontweight="bold",color=INK,pad=8)
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.png",dpi=300,facecolor=BG,bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf",facecolor=BG,bbox_inches="tight"); plt.close(fig); print("  ",name)

# =====================================================================================
# S1: validation of the coalescent model choice (old S1 + old S3)
# =====================================================================================
# 18 Nov 2026: S1 A-E read the literature-mu clock ABC (sim9), the primary demographic analysis;
# 12 shape statistics plus 3 log clock statistics. The shape-only inputs (abc_for_fig/abc_supp) are retired.
# abcrf 2.0 results (abcrf/abcrf_analysis.R); the confusion matrices and importances come from the first
# of ten replicate forests, and the LDA axes abcrf adds by default appear in the importance panel
AB=json.load(open("figdata/abcrf_results.json"))
RT=pd.read_csv("abcrf/reftable_sim9.csv"); OB=pd.read_csv("abcrf/observed_sim9.csv",index_col=0)
STAT=list(OB.columns); X=RT[STAT].values; obs=OB.loc["primary",STAT].values
cmat=lambda k,n: np.array([r[:n] for r in AB[k]["confusion"]]).astype(int)
A=dict(cm4=cmat("four",4),cm1=cmat("stage1",2),cm2=cmat("stage2",2))
IMPn=list(AB["four"]["importance"].keys()); IMPv=np.array(list(AB["four"]["importance"].values()))
G9=dict(acc4=1-AB["four"]["prior_error_mean"],acc1=1-AB["stage1"]["prior_error_mean"],acc2=1-AB["stage2"]["prior_error_mean"])
fig=plt.figure(figsize=(12.6,14.6)); fig.patch.set_facecolor(BG)
gs=GridSpec(4,2,figure=fig,left=0.085,right=0.965,top=0.975,bottom=0.05,hspace=0.62,wspace=0.30,
            height_ratios=[1.0,0.95,0.95,0.95])

# A: prior-predictive check
ax=fig.add_subplot(gs[0,:]); style(ax)
for j in range(len(STAT)):
    v=X[:,j]; v=v[np.isfinite(v)]
    # plot the PERCENTILE of the observed value within the simulated distribution, which is what
    # the tick labels assert; linear value space put right-skewed statistics at the bottom edge
    ax.plot([j,j],[0.025,0.975],color=RULE,lw=6,solid_capstyle="round")
    ax.plot([j,j],[0.25,0.75],color=MUTED,lw=6,alpha=.55,solid_capstyle="round")
    ax.plot(j,0.5,"o",ms=4.5,color=MUTED)
    ax.plot(j,np.mean(v<=obs[j]),"*",ms=17,color=GREEN,markeredgecolor="white",markeredgewidth=0.8,zorder=5)
ax.set_xticks(range(len(STAT))); ax.set_xticklabels([s.replace("_"," ") for s in STAT],rotation=34,ha="right",fontsize=8)
ax.set_ylim(-0.2,1.2); ax.set_yticks([0,.5,1]); ax.set_yticklabels(["2.5%","50%","97.5%"])
ax.set_ylabel("percentile of the observed value\nwithin the simulated distribution")
tag(ax,"A")

def cmplot(ax,cm,labs,t,acc):
    style(ax); n=cm/cm.sum(1,keepdims=True)
    ax.imshow(n,cmap=mpl.colors.LinearSegmentedColormap.from_list("g",[PANEL,GREEN]),vmin=0,vmax=1)
    for i in range(len(labs)):
        for j in range(len(labs)):
            ax.text(j,i,f"{cm[i,j]:,}",ha="center",va="center",fontsize=9.5,
                    color="white" if n[i,j]>0.5 else INK)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs,fontsize=9.5)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=9.5)
    ax.set_xlabel("classified as"); ax.set_ylabel("true model")
    tag(ax,t,f"OOB {acc:.3f}")
cmplot(fig.add_subplot(gs[1,0]),A["cm4"],["H1′","H1b′","H2′","H3′"],"B",G9["acc4"])
ax=fig.add_subplot(gs[1,1]); style(ax)
o=np.argsort(IMPv)[::-1]
ax.barh(range(len(o)),IMPv[o]/IMPv.sum(),color=GREEN,alpha=.85)
ax.set_yticks(range(len(o))); ax.set_yticklabels([IMPn[i].replace("_"," ") for i in o],fontsize=8.5)
ax.invert_yaxis(); ax.set_xlabel("variable importance (share of total)")
tag(ax,"C")
cmplot(fig.add_subplot(gs[2,0]),A["cm1"],["one\nrefugium","two\nrefugia"],"D",G9["acc1"])
cmplot(fig.add_subplot(gs[2,1]),A["cm2"],["isolation","contact"],"E",G9["acc2"])

# F: a single axis captures the stand layout within each range
S=np.load("figdata/stand_coords.npz",allow_pickle=True); SC,CO=S["coords"],S["coastal"]
def km(a,b):
    R=6371.0;p1,p2=np.radians(a[1]),np.radians(b[1])
    h=np.sin((p2-p1)/2)**2+np.cos(p1)*np.cos(p2)*np.sin(np.radians(b[0]-a[0])/2)**2
    return 2*R*np.arcsin(np.sqrt(h))
ax=fig.add_subplot(gs[3,0]); style(ax)
for lab,idx,c in (("Andean",[i for i in range(12) if not CO[i]],AND),("Coastal",[i for i in range(12) if CO[i]],COAST)):
    n=len(idx); Dm=np.array([[km(SC[a],SC[b]) for b in idx] for a in idx])
    J=np.eye(n)-np.ones((n,n))/n; B=-0.5*J@(Dm**2)@J
    w,v=np.linalg.eigh(B); oo=np.argsort(w)[::-1]; x=v[:,oo[0]]*np.sqrt(max(w[oo[0]],0))
    D1=np.abs(x[:,None]-x[None,:]); iu=np.triu_indices(n,1)
    r=np.corrcoef(Dm[iu],D1[iu])[0,1]
    ax.scatter(Dm[iu],D1[iu],s=52,color=c,edgecolor="white",lw=0.7,label=f"{lab}  r = {r:+.3f}")
ax.plot([0,200],[0,200],ls="--",color=MUTED,lw=1.2)
ax.set_xlabel("true pairwise distance (km)"); ax.set_ylabel("distance along the 1D chain (km)")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9.5,loc="upper left")
tag(ax,"F")

# G: isolation by distance produced by each simulation architecture
ax=fig.add_subplot(gs[3,1]); style(ax)
arch=["12 demes,\none per stand","343-node\n2D lattice","1D chains\n(adopted)"]; ibd=[0.06,0.005,0.95]
ax.bar(range(3),ibd,color=[MUTED,MUTED,GREEN],alpha=.9)
ax.axhline(0.647,color=AND,lw=2.2,ls="--")
ax.text(-0.38,0.672,"observed  +0.647",fontsize=9.5,color=AND,fontweight="bold",ha="left")  # left, clear of the bars
for i,v in enumerate(ibd): ax.text(i,v+0.03,f"{v:.3f}",ha="center",fontsize=10,fontweight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels(arch,fontsize=9.5)
ax.set_ylabel("isolation by distance produced\n(correlation of F$_{ST}$ with log distance)")
ax.set_ylim(0,1.08)
tag(ax,"G")
save(fig,"FigS1_abc_validation")

# =====================================================================================
# S2: the heterozygote deficit (old S2A + old S8 A, B, D + old S6B)
# =====================================================================================
H=json.load(open("figdata/het_diagnosis.json"))
fig=plt.figure(figsize=(13.6,11.0)); fig.patch.set_facecolor(BG)
gsT=GridSpec(1,3,figure=fig,left=0.065,right=0.975,top=0.965,bottom=0.625,wspace=0.34)
gsB=GridSpec(1,2,figure=fig,left=0.065,right=0.975,top=0.505,bottom=0.075,wspace=0.62,width_ratios=[1,1.45])

# A: allele balance at >=30x (old S2A)
txt=open("figdata/allele_balance.txt").read()
bins=re.findall(r"(\d\.\d\d)-(\d\.\d\d)\s+(\d+)\s+\(\s*([\d.]+)%\)",txt)
lo=np.array([float(b[0]) for b in bins]); pc=np.array([float(b[3]) for b in bins])
ax=fig.add_subplot(gsT[0,0]); style(ax)
ax.bar(lo+0.025,pc,width=0.045,color=[AND if (l<0.15 or l>=0.85) else COAST for l in lo],
       edgecolor="white",linewidth=0.5)
ax.axvspan(0.15,0.85,color=COAST,alpha=.07,zorder=0)
ax.set_yscale("log"); ax.set_xlabel("allele balance (fraction ALT reads)")
ax.set_ylabel("% of genotypes at ≥30× (log scale)")
tag(ax,"A")

# B: F within stand against depth (old S8A)
ax=fig.add_subplot(gsT[0,1]); style(ax)
dd=H["F_by_genotype_depth"]; x=np.arange(len(dd["bins"]))
ax.plot(x,dd["F_within_stand"],color=AND,lw=2.6,marker="o",ms=7,mfc=PANEL,mew=1.7,zorder=4)
ax.axhline(0,color=INK,lw=1.2,zorder=3)
for xi,yi,ni in zip(x,dd["F_within_stand"],dd["n"]):
    # counts sit up and to the right of each point so the steep descending segment cannot cross them
    ax.text(xi+0.12,yi+0.035,f"{ni/1000:.0f}k",fontsize=7.4,color=MUTED,ha="left",va="bottom")
ax.set_xticks(x); ax.set_xticklabels(dd["bins"],rotation=45,ha="right",fontsize=8.4)
ax.set_ylim(-0.08,1.02)
ax.set_xlabel("read depth per genotype"); ax.set_ylabel(r"$F_{IS}$ within stand")
tag(ax,"B")

# C: rarefying to 30 reads does not remove it (old S8B)
ax=fig.add_subplot(gsT[0,2]); style(ax)
rr=H["rarefaction_to_30_reads"]; xb=np.arange(len(rr["bins"]))
ax.axhline(0,color=INK,lw=1.2,zorder=3)
for i,(nv,rv) in enumerate(zip(rr["native"],rr["rarefied"])):
    ax.plot([i-0.16,i+0.16],[nv,rv],color=GRAY,lw=1.4,zorder=2)
ax.scatter(xb-0.16,rr["native"],s=62,c=MUTED,edgecolor="white",lw=0.9,zorder=4,label="as sequenced")
ax.scatter(xb+0.16,rr["rarefied"],s=62,c=COAST,edgecolor="white",lw=0.9,zorder=5,label="rarefied to 30 reads")
ax.set_xticks(xb); ax.set_xticklabels(rr["bins"],rotation=45,ha="right",fontsize=8.4)
ax.set_ylim(-0.78,0.78)
ax.set_xlabel("original read depth"); ax.set_ylabel(r"$F$")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=8.2,loc="upper right")
tag(ax,"C")

# D: three independent pipelines, all seven depth bins (old S6B, which superseded the 4-bin S8C)
ax=fig.add_subplot(gsB[0,0]); style(ax)
NM=["6-10","10-15","15-20","20-30","30-50","50-100","100+"]; MID=[8,12.5,17.5,25,40,75,150]
fb={q["depth"]:q for q in csv.DictReader(open("figdata/fb_F_by_depth.csv"))}
pr={q["depth"]:q for q in csv.DictReader(open("figdata/polyrad_F_by_depth.csv"))}
# 19 Nov 2026: the binomial series is now READ from the reproduced output (s2e_polyrad_route.R; binomial
# calls on the same 34,674 polyRAD loci and reads); the hard-coded series it replaces had no producer
IB={q["depth"]:q for q in csv.DictReader(open("figdata/ipyrad_binomial_F.tsv"),delimiter="\t")}
# polyRAD de novo line: same estimator as the FreeBayes line and the S2E baseline (s2e_mechanisms.tsv row 9);
# figdata/polyrad_F_by_depth.csv used a numerator over all depths, giving 1.175 instead of 1.207
import pandas as _pd0; _S2E=_pd0.read_csv("figdata/s2e_mechanisms.tsv",sep="\t")
series=[("ipyrad de novo, binomial calls",[float(IB[n]["F_valid"]) for n in NM],MUTED,(0,(5,3)),1.8),
        ("ipyrad de novo, polyRAD",[float(_S2E.loc[9,"F_"+n]) for n in NM],COAST,(0,(2,2)),2.0),
        ("bowtie2 + FreeBayes, polyRAD",[float(fb[n]["F_IS"]) for n in NM],AND,"-",2.6)]
ax.set_xscale("log"); ax.axhline(0,color=INK,lw=1.0,zorder=3)
for lab,v,c,ls,lw in series:
    ax.plot(MID,v,color=c,lw=lw,ls=ls,marker="o",ms=5,mfc=PANEL,mew=1.2,zorder=4,
            label=f"{lab}\nspread {max(v)-min(v):.3f}")
ax.set_xlabel("read depth per genotype (log scale)"); ax.set_ylabel(r"$F_{IS}$ within stand")
ax.set_ylim(-0.62,1.0)
ax.set_xticks([10,20,50,100]); ax.set_xticklabels(["10","20","50","100"])
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=7.6,loc="lower left",labelspacing=0.8)
tag(ax,"D")

# E: no filter removes it. 19 Nov 2026: read from figdata/s2e_mechanisms.tsv (s2e_*.py/.R); every row
# is compared with the unfiltered spread of its OWN genotype route (tick), on the same seven depth bins
import pandas as _pd
T=_pd.read_csv("figdata/s2e_mechanisms.tsv",sep="\t")
NAMES=["binomial calls, unfiltered (baseline)","rarefy every genotype to 30 reads","locus mean depth \u2264 1.25\u00d7 median",
       "depth-CV filter","locus-relative depth (CNV) \u2264 1.5","heterozygosity ceiling 0.20","structural locus filter",
       "single-copy loci only","exclude segmental duplications",
       "polyRAD, unfiltered (baseline)","polyRAD + heterozygosity ceiling 0.20","polyRAD + structural locus filter",
       "polyRAD model vs binomial calls (same loci)","bowtie2 + FreeBayes + polyRAD (replication)"]
assert len(T)==len(NAMES)
ax=fig.add_subplot(gsB[0,1]); style(ax)
gap=[0]*9+[1]*5; y=-(np.arange(len(T))+np.array(gap)*0.8)
ax.axvspan(0,0.15,color=GREEN,alpha=0.13,zorder=1)
for (i,r),yy,nm in zip(T.iterrows(),y,NAMES):
    b,v=r["baseline_spread"],r["filtered_spread"]; isb="baseline" in nm
    col=INK if isb else (COAST if v<b-0.005 else (AND if v>b+0.005 else MUTED))
    ax.barh(yy,v,height=0.56,color=col,alpha=0.30 if isb else 0.85,edgecolor=col,linewidth=1.2,zorder=4)
    if not isb: ax.plot([b,b],[yy-0.36,yy+0.36],color=INK,lw=1.6,zorder=5)
    ax.text(max(v,b)+0.03,yy,"%.3f"%v+("" if isb else "  (%+.0f%%)"%(100*(v-b)/b)),va="center",fontsize=8.6,color=INK,zorder=6)
ax.set_yticks(y); ax.set_yticklabels(NAMES,fontsize=8.6)
ax.set_xlim(0,2.05); ax.set_ylim(y.min()-0.6,0.6); ax.tick_params(axis="y",length=0)
ax.set_xlabel(r"$F$ spread across seven depth bins  (max $-$ min)")
tag(ax,"E")
save(fig,"FigS2_heterozygote_deficit")
print(f"     S2E rows {len(T)}, best filtered spread {T.filtered_spread.min():.3f}")

# =====================================================================================
# S3: de novo versus reference-based F_ST (old S6A)
# =====================================================================================
def load(f):
    return {frozenset((r["stand_a"],r["stand_b"])):(float(r["fst"]),int(r["between"]))
            for r in csv.DictReader(open(f),delimiter="\t")}
den=load("figdata/host_fst_pairs.tsv"); ref=load("figdata/fb_fst_indiv.tsv")
k=sorted(set(den)&set(ref),key=lambda s:sorted(s))
x=np.array([den[i][0] for i in k]); y=np.array([ref[i][0] for i in k]); bt=np.array([den[i][1] for i in k],bool)
r=np.corrcoef(x,y)[0,1]; rs=np.corrcoef(np.argsort(np.argsort(x)),np.argsort(np.argsort(y)))[0,1]
fig=plt.figure(figsize=(6.8,6.4)); fig.patch.set_facecolor(BG)
ax=fig.add_axes([0.14,0.12,0.82,0.84]); style(ax)
lim=[0.06,0.27]
ax.plot(lim,lim,color=MUTED,lw=1.0,ls="--",zorder=1,label="1:1")
ax.scatter(x[~bt],y[~bt],s=52,c=COAST,edgecolor="white",lw=0.7,zorder=3,label=f"within range (n = {(~bt).sum()})")
ax.scatter(x[bt],y[bt],s=52,c=AND,edgecolor="white",lw=0.7,marker="D",zorder=3,label=f"between ranges (n = {bt.sum()})")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("$F_{ST}$, ipyrad de novo"); ax.set_ylabel("$F_{ST}$, bowtie2 + FreeBayes on AraucariaDB")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9,loc="upper left")
ax.text(0.97,0.06,f"Pearson r = {r:.3f}\nSpearman r = {rs:.3f}\n{len(k)} stand pairs",transform=ax.transAxes,
        ha="right",fontsize=10.5,color=INK,linespacing=1.5,bbox=dict(fc=PANEL,ec=RULE,boxstyle="round,pad=0.4"))
save(fig,"FigS3_pipeline_concordance")
print(f"     Pearson {r:.4f} | Spearman {rs:.4f} | pairs {len(k)}")

# =====================================================================================
# S4: no core fungal microbiome at equal depth. 19 Nov 2026: forward-read (R1) table on the Figure 5
# sample set; panel A from its_figdata.R (ITS_TABLE=r1), panel B from core_microbiome.R (microbiome 1.26)
# =====================================================================================
prevA=pd.read_csv("figdata/its_prevalence_r1.csv")["prevalence"].values; NLIB=239
CG=pd.read_csv("figdata/core/r1_core_grid.tsv",sep="\t")
fig=plt.figure(figsize=(13.2,5.6)); fig.patch.set_facecolor(BG)
gs=GridSpec(1,2,figure=fig,left=0.070,right=0.972,top=0.935,bottom=0.165,wspace=0.26)
ax=fig.add_subplot(gs[0,0]); style(ax)
ax.hist(prevA,bins=np.arange(0,1.02,0.025),color=GREEN,alpha=.85,edgecolor="white",linewidth=0.5)
ax.set_yscale("log"); ax.axvline(0.5,color=AND,lw=2.4,ls="--")
ax.set_xlabel(f"prevalence (fraction of the {NLIB} libraries, any read)"); ax.set_ylabel("number of ASVs (log scale)")
ax.text(0.97,0.97,f"{100*np.mean(prevA<=1/NLIB+1e-9):.0f}% occur in exactly one library\n"
        f"most widespread reaches {100*prevA.max():.0f}%",transform=ax.transAxes,ha="right",va="top",
        fontsize=9.5,color=INK,bbox=dict(fc=PANEL,ec=RULE,lw=0.8,boxstyle="round,pad=0.4"))
tag(ax,"A")
ax=fig.add_subplot(gs[0,1]); style(ax)
for lev,col,mk in (("ASV",AND,"^"),("Family",GREEN,"o")):
    for dat,ls,lab in (("relative abundance","-","all reads"),("rarefied to 1062","--","rarefied to 1,062 reads")):
        q=CG[(CG.level==lev)&(CG.data==dat)&(CG.detection==0.001)].sort_values("prevalence")
        ax.plot(q.prevalence,q.n_core,ls=ls,marker=mk,color=col,lw=2.1,ms=7,mfc=PANEL if ls=="--" else col,
                mec=col,mew=1.4,label=f"{'ASVs' if lev=='ASV' else 'families'}, {lab}")
ax.axvline(0.5,color=AND,lw=1.4,ls=":")
ax.set_xlabel("prevalence threshold (fraction of libraries, detection 0.1% of reads)")
ax.set_ylabel("taxa in the core")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9.0)
tag(ax,"B")
save(fig,"FigS4_no_core_microbiome")
print(f"     prevalence: {len(prevA)} ASVs, singletons {100*np.mean(prevA<=1/164+1e-9):.1f}%, max {prevA.max():.3f}")
print("\nall four supplementary figures written")
