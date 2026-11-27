#!/usr/bin/env python3
"""F2: host genetic structure: PCA, eigenvalues, pairwise Hudson F_ST on distance, sNMF K = 2.

JOURNAL STYLE (4 Sep 2026): no title, subtitle, panel headlines or on-canvas footnote.
Interpretation lives in ARAUCARIA_figure_legends.docx."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BG,PANEL="#ffffff","#ffffff";  # journal: plain white (was cream #f6f4ec/#fffdf7, 17 Nov 2026)
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,
    "axes.labelcolor":INK,"xtick.color":MUTED,"ytick.color":MUTED,
    "xtick.labelsize":9,"ytick.labelsize":9,"axes.labelsize":11,"pdf.fonttype":42})
COASTAL={"VA","TR","PN"}

pop=dict(l.split() for l in open("figdata/popmap_177.tsv") if l.strip())
# 26 Nov 2026: PCA and sNMF on the 163 analyzed trees (was all 177 libraries, relatives and
# second libraries included). Sources: scripts/revision/structure.sbatch (plink 1.9 PCA on
# LD-pruned SNPs, --indep-pairwise 50 5 0.2, 66,131 SNPs) and scripts/revision/snmf.R
# (LEA sNMF K = 2, 10 repetitions, alpha 100, best run; column "coastal" oriented to Coastal).
ev=pd.read_csv("pca_data/q1_pca_plink163.eigenvec",sep=r"\s+",header=None)
ev["id"]=ev[0].str.replace(r"\.\d+$","",regex=True)
ev["stand"]=ev["id"].map(pop)
ev=ev.dropna(subset=["stand"])
assert len(ev)==163, len(ev)
val=np.loadtxt("pca_data/q1_pca_plink163.eigenval"); pct=100*val/len(ev)
qs=pd.read_csv("pca_data/snmf_K2_163.tsv",sep="\t")[["id","stand","coastal"]]
assert len(qs)==163 and qs["stand"].notna().all()
fst=pd.read_csv("figdata/host_fst_pairs.tsv",sep="\t")

fig=plt.figure(figsize=(13.6,11.6)); fig.patch.set_facecolor(BG)
gs=GridSpec(3,6,figure=fig,height_ratios=[1.06,0.94,0.72],
            hspace=0.40,wspace=0.72,left=0.075,right=0.975,top=0.955,bottom=0.055)

def style(ax):
    ax.set_facecolor(PANEL)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3,width=0.9)
def tag(ax,t,dx=None,dy=None):
    # letter set ABOVE the axes, journal style -- inside the frame it collided with
    # the leading eigenvalue in B and with the sNMF bars in E
    ax.set_title(t,loc="left",fontsize=14,fontweight="bold",color=INK,pad=8)

# ---- A: PCA ----
axA=fig.add_subplot(gs[0,0:4]); style(axA)
for grp,col,lab in ((True,COAST,"Coastal (VA, TR, PN)"),(False,AND,"Andean (9 stands)")):
    m=ev["stand"].isin(COASTAL)==grp
    axA.scatter(ev.loc[m,2],ev.loc[m,3],s=52,c=col,edgecolor="white",linewidth=0.7,
                alpha=.92,label=lab,zorder=3)
axA.axhline(0,color=RULE,lw=0.9,zorder=1); axA.axvline(0,color=RULE,lw=0.9,zorder=1)
axA.set_xlabel(f"PC1  ({pct[0]:.1f}%)"); axA.set_ylabel(f"PC2  ({pct[1]:.1f}%)")
axA.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9.5,loc="upper right")
for s in ("PN","TR","VA"):
    sub=ev[ev["stand"]==s]
    axA.annotate(s,(sub[2].mean(),sub[3].mean()),fontsize=10,fontweight="bold",
                 ha="center",va="center",color=INK,
                 bbox=dict(fc=PANEL,ec=RULE,lw=0.7,boxstyle="round,pad=0.22"),zorder=5)
tag(axA,"A")

# ---- C: eigenvalues ----
axC=fig.add_subplot(gs[0,4:6]); style(axC)
cols=[GREEN]+[ "#dcd8ca"]*9
axC.bar(range(1,11),val[:10],color=cols,edgecolor=MUTED,linewidth=0.6)
axC.set_xticks(range(1,11)); axC.set_xlabel("component"); axC.set_ylabel("eigenvalue")
axC.text(1,val[0]+0.14,f"{val[0]:.2f}",ha="center",fontsize=9.5,fontweight="bold",color=GREEN)
tag(axC,"B")

# ---- D: F_ST vs log distance ----
axD=fig.add_subplot(gs[1,0:4]); style(axD)
w=fst[fst.between_==0] if "between_" in fst else fst[fst["between"]==0]
b=fst[fst["between"]==1]
x=np.log(fst["km"]); y=fst["fst"]
sl,ic=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),100)
axD.plot(xs,ic+sl*xs,color=MUTED,lw=1.4,ls="--",zorder=2,
         label=f"all pairs  r = +{np.corrcoef(x,y)[0,1]:.3f}")
axD.scatter(np.log(w["km"]),w["fst"],s=58,c=AND,edgecolor="white",lw=0.7,zorder=3,
            label=f"within range (n = {len(w)})")
axD.scatter(np.log(b["km"]),b["fst"],s=58,c=COAST,edgecolor="white",lw=0.7,marker="D",
            zorder=3,label=f"between ranges (n = {len(b)})")
axD.axhspan(w["fst"].max(),b["fst"].min(),color=GREEN,alpha=0.10,zorder=1)
tk=[10,20,50,100,200]; axD.set_xticks(np.log(tk)); axD.set_xticklabels(tk)
axD.set_xlabel("geographic distance (km, log scale)"); axD.set_ylabel(r"pairwise Hudson $F_{ST}$")
axD.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9,loc="upper left")
tag(axD,"C")

# ---- E: marginal strip ----
axE=fig.add_subplot(gs[1,4:6]); style(axE)
rng=np.random.default_rng(3)
for i,(d,col,lab) in enumerate(((w["fst"],AND,"within"),(b["fst"],COAST,"between"))):
    axE.scatter(np.full(len(d),i)+rng.uniform(-.12,.12,len(d)),d,s=46,c=col,
                edgecolor="white",lw=0.6,zorder=3)
    axE.hlines(d.median(),i-.28,i+.28,color=INK,lw=2,zorder=4)
axE.axhspan(w["fst"].max(),b["fst"].min(),color=GREEN,alpha=0.10,zorder=1)
axE.set_xticks([0,1]); axE.set_xticklabels(["within","between"])
axE.set_xlim(-.6,1.6); axE.set_ylabel(r"$F_{ST}$")
tag(axE,"D")

# ---- F: sNMF K = 2 ----
axF=fig.add_subplot(gs[2,0:6]); style(axF)
order=[s for s in sorted(qs["stand"].unique()) if s not in COASTAL]+["PN","TR","VA"]
pos=0; ticks=[]
for s in order:
    sub=qs[qs["stand"]==s].sort_values("coastal",ascending=False)
    n=len(sub); xs2=np.arange(pos,pos+n)
    axF.bar(xs2,sub["coastal"],width=1.0,color=COAST,linewidth=0)
    axF.bar(xs2,1-sub["coastal"],width=1.0,bottom=sub["coastal"],color=AND,linewidth=0)
    ticks.append((pos+n/2,s)); pos+=n
    axF.axvline(pos-0.5,color=PANEL,lw=2.2,zorder=4)
axF.set_xticks([t[0] for t in ticks]); axF.set_xticklabels([t[1] for t in ticks],fontweight="bold",fontsize=10)
axF.set_xlim(-0.5,pos-0.5); axF.set_ylim(0,1); axF.set_ylabel("sNMF coefficient (K = 2)")
tag(axF,"E")

fig.savefig("figures/journal/F2_host_structure.png",dpi=300,facecolor=BG,bbox_inches="tight")
fig.savefig("figures/journal/F2_host_structure.pdf",facecolor=BG,bbox_inches="tight")
print("F2 written")
print(f"  within  n={len(w)}  range {w['fst'].min():.4f}-{w['fst'].max():.4f}")
print(f"  between n={len(b)}  range {b['fst'].min():.4f}-{b['fst'].max():.4f}")
print(f"  gap = {b['fst'].min()-w['fst'].max():+.4f}   ibd_r = {np.corrcoef(x,y)[0,1]:+.3f}")
