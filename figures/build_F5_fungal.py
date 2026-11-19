#!/usr/bin/env python3
"""F5 — the fungal endophyte community.

A  dbRDA constrained on range + compartment (rank 2, so dbRDA1 and dbRDA2 both exist;
   range alone is a 2-level factor and yields a single axis), Bray-Curtis
B  ORDER composition per site from TREE means, scaled to 100%, Coastal stands first
C  Shannon diversity by range and compartment, Wilcoxon
D  host genetic effect tested from the SNP matrix directly, not from kinship estimates

Journal style: no title, no subtitle, no on-canvas footnote; panels lettered only.
All values read from figdata/fungfig_r1/, produced by scripts/fungal_fig.R with ITS_TABLE=r1 (McCartney job 45373).
"""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

BG,PANEL="#ffffff","#ffffff";  # journal: plain white (was cream #f6f4ec/#fffdf7, 17 Nov 2026)
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
GRAY="#b9b5a8"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,"axes.labelcolor":INK,
    "xtick.color":MUTED,"ytick.color":MUTED,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.labelsize":10.5,"pdf.fonttype":42})
Q="figdata/fungfig_r1/"   # 19 Nov 2026: forward-read (R1) table, per the 4 Aug pre-registered rule
A=pd.read_csv(Q+"panelA_scores.csv"); AS=pd.read_csv(Q+"panelA_stats.csv").iloc[0]
B=pd.read_csv(Q+"panelB_phylum.csv"); C=pd.read_csv(Q+"panelC_shannon.csv")
CT=pd.read_csv(Q+"panelC_tests.csv")
# BH across the tests actually REPORTED in panel C -- the three tree-level contrasts. The
# stored p_BH corrects over all six rows, which double-counts the same contrast at library and
# tree level and is the wrong family.
_REP=["range within needle (tree means)","range within branch (tree means)","compartment (paired)"]
_sub=CT[CT["test"].isin(_REP)].copy().sort_values("p").reset_index(drop=True)
_n=len(_sub); _adj=[]; _prev=1.0
for _i in range(_n-1,-1,-1):
    _v=min(_prev,_sub.loc[_i,"p"]*_n/(_i+1)); _adj.insert(0,_v); _prev=_v
_sub["p_BH_reported"]=_adj
BH=dict(zip(_sub["test"],_sub["p_BH_reported"])); D=pd.read_csv(Q+"panelD_results.csv")

fig=plt.figure(figsize=(13.8,11.2)); fig.patch.set_facecolor(BG)
gsT=GridSpec(1,2,figure=fig,left=0.070,right=0.975,top=0.955,bottom=0.600,wspace=0.26,
             width_ratios=[1.0,1.12])
gsB=GridSpec(1,2,figure=fig,left=0.070,right=0.975,top=0.470,bottom=0.070,wspace=0.30,
             width_ratios=[1.0,1.15])
def style(ax):
    ax.set_facecolor(PANEL)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3,width=0.9)
def tag(ax,t): ax.set_title(t,loc="left",fontsize=14,fontweight="bold",color=INK,pad=9)

# ---------------- A: dbRDA on range ----------------
ax=fig.add_subplot(gsT[0,0]); style(ax)
for rg,col in (("Coastal",COAST),("Andean",AND)):
    for cp,mk in (("needle","o"),("branch","D")):
        m=(A["range"]==rg)&(A["compartment"]==cp)
        ax.scatter(A.loc[m,"CAP1"],A.loc[m,"CAP2"],s=42,c=col,marker=mk,
                   edgecolor="white",lw=0.6,alpha=.88,zorder=3,label=f"{rg}, {cp}")
ax.axhline(0,color=RULE,lw=0.8,zorder=1); ax.axvline(0,color=RULE,lw=0.8,zorder=1)
ax.set_xlabel(f"dbRDA1 ({AS['cap1_pct']:.2f}% of total inertia)")
ax.set_ylabel(f"dbRDA2 ({AS['cap2_pct']:.2f}%)")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=8.2,loc="upper right",ncol=2,
          columnspacing=0.9,handletextpad=0.4)
# the ordination is drawn at library level so compartment is visible; the RANGE statistic
# quoted is the tree-level one, since tree is the unit
ax.text(0.015,0.035,"range  $R^2$ = %.3f, $p$ = %.4f   (%d trees)\ncompartment  $R^2$ = %.3f, $p$ = %.4f   (within tree)"
        %(AS["tree_range_R2"],AS["tree_range_p"],int(AS["n_tree"]),AS["comp_R2"],AS["comp_p"]),
        transform=ax.transAxes,fontsize=9,color=INK,va="bottom",
        bbox=dict(fc=PANEL,ec=RULE,lw=0.8,boxstyle="round,pad=0.35"))
tag(ax,"A")

# ---------------- B: phylum per site ----------------
ax=fig.add_subplot(gsT[0,1]); style(ax)
phy=[c for c in B.columns if c not in ("site","range")]
# 10 most abundant families plus the two catch-alls; catch-alls kept neutral so the named
# families carry the color
FAMPAL=["#2f6f4e","#a8543a","#35708c","#d9a441","#7fa9c4","#8c6d9c",
        "#5f8c5a","#c47f5a","#6b8fa8","#b07f9c"]
PAL={}
_named=[c for c in B.columns if c not in ("site","range","other assigned","unassigned")]
for _i,_f in enumerate(_named): PAL[_f]=FAMPAL[_i%len(FAMPAL)]
PAL["other assigned"]="#cfcabc"; PAL["unassigned"]=GRAY
x=np.arange(len(B)); bot=np.zeros(len(B))
for p_ in phy:
    v=100*B[p_].values
    ax.bar(x,v,bottom=bot,width=0.78,color=PAL.get(p_,"#cfcabc"),
           edgecolor="white",linewidth=0.6,label=p_,zorder=3)
    bot+=v
ax.set_xticks(x)
ax.set_xticklabels(B["site"],fontsize=9.5,fontweight="bold")
for t,rg in zip(ax.get_xticklabels(),B["range"]):
    t.set_color(COAST if rg=="Coastal" else AND)
ax.axvline(1.5,color=INK,lw=1.3,ls="--",zorder=5)
ax.set_ylim(0,100); ax.set_ylabel("% of reads")
ax.set_xlim(-0.6,len(B)-0.4)
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=6.8,loc="upper center",
          bbox_to_anchor=(0.5,-0.09),ncol=3,columnspacing=0.9,handlelength=1.1,labelspacing=0.35)
tag(ax,"B")

# ---------------- C: Shannon ----------------
ax=fig.add_subplot(gsB[0,0]); style(ax)
groups=[("needle","Coastal"),("needle","Andean"),("branch","Coastal"),("branch","Andean")]
rng=np.random.default_rng(4)
for i,(cp,rg) in enumerate(groups):
    v=C.loc[(C["compartment"]==cp)&(C["range"]==rg),"shannon"].values
    col=COAST if rg=="Coastal" else AND
    bp=ax.boxplot([v],positions=[i],widths=0.52,patch_artist=True,showfliers=False,zorder=2)
    bp["boxes"][0].set(facecolor=col,alpha=0.22,edgecolor=col,linewidth=1.3)
    for k in ("whiskers","caps"): [w.set(color=col,linewidth=1.1) for w in bp[k]]
    bp["medians"][0].set(color=INK,linewidth=2.0)
    ax.scatter(np.full(len(v),i)+rng.uniform(-.15,.15,len(v)),v,s=17,c=col,
               edgecolor="white",lw=0.4,alpha=.75,zorder=3)
ax.set_xticks(range(4))
ax.set_xticklabels([f"{c}\n{r}" for c,r in groups],fontsize=8.8)
ax.set_ylabel("Shannon diversity")
ax.set_xlim(-0.6,3.6)
yt=ax.get_ylim()[1]
def pv(row): return CT.loc[CT["test"]==row,"p"].values[0]
for cp,lab in (("needle",0),("branch",2)):
    key=f"range within {cp} (tree means)"
    ax.plot([lab,lab+1],[yt*0.965]*2,color=MUTED,lw=1.0)
    ax.text(lab+0.5,yt*0.975,"$p_{BH}$ = %.3f"%BH[key],ha="center",fontsize=8.4,
            color=INK if BH[key]<0.05 else MUTED)
ax.set_ylim(top=yt*1.06)
ax.text(0.015,0.030,"needle vs branch, paired within tree:  $p_{BH}$ = %.4f\n"
        "range with both compartments pooled per tree:  $p$ = %.3f"
        %(BH["compartment (paired)"],pv("range (tree level, both compartments)")),
        transform=ax.transAxes,fontsize=8.6,color=INK,va="bottom",
        bbox=dict(fc=PANEL,ec=RULE,lw=0.8,boxstyle="round,pad=0.3"))
tag(ax,"C")

# ---------------- D: host genetics from the SNP matrix ----------------
ax=fig.add_subplot(gsB[0,1]); style(ax)
MODELS=[("genetic PCs, no conditioning",MUTED,"o"),
        ("genetic PCs | stand",COAST,"s"),
        ("genetic PCs | range",AND,"D")]
ks=[2,4,8]
for nm,col,mk in MODELS:
    sub=D[D["model"]==nm].sort_values("k")
    ax.plot(sub["k"],sub["R2"],color=col,lw=2.0,marker=mk,ms=8,mfc=PANEL,mew=1.7,zorder=4,label=nm)
    # the three models converge at k = 2, so their p labels are fanned vertically
    # the two conditioned models nearly coincide, so their labels are offset sideways as well
    off={"genetic PCs, no conditioning":(0,12,"center"),
         "genetic PCs | stand":(-21,-15,"right"),
         "genetic PCs | range":(23,-15,"left")}[nm]
    for _,r in sub.iterrows():
        sig=r["p"]<0.05
        ax.annotate(("%.4f" if r["p"]<0.001 else "%.3f")%r["p"],(r["k"],r["R2"]),textcoords="offset points",
                    xytext=(off[0],off[1]),ha=off[2],fontsize=8.2,
                    color=INK if sig else MUTED,fontweight="bold" if sig else "normal")
for nm,val in (("range",D.loc[D["model"]=="range","R2"].values[0]),
               ("stand",D.loc[D["model"]=="site_code","R2"].values[0])):
    ax.axhline(val,color=GREEN,lw=1.2,ls=":",zorder=2)
    ax.text(2.05,val,"%s, %.3f"%(nm,val),fontsize=8.2,color=GREEN,va="bottom",ha="left",
            bbox=dict(fc=PANEL,ec="none",pad=1.2),zorder=6)
ax.set_xscale("log",base=2); ax.set_xticks(ks); ax.set_xticklabels(ks)
ax.set_xlim(1.80,9.9)
ax.set_xlabel("number of SNP principal components fitted")
ax.set_ylabel(r"variance explained in fungal composition ($R^2$)")
ax.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=8.2,loc="upper left",
          bbox_to_anchor=(0.0,0.90))
tag(ax,"D")

for ext in ("png","pdf"):
    fig.savefig(f"figures/journal/F5_fungal_community.{ext}",
                dpi=300 if ext=="png" else None,facecolor=BG)
print("F5_fungal_community written")
print("  C, BH over the three reported tree-level tests: %s"
      %{k:round(v,4) for k,v in BH.items()})
print("  A: TREE-level range R2 %.4f p %.4f (n=%d: %d Coastal, %d Andean) | compartment R2 %.4f p %.4f"
      %(AS["tree_range_R2"],AS["tree_range_p"],AS["n_tree"],AS["n_coastal"],AS["n_andean"],
        AS["comp_R2"],AS["comp_p"]))
print("     stand, reported not the unit: R2 %.4f p %.4f | betadisper p %.4f"
      %(AS["tree_stand_R2"],AS["tree_stand_p"],AS["disp_p"]))
print("  D: unconditioned p %s | conditioned on stand p %s"
      %(list(D[D["model"]=="genetic PCs, no conditioning"]["p"].round(4)),
        list(D[D["model"]=="genetic PCs | stand"]["p"].round(4))))
