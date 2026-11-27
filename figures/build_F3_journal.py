#!/usr/bin/env python3
"""Figure 3: demographic history. A: the four histories with their random-forest vote shares.
B: posterior probabilities for the two nested contrasts, and the migration detection limit.
C: posterior of the split age in years from the clock ABC (literature mutation-rate prior).

18 Nov 2026: all three panels now read the literature-mu clock ABC (sim9), the primary demographic
analysis. The earlier shape-only analysis is kept as build_F4_journal.bak-shape-20261118.py.

JOURNAL STYLE (4 Sep 2026): no title, subtitle, panel headlines, verdict text or on-canvas
footnote. What each shaded band or reference line MEANS is now carried by a legend entry or by
ARAUCARIA_figure_legends.docx, never by a sentence on the canvas.
"""
import json, numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BG,PANEL="#ffffff","#ffffff";  # journal: plain white (was cream #f6f4ec/#fffdf7, 17 Nov 2026)
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
GRAY="#b9b5a8"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,"axes.labelcolor":INK,
    "xtick.color":MUTED,"ytick.color":MUTED,"xtick.labelsize":9.5,"ytick.labelsize":9.5,
    "axes.labelsize":10.5,"pdf.fonttype":42})
# 18 Nov 2026: every value from abcrf 2.0 (abcrf/abcrf_analysis.R); the sklearn version is
# build_F4_journal.bak-sklearn-20261118.py
import pandas as pd
A=json.load(open("figdata/abcrf_results.json"))
D=dict(votes4=A["four"]["votes_mean"],cm4_recall=A["four"]["recall"],
       stage1=dict(name="two refugia vs one",posterior=A["stage1"]["posterior_mean"],oob=1-A["stage1"]["prior_error_mean"],
                   selected="two refugia"),
       stage2=dict(name="isolation vs secondary contact",posterior=A["stage2"]["posterior_mean"],
                   oob=1-A["stage2"]["prior_error_mean"],selected="unresolved (contact, if any, modest)"),
       migration=dict(Nm_median=A["Nm_contact"]["median"],
                      power=[[p["bin"],p["n"],p["called_contact"]] for p in A["power"]]))

fig=plt.figure(figsize=(13.6,11.8)); fig.patch.set_facecolor(BG)
gsA=GridSpec(1,4,figure=fig,left=0.055,right=0.975,top=0.940,bottom=0.700,wspace=0.16)
gsB=GridSpec(2,1,figure=fig,left=0.075,right=0.470,top=0.585,bottom=0.075,hspace=0.62,
             height_ratios=[1,1])
gsC=GridSpec(1,1,figure=fig,left=0.585,right=0.975,top=0.585,bottom=0.075)
def style(ax):
    ax.set_facecolor(PANEL)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3,width=0.9)

# ================= A: the four histories and their votes =================
NAMES=["one refuge,\nAndean source","one refuge,\nCoastal source",
       "two refuges,\nno reconnection","two refuges,\nlater contact"]
votes=D["votes4"]; recall=D["cm4_recall"]
for j in range(4):
    ax=fig.add_subplot(gsA[0,j]); ax.set_facecolor(PANEL)
    ax.set_xlim(0,1); ax.set_ylim(-0.42,1.02); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_color(RULE)
    two = j>=2
    if not two:
        src,der=(0.68,0.30) if j==0 else (0.30,0.68)
        csrc,cder=(AND,COAST) if j==0 else (COAST,AND)
        ax.plot([src,src],[0.20,0.80],color=csrc,lw=6,solid_capstyle="butt",zorder=3)
        ax.plot([der,der],[0.20,0.52],color=cder,lw=6,solid_capstyle="butt",zorder=3)
        ax.annotate("",xy=(der,0.54),xytext=(src,0.66),
                    arrowprops=dict(arrowstyle="-|>",color=MUTED,lw=1.6,
                                    connectionstyle="arc3,rad=0.35"),zorder=2)
    else:
        ax.plot([0.30,0.30],[0.20,0.62],color=COAST,lw=6,solid_capstyle="butt",zorder=3)
        ax.plot([0.70,0.70],[0.20,0.62],color=AND,lw=6,solid_capstyle="butt",zorder=3)
        ax.plot([0.30,0.50,0.70],[0.62,0.80,0.62],color=MUTED,lw=3,zorder=2)
        if j==3:
            for yy in (0.32,0.42,0.52):
                ax.annotate("",xy=(0.665,yy),xytext=(0.335,yy),
                            arrowprops=dict(arrowstyle="<|-|>",color=COAST,lw=1.3,
                                            mutation_scale=9),zorder=4)
    ax.text(0.5,0.10,"Coastal          Andean",ha="center",fontsize=8,color=MUTED)
    ax.axhline(-0.03,color=RULE,lw=0.9)
    ax.set_title(NAMES[j],loc="center",fontsize=10.5,color=INK,pad=8,linespacing=1.3)
    # the vote share, drawn as a bar under the schematic
    ax.add_patch(plt.Rectangle((0.08,-0.20),0.84,0.085,facecolor=GRAY,alpha=0.30,zorder=2))
    ax.add_patch(plt.Rectangle((0.08,-0.20),0.84*votes[j]/0.5,0.085,
                               facecolor=GREEN if j>=2 else MUTED,alpha=0.95,zorder=3))
    ax.text(0.08,-0.30,"%.1f%% of votes"%(100*votes[j]),fontsize=10,color=INK,
            va="top",fontweight="bold" if j>=2 else "normal")
    ax.text(0.92,-0.30,"recovered %.0f%%"%(100*recall[j]),fontsize=8.4,color=MUTED,
            ha="right",va="top")
fig.text(0.055,0.972,"A",fontsize=15,fontweight="bold",color=INK,va="top")

# ================= B: what is resolved, and what is not =================
axB=fig.add_subplot(gsB[0,0]); style(axB)
rows=[(D["stage1"]["name"],D["stage1"]["posterior"],D["stage1"]["oob"],
       D["stage1"]["selected"],True),
      (D["stage2"]["name"],D["stage2"]["posterior"],D["stage2"]["oob"],
       D["stage2"]["selected"],False)]
yp=[1.40,-0.30]
for (lab,post,oob,selname,ok),y in zip(rows,yp):
    col=GREEN if ok else AND
    axB.barh(y,post,height=0.42,color=col,alpha=0.90 if ok else 0.35,
             edgecolor=col,linewidth=1.6,zorder=3)
    axB.text(post+0.03,y,"%.3f"%post,va="center",fontsize=13.5,fontweight="bold",color=INK)
    axB.text(0.010,y+0.42,lab,fontsize=11,color=INK,fontweight="bold")
    axB.text(0.985,y+0.42,"prior error rate %.1f%%"%(100*(1-oob)),ha="right",
             fontsize=8.8,color=MUTED)
    axB.text(0.010,y-0.44,("favored: %s" if ok else "%s")%selname,fontsize=9.2,color=col,va="top")
axB.axvline(0.5,color=INK,lw=1.2,ls=":",zorder=4)
axB.set_xlim(0,1.0); axB.set_ylim(-1.30,2.20); axB.set_yticks([])
axB.set_xlabel("posterior probability (upper: two refugia; lower: contact)")
axB.spines["left"].set_visible(False)
fig.text(0.075,0.620,"B",fontsize=15,fontweight="bold",color=INK,va="top")

axB2=fig.add_subplot(gsB[1,0]); style(axB2)
pw=D["migration"]["power"]
xc=[0.1,0.5,2.5,8.0]; yv=[p[2]*100 for p in pw]; nn=[p[1] for p in pw]
axB2.plot(xc,yv,color=COAST,lw=2.4,marker="o",ms=7,mfc=PANEL,mew=1.6,zorder=4)
axB2.set_xscale("log")
axB2.axhline(95,color=MUTED,lw=1.0,ls="--",zorder=2)
axB2.text(0.055,96,"95%",fontsize=8.4,color=MUTED)
# posterior median under contact; its 95% interval (0.005 to 26) spans the whole axis, so it is
# stated in the legend document rather than drawn as a band
obsNm=D["migration"]["Nm_median"]
axB2.axvline(obsNm,color=AND,lw=2.2,zorder=5)
axB2.text(obsNm*1.12,30,"posterior median %.2f"%obsNm,
          fontsize=9,color=AND,fontweight="bold",va="top")
for x,y,n in zip(xc,yv,nn):
    axB2.text(x*(1.18 if x<0.2 else 1),y+(-1.5 if x<0.2 else 3.5),"%.0f%%"%y,ha="left" if x<0.2 else "center",fontsize=8.6,color=INK)
axB2.set_xlim(0.04,14); axB2.set_ylim(25,108)
axB2.set_xticks([0.1,1,10]); axB2.set_xticklabels(["0.1","1","10"])
axB2.set_xlabel("true gene flow across the gap (migrants per generation)")
axB2.set_ylabel("% of simulations\ncorrectly called",fontsize=9.5,linespacing=1.3)

# ================= C: how old is the split =================
axC=fig.add_subplot(gsC[0,0]); style(axC)
grid=np.linspace(3,6.6,1200)
def dens(tag):
    W=pd.read_csv(f"figdata/abcrf_weights_{tag}.csv"); lT,w=np.log10(W.value.values),W.weight.values
    d=(w[None,:]*np.exp(-0.5*((grid[:,None]-lT[None,:])/0.12)**2)).sum(1)
    return d/np.trapezoid(d,grid)
LGM=(19e3,26.5e3); LIG=(115e3,130e3)
axC.axvspan(*np.log10(LGM),color=AND,alpha=0.16,zorder=1,lw=0,label="Last Glacial Maximum")
axC.axvspan(*np.log10(LIG),color=GREEN,alpha=0.16,zorder=1,lw=0,label="last interglacial")
da=dens("avg")
axC.fill_between(grid,0,da,color=COAST,alpha=0.22,zorder=2)
axC.plot(grid,da,color=COAST,lw=2.8,zorder=5,label="averaged over both histories")
axC.plot(grid,dens("iso"),color=INK,lw=1.3,zorder=4,label="assuming isolation")
axC.plot(grid,dens("con"),color=INK,lw=1.3,ls=(0,(4,2)),zorder=4,label="assuming contact")
axC.plot(grid,dens("broad"),color=MUTED,lw=1.6,ls=":",zorder=3,label="broad mutation-rate prior")
med=np.log10(A["split_model_averaged"]["median"])
axC.axvline(med,color=COAST,lw=1.2,ls="--",zorder=6)
ticks=[1e3,1e4,1e5,1e6]; axC.set_xticks(np.log10(ticks))
axC.set_xticklabels(["1,000","10,000","100,000","1,000,000"])
axC.set_xticks(np.log10([v*k for v in ticks[:-1] for k in range(2,10)]),minor=True)
axC.set_xlim(3,6.6); axC.set_ylim(0,None); axC.set_yticks([])
axC.spines["left"].set_visible(False)
axC.set_xlabel("age of the Coastal-Andean split (years before present)")
axC.set_ylabel("posterior density")
axC.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=8.6,loc="upper right")
fig.text(0.585,0.620,"C",fontsize=15,fontweight="bold",color=INK,va="top")

for ext in ("png","pdf"):
    fig.savefig(f"figures/journal/F3_demographic_history.{ext}",dpi=200,facecolor=BG)
print("F3 rebuilt from computed values")
print("  votes %s | stage1 %.3f/%.3f | stage2 %.3f/%.3f"
      %(["%.3f"%v for v in votes],D["stage1"]["posterior"],D["stage1"]["oob"],
        D["stage2"]["posterior"],D["stage2"]["oob"]))
print("  split median {:,.0f} yr".format(10**med))
