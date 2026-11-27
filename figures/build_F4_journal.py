#!/usr/bin/env python3
"""Figure 4: Q_ST across heritabilities against F_ST, and per-trait crossover heritability.

Panel B reports each trait's CROSSOVER heritability rather than its Q_ST at some assumed h2:
the crossover is the value at which that trait would stop exceeding F_ST, so nothing has to be
assumed to read the panel.

JOURNAL STYLE (4 Sep 2026): no title, subtitle, panel headlines or on-canvas conclusions.
Interpretation lives in ARAUCARIA_figure_legends.docx.
"""
import numpy as np, json, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BG,PANEL="#ffffff","#ffffff";  # journal: plain white (was cream #f6f4ec/#fffdf7, 17 Nov 2026)
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,"axes.labelcolor":INK,
    "xtick.color":MUTED,"ytick.color":MUTED,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.labelsize":11,"pdf.fonttype":42})

# McIntosh et al. 2025 Table S4: trait, MS_region, MS_pop(region), MS_residual
S4=[("Needle width",26.7,5.27,148.5/192),("Needle area",15.6,0.32,43.6/194),
    ("Needle succulence",11.0,0.86,76.0/194),("Number of whorls",9.5,0.80,66.9/193),
    ("Needle thickness",0.5,0.06,3.7/194),("Needle length",14.9,0.31,81.2/192),
    ("Number of branches",24.7,2.67,237.5/193),("Needle mass/area",0.6,0.21,15.3/194),
    ("Branch length",133.2,22.24,2268.7/192),("Percent germination",0.021,0.025,2.4/192)]
d=json.load(open("figdata/qst_fst.json")); FST=d["fst_regional"]; DIV=d["n0"]*d["k0"]
PLAUS=(0.15,0.35)          # conifer seedling traits, literature range

def qst(a,b,c,h): 
    v=max((a-b)/DIV,0.0); return v/(v+2*h*c)
def crossover(a,b,c):
    v=max((a-b)/DIV,0.0)
    return v*(1-FST)/(2*FST*c) if c>0 and v>0 else 0.0

fig=plt.figure(figsize=(13.4,5.9)); fig.patch.set_facecolor(BG)
gs=GridSpec(1,2,figure=fig,wspace=0.30,left=0.075,right=0.975,top=0.925,bottom=0.135)
def style(ax):
    ax.set_facecolor(PANEL)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3,width=0.9)
def tag(ax,t):
    ax.set_title(t,loc="left",fontsize=14,fontweight="bold",color=INK,pad=9)

# ---------- A: Q_ST across every possible heritability ----------
axA=fig.add_subplot(gs[0,0]); style(axA)
hs=np.linspace(0.03,1.0,300)
mq=np.array([np.mean([qst(a,b,c,h) for _,a,b,c in S4]) for h in hs])
axA.fill_between(hs,mq,FST,where=mq>FST,color=AND,alpha=0.14,zorder=1)
axA.axvspan(*PLAUS,color=GREEN,alpha=0.10,zorder=0,
            label=r"plausible $h^2$, conifer seedlings")
axA.plot(hs,mq,color=AND,lw=2.8,zorder=4,label=r"$Q_{ST}$  Coastal vs Andes traits")
axA.axhline(FST,color=COAST,lw=2.0,ls="--",zorder=3,
            label=rf"$F_{{ST}}$ = {FST:.3f}  (RAD-seq, Hudson)")
xc=d["crossover_h2"]
axA.plot([xc],[FST],"o",ms=8,color=INK,zorder=5)
axA.annotate(f"$h^2$ = {xc:.3f}",xy=(xc,FST),xytext=(0.60,0.30),
             fontsize=10,color=INK,arrowprops=dict(arrowstyle="->",color=INK,lw=1.1))
axA.set_xlim(0,1.0); axA.set_ylim(0,0.60)
axA.set_xlabel(r"assumed heritability  $h^2$"); axA.set_ylabel(r"$Q_{ST}$   /   $F_{ST}$")
axA.legend(frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=9.5,loc="upper right")
tag(axA,"A")

# ---------- B: per-trait crossover, nothing assumed ----------
axB=fig.add_subplot(gs[0,1]); style(axB)
tr=sorted(((crossover(a,b,c),nm) for nm,a,b,c in S4),reverse=True)
XMAX=1.5
axB.axvspan(*PLAUS,color=GREEN,alpha=0.10,zorder=0)
axB.axvline(1.0,color=INK,lw=1.4,ls=":",zorder=3)
for i,(h,nm) in enumerate(tr):
    col=AND if h>1 else (COAST if h>PLAUS[1] else MUTED)
    axB.barh(i,min(h,XMAX),height=0.62,color=col,alpha=0.92,zorder=2)
    if h>XMAX:
        axB.text(XMAX-0.02,i,f"{h:.2f} >",ha="right",va="center",fontsize=9,
                 color="white",fontweight="bold",zorder=4)
    elif h>0:
        if 0.85<h<1.0:   # would collide with the h2 = 1 rule; set it inside the bar
            axB.text(h-0.02,i,f"{h:.2f}",ha="right",va="center",fontsize=9,
                     color="white",fontweight="bold")
        else:
            axB.text(min(h,XMAX)+0.02,i,f"{h:.2f}",va="center",fontsize=9,color=MUTED)
    else:
        axB.text(0.02,i,"no divergence",va="center",fontsize=9,color=MUTED,style="italic")
axB.set_yticks(range(len(tr))); axB.set_yticklabels([n for _,n in tr],fontsize=10)
axB.invert_yaxis(); axB.set_xlim(0,XMAX)
axB.set_xlabel(r"crossover $h^2$")
tag(axB,"B")

OUT="figures/journal/F4_qst_fst_selection"
for ext in ("png","pdf"):
    fig.savefig(f"{OUT}.{ext}",dpi=200,facecolor=BG)
n_any=sum(1 for h,_ in tr if h>1); n_pl=sum(1 for h,_ in tr if h>PLAUS[1])
print(f"traits exceeding F_ST at ANY h2 : {n_any}/10")
print(f"traits with crossover above the plausible range (>{PLAUS[1]}) : {n_pl}/10")
print(f"saved {OUT}.png/.pdf")
