#!/usr/bin/env python3
"""F1: distribution, sampling, and the environmental setting.

JOURNAL STYLE (4 Sep 2026): no figure title, no subtitle, no on-canvas footnote and no
interpretive phrasing. Panels carry bare letters; every word that argued a point now lives in
ARAUCARIA_figure_legends.docx. Text kept on the canvas is restricted to axis labels, units,
legends, map furniture (scale bar, north arrow) and data values.

CHANGED 10 Sep 2026: the fungal variance partition (former panel C) is removed -- it is a
microbiome result in a population-genetics figure and its terms came from three different
analyses. Panel B now maps ALL SEVEN variables of the final literature set rather than one
per "mechanism"; showing four left three model variables invisible to the reader."""
import numpy as np, pickle, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle

BG,PANEL="#ffffff","#ffffff";  # journal: plain white (was cream #f6f4ec/#fffdf7, 17 Nov 2026)
AND,COAST="#a8543a","#35708c"
RULE,MUTED,INK,GREEN="#dcd8ca","#6c6b5e","#14180f","#40663f"
LAND="#e4e0d3"
mpl.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
    "axes.edgecolor":RULE,"axes.linewidth":0.9,"text.color":INK,"axes.labelcolor":INK,
    "xtick.color":MUTED,"ytick.color":MUTED,"xtick.labelsize":9,"ytick.labelsize":9,
    "pdf.fonttype":42})

D=pickle.load(open("figdata/distribution_cache.pkl","rb"))
S=np.load("figdata/stand_coords.npz",allow_pickle=True)
SC,ST,CO=S["coords"],[str(x) for x in S["stands"]],S["coastal"]
NPER=dict(zip(ST,[16,11,16,18,11,12,9,16,16,11,14,13]))

fig=plt.figure(figsize=(11.6,13.0)); fig.patch.set_facecolor(BG)
gs=GridSpec(1,1,figure=fig,left=0.085,right=0.965,top=0.970,bottom=0.500)
ax=fig.add_subplot(gs[0,0]); ax.set_facecolor(PANEL)
ax.set_xlim(-74.35,-70.85); ax.set_ylim(-39.95,-37.05)
ax.imshow(D["grid"],extent=[D["gx"][0],D["gx"][-1],D["gy"][0],D["gy"][-1]],
          origin="lower",cmap=mpl.colors.ListedColormap([PANEL,LAND]),interpolation="nearest",zorder=1)
ax.add_collection(PolyCollection(D["rings_extra"],facecolors="#8fae86",edgecolors="none",zorder=2))
ax.add_collection(PolyCollection(D["rings_dom"],facecolors=GREEN,edgecolors="none",zorder=3))
for j,s in enumerate(ST):
    c=COAST if CO[j] else AND
    ax.plot(SC[j,0],SC[j,1],"o",ms=11,color=c,markeredgecolor="white",markeredgewidth=1.6,zorder=6)
off={"TR":(-0.20,0.06),"PN":(-0.20,-0.06),"VA":(0.16,0.05),"LM":(0.16,0.05),"RA":(0.17,-0.05),
     "TH":(0.16,0.02),"MC":(-0.19,0.06),"LR":(0.17,-0.05),"PC":(-0.19,0.0),"HL":(0.16,0.02),
     "BP":(-0.19,0.02),"MM":(0.16,-0.02)}
for j,s in enumerate(ST):
    dx,dy=off.get(s,(0.15,0.05))
    ax.annotate(s,(SC[j,0],SC[j,1]),xytext=(SC[j,0]+dx,SC[j,1]+dy),fontsize=11,fontweight="bold",
                ha="center",va="center",color=INK,zorder=7,
                bbox=dict(fc=PANEL,ec=RULE,lw=0.7,boxstyle="round,pad=0.22"))
# the narrowest crossing, drawn between the actual closest cross-range stand pair (PN-TH)
iPN,iTH=ST.index("PN"),ST.index("TH")
p,t=SC[iPN],SC[iTH]
def _km(a,b):
    R=6371.0; p1,p2=np.radians(a[1]),np.radians(b[1])
    h=np.sin((p2-p1)/2)**2+np.cos(p1)*np.cos(p2)*np.sin(np.radians(b[0]-a[0])/2)**2
    return 2*R*np.arcsin(np.sqrt(h))
gap=_km(p,t)
ax.annotate("",xy=(t[0],t[1]),xytext=(p[0],p[1]),
            arrowprops=dict(arrowstyle="<|-|>",color=INK,lw=2.0,mutation_scale=15,
                            shrinkA=9,shrinkB=9),zorder=8)
mx,my=(p[0]+t[0])/2,(p[1]+t[1])/2
ax.text(mx,my+0.075,f"{gap:.0f} km",fontsize=11.5,ha="center",color=INK,fontweight="bold",zorder=9,
        bbox=dict(fc=PANEL,ec=INK,lw=1.0,boxstyle="round,pad=0.26"))
# The arrow measures the closest SAMPLED cross-range pair, not the gap between the
# ranges themselves -- those are 83.8 km (Araucaria-dominant polygons) and 54.8 km
# (including component stands), both well under 117 km. Both were computed twice:
# on every polygon vertex (2.26 M Andean) and on a 1-in-6 thinned set (380 k), which
# agree to 0.1 km and return the same closest-approach coordinates. The old label claimed "no
# Araucaria in the Central Valley", which a reader would attach to the 117 km and
# read as the range gap. The distributional numbers belong in the caption, beside
# the class definitions they depend on.
# The arrow slopes down to the right, so the label sits clear below it with the
# land masked behind.
# scale bar
y0=-37.25; x0=-74.15; deg=50/(111.32*np.cos(np.radians(37.3)))
ax.plot([x0,x0+deg],[y0,y0],color=INK,lw=2.6,solid_capstyle="butt",zorder=9)
for xx in (x0,x0+deg): ax.plot([xx,xx],[y0-0.035,y0+0.035],color=INK,lw=2.0,zorder=9)
ax.text(x0+deg/2,y0+0.07,"50 km",fontsize=10,ha="center",fontweight="bold",color=INK,zorder=9)
ax.annotate("",xy=(-71.05,-37.18),xytext=(-71.05,-37.42),
            arrowprops=dict(arrowstyle="-|>",color=INK,lw=2.0,mutation_scale=16),zorder=9)
ax.text(-71.05,-37.50,"N",fontsize=11,ha="center",fontweight="bold",color=INK,zorder=9)
ax.set_xticks([-74,-73,-72,-71]); ax.set_xticklabels(["74°W","73°W","72°W","71°W"])
ax.set_yticks([-39.5,-39,-38.5,-38,-37.5]); ax.set_yticklabels(["39.5°S","39°S","38.5°S","38°S","37.5°S"])
ax.grid(color=RULE,lw=0.6,alpha=.55,zorder=4)
for s in ("top","right"): ax.spines[s].set_visible(False)
h=[plt.Line2D([],[],marker="s",ls="",ms=11,color=GREEN,label="Araucaria-dominant forest"),
   plt.Line2D([],[],marker="s",ls="",ms=11,color="#8fae86",label="Araucaria present as a component"),
   plt.Line2D([],[],marker="o",ls="",ms=10,color=COAST,markeredgecolor="white",label="Coastal stand (3)"),
   plt.Line2D([],[],marker="o",ls="",ms=10,color=AND,markeredgecolor="white",label="Andean stand (9)")]
ax.legend(handles=h,frameon=True,facecolor=PANEL,edgecolor=RULE,fontsize=10,loc="lower left",
          borderpad=0.8,labelspacing=0.7)
# VA inset. imshow forces aspect="equal", so the map axes shrinks and centers
# inside its GridSpec cell -- its drawn left edge is well right of the cell's
# left. Anchoring the inset to figure coordinates put it on top of the latitude
# labels, so take the position from the axes after a draw.
fig.canvas.draw()
_p=ax.get_position()
axi=fig.add_axes([_p.x0+0.022,_p.y0+0.095,0.205,0.083]); axi.set_facecolor(PANEL)
va=SC[ST.index("VA")]
axi.set_xlim(va[0]-0.085,va[0]+0.085); axi.set_ylim(va[1]-0.055,va[1]+0.055)
axi.add_collection(PolyCollection(D["rings_dom"],facecolors=GREEN,edgecolors="none"))
axi.plot(va[0],va[1],"o",ms=10,color=COAST,markeredgecolor="white",markeredgewidth=1.5)
axi.set_xticks([]); axi.set_yticks([])
for s in axi.spines.values(): s.set_edgecolor(AND); s.set_linewidth(1.8)
axi.set_title("VA",fontsize=9,color=MUTED,pad=4)
d2=2/(111.32*np.cos(np.radians(38.5))); xi=va[0]-0.070
axi.plot([xi,xi+d2],[va[1]-0.045,va[1]-0.045],color=INK,lw=2.2)
axi.text(xi+d2/2,va[1]-0.040,"2 km",fontsize=8.5,ha="center",color=INK,fontweight="bold")
ax.add_patch(Rectangle((va[0]-0.085,va[1]-0.055),0.17,0.11,fill=False,ec=AND,lw=1.6,zorder=8))

cov=D["cov"]["dominant"]

# ============================ climate panels ============================
# CHELSA v2.1 windows over the same extent as the map above
# (environment/pull_chelsa_rasters.py). CHELSA carries interpolated values over
# the ocean -- bio7 reads ~5 degC there against 19-31 inland -- which would take
# over the color scale, so cells are masked to the same land grid the main map
# draws and limits are taken from land only.
import csv as _csv
CH=np.load("figdata/chelsa_surfaces.npz")
CL=np.load("figdata/soilgrids_soil.npz")
COASTAL_CODES={"PN","TR","VA"}

def _standmeans(path,key,var):
    """Stand means, d (range difference / SD across the 12 stand means) between ranges, and whether the ranges overlap.
    Computed here rather than hardcoded so the panel labels cannot drift from
    the data."""
    per={}
    for r in _csv.DictReader(open(path),delimiter="\t"):
        try: x=float(r[var])
        except (ValueError,TypeError): continue
        per.setdefault(r[key],[]).append(x)
    m={k:sum(v)/len(v) for k,v in per.items()}
    co=[v for k,v in m.items() if k in COASTAL_CODES]
    an=[v for k,v in m.items() if k not in COASTAL_CODES]
    pooled=np.std(list(m.values()))
    d=(np.mean(co)-np.mean(an))/pooled if pooled else 0.0
    return d, (max(co)<min(an) or max(an)<min(co))

_gx,_gy,_grid=D["gx"],D["gy"],D["grid"]
def _mask(a,ext):
    lon=np.linspace(ext[0],ext[1],a.shape[1]); lat=np.linspace(ext[3],ext[2],a.shape[0])
    ix=np.clip(np.searchsorted(_gx,lon)-1,0,len(_gx)-1)
    iy=np.clip(np.searchsorted(_gy,lat)-1,0,len(_gy)-1)
    out=a.copy(); out[~_grid[np.ix_(iy,ix)]]=np.nan; return out

CHE="environment/chelsa_per_tree.tsv"; SOI="environment/soil_per_stand.tsv"
# ordered by |d| over all 39 variables tested (climate + terrain + soil)
# The FOUR mechanisms of the final 7-variable set, one panel each. bio4, bio7 and
# bio12 were shown in earlier versions but are NOT in the model -- they were
# dropped as redundant (|r| = 0.80, 0.80, 0.96 with bio6/bio18). CEC, pH and
# nitrogen are in the model but are secondary soil terms and are not mapped.
# All seven variables of the final literature set. bio4, bio7 and bio12 appeared in earlier
# versions but are NOT in the model -- dropped as redundant (|r| = 0.80, 0.80, 0.96 with
# bio6/bio18), as was elevation (|r| = 0.91 with bio6).
CLIMATE=[(CH["bio5"], CH["bio5_extent"], "max temperature, warmest month","\u00b0C","OrRd", CHE,"bio5"),
         (CH["bio18"],CH["bio18_extent"],"precipitation, warmest quarter","mm","GnBu",      CHE,"bio18"),
         (CH["bio6"], CH["bio6_extent"], "min temperature, coldest month","\u00b0C","OrRd", CHE,"bio6")]
SOIL=[(CL["clay"], CL["clay_extent"], "clay content","%","YlOrBr",                        SOI,"clay_0-5cm"),
      (CL["cec"],  CL["cec_extent"],  "cation exchange capacity","cmol(c)/kg","PuBuGn",   SOI,"cec_0-5cm"),
      (CL["ph"],   CL["ph_extent"],   "pH (H\u2082O)","pH","RdYlBu",                       SOI,"phh2o_0-5cm"),
      (CL["nitro"],CL["nitro_extent"],"total nitrogen","g/kg","YlGn",                     SOI,"nitrogen_0-5cm")]

def draw_row(gsrow,items):
    for j,(arr,ext,lab,unit,cmap,src,var) in enumerate(items):
        axc=fig.add_subplot(gsrow[0,j]); axc.set_facecolor(PANEL)
        A=_mask(arr,ext)
        cm=mpl.colormaps[cmap].copy(); cm.set_bad(PANEL)
        lo,hi=np.nanpercentile(A,2),np.nanpercentile(A,98)
        im=axc.imshow(A,extent=list(ext),origin="upper",cmap=cm,vmin=lo,vmax=hi,
                      interpolation="nearest",zorder=1)
        axc.set_xlim(-74.35,-70.85); axc.set_ylim(-39.95,-37.05)
        for xy,c in zip(SC,CO):
            axc.plot(xy[0],xy[1],"o" if c else "D",ms=6.0 if c else 5.4,
                     mfc=(COAST if c else AND),mec=INK,mew=1.0,zorder=4)
        for sp in ("top","right"): axc.spines[sp].set_visible(False)
        for sp in ("left","bottom"): axc.spines[sp].set_color(RULE)
        axc.set_xticks([]); axc.set_yticks([])
        d,_sep=_standmeans(src,"site_code",var)
        axc.set_title(lab,loc="left",fontsize=9.0,color=INK,pad=5)
        axc.text(0.0,-0.075,"$d$ = %+.2f"%d,transform=axc.transAxes,fontsize=8.2,color=MUTED)
        cb=fig.colorbar(im,ax=axc,orientation="horizontal",fraction=0.055,pad=0.10,aspect=18)
        cb.outline.set_edgecolor(RULE); cb.ax.tick_params(labelsize=6.8,colors=MUTED,length=2)
        cb.locator=mpl.ticker.MaxNLocator(4); cb.update_ticks()
        if unit: cb.set_label(unit,fontsize=7.2,color=MUTED,labelpad=1)

# climate on the upper row (3 of 4 cells), soil on the lower row (4 of 4), so the two
# groups of the model read as groups and every panel keeps the same width
gsB1=GridSpec(1,4,figure=fig,left=0.085,right=0.965,bottom=0.300,top=0.462,wspace=0.16)
gsB2=GridSpec(1,4,figure=fig,left=0.085,right=0.965,bottom=0.098,top=0.260,wspace=0.16)
draw_row(gsB1,CLIMATE)
draw_row(gsB2,SOIL)

fig.text(0.085,0.978,"A",fontsize=15,fontweight="bold",color=INK,ha="left",va="top")
fig.text(0.085,0.487,"B",fontsize=15,fontweight="bold",color=INK,ha="left",va="top")
fig.savefig("figures/journal/F1_distribution_sampling.png",dpi=300,facecolor=BG)
fig.savefig("figures/journal/F1_distribution_sampling.pdf",facecolor=BG)
print(f"F1 rebuilt: {D['n_dom']:,} dominant polygons, {D['ha_dom']:,.0f} ha, coverage {cov['pct50']:.1f}%")
print(f"  narrowest cross-range gap drawn: PN-TH = {gap:.1f} km")
