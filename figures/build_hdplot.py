#!/usr/bin/env python3
"""HDplot (McKinney et al. 2017) for the ARAUCARIA RAD loci."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'

H=[]; D=[]
for r in csv.DictReader(open('figdata/hdplot.tsv'), delimiter='\t'):
    try: h=float(r['H']); d=float(r['D'])
    except ValueError: continue
    H.append(h); D.append(d)
H=np.array(H); D=np.array(D); n=len(H)
HT, DT = 0.55, 7.0
para=(H>HT)&(np.abs(D)>DT)

fig=plt.figure(figsize=(12.2,5.9))
fig.patch.set_facecolor(PAGE)
ax =fig.add_axes([0.070,0.215,0.520,0.640])
cax=fig.add_axes([0.605,0.215,0.013,0.640])
ax2=fig.add_axes([0.715,0.215,0.265,0.640])
for a in (ax,ax2):
    a.set_facecolor(PANEL)
    for s in ('top','right'): a.spines[s].set_visible(False)
    for s in ('left','bottom'): a.spines[s].set_color(RULE)
    a.tick_params(colors=MUTED, labelsize=9)

ylo,yhi=np.percentile(D,0.2), np.percentile(D,99.8)
hb=ax.hexbin(H,D,gridsize=64,bins='log',mincnt=1,cmap='YlGnBu',linewidths=0,
             extent=(0,1,ylo,yhi))
cb=fig.colorbar(hb,cax=cax)
cb.set_label('loci per bin (log)',color=MUTED,fontsize=8.6)
cb.ax.tick_params(colors=MUTED,labelsize=8); cb.outline.set_edgecolor(RULE)

ax.axhline(0,color=INK,lw=.9,alpha=.5)
ax.axvline(HT,color=RED,lw=1.3,ls='--')
ax.hlines([DT,-DT],HT,1.0,color=RED,lw=1.3,ls='--')
ax.add_patch(plt.Rectangle((HT,DT),1-HT,yhi-DT,facecolor=RED,alpha=.11,ec='none'))
ax.add_patch(plt.Rectangle((HT,ylo),1-HT,-DT-ylo,facecolor=RED,alpha=.11,ec='none'))
ax.text(0.575,yhi*0.72,'paralog rejection\nregion\nH>%.2f, |D|>%g\n\n%d loci'%(HT,DT,para.sum()),
        color=RED,fontsize=8.8,va='top',linespacing=1.45)
ax.set_xlim(0,1); ax.set_ylim(ylo,yhi)
ax.set_xlabel('H  —  proportion of individuals called heterozygous',color=INK,fontsize=10)
ax.set_ylabel('D  —  read-ratio deviation',color=INK,fontsize=10)
ax.set_title('HDplot — ARAUCARIA RAD loci  (%s loci, ct091 cr30)'%format(n,','),
             color=INK,fontsize=11.5,loc='left',pad=9)

cnt,edges,_=ax2.hist(H,bins=60,color=BLUE,alpha=.85,edgecolor=PANEL,lw=.4)
ax2.set_yscale('log'); ax2.set_xlim(0,1)
top=cnt.max()
ax2.axvline(0.5,color=GREEN,lw=1.1,ls=':')
ax2.axvline(HT,color=RED,lw=1.3,ls='--')
ax2.text(0.475,top*0.30,'H = 0.5',color=GREEN,fontsize=8.4,ha='right')
ax2.text(0.575,top*0.030,'McKinney\ncut-off',color=RED,fontsize=8.4,linespacing=1.3)
ax2.set_xlabel('H',color=INK,fontsize=10)
ax2.set_ylabel('loci (log scale)',color=INK,fontsize=10)
ax2.set_title('distribution of H',color=INK,fontsize=10.5,loc='left',pad=9)

fig.text(0.070,0.075,
 'Collapsed paralogs make every individual appear heterozygous, so they sit at HIGH H with LARGE |D|. '
 'Zero loci fall in the rejection region.',color=MUTED,fontsize=8.2)
fig.text(0.070,0.046,
 'Loci with H > 0.40: %d      H > 0.50: %d      H > 0.55: %d      H > 0.60: %d'
 %((H>0.4).sum(),(H>0.5).sum(),(H>0.55).sum(),(H>0.6).sum()),color=MUTED,fontsize=8.2)
fig.text(0.070,0.017,
 'The classic fixed-difference paralog signature is absent. The cut-off is strict, so milder '
 'heterozygote excess is NOT excluded by this test.',color=MUTED,fontsize=8.2)
fig.savefig('HDplot.png',dpi=200,facecolor=PAGE)
fig.savefig('HDplot.pdf',facecolor=PAGE)
print('loci %d | rejection region %d | D range %.1f..%.1f (plotted %.1f..%.1f)'
      %(n,para.sum(),D.min(),D.max(),ylo,yhi))
