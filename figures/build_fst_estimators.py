#!/usr/bin/env python3
"""Why the two F_ST values differed, and why Hudson is the one to report.

Data: figdata/fst_estimator_control.tsv — all three quantities computed on ONE
identical SNP set, ratio-of-averages, 66 stand pairs (cluster job 41051).
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'

r=list(csv.DictReader(open('figdata/fst_estimator_control.tsv'),delimiter='\t'))
hud=np.array([float(x['hudson']) for x in r])
obs=np.array([float(x['wc_obs'])  for x in r])
hwe=np.array([float(x['wc_hwe'])  for x in r])
btw=np.array([int(x['between'])   for x in r],bool)

fig,(axA,axB)=plt.subplots(1,2,figsize=(12.4,5.4),gridspec_kw={'width_ratios':[1.1,1]})
fig.patch.set_facecolor(PAGE)
for a in (axA,axB):
    a.set_facecolor(PANEL)
    for s in ('top','right'): a.spines[s].set_visible(False)
    for s in ('left','bottom'): a.spines[s].set_color(RULE)
    a.tick_params(colors=MUTED,labelsize=9); a.grid(alpha=.25,color=RULE,lw=.7)

# --- A: each W&C variant against Hudson, same pairs, same sites
lim=[0,0.23]
axA.plot(lim,lim,color=MUTED,lw=1,ls='--',zorder=1)
axA.scatter(hud,obs,s=34,color=BLUE,alpha=.85,zorder=3,
            label='W&C using OBSERVED heterozygosity')
axA.scatter(hud,hwe,s=34,color=GREEN,alpha=.85,zorder=3,marker='^',
            label='W&C with $h$ set to its HWE expectation')
axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel('Hudson $F_{ST}$',color=INK,fontsize=10.5)
axA.set_ylabel('Weir–Cockerham $F_{ST}$',color=INK,fontsize=10.5)
axA.set_title('One SNP set, 66 stand pairs',color=INK,fontsize=11.5,loc='left',pad=9)
axA.legend(frameon=False,fontsize=9,labelcolor=MUTED,loc='upper left')
axA.annotate('correcting only the heterozygosity term\nrecovers 83% of the gap',
             xy=(0.155,0.145),xytext=(0.085,0.045),fontsize=9,color=GREEN,
             arrowprops=dict(arrowstyle='->',color=GREEN,lw=1.1))

# --- B: the quantity the manuscript actually leans on
labs=['Hudson','W&C\n(observed $h$)','W&C\n($h$ at HWE)']
bw=[hud[btw].mean(),obs[btw].mean(),hwe[btw].mean()]
wi=[hud[~btw].mean(),obs[~btw].mean(),hwe[~btw].mean()]
x=np.arange(3); w=0.36
axB.bar(x-w/2,bw,w,color=RED,label='between region (n=27)')
axB.bar(x+w/2,wi,w,color=BLUE,label='within region (n=39)')
for i,(b,v) in enumerate(zip(bw,wi)):
    axB.text(i,max(b,v)+0.007,'%.2f×'%(b/v),ha='center',fontsize=11,
             color=INK,fontweight='bold')
axB.set_xticks(x); axB.set_xticklabels(labs,color=MUTED,fontsize=9.5)
axB.set_ylabel('mean pairwise $F_{ST}$',color=INK,fontsize=10.5)
axB.set_ylim(0,0.20)
axB.set_title('Between- vs within-region differentiation',color=INK,fontsize=11.5,
              loc='left',pad=9)
axB.legend(frameon=False,fontsize=9,labelcolor=MUTED)

fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'fst_estimators.{ext}',dpi=200,facecolor=PAGE)
print('hudson %.4f | wc_obs %.4f | wc_hwe %.4f'%(hud.mean(),obs.mean(),hwe.mean()))
print('between/within  hudson %.2fx  wc_obs %.2fx  wc_hwe %.2fx'
      %(bw[0]/wi[0],bw[1]/wi[1],bw[2]/wi[2]))
print('gap recovered by the HWE correction: %.1f%%'
      %(100*(hwe.mean()-obs.mean())/(hud.mean()-obs.mean())))
print('saved fst_estimators.png/.pdf')
