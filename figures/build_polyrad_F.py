#!/usr/bin/env python3
"""F_IS by depth: ipyrad binomial calls vs polyRAD beta-binomial (d=5)."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'

NAMES=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
MID  =[8,12.5,17.5,25,40,75,150]

pr={r['depth']:r for r in csv.DictReader(open('figdata/polyrad_F_by_depth.csv'))}
Fp=[float(pr[n]['F_IS']) for n in NAMES]
Hop=[float(pr[n]['Ho'])  for n in NAMES]
Hep=[float(pr[n]['He'])  for n in NAMES]
Np =[int(pr[n]['n_geno'])for n in NAMES]
# ipyrad binomial baseline: SAME bins, SAME within-stand method,
# SAME locus-level read-based allele frequency as the polyRAD curve
Fi=[0.8694,0.8153,0.7847,0.6512,0.3061,-0.1960,-0.3363]

fig,(ax,ax2)=plt.subplots(2,1,figsize=(8.6,8.4),sharex=True,
                          gridspec_kw={'height_ratios':[2.05,1]})
fig.patch.set_facecolor(PAGE)
for a in (ax,ax2):
    a.set_facecolor(PANEL); a.set_xscale('log')
    for s in ('top','right'): a.spines[s].set_visible(False)
    for s in ('left','bottom'): a.spines[s].set_color(RULE)
    a.tick_params(colors=MUTED, labelsize=9)
    a.grid(alpha=.25, color=RULE, lw=.7)

ax.axhline(0, color=INK, lw=1.1, zorder=3)
ax.axhspan(-0.05,0.05, color=GREEN, alpha=.10, zorder=1)
ax.text(6.4,0.005,'F ≈ 0  (Hardy–Weinberg)', color=GREEN, fontsize=8.4, va='bottom')

ax.plot(MID, Fi, color=MUTED, lw=2.0, ls=(0,(5,3)), marker='o', ms=5,
        mfc=PANEL, zorder=4, label='ipyrad — binomial genotype calls')
ax.plot(MID, Fp, color=RED, lw=2.6, marker='o', ms=7, mfc=PANEL, mew=1.6,
        zorder=5, label='polyRAD — beta-binomial, d = 5, PopStructLD')
ax.fill_between(MID, Fp, Fi, color=RED, alpha=.08, zorder=2)

ax.annotate('polyRAD  +0.74 → −0.43\nipyrad   +0.87 → −0.34', xy=(150,-0.434), xytext=(30,-0.66),
            color=RED, fontsize=9, linespacing=1.35,
            arrowprops=dict(arrowstyle='->', color=RED, lw=1.0))
ax.set_ylabel('$F_{IS}$   (within stand)', color=INK, fontsize=10.5)
ax.set_ylim(-0.72,1.0)
ax.legend(frameon=False, fontsize=9, loc='upper right', labelcolor=MUTED)
ax.set_title('Modeling overdispersion does not remove the depth dependence',
             color=INK, fontsize=12.5, loc='left', pad=12)

ax2.plot(MID, Hep, color=MUTED, lw=1.9, label='H$_{exp}$')
ax2.plot(MID, Hop, color=RED,  lw=2.3, marker='o', ms=5, mfc=PANEL,
         label='H$_{obs}$ (polyRAD posterior calls)')
ax2.fill_between(MID, Hop, Hep, color=RED, alpha=.10)
ax2.set_ylabel('heterozygosity', color=INK, fontsize=10.5)
ax2.set_xlabel('read depth per genotype  (log scale)', color=INK, fontsize=10.5)
ax2.legend(frameon=False, fontsize=9, loc='center left', labelcolor=MUTED,
           bbox_to_anchor=(0.02,0.62))
ax2.set_ylim(0,0.48)
ticks=[6,10,20,30,50,100,150]
ax2.set_xticks(ticks); ax2.set_xticklabels([str(t) for t in ticks])
for x,nv in zip(MID,Np):
    ax2.annotate('%.0fk'%(nv/1000), xy=(x,0.445), ha='center', fontsize=7.2, color=MUTED)
ax2.text(6.2,0.462,'genotypes per bin', fontsize=7.4, color=MUTED)

fig.text(.012,.012,
 'polyRAD 2.0.1: overdispersion d = 5 fitted by TestOverdispersion (default is 9); '
 'IteratePopStructLD with LDdist = 150 bp (one RAD tag); contamination rate 0.001 (default —\n'
 'no blank barcodes exist in these libraries). 956,029 genotypes, 34,674 loci, 177 trees, '
 'within stands so Wahlund cannot contribute. Allele frequency per locus per stand across all\n'
 'depths, so H$_{exp}$ is a locus property and only H$_{obs}$ varies with depth. '
 'Spread across depth: binomial 1.206, polyRAD 1.175.',
 color=MUTED, fontsize=7.5, linespacing=1.5)
fig.tight_layout(rect=[0,.075,1,1])
fig.savefig('polyRAD_F_by_depth.png', dpi=200, facecolor=PAGE)
fig.savefig('polyRAD_F_by_depth.pdf', facecolor=PAGE)
print('binomial spread %.3f | polyRAD spread %.3f'%(max(Fi)-min(Fi), max(Fp)-min(Fp)))
print('saved polyRAD_F_by_depth.png / .pdf')
