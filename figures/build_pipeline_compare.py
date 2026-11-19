#!/usr/bin/env python3
"""Three independent pipelines, same F-vs-depth curve."""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'
NM=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
MID=[8,12.5,17.5,25,40,75,150]
fb={r['depth']:r for r in csv.DictReader(open('figdata/fb_F_by_depth.csv'))}
pr={r['depth']:r for r in csv.DictReader(open('figdata/polyrad_F_by_depth.csv'))}
S={
 'ipyrad de novo · binomial calls':      ([0.869,0.815,0.785,0.651,0.306,-0.196,-0.336], MUTED,(0,(5,3)),1.8),
 'ipyrad de novo · polyRAD (d=5)':       ([float(pr[n]['F_IS']) for n in NM], BLUE,(0,(2,2)),2.0),
 'bowtie2 + FreeBayes · polyRAD (d=7)':  ([float(fb[n]['F_IS']) for n in NM], RED,'-',2.8),
}
fig,(ax,ax2)=plt.subplots(2,1,figsize=(8.8,8.2),sharex=True,gridspec_kw={'height_ratios':[2.05,1]})
fig.patch.set_facecolor(PAGE)
for a in (ax,ax2):
    a.set_facecolor(PANEL); a.set_xscale('log')
    for s in ('top','right'): a.spines[s].set_visible(False)
    for s in ('left','bottom'): a.spines[s].set_color(RULE)
    a.tick_params(colors=MUTED,labelsize=9); a.grid(alpha=.25,color=RULE,lw=.7)
ax.axhline(0,color=INK,lw=1.1,zorder=3)
ax.axhspan(-0.05,0.05,color=GREEN,alpha=.10,zorder=1)
for lab,(v,c,ls,lw) in S.items():
    ax.plot(MID,v,color=c,lw=lw,ls=ls,marker='o',ms=5.5,mfc=PANEL,mew=1.3,zorder=4,label=lab)
ax.text(6.4,0.005,'F ≈ 0  (Hardy–Weinberg)',color=GREEN,fontsize=8.4,va='bottom')
ax.set_ylabel('$F_{IS}$   (within stand)',color=INK,fontsize=10.5)
ax.set_ylim(-0.62,1.0)
ax.legend(frameon=False,fontsize=9,loc='lower left',labelcolor=MUTED)
ax.set_title('The same depth artifact under three independent pipelines',
             color=INK,fontsize=12.5,loc='left',pad=12)
ax.text(46,0.72,'spread across depth\n  ipyrad binomial  1.206\n  polyRAD de novo  1.175\n  FreeBayes ref    1.231',
        color=MUTED,fontsize=8.8,family='monospace',linespacing=1.5,va='top')
Ho=[float(fb[n]['Ho']) for n in NM]; He=[float(fb[n]['He']) for n in NM]
N=[int(fb[n]['n_geno']) for n in NM]
ax2.plot(MID,He,color=MUTED,lw=1.9,label='H$_{exp}$')
ax2.plot(MID,Ho,color=RED,lw=2.3,marker='o',ms=5,mfc=PANEL,label='H$_{obs}$ (FreeBayes route)')
ax2.fill_between(MID,Ho,He,color=RED,alpha=.10)
ax2.set_ylabel('heterozygosity',color=INK,fontsize=10.5)
ax2.set_xlabel('read depth per genotype  (log scale)',color=INK,fontsize=10.5)
ax2.legend(frameon=False,fontsize=9,loc='center left',labelcolor=MUTED,
           bbox_to_anchor=(0.02,0.58))
ax2.set_ylim(0,0.50)
t=[6,10,20,30,50,100,150]; ax2.set_xticks(t); ax2.set_xticklabels([str(x) for x in t])
for x,n in zip(MID,N): ax2.annotate('%.0fk'%(n/1000),xy=(x,0.465),ha='center',fontsize=7.2,color=MUTED)
ax2.text(6.2,0.482,'genotypes per bin',fontsize=7.4,color=MUTED)
fig.text(.012,.048,
 'Reference route: 177 libraries mapped with bowtie2 (median 94.1% aligned), FreeBayes called jointly on all 177 (--pooled-continuous, --use-best-n-alleles 4),',
 color=MUTED,fontsize=7.4)
fig.text(.012,.029,
 '428,396 biallelic SNPs -> polyRAD: 122,051 loci, 2,121,746 genotypes. De novo route: ipyrad, 34,674 loci, 956,029 genotypes.',
 color=MUTED,fontsize=7.4)
fig.text(.012,.010,
 'Overdispersion was fitted independently in each route: d = 7 (reference) and d = 5 (de novo), both far below polyRAD\'s default of 9.',
 color=MUTED,fontsize=7.4)
fig.tight_layout(rect=[0,.075,1,1])
fig.savefig('pipeline_comparison.png',dpi=200,facecolor=PAGE)
fig.savefig('pipeline_comparison.pdf',facecolor=PAGE)
print('saved pipeline_comparison.png / .pdf')
