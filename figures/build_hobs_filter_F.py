#!/usr/bin/env python3
"""F_IS by depth under a species-wide Hobs ceiling (fixed-paralog filter)."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'
BINS=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
MID =[8,12.5,17.5,25,40,75,150]
COL =[INK,'#8a8778',BLUE,GREEN,RED]
LW  =[2.8,1.7,1.8,1.9,2.2]
LS  =['-', (0,(6,3)), (0,(3,2)), (0,(1,1.6)), (0,(5,2,1,2))]

rows=list(csv.DictReader(open('figdata/hobs_filter_F.tsv'), delimiter='\t'))
fig,ax=plt.subplots(figsize=(9.0,6.0))
fig.patch.set_facecolor(PAGE); ax.set_facecolor(PANEL)
ax.set_xscale('log')
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(RULE)
ax.tick_params(colors=MUTED, labelsize=9)
ax.grid(alpha=.25, color=RULE, lw=.7)
ax.axhline(0, color=INK, lw=1.1, zorder=3)
ax.axhspan(-0.05,0.05, color=GREEN, alpha=.10, zorder=1)

for k,r in enumerate(rows):
    ys=[float(r[b]) for b in BINS]
    ax.plot(MID,ys,color=COL[k],lw=LW[k],ls=LS[k],marker='o',ms=5,
            mfc=PANEL,mew=1.3,zorder=4+k,
            label='%s   %s loci (%s%%)   spread %s'%(
                r['filter'].ljust(15), format(int(r['loci']),',').rjust(6), r['pct'].rjust(3), r['spread']))

# annotations placed in the empty upper-right region so they miss the curves
ax.annotate('', xy=(8,0.656), xytext=(8,0.762),
            arrowprops=dict(arrowstyle='<->', color=RED, lw=1.2))
ax.annotate('', xy=(150,-0.256), xytext=(150,-0.445),
            arrowprops=dict(arrowstyle='<->', color=RED, lw=1.2))
ax.text(42, 0.80,
        'the ceiling compresses BOTH ends\n'
        '   6-10x   0.762 → 0.656\n'
        '   100+    −0.445 → −0.256\n'
        'but the slope between them is unchanged:\n'
        'spread 1.207 → 1.146, removing 41% of loci',
        color=MUTED, fontsize=8.6, va='top', linespacing=1.5)
ax.text(6.4,0.005,'F ≈ 0  (Hardy–Weinberg)', color=GREEN, fontsize=8.4, va='bottom')

ax.set_xlabel('read depth per genotype  (log scale)', color=INK, fontsize=10.5)
ax.set_ylabel('$F_{IS}$   (within stand)', color=INK, fontsize=10.5)
ax.set_ylim(-0.62,0.92)
ticks=[6,10,20,30,50,100,150]
ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
ax.legend(frameon=False, fontsize=8.0, loc='lower left', labelcolor=MUTED,
          handlelength=2.8, labelspacing=0.55, prop={'family':'monospace','size':7.6})
ax.set_title('A species-wide heterozygosity ceiling does not remove the depth dependence',
             color=INK, fontsize=12, loc='left', pad=12)
fig.text(.012,.014,
 'Fixed ancestral duplications should read as heterozygous in every tree, so a ceiling on species-wide H$_{obs}$ should remove them.\n'
 'Per-locus H$_{obs}$ from polyRAD weighted genotypes (het window 0.2-0.8): median 0.147, p95 0.576, p99 0.672 — so a 0.35 ceiling\n'
 'reaches well below HDplot\'s 0.55 cut-off and removes 25% of loci that HDplot never flags. polyRAD d = 5, IteratePopStructLD, within stands.',
 color=MUTED, fontsize=7.5, linespacing=1.5)
fig.tight_layout(rect=[0,.085,1,1])
fig.savefig('hobs_filter_F_by_depth.png', dpi=200, facecolor=PAGE)
fig.savefig('hobs_filter_F_by_depth.pdf', facecolor=PAGE)
print('saved hobs_filter_F_by_depth.png / .pdf')
