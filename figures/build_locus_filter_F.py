#!/usr/bin/env python3
"""F_IS by depth under structural locus filters."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'
BINS=['6-10','10-15','15-20','20-30','30-50','50-100','100+']
MID =[8,12.5,17.5,25,40,75,150]
COL =[INK, MUTED, BLUE, GREEN, RED]
LW  =[2.8,1.8,1.8,2.0,2.0]
LS  =['-', (0,(6,3)), (0,(3,2)), (0,(1,1.6)), (0,(5,2,1,2))]

rows=list(csv.DictReader(open('figdata/locus_filter_F.tsv'), delimiter='\t'))
fig,ax=plt.subplots(figsize=(9.0,6.2))
fig.patch.set_facecolor(PAGE); ax.set_facecolor(PANEL)
ax.set_xscale('log')
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(RULE)
ax.tick_params(colors=MUTED, labelsize=9)
ax.grid(alpha=.25, color=RULE, lw=.7)
ax.axhline(0, color=INK, lw=1.1, zorder=3)
ax.axhspan(-0.05,0.05, color=GREEN, alpha=.10, zorder=1)

for k,r in enumerate(rows):
    xs=[]; ys=[]
    for b,m in zip(BINS,MID):
        v=r[b]
        if v=='NA': continue
        xs.append(m); ys.append(float(v))
    lab='%s\n   %s loci (%s%%)  spread %s'%(r['filter'],format(int(r['loci']),','),r['pct'],r['spread'])
    ax.plot(xs,ys,color=COL[k],lw=LW[k],ls=LS[k],marker='o',ms=5,
            mfc=PANEL,mew=1.3,zorder=4+k,label=lab)

# the strict filter's 100+ bin was degenerate (F = 1.000 on a near-empty bin)
ax.plot([150],[1.0],marker='x',ms=9,color=RED,mew=2,zorder=6)
ax.annotate('degenerate bin\n(619 loci, F = 1.000)', xy=(150,1.0), xytext=(58,0.93),
            color=RED, fontsize=8.2, linespacing=1.3,
            arrowprops=dict(arrowstyle='-', color=RED, lw=.8))

ax.text(6.4,0.83,'every filter overlaps here', color=MUTED, fontsize=9)
ax.annotate('', xy=(6.6,0.80), xytext=(20,0.80),
            arrowprops=dict(arrowstyle='<->', color=MUTED, lw=.9))
ax.text(6.4,0.005,'F ≈ 0  (Hardy–Weinberg)', color=GREEN, fontsize=8.4, va='bottom')

ax.set_xlabel('read depth per genotype  (log scale)', color=INK, fontsize=10.5)
ax.set_ylabel('$F_{IS}$   (within stand)', color=INK, fontsize=10.5)
ax.set_ylim(-0.95,1.05)
ticks=[6,10,20,30,50,100,150]
ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
ax.legend(frameon=False, fontsize=7.6, loc='lower left', labelcolor=MUTED,
          handlelength=2.6, labelspacing=0.9)
ax.set_title('Structural locus filters do not remove the depth dependence',
             color=INK, fontsize=12.5, loc='left', pad=12)
fig.text(.012,.014,
 'polyRAD posterior genotypes (d = 5, IteratePopStructLD), within stands. Locus metrics computed over PRESENT genotypes only\n'
 '(zeros are missing data, not depth 0). Filters target segregating duplicates (depth CV), the multi-copy tail (mean depth)\n'
 'and null alleles (missingness). Two moderate filters make the spread WORSE; the strict one reduces it only by discarding 98% of loci.',
 color=MUTED, fontsize=7.5, linespacing=1.5)
fig.tight_layout(rect=[0,.085,1,1])
fig.savefig('locus_filter_F_by_depth.png', dpi=200, facecolor=PAGE)
fig.savefig('locus_filter_F_by_depth.pdf', facecolor=PAGE)
print('saved locus_filter_F_by_depth.png / .pdf')
