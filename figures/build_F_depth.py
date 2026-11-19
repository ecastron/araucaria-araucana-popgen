#!/usr/bin/env python3
"""F_IS as a function of read depth, with the HWE null."""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'

rows=list(csv.DictReader(open('figdata/F_vs_depth.tsv'), delimiter='\t'))
mid=[float(r['mid']) for r in rows]
Fo =[float(r['F_obs']) for r in rows]
Fn =[float(r['F_null']) for r in rows]
Ho =[float(r['Ho']) for r in rows]
He =[float(r['He']) for r in rows]
n  =[int(r['n_geno']) for r in rows]

# where the observed curve crosses zero (linear interpolation)
cross=None
for i in range(len(Fo)-1):
    if Fo[i]>0>=Fo[i+1]:
        t=Fo[i]/(Fo[i]-Fo[i+1]); cross=mid[i]+t*(mid[i+1]-mid[i]); break

fig,(ax,ax2)=plt.subplots(2,1,figsize=(8.4,8.2),sharex=True,
                          gridspec_kw={'height_ratios':[2.05,1]})
fig.patch.set_facecolor(PAGE)
for a in (ax,ax2):
    a.set_facecolor(PANEL)
    a.set_xscale('log')
    for s in ('top','right'): a.spines[s].set_visible(False)
    for s in ('left','bottom'): a.spines[s].set_color(RULE)
    a.tick_params(colors=MUTED, labelsize=9)
    a.grid(alpha=.25, color=RULE, lw=.7)

ax.axhline(0, color=INK, lw=1.1, zorder=2)
ax.fill_between(mid, Fo, 0, where=[f>0 for f in Fo], color=RED,   alpha=.13, zorder=1)
ax.fill_between(mid, Fo, 0, where=[f<0 for f in Fo], color=BLUE,  alpha=.13, zorder=1)
ax.plot(mid, Fn, color=MUTED, lw=1.6, ls=(0,(5,3)), zorder=3,
        label='HWE null (same depths, same caller)')
ax.plot(mid, Fo, color=RED, lw=2.4, zorder=4)
ax.scatter(mid, Fo, s=[max(14,min(95,x/900)) for x in n], color=RED,
           edgecolor=PANEL, lw=.9, zorder=5, label='observed $F_{IS}$ (within stand)')
if cross:
    ax.axvline(cross, color=GREEN, lw=1.2, ls=':', zorder=3)
    ax.annotate('crosses zero\nat ~%.0f×'%cross, xy=(cross,0), xytext=(cross*1.13,0.30),
                color=GREEN, fontsize=9,
                arrowprops=dict(arrowstyle='-', color=GREEN, lw=.9))
ax.axvline(32.6, color=MUTED, lw=1.0, ls='-.', alpha=.8, zorder=2)
ax.annotate('mean coverage 32.6×', xy=(32.6,-0.52), xytext=(11.5,-0.545),
            color=MUTED, fontsize=8.5)
ax.text(7.4, .70, 'heterozygotes MISSED\n(called homozygous)', color=RED,
        fontsize=9.5, va='top', linespacing=1.35)
ax.text(96, -.30, 'heterozygotes\nin EXCESS', color=BLUE, fontsize=9.5,
        va='top', ha='left', linespacing=1.35)
ax.set_ylabel('$F_{IS}$   (1 − H$_{obs}$/H$_{exp}$)', color=INK, fontsize=10.5)
ax.set_ylim(-0.72,1.0)
ax.legend(frameon=False, fontsize=9, loc='upper right', labelcolor=MUTED)
ax.set_title('Heterozygote calling in the ARAUCARIA RAD data depends on read depth',
             color=INK, fontsize=12.5, pad=12, loc='left')
ax.text(0,1.005,'', transform=ax.transAxes)

ax2.plot(mid, He, color=MUTED, lw=1.9, label='H$_{exp}$ (from read-based allele freq.)')
ax2.plot(mid, Ho, color=RED,  lw=2.2, label='H$_{obs}$ (called heterozygotes)')
ax2.fill_between(mid, Ho, He, color=RED, alpha=.10)
ax2.set_ylabel('heterozygosity', color=INK, fontsize=10.5)
ax2.set_xlabel('read depth per genotype  (log scale)', color=INK, fontsize=10.5)
ax2.legend(frameon=False, fontsize=9, loc='lower right', labelcolor=MUTED,
           bbox_to_anchor=(1.0,0.02))
ticks=[6,10,20,30,50,75,100,200]
ax2.xaxis.set_major_locator(FixedLocator(ticks))
ax2.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
ax2.set_ylim(0,0.28)

fig.text(.012,.012,
 'Within-stand $F_{IS}$, %s genotypes, ct091 cr30. Allele frequency estimated per locus per stand from READ counts\n'
 'so it is independent of genotype calling. Null: simulated Hardy–Weinberg population, identical depth distribution,\n'
 'identical calling rule — it stays flat at ~0, so the observed gradient is not a mechanical artifact of binning.'
 % format(sum(n),','), color=MUTED, fontsize=7.6, linespacing=1.5)
fig.tight_layout(rect=[0,.058,1,1])
fig.savefig('F_vs_depth.png', dpi=200, facecolor=PAGE)
fig.savefig('F_vs_depth.pdf', facecolor=PAGE)
print('crossing point: %.1fx'%(cross or float('nan')))
print('saved F_vs_depth.png / .pdf')
