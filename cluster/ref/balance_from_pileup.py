#!/usr/bin/env python3
"""Depth-STRATIFIED allele balance from a samtools mpileup stream.

Why stratified: a fair binomial tightens as depth rises (SD = sqrt(p(1-p)/n),
so 0.091 at 30x but 0.050 at 100x). Real diploid heterozygotes must therefore
sharpen at 0.5 with depth. Mixed templates (collapsed paralogs, contamination,
index hopping) stay ~uniform (SD 0.202) at every depth. The TREND of the
overdispersion ratio across depth is the diagnostic; any single threshold is not.

Bins match the strata used in the original de novo characterisation so the two
are directly comparable, and start at 6x because that is ipyrad's mindepth --
the actual floor of the working dataset (cr30), not the 20x/30x diagnostic cuts.
"""
import sys, argparse, re, math
ap = argparse.ArgumentParser()
ap.add_argument('--label', default='')
ap.add_argument('--depths', help='write per-site depth+balance TSV here')
a = ap.parse_args()

BINS = [(6,10),(10,15),(15,20),(20,30),(30,50),(50,100),(100,10**9)]
NAMES = ['6-10','10-15','15-20','20-30','30-50','50-100','100+']
NB = 20
INDEL = re.compile(r'[+-](\d+)')

def parse(bases):
    i, n = 0, len(bases); ref = 0; alt = {}
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            m = INDEL.match(bases, i)
            if m: i = m.end() + int(m.group(1)); continue
            i += 1; continue
        if c in '.,': ref += 1; i += 1; continue
        u = c.upper()
        if u in 'ACGT': alt[u] = alt.get(u, 0) + 1; i += 1; continue
        i += 1
    return ref, alt

nb = len(BINS)
tot=[0]*nb; hist=[[0]*NB for _ in range(nb)]
hs=[0]*nb; hsum=[0.0]*nb; hsq=[0.0]*nb; dsum=[0.0]*nb
dh = open(a.depths,'w') if a.depths else None
if dh: dh.write('contig\tpos\tdepth\tbalance\n')

for line in sys.stdin:
    f = line.rstrip('\n').split('\t')
    if len(f) < 5: continue
    rc, ad = parse(f[4])
    if ad:
        ab, ac = max(ad.items(), key=lambda x: x[1])
        dp = rc + ac; bal = ac/dp if dp else 0.0
    else:
        dp = rc; bal = 0.0
    if dp < 6: continue
    k = next((j for j,(lo,hi) in enumerate(BINS) if lo <= dp < hi), None)
    if k is None: continue
    tot[k]+=1; hist[k][min(int(bal*NB),NB-1)]+=1
    if 0.15 < bal < 0.85:
        hs[k]+=1; hsum[k]+=bal; hsq[k]+=bal*bal; dsum[k]+=dp
    if dh: dh.write('%s\t%s\t%d\t%.4f\n'%(f[0],f[1],dp,bal))
if dh: dh.close()

print('label: %s'%a.label)
print('\n%-7s %10s %10s %7s %7s %7s %8s %8s %7s'%(
    'depth','sites','het band','het%','mean','SD','binom SD','overdisp','c0.5%'))
print('-'*82)
for k,nm in enumerate(NAMES):
    if not tot[k]: 
        print('%-7s %10d %10s'%(nm,0,'-')); continue
    row = '%-7s %10d %10d %6.2f%%'%(nm,tot[k],hs[k],100*hs[k]/tot[k])
    if hs[k] >= 50:
        m = hsum[k]/hs[k]
        sd = math.sqrt(max(hsq[k]/hs[k]-m*m,0))
        dbar = dsum[k]/hs[k]
        bsd = math.sqrt(0.25/dbar)          # fair binomial SD at this bin's mean depth
        c = 100*sum(hist[k][9:11])/hs[k]
        row += ' %7.4f %7.4f %8.4f %7.2fx %6.2f%%'%(m,sd,bsd,sd/bsd,c)
    print(row)
print('\nREAD: overdisp ~1x and FALLING toward 1 with depth => real heterozygotes.')
print('      overdisp RISING with depth (SD pinned near 0.202) => mixed templates.')
print('      uniform-distribution SD = 0.2021; c0.5%% under uniform ~14.3%%')

print('\n=== full histograms per depth bin (fraction ALT, width 0.05) ===')
for k,nm in enumerate(NAMES):
    if not tot[k]: continue
    print('\n-- depth %s  (n=%d) --'%(nm,tot[k]))
    mx = max(hist[k]) or 1
    for i,h in enumerate(hist[k]):
        print('  %.2f-%.2f %9d (%5.2f%%) %s'%(i/NB,(i+1)/NB,h,100*h/tot[k],'#'*int(50*h/mx)))
