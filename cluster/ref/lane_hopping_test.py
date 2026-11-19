#!/usr/bin/env python3
"""Does index hopping / demultiplexing explain the flat allele balance?

Index hopping happens WITHIN a lane, never across. So if a library's spurious
minor alleles come from lane-mates, they should be carried by same-lane
libraries more often than by different-lane libraries. Cross-lane sharing is
the built-in null.

CONFOUND CHECK FIRST: if libraries were loaded onto lanes by geographic batch,
lane-mates are also population-mates and elevated sharing would be BIOLOGY, not
hopping. The script reports the lane x stand table and, if they are associated,
restricts the comparison to WITHIN-STAND pairs so population is held constant.
"""
import sys, argparse, collections
ap=argparse.ArgumentParser()
ap.add_argument('--vcf',required=True); ap.add_argument('--lanes',required=True)
ap.add_argument('--popmap',required=True); ap.add_argument('--mindp',type=int,default=20)
ap.add_argument('--maxsites',type=int,default=40000)
a=ap.parse_args()

lane={}
for ln in open(a.lanes):
    p=ln.rstrip('\n').split('\t')
    if len(p)>=2: lane[p[0]]=p[1]
stand={}
for ln in open(a.popmap):
    p=ln.split()
    if len(p)>=2: stand[p[0]]=p[1]

print('=== lane x stand table (confound check) ===')
tab=collections.Counter((lane.get(s,'?'),stand[s]) for s in stand if s in lane)
lanes=sorted({k[0] for k in tab}); stands=sorted({k[1] for k in tab})
print('%-8s'%'lane', ' '.join('%5s'%x for x in stands))
for L in lanes:
    print('%-8s'%L, ' '.join('%5d'%tab[(L,s)] for s in stands))
conf = any(all(tab[(L,s)]==0 for L in lanes if L!=L0) for s in stands for L0 in lanes if tab[(L0,s)]>0)
print('\nlane and stand CONFOUNDED (some stand entirely on one lane): %s'%conf)
print('-> restricting to WITHIN-STAND pairs regardless, so population is held constant.\n')

IDX={'C':0,'A':1,'T':2,'G':3}
same_num=same_den=diff_num=diff_den=0
nsites=0
with open(a.vcf) as fh:
    for line in fh:
        if line.startswith('##'): continue
        f=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            cols=[(i,s.split('.')[0]) for i,s in enumerate(f[9:],9)]
            cols=[(i,s) for i,s in cols if s in lane and s in stand]
            print('libraries usable (lane + stand known): %d'%len(cols))
            continue
        ref,alt=f[3],f[4]
        if len(alt)!=1 or ref not in IDX or alt not in IDX: continue
        ri,ai=IDX[ref],IDX[alt]
        fmt=f[8].split(':')
        try: ci=fmt.index('CATG')
        except ValueError: continue
        # per-library allele counts at this site
        cnt={}
        for i,s in cols:
            p=f[i].split(':')
            if len(p)<=ci: continue
            try: catg=[int(x) for x in p[ci].split(',')]
            except ValueError: continue
            rc,ac=catg[ri],catg[ai]
            if rc+ac>0: cnt[s]=(rc,ac)
        if len(cnt)<20: continue
        nsites+=1
        if nsites>a.maxsites: break
        for s,(rc,ac) in cnt.items():
            dp=rc+ac
            if dp<a.mindp: continue
            bal=ac/dp
            if not (0.15<bal<0.85): continue          # a flat-balance genotype
            minor = 'alt' if ac<rc else 'ref'         # the suspicious allele
            for s2,(rc2,ac2) in cnt.items():
                if s2==s: continue
                if stand[s2]!=stand[s]: continue      # hold population constant
                carries = (ac2>0) if minor=='alt' else (rc2>0)
                if lane[s2]==lane[s]: same_den+=1; same_num+=carries
                else:                 diff_den+=1; diff_num+=carries

print('\n=== does the minor allele come from lane-mates? (within-stand pairs only) ===')
print('sites scanned: %d'%nsites)
if same_den and diff_den:
    ps=same_num/same_den; pd=diff_num/diff_den
    print('  same-lane  : %d/%d carry it = %.4f'%(same_num,same_den,ps))
    print('  diff-lane  : %d/%d carry it = %.4f'%(diff_num,diff_den,pd))
    print('  ratio same/diff = %.3f'%(ps/pd if pd else float('nan')))
    print('\n  >1.05 and significant -> consistent with index hopping')
    print('  ~1.0                  -> hopping does NOT explain these alleles')
else:
    print('  insufficient pairs (same_den=%d diff_den=%d)'%(same_den,diff_den))
