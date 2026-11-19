#!/usr/bin/env python3
"""POSITIVE CONTROL for the overdispersion statistic.

Simulates a diploid with KNOWN heterozygous sites, RAD-like read structure,
and our read lengths/error rate. If balance_from_pileup.py does not return
~1.0x overdispersion on this, the statistic is broken and every conclusion
drawn from it today is an artifact of the instrument rather than the data.

Also emits a NEGATIVE control arm (3 haplotypes mixed at random ratios per
locus) which MUST return high overdispersion -- otherwise the statistic has
no power and 1.0x would be meaningless too.
"""
import random, sys, argparse, gzip
ap = argparse.ArgumentParser()
ap.add_argument('--ref', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--mode', choices=['diploid','mixture'], default='diploid')
ap.add_argument('--het', type=float, default=0.005)   # per-bp heterozygosity
ap.add_argument('--loci', type=int, default=20000)
ap.add_argument('--err', type=float, default=0.00016) # Q38
ap.add_argument('--seed', type=int, default=42)
a = ap.parse_args()
rng = random.Random(a.seed)

# load reference contigs
seqs={}; name=None; buf=[]
for ln in open(a.ref):
    if ln.startswith('>'):
        if name: seqs[name]=''.join(buf)
        name=ln[1:].split()[0]; buf=[]
    else: buf.append(ln.strip())
if name: seqs[name]=''.join(buf)
names=[n for n in seqs if len(seqs[n])>50000]
sys.stderr.write('contigs usable: %d\n'%len(names))

COMP=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(COMP)[::-1]
ALT={'A':'CGT','C':'AGT','G':'ACT','T':'ACG'}
def mutate(s,rate):
    l=list(s)
    for i,c in enumerate(l):
        if c in ALT and rng.random()<rate: l[i]=rng.choice(ALT[c])
    return ''.join(l)
def adderr(s,rate):
    l=list(s)
    for i,c in enumerate(l):
        if c in ALT and rng.random()<rate: l[i]=rng.choice(ALT[c])
    return ''.join(l)

R1L,R2L = 101,91
nhap = 2 if a.mode=='diploid' else 3
o1=gzip.open(a.out+'_1.fq.gz','wt'); o2=gzip.open(a.out+'_2.fq.gz','wt')
q1='I'*R1L; q2='I'*R2L
nread=0
for li in range(a.loci):
    cn=rng.choice(names); L=len(seqs[cn])
    pos=rng.randrange(1000, L-1200)
    region=seqs[cn][pos:pos+900]
    if 'N' in region[:400]: continue
    haps=[region]+[mutate(region,a.het) for _ in range(nhap-1)]
    depth=rng.choice([7,9,12,17,25,35,45,70,120,200])
    if a.mode=='diploid':
        w=[0.5,0.5]
    else:
        # mixture arm: arbitrary per-locus ratios -> the artifact we think we see
        r=[rng.random() for _ in range(nhap)]; s=sum(r); w=[x/s for x in r]
    for d in range(depth):
        h=rng.choices(range(nhap),weights=w)[0]
        ins=rng.randrange(180,420)
        frag=haps[h][:ins]
        if len(frag)<max(R1L,R2L): continue
        r1=adderr(frag[:R1L],a.err); r2=adderr(rc(frag)[:R2L],a.err)
        nread+=1
        o1.write('@sim%d_%d/1\n%s\n+\n%s\n'%(li,d,r1,q1))
        o2.write('@sim%d_%d/2\n%s\n+\n%s\n'%(li,d,r2,q2))
o1.close(); o2.close()
sys.stderr.write('mode=%s loci=%d pairs=%d\n'%(a.mode,a.loci,nread))
