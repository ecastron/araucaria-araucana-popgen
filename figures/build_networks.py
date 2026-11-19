#!/usr/bin/env python3
"""Two network panels: NeighborNet of stands (FST) and MSN of individuals (read dosage)."""
import csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.sparse.csgraph import minimum_spanning_tree

PAGE='#f6f4ec'; PANEL='#fffdf7'; RULE='#dcd8ca'
INK='#14180f'; MUTED='#6c6b5e'; RED='#a8543a'; BLUE='#35708c'; GREEN='#40663f'
COASTAL={'PN','TR','VA'}

# ---------- Panel A : NeighborNet of the 12 stands ----------
V={int(r['v']):(float(r['x']),float(r['y'])) for r in csv.DictReader(open('figdata/nnet_vertices.tsv'),delimiter='\t')}
E=[(int(r['v1']),int(r['v2'])) for r in csv.DictReader(open('figdata/nnet_edges.tsv'),delimiter='\t')]
T={r['label']:int(r['v']) for r in csv.DictReader(open('figdata/nnet_tips.tsv'),delimiter='\t')}

# ---------- Panel B : individuals ----------
L=open('figdata/indiv_dist.tsv').read().rstrip('\n').split('\n')
names=L[0].split('\t')[1:]
D=np.array([[float(x) if x!='NA' else np.nan for x in l.split('\t')[1:]] for l in L[1:]])
D=(D+D.T)/2
D[np.isnan(D)]=np.nanmax(D)                     # 1.9% of pairs lack shared loci
np.fill_diagonal(D,0)
pop={}
for ln in open('figdata/popmap_177.tsv'):
    p=ln.split()
    if len(p)>=2: pop[p[0]]=p[1]
stands=[pop.get(n,'?') for n in names]
# classical MDS for layout
n=len(D); J=np.eye(n)-np.ones((n,n))/n
B=-0.5*J@(D**2)@J
w,v=np.linalg.eigh(B); order=np.argsort(w)[::-1]
XY=v[:,order[:2]]*np.sqrt(np.maximum(w[order[:2]],0))
mst=minimum_spanning_tree(D).tocoo()

ORDER=['VA','TR','PN','BP','HL','LM','LR','MC','MM','PC','RA','TH']
COL={s:(BLUE if s in COASTAL else RED) for s in ORDER}
# distinct fills so 12 stands are separable; ring color still encodes cordillera
SHADE={'VA':'#9ecfe0','TR':'#4f9bb8','PN':'#1f5f7a',
       'BP':'#e0a48f','HL':'#d18a6e','LM':'#c47150','LR':'#a8543a',
       'MC':'#8d4430','MM':'#f0c9a8','PC':'#c98b5e','RA':'#7a3a28','TH':'#e8b48c'}

fig,(axA,axB)=plt.subplots(1,2,figsize=(13.4,6.4))
fig.patch.set_facecolor(PAGE)
for a in (axA,axB):
    a.set_facecolor(PANEL); a.set_aspect('equal'); a.axis('off')

for v1,v2 in E:
    if v1 in V and v2 in V:
        (x1,y1),(x2,y2)=V[v1],V[v2]
        axA.plot([x1,x2],[y1,y2],color=MUTED,lw=1.0,alpha=.55,zorder=1)
for lab,vi in T.items():
    if vi not in V: continue
    x,y=V[vi]; c=BLUE if lab in COASTAL else RED
    axA.scatter([x],[y],s=180,color=c,edgecolor=PANEL,lw=1.6,zorder=3)
    axA.annotate(lab,(x,y),xytext=(0,-16),textcoords='offset points',
                 ha='center',fontsize=9.5,color=INK,weight='bold',zorder=4)
axA.set_title('A   NeighborNet of stands  (pairwise $F_{ST}$)',
              color=INK,fontsize=11.5,loc='left',pad=6)
axA.legend(handles=[Line2D([],[],marker='o',ls='',ms=9,mfc=BLUE,mec=PANEL,label='Coastal'),
                    Line2D([],[],marker='o',ls='',ms=9,mfc=RED,mec=PANEL,label='Andean')],
           frameon=False,fontsize=9,loc='lower left',labelcolor=MUTED)

for i,j,_ in zip(mst.row,mst.col,mst.data):
    axB.plot([XY[i,0],XY[j,0]],[XY[i,1],XY[j,1]],color=RULE,lw=.55,alpha=.45,zorder=1)
for s in ORDER:
    k=[i for i,x in enumerate(stands) if x==s]
    if not k: continue
    axB.scatter(XY[k,0],XY[k,1],s=46,color=SHADE[s],edgecolor=COL[s],lw=1.2,
                zorder=3,label=s)
axB.set_title('B   Minimum spanning network of 177 trees  (read-dosage distance)',
              color=INK,fontsize=11.5,loc='left',pad=6)
axB.legend(frameon=False,fontsize=7.6,loc='center left',bbox_to_anchor=(1.0,0.5),
           labelcolor=MUTED,ncol=1,handletextpad=.3)

fig.text(.012,.052,
 'A  NeighborNet (phangorn) on the 66 pairwise $F_{ST}$ values, 12 stands. Box-like regions indicate conflicting signal — relationships that are not strictly tree-like.',
 color=MUTED,fontsize=7.6)
fig.text(.012,.030,
 'B  Distance from per-locus ALT read fraction (dosage), never passing through a genotype call, so it is unaffected by the heterozygote-calling artifact.',
 color=MUTED,fontsize=7.6)
fig.text(.012,.009,
 '    Layout by classical MDS; edges are the minimum spanning tree — long crossing edges are a projection artifact of the 2-D layout, not real adjacency. Ring color = cordillera, fill = stand.',
 color=MUTED,fontsize=7.6)
fig.tight_layout(rect=[0,.085,1,1])
fig.savefig('networks.png',dpi=200,facecolor=PAGE)
fig.savefig('networks.pdf',facecolor=PAGE)
print('stands %d | individuals %d | MST edges %d'%(len(T),n,mst.nnz))
print('saved networks.png / .pdf')
