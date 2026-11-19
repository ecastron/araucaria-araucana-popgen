#!/usr/bin/env python
"""Stepping-stone lattice over the Chilean Araucaria range, used by all three ABC models.

The 12 stands are 12 points sampled from a CONTINUOUS population; a model with only 12 demes
exchanging migrants directly cannot generate isolation by distance, because the intervening
habitat that carries the stepping-stone chain is absent (verified: ibd_r = +0.06 vs +0.65 observed
even with no missing data). The lattice restores it.

Nodes come from the triangular grid already built for feems, masked to within MASK_KM of a stand.
That mask reproduces the real disjunct distribution: the Coastal-Andean gap is 117 km, wider than
2 x MASK_KM, so the two ranges form separate components, while VA (78 km from PN) stays attached
to the Coastal one. Araucaria occupies only the places we sampled -- the Central Valley holds no
trees -- so a bounding-box grid would let gene flow through habitat that does not exist.
"""
import numpy as np, pickle

P = "/home/ecastron/araucaria"
MASK_KM = 50.0

def _km(a, b):
    R = 6371.0
    la1, lo1 = np.radians(a[1]), np.radians(a[0]); la2, lo2 = np.radians(b[1]), np.radians(b[0])
    h = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*R*np.arcsin(np.sqrt(h))

def build_lattice(tag="tri12", mask_km=MASK_KM):
    D = pickle.load(open(f"{P}/results/Q1/abc/observed.pkl", "rb"))
    S = D["coords"]; coastal_stand = D["coastal"]; stands = D["stands"]
    g = np.load(f"{P}/results/Q1/feems/grid_{tag}.npy")
    e = np.load(f"{P}/results/Q1/feems/edges_{tag}.npy")
    off = 1 if e.min() == 1 else 0
    e = e - off
    # distance from every grid node to every stand
    D2 = np.array([[_km(g[i], S[j]) for j in range(len(S))] for i in range(g.shape[0])])
    keep = np.where(D2.min(axis=1) <= mask_km)[0]
    remap = -np.ones(g.shape[0], dtype=int); remap[keep] = np.arange(len(keep))
    edges = np.array([[remap[a], remap[b]] for a, b in e if remap[a] >= 0 and remap[b] >= 0])
    gk = g[keep]; D2k = D2[keep]
    # each node inherits the range of its nearest stand
    nearest = D2k.argmin(axis=1)
    node_coastal = coastal_stand[nearest]
    # each stand sits on its nearest kept node (bijection: 12 stands -> 12 distinct nodes)
    stand_node = np.array([np.argmin([_km(gk[i], S[j]) for i in range(len(gk))]) for j in range(len(S))])
    # connected components of the masked graph
    n = len(keep); adj = [[] for _ in range(n)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    comp = -np.ones(n, dtype=int); nc = 0
    for s in range(n):
        if comp[s] >= 0: continue
        st = [s]; comp[s] = nc
        while st:
            u = st.pop()
            for v in adj[u]:
                if comp[v] < 0: comp[v] = nc; st.append(v)
        nc += 1
    return dict(coords=gk, edges=edges, node_coastal=node_coastal, stand_node=stand_node,
                comp=comp, ncomp=nc, stands=stands, stand_coords=S, coastal_stand=coastal_stand)

if __name__ == "__main__":
    for tag in ("tri20", "tri12", "tri8"):
        L = build_lattice(tag)
        gk, comp = L["coords"], L["comp"]
        nc_ = L["node_coastal"]
        print(f"\n=== {tag} masked to {MASK_KM:.0f} km ===")
        print(f"  nodes kept : {gk.shape[0]}   edges kept : {L['edges'].shape[0]}")
        print(f"  components : {L['ncomp']}  sizes={[int((comp==k).sum()) for k in range(L['ncomp'])]}")
        print(f"  Coastal nodes {int(nc_.sum())} | Andean nodes {int((~nc_).sum())}")
        sn = L["stand_node"]
        print(f"  stands map to {len(set(sn.tolist()))} distinct nodes (need 12)")
        for k in range(L["ncomp"]):
            mem = np.where(comp == k)[0]
            which = sorted({L['stands'][j] for j in range(12) if comp[sn[j]] == k})
            print(f"    comp{k}: n={len(mem):3d} coastal={int(nc_[mem].sum()):3d} stands={which}")
