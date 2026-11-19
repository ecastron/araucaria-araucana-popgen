#!/usr/bin/env python
"""1D stepping-stone geometry for the Araucaria ranges.

The 343-node 2D lattice produced NO isolation by distance (ibd_r=+0.005 even with the ancestral
collapse set ~3x beyond the lattice coalescence time) and cost 3260 s per simulation. A 1D chain
does produce IBD (verified: ibd_r=+0.95, invariant to mutation rate) and runs in seconds.

This is not a computational compromise: the Andean distribution is a narrow band along the
cordillera, and classical 1D MDS of the stand coordinates recovers 95.6% of the variance on a
single axis with r=+0.997 against true pairwise distance (Coastal: 98.7%, r=+1.000).
"""
import numpy as np, pickle
P = "/home/ecastron/araucaria"
SPACING = 4.0          # km between adjacent demes

def _km(a, b):
    R = 6371.0
    la1, lo1 = np.radians(a[1]), np.radians(a[0]); la2, lo2 = np.radians(b[1]), np.radians(b[0])
    h = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*R*np.arcsin(np.sqrt(h))

def _mds1(C, idx):
    n = len(idx)
    Dm = np.array([[_km(C[a], C[b]) for b in idx] for a in idx])
    J = np.eye(n) - np.ones((n, n))/n
    B = -0.5*J@(Dm**2)@J
    w, v = np.linalg.eigh(B); o = np.argsort(w)[::-1]
    x = v[:, o[0]]*np.sqrt(max(w[o[0]], 0))
    x = x - x.min()
    return x, Dm

def build_chains(spacing=SPACING):
    D = pickle.load(open(f"{P}/results/Q1/abc/observed.pkl", "rb"))
    C, CO, stands = D["coords"], D["coastal"], D["stands"]
    out = {}
    an_idx = [i for i in range(12) if not CO[i]]
    co_idx = [i for i in range(12) if CO[i]]
    info = {}
    for tag, idx in (("A", an_idx), ("C", co_idx)):
        x, Dm = _mds1(C, idx)
        ndeme = int(np.ceil(x.max()/spacing)) + 1
        # assign each stand to its nearest deme, shifting on collision so every stand is distinct
        order = np.argsort(x)
        taken = {}
        assign = np.zeros(len(idx), dtype=int)
        for k in order:
            d = int(round(x[k]/spacing))
            while d in taken: d += 1
            taken[d] = k; assign[k] = d
        ndeme = max(ndeme, max(taken)+1)
        info[tag] = dict(idx=idx, x=x, assign=assign, ndeme=ndeme, Dm=Dm)
    nA, nC = info["A"]["ndeme"], info["C"]["ndeme"]
    # deme ids: Andean 0..nA-1, Coastal nA..nA+nC-1
    K = nA + nC
    stand_deme = np.zeros(12, dtype=int)
    for k, j in enumerate(info["A"]["idx"]): stand_deme[j] = info["A"]["assign"][k]
    for k, j in enumerate(info["C"]["idx"]): stand_deme[j] = nA + info["C"]["assign"][k]
    edges = [(i, i+1) for i in range(nA-1)] + [(nA+i, nA+i+1) for i in range(nC-1)]
    node_coastal = np.array([False]*nA + [True]*nC)
    # the H1' colonisation source / H3' contact corridor: the closest true stand pair across the gap
    best = min(((_km(C[a], C[b]), a, b) for a in an_idx for b in co_idx), key=lambda t: t[0])
    out = dict(nA=nA, nC=nC, K=K, edges=np.array(edges), stand_deme=stand_deme,
               node_coastal=node_coastal, stands=stands, coastal=CO, coords=C,
               gap_km=best[0], gap_andean_deme=int(stand_deme[best[1]]),
               gap_coastal_deme=int(stand_deme[best[2]]), spacing=spacing, info=info)
    return out

if __name__ == "__main__":
    L = build_chains()
    st = L["stands"]
    print(f"Andean chain: {L['nA']} demes ({(L['nA']-1)*L['spacing']:.0f} km)   "
          f"Coastal chain: {L['nC']} demes ({(L['nC']-1)*L['spacing']:.0f} km)   total K={L['K']}")
    print(f"edges: {L['edges'].shape[0]}   spacing {L['spacing']} km")
    print(f"gap between ranges: {L['gap_km']:.1f} km  "
          f"(Andean deme {L['gap_andean_deme']} <-> Coastal deme {L['gap_coastal_deme']})")
    print("\nstand -> deme:")
    for j in range(12):
        print(f"  {st[j]:>3} {'Coastal' if L['coastal'][j] else 'Andean ':>7} -> deme {L['stand_deme'][j]:>3}")
    assert len(set(L["stand_deme"].tolist())) == 12, "stand->deme is not injective"
    print("\nall 12 stands on distinct demes: OK")
    # fidelity of chain distance vs true distance
    import itertools as it
    tru, chn = [], []
    for a, b in it.combinations(range(12), 2):
        if L["coastal"][a] != L["coastal"][b]: continue
        tru.append(_km(L["coords"][a], L["coords"][b]))
        chn.append(abs(L["stand_deme"][a]-L["stand_deme"][b])*L["spacing"])
    print(f"within-range chain distance vs true: r={np.corrcoef(tru,chn)[0,1]:+.4f} (n={len(tru)})")
