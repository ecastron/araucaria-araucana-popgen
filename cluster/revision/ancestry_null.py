#!/usr/bin/env python
"""Revision (26 Nov 2026): null distribution of sNMF K = 2 ancestry coefficients
under the two-refugia history with zero gene flow (H2'), using the Figure 2E method.

Replaces the 26 Aug check (plain NMF, mean-imputed, 3 replicates; ANALYSIS_SPEC.md).
Parameter sets: the 50 H2 rows of reftable_sim9.csv (the published ABC reference table) nearest
the observed summary statistics (Euclidean distance on the 15 statistics scaled by their MAD
across the whole table; rejection ABC, Beaumont et al. 2002, Genetics 162:2025-2035), so the
null is drawn from histories that reproduce the observed data. As a contrast, the 50 nearest H3
rows (two refugia with secondary contact) are simulated the same way.
Genotypes: sim9.one() pipeline exactly (163 diploids in the observed stand sizes, 1 Mb,
r = 1e-7, mutation rate tuned to 20,000 SNPs, evenly spaced subsample, missingness copied
column-wise from the observed 163-tree mask). New seeds: 31,000,000 + 1,000 * model + rank.
Each replicate is written as a LEA .geno file (9 = missing) for ancestry_null.R.
Usage: ancestry_null.py TASK NTASKS   (replicates rank % NTASKS == TASK - 1)
Output: results/revision/ancnull/<model>_<rank>.geno, results/revision/ancnull/selected.tsv
"""
import sys, os, numpy as np, pandas as pd, msprime
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim5 as B

P = "/home/ecastron/araucaria"; OUT = f"{P}/results/revision/ancnull"; os.makedirs(OUT, exist_ok=True)
REF = f"{P}/results/Q1/abc_clock_lit/reftable_sim9.csv"
OBS = f"{P}/results/Q1/abc_clock_lit/observed_sim9.csv"
NSEL = 50
task, ntask = int(sys.argv[1]), int(sys.argv[2])

ref = pd.read_csv(REF); obs = pd.read_csv(OBS, index_col=0).loc["primary"]
stats = list(obs.index)
X = ref[stats].to_numpy(float)
mad = np.median(np.abs(X - np.median(X, 0)), 0)
dist = np.sqrt((((X - obs.to_numpy(float)) / mad) ** 2).sum(1))
ref["dist"] = dist
sel = []
for mi, mod in enumerate(("H2", "H3")):
    s = ref[ref.model == mod].nsmallest(NSEL, "dist").reset_index(drop=True)
    s["rank"] = np.arange(NSEL); s["seed"] = 31_000_000 + 1_000 * mi + s["rank"]
    sel.append(s)
sel = pd.concat(sel, ignore_index=True)
if task == 1:
    sel[["model", "rank", "seed", "dist", "p_Ne", "p_m", "p_T1", "p_Nanc", "p_Tsep", "p_mcross",
         "p_Tcontact", "p_T_years"]].to_csv(f"{OUT}/selected.tsv", sep="\t", index=False)

for _, row in sel.iterrows():
    if row["rank"] % ntask != task - 1: continue
    mod, seed = row["model"], int(row["seed"])
    fn = f"{OUT}/{mod}_{int(row['rank']):02d}.geno"
    if os.path.exists(fn): continue
    p = dict(Ne=row.p_Ne, m=row.p_m, T1=row.p_T1, Nanc=row.p_Nanc, Tsep=row.p_Tsep)
    if mod == "H3": p.update(mcross=row.p_mcross, Tcontact=row.p_Tcontact)
    dem = B.build(mod, p)
    samples = {f"d{B.STAND_DEME[j]}": int(B.NPER[j]) for j in range(B.ND)}
    anc = msprime.sim_ancestry(samples=samples, demography=dem, sequence_length=1_000_000,
                               recombination_rate=1e-7, random_seed=seed, ploidy=2)
    ts = None; mu = 1e-8
    for k in range(8):
        ts = msprime.sim_mutations(anc, rate=mu, random_seed=seed + 1 + k,
                                   model=msprime.BinaryMutationModel())
        if ts.num_sites >= B.NSNP or mu >= B.MU_CAP: break
        mu = min(mu * max(2.0, 1.3 * B.NSNP / max(ts.num_sites, 1)), B.MU_CAP)
    Gm = ts.genotype_matrix()
    Gd = (Gm[:, 0::2] + Gm[:, 1::2]).T.astype(float)
    if Gd.shape[1] > B.NSNP:
        Gd = Gd[:, np.linspace(0, Gd.shape[1] - 1, B.NSNP).astype(int)]
    cols = np.random.default_rng(seed + 2).choice(B.MASK.shape[1], size=Gd.shape[1], replace=True)
    Gd[B.MASK[:, cols]] = np.nan
    G = np.where(np.isnan(Gd), 9, Gd).astype(int).T          # SNPs x individuals
    keep = [(r[r != 9].size > 0) and (np.unique(r[r != 9]).size > 1) for r in G]
    G = G[np.array(keep)]
    with open(fn, "w") as fh:
        for r in G: fh.write("".join(map(str, r)) + "\n")
    print(f"{mod} rank {int(row['rank'])}: {G.shape[0]} SNPs, {G.shape[1]} individuals", flush=True)

# observed 163-tree matrix subsampled to the same SNP count, evenly spaced (task 1 only)
if task == 1 and not os.path.exists(f"{OUT}/obs20k.geno"):
    L = open(f"{P}/results/revision/all163.geno").read().split()
    idx = np.linspace(0, len(L) - 1, B.NSNP).astype(int)
    open(f"{OUT}/obs20k.geno", "w").write("\n".join(L[i] for i in idx) + "\n")
