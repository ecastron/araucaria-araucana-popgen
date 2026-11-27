# Population history and adaptive divergence in *Araucaria araucana*

Code and derived data for the article "Population history and adaptive divergence in *Araucaria
araucana*" (Castro-Nallar, in preparation for *Molecular Ecology*).

The study combines RAD-seq genotypes of 163 trees from 12 stands spanning the Chilean range
(78,026 SNPs), published common-garden traits and ITS1 sequencing of foliar endophytes to ask
how the Coastal and Andean ranges are structured, how and when they separated, whether traits
diverged more than neutral markers, and whether the endophyte community tracks host lineage.

## Data

Raw reads are in the NCBI Sequence Read Archive:

| Data | BioProject | Runs |
|---|---|---|
| RAD-seq (host genotypes) | PRJNA634877 | 188 runs, 173 trees |
| ITS1 and 16S amplicons (endophytes) | PRJNA517193 | January 2018 campaign |

Run-level metadata, including stand coordinates and collection dates, are in `metadata/sra_runs_radseq.csv` and
`metadata/sra_runs_amplicon.csv`. The reference genome is AraucariaDB v1. Climate is CHELSA v2.1,
soil is SoilGrids 2.0, and the forest distribution is the CONAF forest cadastre (Catastro);
`environment/` holds the extraction scripts and the extracted values.

## Repository layout

| Folder | Contents |
|---|---|
| `cluster/` | The scripts on the path from raw reads to the data behind each figure, as run on the McCartney HPC cluster (Slurm): assembly, SNP calling and filtering, structure, F_ST, coalescent simulations, fungal community analyses. Development, exploratory and superseded scripts are not included. Scripts keep their cluster paths (`/home/ecastron/araucaria/...`); adjust the project root to rerun them. |
| `cluster/reconstructed/` | Producers rewritten on 18 to 19 Nov 2026 for derived files whose original code was run interactively and not saved. Each header states the verification (most outputs are byte-identical to the originals). |
| `abcrf/` | Model choice and parameter estimation with the R package abcrf 2.0, with the reference tables (`reftable_sim9.csv`, literature mutation-rate prior; `reftable_sim8.csv`, broad prior) and observed statistics. |
| `figures/` | Figure builders (Python, matplotlib). Run from the repository root; they read `figdata/` and `pca_data/` and write `figures/journal/` and `figures/supplementary/journal/`. |
| `figdata/`, `pca_data/` | The derived data the figure builders read. |
| `analysis/` | Local analyses added in the November 2026 revision (see `SCRIPT_MAP.md`, section 8). Run from the repository root; they write `figdata/revision/`. |
| `data/` | Genotypes of the 163 analyzed trees at 78,026 SNPs (`araucaria_163trees_78026snps.vcf.gz`, samples named by SRA run accession) and the ITS1 ASV table with taxonomy (`its1_asv_counts.csv.gz`, `its1_asv_taxonomy.csv.gz`). |
| `environment/` | Climate, soil and terrain extraction per tree and per stand. |
| `metadata/` | SRA run manifests. |

`SCRIPT_MAP.md` maps each analysis in the paper to its scripts in run order, with a status for
each script (final, diagnostic, superseded or exploratory) and the files it writes.

## Reproducing the figures

```bash
python figures/build_F1_journal.py      # Figure 1
python figures/build_F2_journal.py      # Figure 2
python figures/build_F4_journal.py      # Figure 3, demography (file F4_demographic_history; reads figdata/abcrf_results.json)
python figures/build_F3_journal.py      # Figure 4, traits (file F3_qst_fst_selection)
python figures/build_F5_fungal.py       # Figure 5
python figures/build_supp_journal.py    # Figures S1 to S4
python figures/range_gap.py             # distance between the Coastal and Andean forests (figdata/range_gap.json)
python figures/split_sensitivity_table.py  # Table S5, sensitivity of the split age (figdata/split_sensitivity.tsv)
Rscript abcrf/abcrf_analysis.R          # ABC-RF model choice and estimation (writes figdata/abcrf_*)
```

Python 3 with numpy, pandas, matplotlib and pyshp/pyproj (Figure 1 cache). R 4.4 or later with
abcrf 2.0, ranger, vegan, phyloseq, microbiome, LDM and genio for the cluster analyses.

## Published packages and bespoke code

Standard steps use published packages: ipyrad, bowtie2, FreeBayes, polyRAD, PLINK, sNMF (LEA),
conStruct, ngsRelate, msprime, abcrf, DADA2, vegan, LDM, phyloseq and microbiome. The following
parts are bespoke and are described in the Methods:

- held-out coverage of the ABC-RF credible intervals, posterior tail areas and the mode, the
  contact-detection curve and the zero-migration calibration (no abcrf function exists);
- the posterior-predictive check of the coalescent models;
- the clock statistics (branch-length diversity and divergence on the simulated genealogies);
- Hudson F_ST as a ratio of averages (`hudson_fst` in `cluster/abc/sumstats.py`);
- the heterozygote-deficit diagnostics and the depth-spread table behind Figure S2E;
- the low-count check on the most prevalent fungal ASVs (no negative controls were sequenced,
  so decontam cannot be applied).

## License and citation

Code is released under the MIT license (`LICENSE`); derived data under CC BY 4.0. Citation
metadata are in `CITATION.cff`. Archived versions carry a Zenodo DOI; cite all versions with the concept DOI [10.5281/zenodo.24357836](https://doi.org/10.5281/zenodo.24357836).
