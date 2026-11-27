# Script map

Maps each Results/Methods unit of the paper to the scripts in this repository, in run order. Only the scripts on the path
to the reported results and figures are included; development, exploratory and superseded scripts are not distributed.
Sources: script header comments and the project log (`ANALYSIS_SPEC.md`). Where neither states a
purpose, the script is listed under "Unclassified".

Status labels:

- **FINAL**: produces a number, table or figure panel reported in the paper, or is a required step
  on the path to one.
- **DIAGNOSTIC**: a control or mechanism test whose outcome is recorded in the project log.
  Whether each one is cited in the manuscript has not been checked here.

Path conventions:

- `cluster/` holds copies of the McCartney `scripts/` directory. The scripts keep their absolute
  cluster paths (`/home/ecastron/araucaria/...`), so `scripts/abc/x.py` in a header or in the log
  is `cluster/abc/x.py` here. Most `.sh`/`.sbatch` files are Slurm wrappers for the script with the
  same stem.
- Cluster outputs are written under `results/` on McCartney. The `figdata/` files listed below are
  copies of those outputs.
- `figures/` builders run locally from the repository root. Journal figures go to
  `figures/journal/` and supplementary figures to `figures/supplementary/journal/`.

---

## 1. Sampling coverage and the environmental contrast between ranges (Figure 1)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `figures/cache_distribution.py` | Caches the CONAF Catastro *Araucaria* distribution (dominant-type and any-component definitions) and computes sampling coverage under each | `figdata/distribution_cache.pkl`, `figdata/stand_coords.npz` | FINAL |
| 1b | `figures/range_gap.py` | Shortest distance between Coastal and Andean *Araucaria* forest polygons (dominant type and any component) and between the closest sampled stands | `figdata/range_gap.json` | FINAL |
| 2 | `environment/pull_chelsa.py` | CHELSA v2.1 bioclim variables per tree, read remotely | `environment/chelsa_per_tree.tsv` | FINAL |
| 3 | `environment/pull_soil.py` | SoilGrids 2.0 covariates per stand (REST point API) | `environment/soil_per_stand.tsv` | FINAL |
| 4 | `environment/pull_terrain.py` | Per-tree terrain covariates from SRTM 30 m | `environment/terrain_per_tree.tsv` | FINAL |
| 5 | `environment/pull_chelsa_rasters.py` | Windowed CHELSA surfaces over the F1 map extent | `figdata/chelsa_surfaces.npz` | FINAL |
| 6 | `environment/pull_soil_raster.py` | SoilGrids surfaces for the four final soil variables | `figdata/soilgrids_soil.npz` | FINAL |
| 7 | `figures/build_F1_journal.py` | Builds Figure 1 | `figures/journal/F1_distribution_sampling.{pdf,png}` | FINAL |

**Producer:** `figures/env_contrast_all.py` writes `figdata/env_contrast_all.tsv` from the files in `environment/` (stand means, Coastal = VA, TR, PN; d = Coastal minus Andean mean divided by the SD across the 12 stand means; exact permutation p over all 220 labelings). Means and d reproduce exactly; 3 of 39 permutation p-values (pH at two depths, sand), none used in the paper, differ by tie handling. `figdata/soilgrids_clay.npz` comes from `environment/pull_soil_raster.py`.

---

## 2. Read retrieval, assembly, reference-based calling and SNP filtering

### 2a. Retrieval from SRA

| order | script | purpose | status |
|---|---|---|---|
| 1 | `cluster/sra_download_array.sh` | Downloads one SRA run per array task from NCBI and converts it to gzipped FASTQ | FINAL |
| 2 | `cluster/verify_downloads.sh` | Checks FASTQ completeness (gzip integrity and other checks) | FINAL |

### 2b. De novo assembly (ipyrad)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/ipyrad_sweep.sh` | ipyrad sweep over clustering threshold on the 30 replicate-pair libraries | sweep assemblies | FINAL (threshold choice) |
| 3 | `cluster/r80_and_nulls.py` | Re-derives the clustering threshold without replicates (r80 plus null expectations) | | FINAL (replaces step 2 as the justification) |
| 4 | `cluster/ipyrad_full.sh` | Full ipyrad assembly at the chosen threshold (ct 0.91) | `full_ct091.vcf`, `.loci` | FINAL |
| 5 | `cluster/filter_vcf_callrate.py` | Applies per-analysis call-rate filters to the permissive assembly | `full_ct091.cr30.vcf` (78,802 SNPs, 177 libraries) | FINAL |
| 6 | `cluster/tstv_recovery.py` | Ts/Tv under filtering (0.48 unfiltered, 1.52 at cr30) used to validate the 30% floor | | FINAL (validation) |
| 7 | `cluster/run_q2.sh` with `cluster/catg_to_beagle.py` | CATG read counts to genotype likelihoods, then ngsRelate kinship | `results/Q2/ngsrelate.tsv` | FINAL (kinship used for relative pruning) |

**Producer:** `cluster/make_keep163.py` writes `keep.txt`: one member of each within-stand pair at ngsRelate theta >= 0.25 sharing at least 1,000 sites is removed, keeping the tree with more called genotypes (14 pairs, 177 to 163 trees). The 78,026 SNPs are the biallelic subset of the 78,802-record cr30 VCF, as filtered in `cluster/abc/load_observed.py`.

### 2c. Reference-based calling (bowtie2 and FreeBayes, AraucariaDB v1 genome)

| order | script | purpose | status |
|---|---|---|---|
| 1 | `cluster/ref/dl_genome.sh` | Downloads the AraucariaDB reference genome | FINAL |
| 2 | `cluster/ref/prep_ref.sh` | Decompresses and indexes the FASTA (samtools faidx) | FINAL |
| 3 | `cluster/ref/bt2_build.sh` | bowtie2 index | FINAL |
| 4 | `cluster/ref/bt2_map_all.sh` | bowtie2 `--very-sensitive` mapping of all 177 libraries | FINAL |
| 5 | `cluster/ref/freebayes_call.sh` | FreeBayes calling over 200 region chunks on the 177 BAMs | FINAL |
| 6 | `cluster/ref/fb_to_polyrad.sh` | Concatenates chunks and prepares polyRAD input | FINAL |
| 6a | `cluster/ref/fb_add_AD.py` | Keeps biallelic SNPs and adds an AD field from RO/AO | FINAL |
| 6b | `cluster/ref/fb_convert_only.sh`, `cluster/ref/fb_bgzip.sh` | Re-run helpers for the AD conversion and bgzip/tabix | FINAL (helpers) |
| 7 | `cluster/ref/run_polyrad_fb.sh` with `cluster/ref/polyrad_fb.R` | polyRAD genotyping of the reference calls | FINAL |

**Producer:** `cluster/make_fb_regions.py` writes `rad_regions.bed`, `rad_regions_clamped.bed` and the 200 `fb_chunks/chunk_NNN` files from `loci_vs_genome.paf` (primary alignments, query coverage >= 0.8, identity >= 0.9, 50 bp padding, merged, clamped to contig length, then GNU `split -n l/200`).

## 3. Genotype-quality diagnosis (Figures S2, S3)

### 3a. Figure data and builders (FINAL)

| script | purpose | output |
|---|---|---|
| `cluster/ref/catg_to_AD.py` | Adds AD to the ipyrad VCF for polyRAD | `ct091_cr30_AD.vcf` |
| `cluster/ref/run_polyrad_step1.sh`, `cluster/ref/run_polyrad.R` | polyRAD on the de novo calls, step 1 | `rad_raw.rds` |
| `cluster/ref/run_polyrad_step2.sh`, `cluster/ref/run_polyrad_step2.R` | polyRAD genotyping with population-structure and LD priors, then F by depth | `figdata/polyrad_F_by_depth.csv` |
| `cluster/ref/polyrad_F_fixed.R` | Corrected F-by-depth on the polyRAD object (also writes `polyrad_F_by_depth`) | `figdata/polyrad_F_by_depth.csv` |
| `cluster/ref/run_fb_F.sh`, `cluster/ref/fb_F_by_depth.R` | F by depth on the reference (FreeBayes plus polyRAD) calls, the third pipeline | `figdata/fb_F_by_depth.csv` |
| `cluster/ref/fb_fst_indiv.sh`, `cluster/ref/fb_fst_indiv.py` | Hudson FST per stand pair from reference calls, trees weighted equally (S3: de novo vs reference FST) | `figdata/fb_fst_indiv.tsv` |
| `figures/build_supp_journal.py` | Builds S1 to S4. S2 = heterozygote deficit (A to E), S3 = de novo vs reference FST | `figures/supplementary/journal/` |

**Producers:** `cluster/allele_balance.py` writes `allele_balance.txt`. `het_diagnosis.json` is assembled from the diagnostic outputs of `cluster/het_diagnosis_inputs.sbatch` and `het_diagnosis_assemble.py`, and corrected by `figures/het_diagnosis_update.py`. Figure S2E reads `figdata/s2e_mechanisms.tsv`, written by `cluster/s2e_mechanisms.sbatch` (`s2e_binomial_route.py`, `s2e_polyrad_route.R`, `s2e_assemble.py`), which also writes the binomial series of S2D (`figdata/ipyrad_binomial_F.tsv`). `figdata/polyrad_F_by_depth.csv` was written by `cluster/ref/polyrad_F_fixed.R` (identical to `results/polyrad/polyrad_F_by_depth_fixed.csv`).

### 3b. Mechanism tests and controls (DIAGNOSTIC)

Grouped by question. All read the cr30 VCF, CATG counts or reference BAMs.

The individual mechanism tests summarized in Supplementary Note 1 (null models of F versus depth, dosage, copy number, segmental duplication, lane hopping and others) are development diagnostics and are not distributed; the six whose outputs feed `figdata/het_diagnosis.json` are rerun by `cluster/het_diagnosis_inputs.sbatch`.

---

## 4. Population structure (Figure 2)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/q1_popgen.py` | Per-stand diversity, genetic distance matrix, stand table (its FST and PCA were replaced, see `q1_fixed.sh`) | `results/Q1/q1_genetic_distance.tsv`, `q1_stands.tsv` | FINAL (inputs to AMOVA) |
| 2 | `cluster/q1_amova.R` | Distance-based AMOVA (adonis2) and beta dispersion | AMOVA 15.8% among stands | FINAL |
| 3 | `cluster/q1_fixed.sh` | PLINK LD-pruned PCA, ADMIXTURE CV, and vcftools W&C FST | `pca_data/q1_pca_plink.eigen{vec,val}`, `results/Q1/q1_admixture_cv.txt`, `figdata/q1_fst_vcftools.tsv` | SUPERSEDED for Figure 2 (all 177 libraries; replaced by the 163-tree PCA in section 8); its W&C FST is SUPERSEDED by Hudson |
| 4 | `cluster/vcf2geno.py` | VCF to LEA `.geno` | `pca_data/ours91.geno.samples` | FINAL |
| 5 | `cluster/run_snmf.R` | sNMF (LEA) on the working set and on a published-style set | `pca_data/snmf_K2.Q` | SUPERSEDED for Figure 2 (all 177 libraries; replaced by `cluster/revision/snmf.R` on the 163 trees) |
| 6 | `cluster/run_construct.R` | conStruct spatial vs nonspatial models | | FINAL |
| 7 | `cluster/run_construct_cv.R` | conStruct cross-validation at 16 replicates (K = 2 vs 3) | | FINAL |
| 8 | `cluster/ref/fst_estimator_control.py` with `cluster/ref/fst_ctl.sh` | Hudson vs W&C (observed and HWE heterozygosity) on one SNP set | `figdata/fst_estimator_control.tsv` | FINAL (estimator choice) |
| 9 | `cluster/abc/mc1_fst_ci.py` | Block-jackknife CIs on every pairwise Hudson FST | `host_fst_pairs_ci.tsv` | FINAL |
| 10 | `cluster/ref/regional_fst.py` | Two-group Hudson FST, Coastal vs Andes, with jackknife SE | stdout | FINAL |
| 11 | `cluster/q4_ibd_residual.R` | IBD and the Coastal excess beyond distance (log: +0.0438, p = 0.0042) | | FINAL |
| 12 | `cluster/f3_corrected.py` | Bias-corrected f3 over all 660 triples | | FINAL |
| 13 | `cluster/abc/f3_targeted.py` | f3 restricted to the 270 triples that can detect cross-gap admixture, 50-block jackknife | | FINAL |
| 14 | `figures/build_F2_journal.py` | Builds Figure 2 (PCA, eigenvalues, Hudson FST on distance, sNMF K = 2) from the 163-tree PCA and sNMF files of section 8 (`pca_data/q1_pca_plink163.eigen{vec,val}`, `pca_data/snmf_K2_163.tsv`) | `figures/journal/F2_host_structure.{pdf,png}` | FINAL |

Reference-genome control on individual-level structure: `cluster/ref/fb_fst_indiv.py` with `cluster/ref/fb_fst_indiv.sh` (`figdata/fb_fst_indiv.tsv`, Figure S3).

**Producer:** `cluster/host_fst_pairs.py` writes `figdata/host_fst_pairs.tsv` and `stand_fst_matrix.tsv` with the project's Hudson estimator (`hudson_fst` in `cluster/abc/sumstats.py`; ratio of averages, at least 6 called chromosomes per stand) and haversine distances between stand centroids.

**SUPERSEDED or EXPLORATORY:**

---

## 5. QST-FST comparison with common-garden traits (Figure 4)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `figures/qst_fst_test.py` | QST from McIntosh et al. (2025) Table S4 nested ANOVA (values hardcoded) across a range of h2, against the regional Hudson FST | `figdata/qst_fst.json` | FINAL |
| 2 | `figures/build_F4_journal.py` | Builds Figure 4 (QST across h2, per-trait crossover heritability) | `figures/journal/F4_qst_fst_selection.{pdf,png}` | FINAL |
| | `figures/qst_fst.py` | Earlier QST vs FST using `figdata/host_fst_pairs.tsv` (QST as lower bounds) | | SUPERSEDED |

Note: `qst_fst_test.py` hardcodes the regional FST (0.1059, SE 0.0011). `cluster/ref/regional_fst.py`
computes that quantity, but no file records the link.

---

## 6. Demographic inference (Figure 3, Figure S1, Tables S3 to S5)

### 6a. FINAL chain

| order | script | purpose | output |
|---|---|---|---|
| 1 | `cluster/abc/sumstats2.py` | Scale-free summary statistics (no SFS, no heterozygosity), identical on observed and simulated data | module |
| 2 | `cluster/abc/chain.py` | 1D stepping-stone chain geometry for the two ranges | module |
| 3 | `cluster/abc/sim5.py` | Base chain simulator for the four histories. `sim9.py` imports it | module |
| 4 | `cluster/abc/observed2.py` | Observed shape-statistic vector and missingness mask | `results/Q1/abc/observed2.pkl` |
| 5 | `cluster/clock_da.py` with `cluster/clock_da.sbatch` | Model-free net divergence d_a and the three clock statistics (pi Coastal, pi Andean, d_xy), on 163 trees | `figdata/clock_da_dates.json`. The clock values are hardcoded into `f4_numbers8.py`/`f4_numbers9.py` |
| 6 | `cluster/abc/sim9.py` | sim8 with the mutation-rate prior narrowed to the gymnosperm range (primary) | module |
| 7 | `cluster/abc/run_sims9.py` with `cluster/abc/prod9.sbatch` | Slurm array: prior draw, simulate, summarize | `results/Q1/abc_clock_lit/s9_*.npz` |
| 8 | `cluster/abc/sim8.py`, `cluster/abc/run_sims8.py`, `cluster/abc/prod8.sbatch` | Same with the broad mutation-rate prior (sensitivity) | `s8_*.npz` |
| 9 | `cluster/abc/f4_numbers9.py`, `cluster/abc/f4_numbers8.py` | REQUIRED as data loaders: `export_reftables.py` executes their top section (loading, seed-0 balancing, prior-predictive statistic filter). Their scikit-learn analysis below that point is SUPERSEDED | |
| 10 | `cluster/abc/export_reftables.py` | Exports the balanced sim9 and sim8 reference tables as CSV | `abcrf/reftable_sim9.csv`, `observed_sim9.csv`, `reftable_sim8.csv`, `observed_sim8.csv` |
| 11 | `abcrf/abcrf_analysis.R` | abcrf 2.0 model choice (grouped and four-way, 10 replicate forests) and regAbcrf parameter estimation | `figdata/abcrf_results.json`, `figdata/abcrf_weights_{avg,iso,con,broad}.csv` |
| 12 | `abcrf/abcrf_logscale.R` | Control: regAbcrf on log10 scale parameters | `figdata/abcrf_logscale.json` |
| 13 | `figures/build_F3_journal.py` | Builds Figure 3 from `abcrf_results.json` and the weights | `figures/journal/F3_demographic_history.{pdf,png}` |
| 14 | `figures/build_supp_journal.py` | Builds S1 (model choice and chain architecture) from `abcrf_results.json` and the sim9 tables | `figures/supplementary/journal/` |
| 15 | `figures/split_sensitivity_table.py` | Table S5: split age under each sensitivity analysis, from `abcrf_results.json` and `abcrf_logscale.json` | `figdata/split_sensitivity.tsv` |

**Tables S3 to S5** are written by the supplementary-document builder, which is not distributed here because it holds manuscript text. Table S3 lists the priors coded in `cluster/abc/sim9.py`; Table S4 reads `figdata/abcrf_results.json` (written by `abcrf/abcrf_analysis.R`); Table S5 reads `figdata/split_sensitivity.tsv`.

## 7. Fungal endophyte community from ITS1 (Figure 5, Figure S4)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/run_dada2.sh` with `cluster/dada2_step1.R` | DADA2 ASV inference per amplicon (16S and ITS, separate parameters) | `microbiome/ITS/seqtab_nochim.rds`, `taxonomy.rds` | FINAL (this is the table `fungal_fig.R` reads) |
| 2 | `cluster/run_its_r1.sh` with `cluster/dada2_its_r1.R` | ITS forward-reads-only run, declared the production run on 4 Aug 2026 | `seqtab_nochim_r1.rds`, `taxonomy_r1.rds` | see note |
| 4 | `cluster/final_protocol.R` | Final endophyte protocol: MRA 2020 settings plus organelle removal, sample floor, prevalence and abundance filter | `results/Q5/final_protocol.txt` | FINAL (protocol) |
| 5 | `cluster/fungal_fig.R` with `cluster/fungal_fig.sh` | F5 statistics at tree level: dbRDA on Bray-Curtis, PERMANOVA, PERMDISP, order composition, Shannon, host-genetic test on SNP PCs (genio `read_plink`) | `figdata/fungfig_r1/panel*.csv` (paper; `figdata/fungfig/` = merged table), `fungal_fig.txt` | FINAL |
| 6 | `cluster/disp_dir.R` | PERMDISP direction at tree level (which group is less dispersed) | stdout | FINAL (dispersion caveat) |
| 7 | `cluster/ord_check.R` | Sequential vs marginal tests for range and compartment | stdout | FINAL (term orthogonality) |
| 8 | `cluster/ldm_run.R` with `cluster/ldm_run.sh` | LDM range test at tree and library level | `figdata/fungfig/ldm.txt` | FINAL |
| 9 | `cluster/ldm_taxa.R` with `cluster/ldm_taxa.sh` | LDM taxon-level FDR and what the VE field means | `figdata/fungfig/ldm_taxa.{csv,txt}` | FINAL |
| 10 | `figures/build_F5_fungal.py` | Builds Figure 5 | `figures/journal/F5_fungal_community.{pdf,png}` | FINAL |
| 11 | `cluster/core_microbiome.R` | Core community with microbiome 1.26 and phyloseq 1.48 (ASV, family; raw and rarefied) | `figdata/core/` | FINAL |
| 12 | `figures/build_supp_journal.py` | Builds S4 from `figdata/its_prevalence_r1.csv` and `figdata/core/r1_core_grid.tsv` | `figures/supplementary/journal/` | FINAL |

**ITS table (resolved 19 Nov 2026):** the paper uses the forward-read table (`seqtab_nochim_r1.rds`, `dada2_its_r1.R`), as decided on 4 August 2026 by a rule fixed before the comparison. The fungal scripts take the table as an argument or `ITS_TABLE` (`merged` by default, which reproduces the earlier merged-pair outputs byte for byte; `r1` for the paper). Paper outputs: `figdata/fungfig_r1/`, `figdata/its_prevalence_r1.csv`, `figdata/core/`.

**Producers:** `cluster/its_figdata.R` writes the `its_*` tables and the S4 prevalence. The core community (Figure S4B) is computed with the R package microbiome on phyloseq objects by `cluster/core_microbiome.R` (ASV and family, all reads and rarefied to equal depth). ITS needed no decontamination step (all ASVs assign to Fungi); no negative controls were sequenced, so decontam does not apply.

**EXPLORATORY or SUPERSEDED (microbiome):**

---


## 8. Revision analyses (November 2026)

Robustness checks and additions made in the November 2026 revision. Cluster scripts live in
`cluster/revision/` (on McCartney `scripts/revision/`, outputs in `results/revision/`); local scripts
in `analysis/` run from the repository root and write `figdata/revision/`.

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/revision/structure.py` with `cluster/revision/structure.sbatch` | 163-tree genotype files for sNMF, LD-pruned PLINK PCA (`--indep-pairwise 50 5 0.2`, 66,131 SNPs), one-SNP-per-locus set, ADMIXTURE CV | `results/revision/pca_ld.eigen{vec,val}` (copied to `pca_data/q1_pca_plink163.eigen{vec,val}`), `structure.json` | FINAL (Figure 2A, B) |
| 2 | `cluster/revision/snmf.R` (called by `structure.sbatch`) | sNMF K = 1 to 6 on the 163 trees (10 repetitions, alpha 100) | `snmf_all163_K2.Q` (copied to `pca_data/snmf_K2_163.tsv`), `figdata/revision/snmf_all163.tsv` | FINAL (Figure 2E) |
| 3 | `cluster/revision/ancestry_null.py`, `ancestry_null.R` with `ancestry_null.sbatch` | sNMF K = 2 on 50 two-refugia (no migration) and 50 secondary-contact simulations at the parameter sets of the sim9 reference table nearest the observed statistics, and on the observed data | `figdata/revision/ancnull/` | FINAL |
| 4 | `analysis/ancestry_null_summary.py` | Observed sNMF coefficients against the simulated distributions | `figdata/revision/ancestry_null_summary.json` | FINAL (Results, ancestry) |
| 5 | `analysis/ibd_mrm.R` | Multiple regression on distance matrices (`ecodist::MRM`): FST/(1 - FST) on log distance and a binary different-range matrix, with and without VA | `figdata/revision/ibd_mrm.json` | FINAL (Results, isolation by distance) |
| 6 | `analysis/abcrf_diagnostics.R` | ABC-RF variable importance, posterior calibration on held-out simulations, prior probabilities and Bayes factors of the split-age events | `figdata/revision/abcrf_*` | FINAL (Results, split age priors) |
| 7 | `analysis/env_stand_values.py` | Stand values of the environmental variables and leave-one-stand-out d | `figdata/revision/env_stand_values.tsv`, `env_loo.tsv` | FINAL (Table S6) |
| 8 | `cluster/revision/qc_table.sh` | Per-tree read counts, loci, call rate and bowtie2 mapping statistics | `figdata/revision/qc_per_tree.tsv` | FINAL (Table S7) |
| 9 | `cluster/revision/its.R` with `cluster/revision/its.sbatch` | Endophyte range effect with the stand as the unit (all 36 assignments of two of nine stands to the Coastal range); Shannon tests at four rarefaction depths | `figdata/revision/its_stand_perm.tsv`, `its_rarefaction_sensitivity.tsv` | FINAL (Results and Figure 5 legend) |
| 10 | `cluster/revision/export_genotypes.py`, `export_asv_table.R` with `export_data.sbatch` | Genotypes of the 163 analyzed trees at the 78,026 biallelic SNPs (all FORMAT fields; samples named by SRA run) and the ITS1 ASV table with UNITE taxonomy, for archiving | `data/araucaria_163trees_78026snps.vcf.gz`, `data/its1_asv_counts.csv.gz`, `data/its1_asv_taxonomy.csv.gz` | FINAL (Data Accessibility) |

---

**Input metadata (not derived):** `figdata/popmap_177.tsv` (library to stand map) and `metadata/` (SRA run manifests; the `collection_date` column is added by `metadata/add_collection_dates.py` from the BioSample records). `figdata/stand_coords.npz` is rebuilt by `figures/stand_coords.py` from `metadata/sra_runs_radseq.csv` and `figdata/keep163.txt`.
