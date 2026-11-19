# Script map

Maps each Results/Methods unit of the paper to the scripts in this repository, in run order.
Sources: script header comments and the project log (`ANALYSIS_SPEC.md`). Where neither states a
purpose, the script is listed under "Unclassified".

Status labels:

- **FINAL**: produces a number, table or figure panel reported in the paper, or is a required step
  on the path to one.
- **DIAGNOSTIC**: a control or mechanism test whose outcome is recorded in the project log.
  Whether each one is cited in the manuscript has not been checked here.
- **SUPERSEDED**: replaced by a later script or a later decision; kept as history.
- **EXPLORATORY**: a side analysis that did not enter the paper.

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
| 2 | `environment/pull_chelsa.py` | CHELSA v2.1 bioclim variables per tree, read remotely | `environment/chelsa_per_tree.tsv` | FINAL |
| 3 | `environment/pull_soil.py` | SoilGrids 2.0 covariates per stand (REST point API) | `environment/soil_per_stand.tsv` | FINAL |
| 4 | `environment/pull_terrain.py` | Per-tree terrain covariates from SRTM 30 m | `environment/terrain_per_tree.tsv` | FINAL |
| 5 | `environment/pull_chelsa_rasters.py` | Windowed CHELSA surfaces over the F1 map extent | `figdata/chelsa_surfaces.npz` | FINAL |
| 6 | `environment/pull_soil_raster.py` | SoilGrids surfaces for the four final soil variables | `figdata/soilgrids_soil.npz` | FINAL |
| 7 | `figures/build_F1_journal.py` | Builds Figure 1 | `figures/journal/F1_distribution_sampling.{pdf,png}` | FINAL |

`cluster/pull_chelsa.py`, `cluster/pull_chelsa_rasters.py` and `cluster/pull_soil_raster.py` are
byte-identical to the `environment/` copies.

**Producer (reconstructed 19 Nov 2026):** `figures/reconstructed/env_contrast_all.py` writes `figdata/env_contrast_all.tsv` from the files in `environment/` (stand means, Coastal = VA, TR, PN; d = Coastal minus Andean mean divided by the SD across the 12 stand means; exact permutation p over all 220 labelings). Means and d reproduce exactly; 3 of 39 permutation p-values (pH at two depths, sand), none used in the paper, differ by tie handling. `figdata/soilgrids_clay.npz` comes from `cluster/pull_soil_raster.py`.

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
| 2 | `cluster/analyse_ipyrad.py` | Replicate concordance across the ipyrad sweep | | SUPERSEDED (the "replicate pairs" turned out to be relatives) |
| 3 | `cluster/r80_and_nulls.py` | Re-derives the clustering threshold without replicates (r80 plus null expectations) | | FINAL (replaces step 2 as the justification) |
| 4 | `cluster/ipyrad_full.sh` | Full ipyrad assembly at the chosen threshold (ct 0.91) | `full_ct091.vcf`, `.loci` | FINAL |
| 5 | `cluster/filter_vcf_callrate.py` | Applies per-analysis call-rate filters to the permissive assembly | `full_ct091.cr30.vcf` (78,802 SNPs, 177 libraries) | FINAL |
| 6 | `cluster/tstv_recovery.py` | Ts/Tv under filtering (0.48 unfiltered, 1.52 at cr30) used to validate the 30% floor | | FINAL (validation) |
| 7 | `cluster/run_q2.sh` with `cluster/catg_to_beagle.py` | CATG read counts to genotype likelihoods, then ngsRelate kinship | `results/Q2/ngsrelate.tsv` | FINAL (kinship used for relative pruning) |

**Producer (reconstructed 19 Nov 2026):** `cluster/reconstructed/make_keep163.py` writes `keep.txt` (byte-identical): one member of each within-stand pair at ngsRelate theta >= 0.25 sharing at least 1,000 sites is removed, keeping the tree with more called genotypes (14 pairs, 177 to 163 trees). The 78,026 SNPs are the biallelic subset of the 78,802-record cr30 VCF, as filtered in `cluster/abc/load_observed.py`.

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

**Producer (reconstructed 19 Nov 2026):** `cluster/reconstructed/make_fb_regions.py` writes `rad_regions.bed`, `rad_regions_clamped.bed` and the 200 `fb_chunks/chunk_NNN` files from `loci_vs_genome.paf` (primary alignments, query coverage >= 0.8, identity >= 0.9, 50 bp padding, merged, clamped to contig length, then GNU `split -n l/200`); all byte-identical to the originals.

Earlier reference-mapping work on an 8-library subset (minimap2): `cluster/ref/index_genome.sh`,
`map_balance.sh`, `bt2_map.sh`, `bwa_check.sh`, `dedup_all.sh`, `feasibility.sh`, `callsnps.sh`
(bcftools on callable positions). These are DIAGNOSTIC and feed unit 3.

### 2d. SUPERSEDED: Stacks assembly (replaced by ipyrad)

`cluster/fix_deflines.sh` (Stacks-compatible read names), `cluster/rad_clean_array.sh`
(process_radtags), `cluster/rad_sweep_array.sh` (M sweep), `cluster/minitest_denovo.sh`,
`cluster/test_dedup.sh`, `cluster/rerun_populations.sh`, `cluster/rerun_populations_r30.sh`,
`cluster/analyse_sweep.py`, `cluster/pick_and_assemble.sh`, `cluster/snp_filter_grid.py` with
`cluster/run_filter_grid.sh` (filter grid scored by replicate concordance),
`cluster/run_filter_and_retest.sh`.

---

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

**Producers (reconstructed 19 Nov 2026):** `cluster/reconstructed/allele_balance.py` writes `allele_balance.txt` (byte-identical). `het_diagnosis.json` is assembled from rerun diagnostics by `cluster/reconstructed/het_diagnosis_inputs.sbatch` and `het_diagnosis_assemble.py`, and corrected by `figures/reconstructed/het_diagnosis_update.py`. Figure S2E reads `figdata/s2e_mechanisms.tsv`, written by `cluster/reconstructed/s2e_mechanisms.sbatch` (`s2e_binomial_route.py`, `s2e_polyrad_route.R`, `s2e_assemble.py`), which also writes the binomial series of S2D (`figdata/ipyrad_binomial_F.tsv`). `figdata/polyrad_F_by_depth.csv` was written by `cluster/polyrad_F_fixed.R` (identical to `results/polyrad/polyrad_F_by_depth_fixed.csv`).

### 3b. Mechanism tests and controls (DIAGNOSTIC)

Grouped by question. All read the cr30 VCF, CATG counts or reference BAMs.

- **F versus depth, with nulls:** `cluster/ref/F_vs_depth_fine.py` (fine-resolution F_IS vs depth
  with matched HWE null), `cluster/ref/null_F_by_depth.py` (HWE null with observed depth and calling
  rule), `cluster/simulate_HWE_F.py` with `cluster/run_sim_F.sh` and `cluster/run_sim_F2.sh`,
  `cluster/ref/paired_depth_F.py` (same loci at shallow and deep depth), `cluster/F_diagnosis.R`
  (reconciles the F-like quantities), `cluster/F_extrapolate.R` (F remaining at high library
  quality), `cluster/F_by_depth.sh`, `cluster/F_rand_control.sh`, `cluster/gl_highdepth.py`,
  `cluster/catg_to_beagle_mindp.py`, `cluster/catg_to_beagle_rand.py` (depth-floor and
  random-thinning controls for ngsRelate F).
- **Rarefaction:** `cluster/ref/rarefy_and_cv.py` (subsample reads to a fixed ceiling, then
  stratify F by original depth).
- **Allele balance:** `cluster/ref/balance_from_pileup.py` (depth-stratified balance from
  mpileup), `cluster/ref/balance_denovo.py` with `cluster/ref/run_denovo.sh` (de novo control),
  `cluster/ref/simulate_reads.py` with `cluster/ref/simcontrol.sh` (positive control for the
  overdispersion statistic), `cluster/ref/two_phenomena.py` (deficit and flat balance, one
  phenomenon or two), `cluster/shape_test.py`, `cluster/shape_test_p.py` (depth-corrected shape
  test).
- **Wahlund and HWE:** `cluster/ref/hwe_by_maf.py`, `cluster/ref/hwe_per_stand.py`,
  `cluster/run_sim_wahlund.sh`, `cluster/run_sim_wahlund50.sh`.
- **Paralogs, copy number, segmental duplications:** `cluster/ref/hdplot.py` with
  `cluster/ref/run_hdplot.sh` (writes `figdata/hdplot.tsv`), `cluster/ref/extract_loci.sh` with
  `cluster/ref/loci_to_fasta.py`, `cluster/ref/map_loci_to_genome.sh`,
  `cluster/ref/F_by_copynumber.py`, `cluster/ref/selfaln_chunk.py`,
  `cluster/ref/selfaln_index.sh`, `cluster/ref/selfaln_map.sh`, `cluster/ref/sd_label_loci.py`
  with `cluster/ref/sd_label.sh`, `cluster/ref/F_by_sd.py` with `cluster/ref/f_sd.sh`,
  `cluster/ref/cnv_filter.py`, `cluster/ref/locus_depth_filter.py`,
  `cluster/ref/locus_structural_filter.R`, `cluster/ref/hobs_paralog_filter.R` with
  `cluster/ref/run_hobs_filter.sh`, `cluster/ref/capped_depth_HF.py`,
  `cluster/ref/cap_and_balance_filter.py`, `cluster/ref/nupt_scan.sh`.
- **Hemizygosity and dosage:** `cluster/ref/dosage_bed.py`, `cluster/ref/dosage_cov.sh`,
  `cluster/ref/dosage_test.py` with `cluster/ref/dosage_run.sh`, `cluster/ref/dosage_control.py`
  with `cluster/ref/dosage_ctl.sh`.
- **Null alleles:** `cluster/ref/null_allele_test.py`.
- **Index hopping:** `cluster/ref/get_lanes.sh`, `cluster/ref/lane_hopping_test.py` with
  `cluster/ref/run_hop.sh`.
- **PCR duplicates:** `cluster/ref/clone_filter.sh`, `cluster/ref/clonemap.sh`,
  `cluster/ref/dupe_test.sh`, `cluster/ref/downsample_control.sh`.
- **Genotype calling rules:** `cluster/ref/recall_threshold.py`, `cluster/ref/betabinom_recall.py`,
  `cluster/ref/maf_from_reads.py`.
- **Organellar haplotypes:** `cluster/ref/cp_haplotypes.sh`.
- **polyRAD install:** `cluster/ref/install_polyrad.sh`.
- **Earlier heterozygote-deficit diagnosis:** `cluster/diagnose_het_deficit.py`.

Local diagnostic figure builders (presentation figures, not in the journal set; write to the
repository root): `figures/build_F_depth.py` (`figdata/F_vs_depth.tsv`),
`figures/build_hdplot.py` (`figdata/hdplot.tsv`), `figures/build_hobs_filter_F.py`
(`figdata/hobs_filter_F.tsv`), `figures/build_locus_filter_F.py` (`figdata/locus_filter_F.tsv`),
`figures/build_polyrad_F.py`, `figures/build_pipeline_compare.py` (three pipelines, same F vs depth
curve). Producers of `F_vs_depth.tsv`, `hobs_filter_F.tsv` and `locus_filter_F.tsv` are not
identified.

---

## 4. Population structure (Figure 2)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/q1_popgen.py` | Per-stand diversity, genetic distance matrix, stand table (its FST and PCA were replaced, see `q1_fixed.sh`) | `results/Q1/q1_genetic_distance.tsv`, `q1_stands.tsv` | FINAL (inputs to AMOVA) |
| 2 | `cluster/q1_amova.R` | Distance-based AMOVA (adonis2) and beta dispersion | AMOVA 15.8% among stands | FINAL |
| 3 | `cluster/q1_fixed.sh` | PLINK LD-pruned PCA, ADMIXTURE CV, and vcftools W&C FST | `pca_data/q1_pca_plink.eigen{vec,val}`, `results/Q1/q1_admixture_cv.txt`, `figdata/q1_fst_vcftools.tsv` | FINAL for PCA and ADMIXTURE; its W&C FST is SUPERSEDED by Hudson |
| 4 | `cluster/vcf2geno.py` | VCF to LEA `.geno` | `pca_data/ours91.geno.samples` | FINAL |
| 5 | `cluster/run_snmf.R` | sNMF (LEA) on the working set and on a published-style set | `pca_data/snmf_K2.Q` | FINAL |
| 6 | `cluster/run_construct.R` | conStruct spatial vs nonspatial models | | FINAL |
| 7 | `cluster/run_construct_cv.R` | conStruct cross-validation at 16 replicates (K = 2 vs 3) | | FINAL |
| 8 | `cluster/ref/fst_estimator_control.py` with `cluster/ref/fst_ctl.sh` | Hudson vs W&C (observed and HWE heterozygosity) on one SNP set | `figdata/fst_estimator_control.tsv` | FINAL (estimator choice) |
| 9 | `cluster/abc/mc1_fst_ci.py` | Block-jackknife CIs on every pairwise Hudson FST | `host_fst_pairs_ci.tsv` | FINAL |
| 10 | `cluster/ref/regional_fst.py` | Two-group Hudson FST, Coastal vs Andes, with jackknife SE | stdout | FINAL |
| 11 | `cluster/q4_ibd_residual.R` | IBD and the Coastal excess beyond distance (log: +0.0438, p = 0.0042) | | FINAL |
| 12 | `cluster/f3_corrected.py` | Bias-corrected f3 over all 660 triples | | FINAL |
| 13 | `cluster/abc/f3_targeted.py` | f3 restricted to the 270 triples that can detect cross-gap admixture, 50-block jackknife | | FINAL |
| 14 | `figures/build_F2_journal.py` | Builds Figure 2 (PCA, eigenvalues, Hudson FST on distance, sNMF K = 2) | `figures/journal/F2_host_structure.{pdf,png}` | FINAL |

Reference-genome controls on structure (FINAL, cited in the log as replicating structure):
`cluster/ref/fb_fst_pairs.py` with `cluster/ref/fb_fst.sh` (`figdata/fb_fst_pairs.tsv`),
`cluster/ref/fb_fst_indiv.py` (unit 3), `cluster/ref/fst_by_copynumber.py` (FST on single-copy vs
multi-copy loci).

**Producer (reconstructed 19 Nov 2026):** `cluster/reconstructed/host_fst_pairs.py` writes `figdata/host_fst_pairs.tsv` and `stand_fst_matrix.tsv` (both byte-identical) with the project's Hudson estimator (`hudson_fst` in `cluster/abc/sumstats.py`; ratio of averages, at least 6 called chromosomes per stand) and haversine distances between stand centroids.

**Sensitivity (FINAL, dp20 set):** `cluster/dp20/mask_depth.py` (genotypes with DP < 20 set to
missing), `cluster/dp20/run_dp20.sh`, then path-edited reruns `cluster/dp20/q1_popgen.py`,
`cluster/dp20/q1_amova.R`, `cluster/dp20/q2_sgs.R`, `cluster/dp20/q4_ibd_residual.R`.

**SUPERSEDED or EXPLORATORY:**

- `cluster/run_q1.sh`: first Q1 run, replaced by `q1_fixed.sh` (negative FST mean).
- `cluster/f3_admixture.py`: f3 without bias correction, replaced by `f3_corrected.py`.
- `cluster/isolation_by_distance.py`: early check that the stand effect is real structure.
- `cluster/pop_structure_prelim.py`, `cluster/pop_structure_full.py`: early checks that stand-level
  structure survives unreliable genotypes.
- `cluster/pub85_faithful.py`, `cluster/pub85_replicate.sh`: replication of the Varas-Myrik et al.
  (2022) filter chain.
- `cluster/run_feems.py`, `cluster/run_feems2.py`, `cluster/make_tri_grid.py`: feems, judged not
  applicable (disjoint ranges).
- `cluster/vcf_to_treemix.py`, `cluster/run_treemix.sh`, `cluster/treemix_summary.R`: TreeMix.
- Presentation network panels: `figures/build_neighbornet.R` (`figdata/stand_fst_matrix.tsv` to
  `figdata/nnet_*.tsv`), `figures/build_networks.py` (NeighborNet of stands plus MSN of
  individuals). `cluster/ref/indiv_distance.py` with `cluster/ref/run_indiv_dist.sh` writes
  `figdata/indiv_dist.tsv`. The producer of `stand_fst_matrix.tsv` is not identified.

---

## 5. QST-FST comparison with common-garden traits (Figure 3)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `figures/qst_fst_test.py` | QST from McIntosh et al. (2025) Table S4 nested ANOVA (values hardcoded) across a range of h2, against the regional Hudson FST | `figdata/qst_fst.json` | FINAL |
| 2 | `figures/build_F3_journal.py` | Builds Figure 3 (QST across h2, per-trait crossover heritability) | `figures/journal/F3_qst_fst_selection.{pdf,png}` | FINAL |
| | `figures/qst_fst.py` | Earlier QST vs FST using `figdata/host_fst_pairs.tsv` (QST as lower bounds) | | SUPERSEDED |

Note: `qst_fst_test.py` hardcodes the regional FST (0.1059, SE 0.0011). `cluster/ref/regional_fst.py`
computes that quantity, but no file records the link.

---

## 6. Demographic inference (Figure 4, Figure S1, Tables S3 and S4)

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
| 13 | `figures/build_F4_journal.py` | Builds Figure 4 from `abcrf_results.json` and the weights | `figures/journal/F4_demographic_history.{pdf,png}` |
| 14 | `figures/build_supp_journal.py` | Builds S1 (model choice and chain architecture) from `abcrf_results.json` and the sim9 tables | `figures/supplementary/journal/` |

**Tables S3 and S4** are written by the supplementary-document builder, which is not distributed here because it holds manuscript text. Table S3 lists the priors coded in `cluster/abc/sim9.py`; Table S4 reads `figdata/abcrf_results.json` (written by `abcrf/abcrf_analysis.R`).

### 6b. SUPERSEDED: scikit-learn ABC-RF analyses (replaced by abcrf, 18 Nov 2026)

`cluster/abc/train.py`, `cluster/abc/final_abc.py` (with `rerun_final.sh`),
`cluster/abc/controls.py` (with `rerun_controls.sh`), `cluster/abc/power.py`,
`cluster/abc/contrasts.py`, `cluster/abc/contrasts4.py`, `cluster/abc/confound6.py`,
`cluster/abc/f4_numbers.py` (with `f4_num.sh`), `cluster/abc/f4_numbers6.py` (with
`f4_num6.sbatch`), `cluster/abc/f4_numbers7.py` (with `f4_num7.sbatch`), the analysis sections of
`f4_numbers8.py` and `f4_numbers9.py` (with `f4_num8.sbatch`, `f4_num9.sbatch`),
`cluster/abc/f4_fig9.py` (with `fig9.sbatch`), `cluster/abc/f4_conv9.py` (with `conv9.sbatch`),
`cluster/abc/f4_contact9.py` (with `con9.sbatch`), `cluster/abc/f4_lgm8.py` (with
`f4_lgm8.sbatch`), `cluster/abc/f4_decomp8.py` (with `f4_decomp8.sbatch`).

### 6c. SUPERSEDED: simulator development history

- Simulators: `cluster/abc/simulate.py` (first three-hypothesis model), `cluster/abc/lattice.py`
  and `cluster/abc/sim2.py` (2D lattice, no IBD), `cluster/abc/sim4.py`, `cluster/abc/sim5.py`
  (still imported by sim9), `cluster/abc/sim6.py` (T1 floor 100 generations),
  `cluster/abc/sim7.py` (floor 10), `cluster/abc/sim8_unlinked.py`.
- Production wrappers: `cluster/abc/run_sims.py`, `run_sims5.py`, `run_sims6.py`, `run_sims7.py`,
  `prod.sbatch`, `prod5.sbatch`, `prod6.sbatch`, `prod7.sbatch`, `prodh1.sbatch`,
  `prodh1b.sbatch`.
- Earlier statistics and observed vectors: `cluster/abc/sumstats.py`, `cluster/abc/observed.py`,
  `cluster/abc/load_observed.py`.
- Feasibility and diagnosis runs: `cluster/abc/feas.py`, `feas2.py`, `feas3.py`, `feas4.py`,
  `feas8.py`, `diag.py`, `diag2.py`, `diag3.py`, `diag4.py`, `rectest.py` (linkage as the source of
  FST scatter), `bench8.py` (timing and noise of clock statistics), each with its `.sbatch`.
- `cluster/abc/ancestry_check.py`: whether a no-gene-flow model reproduces the K = 2 barplot
  (sklearn-era control).
- MIGRATE-N, excluded from the manuscript: `cluster/make_migrate_infile.py`,
  `cluster/run_migrate_long.sbatch`.

---

## 7. Fungal endophyte community from ITS1 (Figure 5, Figure S4)

| order | script | purpose | output | status |
|---|---|---|---|---|
| 1 | `cluster/run_dada2.sh` with `cluster/dada2_step1.R` | DADA2 ASV inference per amplicon (16S and ITS, separate parameters) | `microbiome/ITS/seqtab_nochim.rds`, `taxonomy.rds` | FINAL (this is the table `fungal_fig.R` reads) |
| 2 | `cluster/run_its_r1.sh` with `cluster/dada2_its_r1.R` | ITS forward-reads-only run, declared the production run on 4 Aug 2026 | `seqtab_nochim_r1.rds`, `taxonomy_r1.rds` | see note |
| 3 | `cluster/decontam_16S.R` | Removes host chloroplast and mitochondrial ASVs from 16S | | SUPERSEDED for the paper (16S dropped) |
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

**Producers:** `cluster/reconstructed/its_figdata.R` writes the `its_*` tables and the S4 prevalence (reconstructed; the merged-table outputs are byte-identical or match to 1e-13). The core community (Figure S4B) is computed with the R package microbiome on phyloseq objects by `cluster/core_microbiome.R` (ASV and family, all reads and rarefied to equal depth). ITS needed no decontamination step (all ASVs assign to Fungi); no negative controls were sequenced, so decontam does not apply.

**EXPLORATORY or SUPERSEDED (microbiome):**

- ITS processing choices: `cluster/its_filter_diag.R`, `cluster/its_maxee_test.R`,
  `cluster/its_r1_vs_merged.R`, `cluster/its_r1_taxcomm.R` with `cluster/its_r1_taxcomm.sh`,
  `cluster/export_asvs.R`.
- 16S (dropped from the paper): `cluster/phylo_16S.sh`, `cluster/phylo_16S_gg2.sh`,
  `cluster/compare_three_trees.R`, `cluster/prune_and_compare_trees.R`,
  `cluster/threshold_table.R`, `cluster/depth_distribution.R`, `cluster/organelle_haplotypes.R`,
  `cluster/q6_organelle_control.R`, `cluster/q6_organelle_control2.R`.
- Above-ground pooling and filter settings: `cluster/aboveground_pool.R`,
  `cluster/ag_250_filtered.R`, `cluster/ag_250_nofilter.R`, `cluster/q6_filter_gradient.R`,
  `cluster/q6_gradient_full.R`.
- Q5 community partition (earlier framing): `cluster/q5_community.R` with `cluster/run_q5.sh`,
  `cluster/q5_rerun_clean.R`, `cluster/q9_microbiome_ibd.R` with `cluster/run_ibd.sh`.
- Q6 kinship-based host-genetic tests, replaced by the SNP-PC test in `fungal_fig.R` panel D:
  `cluster/q6_extended_phenotype.R`, `cluster/q6_rerun_clean.R`, `cluster/q6_varpart.R`,
  `cluster/q6_aboveground_varpart.R`, `cluster/q6_ksens_mrm.R` with `cluster/run_ksens.sh`,
  `cluster/q6_power.R` with `cluster/run_power.sh`, `cluster/q6_heritability.R` with
  `cluster/run_herit.sh`, `cluster/q6_family_level.R` with `.sh`, `cluster/q6_family2.R` with
  `.sh`, `cluster/q6_null_values.R` with `.sh`, `cluster/q6_core_sweep.R`, and the dp20 and
  no-floor variants `cluster/dp20/q6_ksens_mrm.R`, `cluster/dp20/q6_power.R`,
  `cluster/dp20/q6_rerun_clean.R`, `cluster/dp20/run_q6_dp20.sh`, `cluster/nofloor/q6_baseline.R`,
  `cluster/nofloor/q6_dp20.R`.
- Cordillera tests with environmental covariates: `cluster/q6_paired_cordillera.R`,
  `cluster/q6_paired_cordillera_lit.R` with `cluster/q6_lit.sh`,
  `cluster/q6_paired_cordillera_soil.R` with `cluster/q6_soil.sh`,
  `cluster/q6_paired_cordillera_pc.R` with `cluster/q6_pc.sh`.
- Old Figure 1C (removed panel): `cluster/f1c_family.R`, `cluster/f1c_robust.R`,
  `cluster/f1c_varpart.R`, each with its `.sh`.
- Functional lifestyle (curiosity test, not in the paper): `cluster/ft_probe.R`,
  `cluster/q6_lifestyle.R` with `cluster/q6_lifestyle.sh` (`figdata/q6_lifestyle.txt`,
  `null_strict.tsv`, `null_propagated.tsv`).
- Q7 and Q8 (folded in or withdrawn): `cluster/q7_q8.R` with `cluster/run_q7q8.sh`.
- Model-based GLM (retracted for non-convergence): `cluster/manyglm_q5q6.R`,
  `cluster/manyglm_decomp.R`, `cluster/run_mglm.sh`.
- Compartment identity of a disputed 2018 sample: `cluster/compartment_test.sh`.

---

## Other scripts

### Exploratory diagnostics: replicate pairs and genotyping error (pre-ipyrad-final phase)

The 15 "replicate pairs" were later shown to be related trees, not replicates. Scripts built on
the replicate criterion are SUPERSEDED.

`cluster/validate_replicate_criterion.py`, `cluster/are_pairs_identical.py`,
`cluster/write_validated_pairs.py`, `cluster/find_best_matches.py`,
`cluster/profile_anomalous_libs.py`, `cluster/depth_confound_test.py`,
`cluster/depth_vs_stand.py`, `cluster/diagnose_p1_failure.py`,
`cluster/test_dropout_hypothesis.py`, `cluster/test_locus_concentration.py`,
`cluster/dropout_vs_nonhomology_cr30.py`, `cluster/verify_calling_mechanism.py`,
`cluster/crossval_locus_filter.py`, `cluster/clean_loci.py`, `cluster/excess_het_filter.py`,
`cluster/simulate_relatedness.py`, `cluster/calibrate_error_model.py`.

### Analyses outside the seven units

- Fine-scale spatial genetic structure (Q2) and dispersal agents (Q3): `cluster/q2_sgs.R`,
  `cluster/check_depth_distance.R` (does depth track within-stand distance),
  `cluster/q3_dispersal.R`.

### Housekeeping

`cluster/status.sh` (queue and result status), `cluster/overnight_summary.sh` (collects overnight
outputs).

### Unclassified

- `cluster/abc/variant_counts.py`: counts of range-private vs shared variants for "main Figure 3B".
  The current Figure 3 is QST-FST, so this refers to an older figure numbering. The output is
  probably `figdata/variant_counts.txt`, which is not traced in the log.
- `cluster/run_pop_test.sh`, `cluster/run_f_test.sh`: headers give only environment setup and a
  VCF path.
- `cluster/ref/dl_kraken.sh`: downloads a Kraken2 database; the header states no purpose.
- `cluster/ref/.jobids`: list of Slurm job IDs, not a script.
