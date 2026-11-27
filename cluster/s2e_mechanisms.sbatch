#!/bin/bash
# Figure S2E: compute every S2E row against a baseline from the same genotype route, and the
# S2-related values figdata/locus_filter_F.tsv, the ipyrad binomial series and the copy-class spreads.
# Steps run side by side; s2e_assemble.py then writes results/derived/s2e_mechanisms.tsv.
#SBATCH -J s2e
#SBATCH -p basic
#SBATCH -c 4
#SBATCH --mem 96G
#SBATCH -t 06:00:00
#SBATCH -o /home/ecastron/araucaria/logs/s2e_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/s2e_%j.err
set -euo pipefail
P=/home/ecastron/araucaria
S=$P/scripts
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
POP=$P/metadata/popmap_177.tsv
O=$P/results/derived/s2e
mkdir -p $O
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
# copy-class spreads
python $P/scripts/ref/F_by_copynumber.py $V $POP $P/results/sv/loci_copies.tsv > $O/F_by_copynumber.txt &
python $S/s2e_binomial_route.py $V $POP $P/results/sv/loci_copies.tsv $P/results/sv/loci_sd.tsv \
       $O/binomial_route.tsv > $O/binomial_route.txt &
(
  module load R/R-4.4.0
  export R_LIBS_USER=/home/ecastron/R/libs
  export LD_LIBRARY_PATH=/usr/lib64:/lib64:${LD_LIBRARY_PATH:-}
  Rscript --vanilla $S/s2e_polyrad_route.R > $O/polyrad_route.txt 2>&1
) &
wait
python $S/s2e_assemble.py > $O/s2e_assemble.txt
echo done
