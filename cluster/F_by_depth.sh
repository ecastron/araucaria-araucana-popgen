#!/bin/bash
#SBATCH -J Fdepth
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -o /home/ecastron/araucaria/logs/Fdepth_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/Fdepth_%j.err
# Does F survive when only confidently-callable genotypes are used?
#
# The GL model assumes a heterozygote's reads are a fair binomial at p=0.5. At depth 6 a 6:0
# split gives ~60:1 for homozygote, so systematic allelic imbalance is indistinguishable from a
# real homozygote. Raising the per-genotype depth floor makes that confusion progressively
# impossible. Decision rule, fixed BEFORE the run:
#   F -> ~0 and flattens        : depth-dependent het undercalling; GLs not absorbing it
#   F stays ~0.665 at all depths: undercalling refuted; cause is elsewhere
#   F < 0 at high depth         : high-depth tail is collapsed repeats; use the banded arm
set -uo pipefail
P=/home/ecastron/araucaria
R=$P/results/Q1/F_by_depth; mkdir -p "$R"
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
DBG=$R/F_by_depth_debug.log; : > "$DBG"
exec > "$R/F_by_depth.txt" 2>&1

echo "=========== F as a function of per-genotype depth floor ==========="
echo "generated $(date)"
echo "input: $V  (ct 0.91 cr30 working set, 78802 SNPs, 177 libraries)"
echo "baseline for comparison: median F = 0.6652"
echo
printf "%-14s %12s %12s %10s %10s\n" arm informativeGT medianSites indsUsed medianF

run_arm () {
  NAME=$1; MIN=$2; MAX=$3
  B=$R/gl_${NAME}.beagle
  echo "### arm $NAME min=$MIN max=$MAX $(date)" >> "$DBG"
  if ! python3 $P/scripts/catg_to_beagle_mindp.py "$V" "$B" 0.01 "$MIN" "$MAX" >> "$DBG" 2>&1; then
    printf "%-14s %12s\n" "$NAME" "!!BEAGLE_FAILED"; return
  fi
  GT=$(grep "informative genotypes kept" "$DBG" | tail -1 | sed 's/.*sites, \([0-9]*\) informative.*/\1/')
  MS=$(grep "per-individual informative" "$DBG" | tail -1 | sed 's/.*median=\([0-9]*\) .*/\1/')
  gzip -f "$B"
  NS=$(zcat "$B.gz" | head -1 | awk '{print (NF-3)/3}')
  NL=$(( $(zcat "$B.gz" | wc -l) - 1 ))
  zcat "$B.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$R/samples_${NAME}.txt"
  if ! ngsRelate -G "$B.gz" -n "$NS" -L "$NL" -z "$R/samples_${NAME}.txt" \
       -O "$R/ngsrelate_${NAME}.tsv" >> "$DBG" 2>&1; then
    printf "%-14s %12s\n" "$NAME" "!!NGSRELATE_FAILED"; return
  fi
  read USED MF < <(python3 - "$R/ngsrelate_${NAME}.tsv" <<'PY'
import sys, statistics as st
hdr=None; F={}
for line in open(sys.argv[1]):
    p=line.rstrip("\n").split("\t")
    if hdr is None: hdr=p; continue
    d=dict(zip(hdr,p))
    if int(d["nSites"])<2000: continue
    for k,v in (("ida","Fa"),("idb","Fb")):
        x=float(d[v])
        if x==x: F.setdefault(d[k],[]).append(x)
m=[st.median(v) for v in F.values()]
print(len(m), f"{st.median(m):.4f}" if m else "!!NO_F")
PY
)
  printf "%-14s %12s %12s %10s %10s\n" "$NAME" "${GT:--}" "${MS:--}" "${USED:--}" "${MF:--}"
}

run_arm baseline     0   1000000000
run_arm dp10        10   1000000000
run_arm dp15        15   1000000000
run_arm dp20        20   1000000000
run_arm dp30        30   1000000000
run_arm dp50        50   1000000000
run_arm dp30_100    30          100

echo
echo "READING THIS:"
echo "  baseline must reproduce ~0.6652 -- if it does not, nothing below is trustworthy"
echo "  F falling toward 0 as the floor rises  -> het undercalling / GL miscalibration"
echo "  F flat near 0.665 across all arms      -> undercalling refuted"
echo "  dp30_100 vs dp30: if dp30 is much lower, the high-depth tail is collapsed repeats"
echo "  Any '!!' is a failure, not a missing value. See $DBG"
