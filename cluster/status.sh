#!/bin/bash
# One-shot status for the ARAUCARIA Q1/Q5 runs. Safe to run any time.
P=/home/ecastron/araucaria
echo "===== ARAUCARIA status ====="; date; echo
echo "--- queue ---"; squeue -u ecastron -o "%.8i %.12j %.9T %.10M %R" | grep -E "arau|JOBID" || echo "(nothing queued)"
echo; echo "--- recent job outcomes ---"
sacct -S $(date -d "-12 days" +%F) --format=JobID,JobName%14,State,Elapsed,ExitCode -P 2>/dev/null \
  | grep -Ev "\.batch|\.extern" | grep -i arau | tail -12
echo; echo "===== Q1: population genomics ====="
[ -f $P/results/Q1/q1_per_stand_diversity.tsv ] && { echo "-- per-stand diversity --"; column -t $P/results/Q1/q1_per_stand_diversity.tsv; }
[ -f $P/results/Q1/q1_fst_vcftools.tsv ] && { echo; echo "-- FST (vcftools, weighted) --"; sort -k3 -g -r $P/results/Q1/q1_fst_vcftools.tsv | head -8; echo "   ... $(wc -l < $P/results/Q1/q1_fst_vcftools.tsv) stand pairs total"; }
[ -f $P/results/Q1/q1_admixture_cv.txt ] && { echo; echo "-- ADMIXTURE cross-validation (lowest CV error = best K) --"; cat $P/results/Q1/q1_admixture_cv.txt; }
[ -f $P/results/Q1/q1_amova.txt ] && { echo; echo "-- AMOVA --"; head -20 $P/results/Q1/q1_amova.txt; }
echo; echo "===== Q5: community analysis ====="
[ -f $P/results/Q5/q5_summary.txt ] && grep -E "^----------|libraries x|stand identity|climate PC|CONTROL|constrained:" $P/results/Q5/q5_summary.txt
echo; echo "===== files produced ====="
ls -la $P/results/Q1/ $P/results/Q5/ 2>/dev/null | grep -vE "^total|^d"
