#!/bin/bash
# Collect every overnight result into one place, so the morning report is read off disk rather
# than reconstructed. Each block prints its own heading and tolerates missing files: a job that
# failed should show as absent, not crash the summary.
P=/home/ecastron/araucaria
echo "===== ARAUCARIA overnight summary ====="; date; echo

echo "--- job outcomes ---"
sacct -S $(date -d yesterday +%F) -u ecastron --format=JobID,JobName%18,State,Elapsed,ExitCode -P \
  | grep -Ev "\.batch|\.extern" | grep -Ei "arau|ipy|its|pop|phylo" || echo "(none)"

echo; echo "--- ipyrad clustering sweep (replicate concordance) ---"
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate ipyrad 2>/dev/null
python $P/scripts/analyse_ipyrad.py 2>/dev/null | head -14 || echo "(unavailable)"

echo; echo "--- full 177-library assembly ---"
grep -hE "^DONE|^FAILED|libraries entering" $P/logs/ipyfull_*.out 2>/dev/null | tail -5 || echo "(not finished)"

echo; echo "--- POPULATION STRUCTURE TEST (the decisive one) ---"
tail -40 $P/logs/poptest_*.out 2>/dev/null || echo "(not run yet)"

echo; echo "--- ITS R1-only production run ---"
grep -hE "^\[ITS-R1\]" $P/logs/its_r1_*.out 2>/dev/null | tail -10 || echo "(not finished)"
