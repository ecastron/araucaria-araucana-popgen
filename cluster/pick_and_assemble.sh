#!/bin/bash
#SBATCH --job-name=arau-assemble
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/assemble_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/assemble_%j.err

# Choose M and n from the sweep by the rule the spec fixed in advance, then run the full
# 173-tree assembly.
#
# The choice is made mechanically from a declared criterion, not by inspection: primary is
# genotyping error measured on the 15 replicate library pairs (Mastretta-Yanes et al. 2015),
# with r80 loci recovered as the tie-break (Paris et al. 2017). Selecting by rule rather than
# by eye is the point - a human picking after seeing the curves reintroduces the discretion
# the pre-registration logic exists to remove.
#
# Gates below stop the job rather than proceed on anomalous input. If any fires, nothing is
# assembled and the reason is logged for review.
#
# The popmap is per LIBRARY (188 entries), not per tree, because the spec requires genotyping
# each library separately first so the replicate pairs can be compared. Merging to one sample
# per tree happens afterwards.

set -uo pipefail

PROJ=/home/ecastron/araucaria
RUNLOG=$PROJ/logs/run_log.md
module load Stacks/stacks-2.68

echo "=== sweep analysis ==="
python3 $PROJ/scripts/analyse_sweep.py > $PROJ/logs/sweep_analysis.txt 2>&1
cat $PROJ/logs/sweep_analysis.txt

SUM=$PROJ/logs/sweep_summary.tsv
if [[ ! -s "$SUM" ]]; then
  echo "GATE FAILED: no sweep summary produced. Nothing assembled." >&2
  exit 1
fi

# --- apply the declared rule, with gates -------------------------------------------------
read -r CHOSEN GATEMSG < <(python3 - "$SUM" <<'PY'
import csv, sys
rows=[r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')]
def f(r,k):
    v=r.get(k,'')
    return None if v in ('','None') else float(v)

usable=[r for r in rows if f(r,'err_mean') is not None]

# Gate 1: enough M values completed
if len(usable) < 6:
    print(f"NONE only {len(usable)}/8 M values produced usable results"); raise SystemExit

# Gate 2: enough replicate pairs contributing
minpairs=min(int(r['pairs_used']) for r in usable)
if minpairs < 10:
    print(f"NONE only {minpairs} replicate pairs usable at some M"); raise SystemExit

# Gate 3: error rates in a sane range
if any(f(r,'err_mean') > 0.5 for r in usable):
    print("NONE implausible error rate (>0.5) at some M"); raise SystemExit

# Gate 4: r80 non-trivial
if max(int(r['r80_snps']) for r in usable) < 500:
    print("NONE fewer than 500 r80 SNPs at every M"); raise SystemExit

best=min(usable, key=lambda r: f(r,'err_mean'))
se=f(best,'err_se') or 0.0
thr=f(best,'err_mean')+se
within=[r for r in usable if f(r,'err_mean') <= thr]
pick=max(within, key=lambda r: int(r['r80_snps']))
print(f"{pick['M']} err={pick['err_mean']} r80={pick['r80_snps']} (min_err_M={best['M']}, thr={thr:.5f})")
PY
)

if [[ "$CHOSEN" == "NONE" || -z "$CHOSEN" ]]; then
  echo "GATE FAILED: $GATEMSG" >&2
  {
    echo ""
    echo "## $(date -u +%Y-%m-%dT%H:%MZ) — assembly NOT run"
    echo "Gate failed: $GATEMSG"
  } >> "$RUNLOG"
  exit 1
fi

M=$CHOSEN
n=$M
m=3
echo "CHOSEN M=$M n=$n m=$m  $GATEMSG"

{
  echo ""
  echo "## $(date -u +%Y-%m-%dT%H:%MZ) — assembly parameters fixed"
  echo ""
  echo "- **M = $M, n = $n, m = $m**"
  echo "- Selected mechanically by the rule declared in the spec: minimum genotyping error on the"
  echo "  15 replicate library pairs, r80 as tie-break within 1 SE of the minimum."
  echo "- \`$GATEMSG\`"
  echo "- Full sweep table: \`logs/sweep_summary.tsv\`"
  echo "- Chosen before any population-genetic or community result was computed."
} >> "$RUNLOG"

# --- full assembly, all 177 libraries ----------------------------------------------------
OUT=$PROJ/rad/final
rm -rf "$OUT"; mkdir -p "$OUT"

echo "=== denovo_map.pl on 177 libraries, M=$M n=$n m=$m ==="
date
if ! denovo_map.pl \
      --samples $PROJ/rad/cleaned \
      --popmap $PROJ/metadata/popmap_177.tsv \
      -o "$OUT" \
      -M "$M" -n "$n" -m "$m" \
      --paired \
      -T "$SLURM_CPUS_PER_TASK" \
      -X "populations: -r 0.80 --vcf --genepop --structure --hwe --fstats" \
      > "$OUT/denovo_map.stdout" 2>&1; then
  echo "ASSEMBLY FAILED" >&2
  tail -30 "$OUT/denovo_map.stdout" >&2
  tail -30 "$OUT/denovo_map.log" >&2 2>/dev/null
  exit 1
fi
date

if [[ ! -s "$OUT/populations.snps.vcf" ]]; then
  echo "ASSEMBLY FAILED: no VCF" >&2
  exit 1
fi

SNPS=$(grep -vc '^#' "$OUT/populations.snps.vcf")
SAMP=$(grep -m1 '^#CHROM' "$OUT/populations.snps.vcf" | awk '{print NF-9}')
echo "ASSEMBLY OK  r80_snps=$SNPS libraries=$SAMP"
{
  echo "- Assembly complete: **$SNPS r80 SNPs across $SAMP libraries** (\`rad/final/\`)"
} >> "$RUNLOG"
