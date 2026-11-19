#!/bin/bash
#SBATCH -J ctcmp
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /home/ecastron/araucaria/logs/ctcmp_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/ctcmp_%j.err
# Computes the full oversplitting-vs-paralogy comparison across ct = 0.85 / 0.88 / 0.91.
#
#   F falls + Ts/Tv stable + shared-het filter stable  -> oversplitting fixed, keep the lower ct
#   F falls + Ts/Tv toward 0.5 + more loci filtered    -> bought with paralogs, reject
#
# 2026-08-24: the 8 Aug run (job 39645) exited 0 and produced a table of dashes. Cause: the
# filter was called as `filter_vcf_callrate.py <vcf> <outpath> 0.30`, but its interface is
# `<vcf> <cut>...` -- it derives output names itself. The output path was parsed as a float,
# raising ValueError, so no .cr30.vcf existed and every downstream step starved. Every failing
# step was wrapped in 2>/dev/null, so the job looked clean. Fixed below: correct arguments,
# all three call-rate cuts in the single streaming pass, and every step now reports failure
# loudly into the table instead of silently emitting "-".
set -uo pipefail
P=/home/ecastron/araucaria
R=$P/results/Q1/ct_compare; mkdir -p "$R"
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
DBG=$R/ct_compare_debug.log
: > "$DBG"
exec > "$R/ct_compare.txt" 2>&1

echo "==================== clustering-threshold comparison ===================="
echo "generated $(date)"
echo "verbose step output: $DBG"
printf "\n%-8s %10s %12s %10s %10s\n" ct ipyradH SNPs_cr30 TsTv medianF

for CT in 085 088 091; do
  D=$P/rad/ipyrad_full
  N=full_ct${CT}
  CL=$D/${N}_clust_0.${CT#0}/s4_joint_estimate.txt
  OUTD=$D/${N}_outfiles
  H=$(awk 'NR>1&&NF>=3{s+=$2;n++}END{if(n)printf "%.5f",s/n}' "$CL" 2>/dev/null)
  VCF=$OUTD/${N}.vcf
  [ -f "$VCF" ] || VCF=$OUTD/${N}.cr30.vcf
  [ -f "$VCF" ] || { printf "%-8s %10s  !! no VCF -- assembly incomplete\n" "0.${CT#0}" "${H:--}"; continue; }

  # Call-rate filters. One streaming pass writes .cr30/.cr50/.cr70 together; cr30 matches the
  # working set used throughout, cr50/cr70 are for methods needing well-genotyped SNPs.
  CR=$OUTD/${N}.cr30.vcf
  if [ ! -f "$CR" ]; then
    echo "### filter ct 0.${CT#0} $(date)" >> "$DBG"
    if ! python3 $P/scripts/filter_vcf_callrate.py "$VCF" 0.30 0.50 0.70 >> "$DBG" 2>&1; then
      printf "%-8s %10s  !! FILTER FAILED -- see %s\n" "0.${CT#0}" "${H:--}" "$DBG"; continue
    fi
  fi
  [ -f "$CR" ] || { printf "%-8s %10s  !! filter wrote no cr30\n" "0.${CT#0}" "${H:--}"; continue; }

  LOCI=$(grep -vc '^#' "$CR")
  TSTV=$(python3 - "$CR" <<'PY'
import sys
ts=tv=0; TS={("A","G"),("G","A"),("C","T"),("T","C")}
for L in open(sys.argv[1]):
    if L[0]=="#": continue
    f=L.split("\t",6)
    r,a=f[3],f[4].split(",")[0]
    if len(r)!=1 or len(a)!=1: continue
    (ts:=ts+1) if (r,a) in TS else (tv:=tv+1)
print(f"{ts/tv:.3f}" if tv else "NA")
PY
)
  # F via the same GL -> ngsRelate path used for the 0.91 working set
  B=$R/gl_${CT}.beagle
  MF="-"
  echo "### beagle ct 0.${CT#0} $(date)" >> "$DBG"
  if ! python3 $P/scripts/catg_to_beagle.py "$CR" "$B" 0.01 >> "$DBG" 2>&1; then
    MF="!!BEAGLE_FAILED"
  else
    gzip -f "$B"
    NS=$(zcat "$B.gz" | head -1 | awk '{print (NF-3)/3}')
    NL=$(( $(zcat "$B.gz" | wc -l) - 1 ))
    zcat "$B.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$R/samples_${CT}.txt"
    echo "### ngsRelate ct 0.${CT#0} n=$NS L=$NL $(date)" >> "$DBG"
    if ! ngsRelate -G "$B.gz" -n "$NS" -L "$NL" -z "$R/samples_${CT}.txt" \
         -O "$R/ngsrelate_${CT}.tsv" >> "$DBG" 2>&1; then
      MF="!!NGSRELATE_FAILED"
    else
      MF=$(python3 - "$R/ngsrelate_${CT}.tsv" <<'PY'
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
print(f"{st.median(m):.4f}" if m else "!!NO_F_PARSED")
PY
) || MF="!!F_PARSE_FAILED"
    fi
  fi
  printf "%-8s %10s %12s %10s %10s\n" "0.${CT#0}" "${H:--}" "${LOCI:--}" "${TSTV:--}" "$MF"
done

echo
echo "---------------- call-rate yield at each cut (for conStruct / feems) ----------------"
printf "%-8s %14s %14s %14s\n" ct cr30 cr50 cr70
for CT in 085 088 091; do
  OUTD=$P/rad/ipyrad_full/full_ct${CT}_outfiles
  row=$(printf "%-8s" "0.${CT#0}")
  for c in 30 50 70; do
    f=$OUTD/full_ct${CT}.cr${c}.vcf
    n=$([ -f "$f" ] && grep -vc '^#' "$f" || echo "-")
    row="$row$(printf "%14s" "$n")"
  done
  echo "$row"
done

echo
echo "---------------- ipyrad's own paralog filters, per assembly ----------------"
for CT in 085 088 091; do
  S=$P/rad/ipyrad_full/full_ct${CT}_outfiles/full_ct${CT}_stats.txt
  echo "--- ct 0.${CT#0} ---"
  [ -f "$S" ] && sed -n '/filtered_by/,/total_filtered/p' "$S" | head -12 || echo "  (stats not written)"
done

echo
echo "READING THIS:"
echo "  F down, Ts/Tv steady (1.5-2.5), filter rates steady -> oversplitting was real; adopt lower ct"
echo "  F down, Ts/Tv falling toward 0.5, more loci filtered -> paralogs collapsing; reject"
echo "  F unchanged                                          -> oversplitting is NOT the residual cause"
echo "  For reference, observed at ct 0.91: median F = 0.6652"
echo "  Simulation ceilings: Wahlund 0.208 | dropout 0.500 | both maxed 0.609"
echo
echo "  Any '!!' above is a real failure, not a missing value. Check $DBG."
