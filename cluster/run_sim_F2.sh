#!/bin/bash
#SBATCH -J simF2
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /home/ecastron/araucaria/logs/simF2_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/simF2_%j.err
# Higher dropout arms, to find where (or whether) independent per-allele dropout reaches the
# observed F = 0.665. Under this model a heterozygote is misread as a homozygote with probability
# 2p(1-p), which PEAKS AT 0.5 when p = 0.5 -- so if the ceiling is real, no rate reaches 0.665 and
# a second mechanism is required. The 0.50 arm sits exactly at the theoretical maximum.
set -uo pipefail
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
O=$P/results/Q1/simF; mkdir -p "$O"
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
for pair in "drop40 0.40" "drop50 0.50"; do
  set -- $pair; tag=$1; drop=$2
  echo "=================== ARM: $tag (dropout=$drop) ==================="
  python3 $P/scripts/simulate_HWE_F.py "$V" "$O/sim_${tag}.beagle" --dropout "$drop" || continue
  gzip -f "$O/sim_${tag}.beagle"
  NS=$(zcat "$O/sim_${tag}.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
  NL=$(( $(zcat "$O/sim_${tag}.beagle.gz" | wc -l) - 1 ))
  zcat "$O/sim_${tag}.beagle.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$O/samples_${tag}.txt"
  ngsRelate -G "$O/sim_${tag}.beagle.gz" -n "$NS" -L "$NL" -z "$O/samples_${tag}.txt" \
            -O "$O/ngsrelate_${tag}.tsv" 2>&1 | tail -1
  python3 - "$O/ngsrelate_${tag}.tsv" "$tag" <<'PY'
import sys, statistics as st
f,tag=sys.argv[1],sys.argv[2]; hdr=None; F={}
for line in open(f):
    p=line.rstrip("\n").split("\t")
    if hdr is None: hdr=p; continue
    d=dict(zip(hdr,p))
    try:
        if int(d["nSites"])<2000: continue
        for k,v in (("ida","Fa"),("idb","Fb")):
            x=float(d[v])
            if x==x: F.setdefault(d[k],[]).append(x)
    except (KeyError,ValueError): continue
m=sorted(st.median(v) for v in F.values())
print(f"  >>> {tag}: {len(m)} individuals | median F = {st.median(m):.4f} | mean {st.mean(m):.4f}")
PY
done
echo "OBSERVED: median F = 0.665"
