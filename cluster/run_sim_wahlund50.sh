#!/bin/bash
#SBATCH -J simW50
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o /home/ecastron/araucaria/logs/simW50_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/simW50_%j.err
set -uo pipefail
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
O=$P/results/Q1/simF
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
for pair in "wahlund_drop50 0.50"; do
  set -- $pair; tag=$1; drop=$2
  echo "=================== ARM: $tag (Wahlund, dropout=$drop) ==================="
  python3 $P/scripts/simulate_HWE_F.py "$V" "$O/sim_${tag}.beagle" --dropout "$drop" \
      --popmap $P/metadata/popmap_177.tsv || continue
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
