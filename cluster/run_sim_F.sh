#!/bin/bash
#SBATCH -J simF
#SBATCH -p basic
#SBATCH --nodelist=n003
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /home/ecastron/araucaria/logs/simF_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/simF_%j.err
set -uo pipefail
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
O=$P/results/Q1/simF; mkdir -p "$O"
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH

run_arm () {
  local tag=$1 drop=$2
  echo "=================== ARM: $tag (dropout=$drop) ==================="
  python3 $P/scripts/simulate_HWE_F.py "$V" "$O/sim_${tag}.beagle" --dropout "$drop" || return 1
  gzip -f "$O/sim_${tag}.beagle"
  NS=$(zcat "$O/sim_${tag}.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
  NL=$(( $(zcat "$O/sim_${tag}.beagle.gz" | wc -l) - 1 ))
  zcat "$O/sim_${tag}.beagle.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$O/samples.txt"
  ngsRelate -G "$O/sim_${tag}.beagle.gz" -n "$NS" -L "$NL" -z "$O/samples.txt" \
            -O "$O/ngsrelate_${tag}.tsv" 2>&1 | tail -2
  python3 - "$O/ngsrelate_${tag}.tsv" "$tag" <<'PY'
import sys, statistics as st
f, tag = sys.argv[1], sys.argv[2]
hdr=None; F={}
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
m=[st.median(v) for v in F.values()]
m.sort()
if m:
    print(f"  >>> {tag}: {len(m)} individuals | median F = {st.median(m):.4f} "
          f"| mean {st.mean(m):.4f} | IQR {m[len(m)//4]:.4f}-{m[3*len(m)//4]:.4f}")
else:
    print(f"  >>> {tag}: no usable F values")
PY
}

echo "OBSERVED for comparison: median F = 0.665, mean 0.741"
run_arm hwe 0.0
run_arm drop10 0.10
run_arm drop25 0.25
