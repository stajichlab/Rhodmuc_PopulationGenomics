#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -J build_mask
#SBATCH -o logs/build_mask.%j.log
# Build repeat + low-complexity mask BEDs for the DH4148 reference.
#   softmask.bed  - lowercase runs already in the GenBank FASTA (provenance: submitter's masking)
#   dust.bed      - low-complexity (NCBI dustmasker, default level 20)
#   trf.bed       - tandem repeats (TRF 2 7 7 80 10 50 500)
#   mask.bed      - union of the three, merged
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load ncbi-blast/2.16.0+ trf/4.11 bedtools/2.30.0

PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
REF=${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=${PROJ}/results/mask
mkdir -p "${OUT}"
cd "${SCRATCH:?}"

# 1. soft-masked (lowercase) intervals, 0-based half-open BED
awk 'BEGIN{OFS="\t"; s=-1}
     /^>/{if(s>=0)print c,s,p; c=substr($1,2); p=0; s=-1; next}
     {n=length($0); for(i=1;i<=n;i++){ch=substr($0,i,1); lc=(ch ~ /[acgtn]/);
        if(lc && s<0) s=p; else if(!lc && s>=0){print c,s,p; s=-1}; p++}}
     END{if(s>=0)print c,s,p}' "${REF}" > softmask.bed

# 2. dustmasker low-complexity (interval output is 0-based inclusive)
dustmasker -in "${REF}" -outfmt interval | awk 'BEGIN{OFS="\t"}
     /^>/{c=substr($1,2); next} {split($0,a," - "); print c,a[1],a[2]+1}' > dust.bed

# 3. Tandem Repeats Finder; -ngs prints "@contig" lines and 1-based inclusive hits
trf "${REF}" 2 7 7 80 10 50 500 -h -ngs | awk 'BEGIN{OFS="\t"}
     /^@/{c=substr($1,2); next} {print c,$1-1,$2}' > trf.bed

for b in softmask dust trf; do sort -k1,1 -k2,2n ${b}.bed | bedtools merge > "${OUT}/${b}.bed"; done
cat "${OUT}"/softmask.bed "${OUT}"/dust.bed "${OUT}"/trf.bed | sort -k1,1 -k2,2n | bedtools merge > "${OUT}/mask.bed"

for b in softmask dust trf mask; do
  awk -v n=$b '{s+=$3-$2} END{printf "%s\t%d intervals\t%d bp\n", n, NR, s}' "${OUT}/${b}.bed"
done
