#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00
#SBATCH -J rm_DH4148
#SBATCH -o logs/repeatmask.%j.log
# De novo repeat library (RepeatModeler 2, with LTR structural search) for the
# DH4148 reference, then RepeatMasker with that library plus the curated fungal
# library (Rmuc_DH4148/lib/fungi_repeat.20170127.lib). Writes:
#   results/mask/repeatmasker.bed  (all RepeatMasker hits: TEs, simple, low complexity)
#   results/mask/mask.bed          (union: repeatmasker + dust + trf)
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load bedtools/2.30.0
# The RepeatModeler modulefile runs `conda activate --stack`, which fails in a
# batch shell without the conda hook. Use the module's env directly; it also
# ships RepeatMasker, rmblastn, LTR_retriever and the other dependencies.
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/RepeatModeler/2.0.6/env/bin:${PATH}
FUNGI_LIB=/bigdata/stajichlab/shared/projects/Rhodotorula/Rmuc_DH4148/lib/fungi_repeat.20170127.lib
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
REF=${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=${PROJ}/results/mask
CPU=${SLURM_CPUS_ON_NODE:-16}
cd "${SCRATCH:?}"
BuildDatabase -name DH4148 "${REF}"
RepeatModeler -database DH4148 -threads "${CPU}" -LTRStruct
cp DH4148-families.fa DH4148-families.stk "${OUT}/"
cat DH4148-families.fa "${FUNGI_LIB}" > DH4148_plus_fungi.lib
cp DH4148_plus_fungi.lib "${OUT}/"
RepeatMasker -engine rmblast -pa $((CPU / 4)) -lib DH4148_plus_fungi.lib -gff -xsmall -dir rm "${REF}"
cp rm/*.out rm/*.tbl rm/*.out.gff "${OUT}/"
# RepeatMasker .out: 3 header lines; cols 5=contig 6=start(1-based) 7=end
awk 'NR>3{OFS="\t"; print $5,$6-1,$7}' rm/*.out | sort -k1,1 -k2,2n | bedtools merge > "${OUT}/repeatmasker.bed"
cat "${OUT}/repeatmasker.bed" "${OUT}/dust.bed" "${OUT}/trf.bed" | sort -k1,1 -k2,2n | bedtools merge > "${OUT}/mask.bed"
for b in repeatmasker dust trf mask; do
  awk -v n=$b '{s+=$3-$2} END{printf "%s\t%d intervals\t%d bp\n", n, NR, s}' "${OUT}/${b}.bed"
done
