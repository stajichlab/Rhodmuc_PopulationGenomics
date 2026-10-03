#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH -J rmask_DH4148
#SBATCH -o logs/repeatmask_only.%j.log
# RepeatMasker step only, using the library repeatmask.sh already built:
# results/mask/DH4148_plus_fungi.lib = RepeatModeler 2 families for DH4148 +
# Rmuc_DH4148/lib/fungi_repeat.20170127.lib. -engine rmblast: the env's
# RepeatMasker defaults to HMMER, which cannot use a FASTA library.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load bedtools/2.30.0
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/RepeatModeler/2.0.6/env/bin:${PATH}
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
REF=${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=${PROJ}/results/mask
CPU=${SLURM_CPUS_ON_NODE:-16}
cd "${SCRATCH:?}"
RepeatMasker -engine rmblast -pa $((CPU / 4)) -lib "${OUT}/DH4148_plus_fungi.lib" -gff -xsmall -dir rm "${REF}"
cp rm/*.out rm/*.tbl rm/*.out.gff "${OUT}/"
# RepeatMasker .out: 3 header lines; cols 5=contig 6=start(1-based) 7=end
awk 'NR>3{OFS="\t"; print $5,$6-1,$7}' rm/*.out | sort -k1,1 -k2,2n | bedtools merge > "${OUT}/repeatmasker.bed"
cat "${OUT}/repeatmasker.bed" "${OUT}/dust.bed" "${OUT}/trf.bed" | sort -k1,1 -k2,2n | bedtools merge > "${OUT}/mask.bed"
for b in repeatmasker dust trf mask; do
  awk -v n=$b '{s+=$3-$2} END{printf "%s\t%d intervals\t%d bp\n", n, NR, s}' "${OUT}/${b}.bed"
done
cat rm/*.tbl
