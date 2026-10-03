#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 32G -t 2:00:00 -J divergence -o logs/divergence.%j.log
# Pairwise strain distances from the all-strain callset: PASS, unmasked,
# biallelic SNPs on the 5 largest contigs, every 4th site; GQ<20 / DP<5 -> missing.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
VCF=${1:-results/all.annotated.vcf.gz}
REG=$(sort -t$'\t' -k2,2nr refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai | head -5 | cut -f1 | paste -sd,)
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -c "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH;
  bcftools view -r $REG -f PASS -m2 -M2 -v snps -T ^results/mask/mask.bed --targets-overlap 1 -Ou $VCF \
   | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GQ<20 | FMT/DP<5' > $SCRATCH/s.bcf
  ( bcftools query -l $SCRATCH/s.bcf | tr '\n' '\t'; echo; bcftools query -f '[%GT\t]\n' $SCRATCH/s.bcf | awk 'NR%4==0' ) > $SCRATCH/gt.tsv"
/usr/bin/python3.12 scripts/variant_qc/pairwise_distance.py results/variant_qc/divergence/all_5contigs < $SCRATCH/gt.tsv
