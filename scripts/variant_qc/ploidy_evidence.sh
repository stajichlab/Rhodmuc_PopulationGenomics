#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 16G -t 2:00:00 -J ploidy_ev -o logs/ploidy_evidence.%j.log
# Read allele balance per strain at PASS, unmasked, biallelic SNPs on the 5 largest contigs.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
REG=$1
VCF=${2:-results/all.annotated.vcf.gz}
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -c "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH;
  bcftools view -r $REG -f PASS -m2 -M2 -v snps -T ^results/mask/mask.bed --targets-overlap 1 -Ou $VCF > $SCRATCH/s.bcf
  echo sites: \$(bcftools view -H $SCRATCH/s.bcf | wc -l)
  ( bcftools query -l $SCRATCH/s.bcf | tr '\n' '\t'; echo; bcftools query -f '[%AD\t]\n' $SCRATCH/s.bcf ) > $SCRATCH/ad.tsv"
/usr/bin/python3.12 scripts/variant_qc/allele_balance.py < $SCRATCH/ad.tsv > results/variant_qc/allele_balance.tsv
echo done
