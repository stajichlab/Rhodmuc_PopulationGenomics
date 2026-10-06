#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 48gb -N 1 -n 1 -t 2:00:00
#SBATCH -J declone --out logs/declone.%j.log
# Near-identical strain groups in rmuc_core and one representative per group.
# Input: the rmuc_core .qc VCF (all sites after VARIANT_QC_FILTER, no MAF filter,
# so private SNPs count), SNPs only.
# Usage: sbatch scripts/variant_qc/declone.sh [CUTOFF=5]
set -euo pipefail
module load bcftools
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
cd "$PROJ"
CUTOFF=${1:-5}
VCF=results/RmucDH4148.rmuc_core.qc.annotated.vcf.gz
OUT=results/variant_qc/declone
mkdir -p "$OUT"
{ bcftools query -l "$VCF" | paste -sd '\t'
  bcftools view -v snps -Ou "$VCF" | bcftools query -f '[%GT\t]\n'
} | /usr/bin/python3.12 scripts/variant_qc/declone.py "$CUTOFF" results/variant_qc/coverage.tsv metadata.txt "${OUT}/rmuc_core"
