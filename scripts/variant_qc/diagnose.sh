#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -J diag_calls
#SBATCH -o logs/diag_calls.%j.log
# Usage: sbatch scripts/variant_qc/diagnose.sh <region|all> <vcf> <outprefix>
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load bcftools
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
REGION=$1; VCF=$2; PREFIX=$3
RARG=(); [[ "$REGION" != all ]] && RARG=(-r "$REGION")
MASK=${PROJ}/results/mask/mask.bed
{ bcftools query -l "$VCF" | paste -sd '\t'
  bcftools query "${RARG[@]}" -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO/DP[\t%GT:%AD:%DP:%GQ]\n' "$VCF"
} | /usr/bin/python3.12 ${PROJ}/scripts/variant_qc/diagnose_calls.py "$PREFIX" $( [[ -s $MASK ]] && echo "$MASK" )
