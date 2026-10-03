#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 8G -t 2:00:00 -J vcf_verify -o logs/vcf_verify.%j.log
# Check the annotated outputs of run_filter.sh + run_snpeff.sh.
# Usage: sbatch scripts/variant_qc/verify_outputs.sh <prefix> [samples.txt]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
P=${1:?prefix}; S=${2:-}
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
singularity exec -B "$PWD" "$IMG" bash -s "$P" "$S" <<'IN'
set -euo pipefail
P=$1; S=$2
Q=results/filtered/$P.qc.annotated.vcf.gz
M=results/filtered/$P.snps.maf.annotated.vcf.gz
for f in $Q $M; do
  echo "== $f"
  echo "line1: $(bcftools view -h $f | head -1)"
  echo "samples: $(bcftools query -l $f | wc -l)"
  echo "records: $(bcftools index -n $f)"
  echo "no_ANN: $(bcftools view -H -e 'INFO/ANN="."' $f | wc -l)"
  echo "nonPASS: $(bcftools view -H -e 'FILTER="PASS"' $f | wc -l)"
  echo "F_MISSING>0.1: $(bcftools view -H -i 'INFO/F_MISSING>0.1' $f | wc -l)"
done
echo "snps multiallelic: $(bcftools view -H -m3 $M | wc -l)"
echo "snps star: $(bcftools view -H -i 'ALT="*"' $M | wc -l)"
echo "snps non-SNP: $(bcftools view -H -V snps $M | wc -l)"
echo "snps MAF<0.05: $(bcftools view -H -i 'INFO/MAF<0.05' $M | wc -l)"
echo "qc types:"; bcftools view -H $Q | awk '{n=split($5,a,","); t=(length($4)==1)?"snp":"other"; for(i=1;i<=n;i++) if(a[i]!="*" && length(a[i])!=1) t="indel/mnp"; print t, (n>1?"multi":"bi")}' | sort | uniq -c
echo "snps ANN impact (first ANN):"; bcftools query -f '%INFO/ANN\n' $M | cut -d'|' -f3 | sort | uniq -c
if [[ -n "$S" ]]; then
  diff <(bcftools query -l $M | sort) <(sort "$S") >/dev/null && echo "sample list matches $S" || echo "SAMPLE LIST DIFFERS from $S"
fi
IN
