#!/usr/bin/env bash
#SBATCH -p short -c 8 --mem 16G -t 2:00:00 -J sp_check -o logs/species_check.%j.log
# 1) SNPs in the curated SNP set that stay polymorphic (MAF>=0.05) in haploids only
#    or in high-het diploids only. 2) Origin of the het allele in high-het diploids.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -s <<'IN'
set -euo pipefail
G=results/variant_qc/groups
M=results/filtered/rmucilaginosa_qc.snps.maf.annotated.vcf.gz
echo "snps.maf total: $(bcftools index -n $M)"
for g in kept_haploid kept_dip_highhet; do
  n=$(bcftools view --threads 8 -S $G/$g.txt -Ou $M | bcftools +fill-tags -Ou -- -t MAF | bcftools view -H -i 'INFO/MAF>=0.05' | wc -l)
  echo "MAF>=0.05 within $g only: $n"
done
# hybrid check on contig CM179498 from the unfiltered all-strain VCF
cat $G/kept_haploid.txt $G/kept_dip_highhet.txt $G/frig_haploid.txt $G/affmuc_haploid.txt > $SCRATCH/s.txt
bcftools view -r CM179498.1 -f PASS -m2 -M2 -v snps -S $SCRATCH/s.txt -Ou results/all.annotated.vcf.gz \
  | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GQ<20 | FMT/DP<5' > $SCRATCH/c.bcf
echo "CM179498 biallelic PASS SNPs: $(bcftools view -H $SCRATCH/c.bcf | wc -l)"
( bcftools query -l $SCRATCH/c.bcf | tr '\n' '\t'; echo; bcftools query -f '[%GT\t]\n' $SCRATCH/c.bcf ) \
  > $SCRATCH/gt.tsv
IN
/usr/bin/python3.12 scripts/variant_qc/hybrid_check.py results/variant_qc/groups < $SCRATCH/gt.tsv > results/variant_qc/hybrid_check.CM179498.tsv
cat results/variant_qc/hybrid_check.CM179498.tsv
