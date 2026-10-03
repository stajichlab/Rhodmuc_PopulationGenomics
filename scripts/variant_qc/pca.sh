#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 24G -t 2:00:00 -J snp_pca -o logs/snp_pca.%j.log
# PCA of the curated SNP set: all kept strains, then without the high-het diploids.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
M=results/filtered/rmucilaginosa_qc.snps.maf.annotated.vcf.gz
D=results/pca/rmucilaginosa_qc.snps.maf.dosage.tsv.gz
G=results/variant_qc/groups
mkdir -p results/pca
singularity exec -B "$PWD" "$IMG" bash -c "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH; bcftools +dosage $M -- -t GT | gzip -c > $D"
/usr/bin/python3.12 scripts/variant_qc/pca_snps.py $D results/variant_qc/strain_qc.tsv metadata.txt $G/pca_groups.tsv results/pca/rmuc_kept267
/usr/bin/python3.12 scripts/variant_qc/pca_snps.py $D results/variant_qc/strain_qc.tsv metadata.txt $G/pca_groups.tsv results/pca/rmuc_no_hybrids $G/kept_dip_highhet.txt
