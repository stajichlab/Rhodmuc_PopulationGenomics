#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 8G -t 2:00:00 -J snp_aln -o logs/snp_aln.%j.log
# SNP alignment for the strain tree, using the pipeline's bin script.
# Usage: sbatch strain_tree.sh <name> <snps.vcf.gz> [max_missing]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
NAME=${1:?name}; VCF=${2:?vcf}; MM=${3:-0}
PIPE=/rhome/jstajich/projects/nf/nf_genotype_population
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
mkdir -p results/strain_tree
cd results/strain_tree
singularity exec -B "$(realpath ../..),${PIPE},${SCRATCH:?},/bigdata/stajichlab/jstajich/projects/nf/nf_genotype_population" "$IMG" bash -c \
  "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH; ${PIPE}/bin/vcf_to_snp_alignment.sh -i ../../${VCF} -o ${NAME} -r RmucDH4148_ref -t 4 --max-missing ${MM}"
cat ${NAME}.alignment_stats.tsv
