#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
#SBATCH -J vcf_qc
#SBATCH -o logs/vcf_qc.%j.log
# One-off: apply the pipeline's VARIANT_QC_FILTER step
# (nf_genotype_population/bin/filter_population_vcf.sh, default thresholds) to
# the already-annotated results/all.annotated.vcf.gz.
# Usage: sbatch scripts/variant_qc/run_filter.sh <out_prefix_name> [samples.txt]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
PIPE=/rhome/jstajich/projects/nf/nf_genotype_population
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
NAME=$1
SAMPLES=${2:-}
OUT=${PROJ}/results/filtered
mkdir -p "${OUT}"
cd "${SCRATCH:?}"
SARG=(); [[ -n "${SAMPLES}" ]] && SARG=(-S "${SAMPLES}")
singularity exec -B "${PROJ},${PIPE},${SCRATCH}" "${IMG}" bash "${PIPE}/bin/filter_population_vcf.sh" \
    -i "${PROJ}/results/all.annotated.vcf.gz" \
    -o "${NAME}" \
    -m "${PROJ}/results/mask/mask.bed" \
    "${SARG[@]}" \
    -t "${SLURM_CPUS_ON_NODE:-4}"
cp -p "${NAME}".* "${OUT}/"
cat "${OUT}/${NAME}.filter_stats.tsv" "${OUT}/${NAME}.depth_threshold.txt"
