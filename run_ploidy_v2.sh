#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -J nf-ploidy-v2
#SBATCH -o logs/nf-ploidy-v2.%j.log
# Ploidy inference with the per-callable-Mb het metric (custom_het_ploidy.py,
# 2026-10-03 fix) on all CRAMs. Separate launch/work dir from the genotyping run.
source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow/26.04.6
module load singularity
export NXF_SYNTAX_PARSER=v1

PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
PIPE=${PIPE:-/rhome/jstajich/projects/nf/nf_genotype_population_hetfix}
OUT=${PROJ}/results/ploidy_v2
mkdir -p "$OUT" "${PROJ}/logs"
cd "$OUT"

nextflow run "${PIPE}/main.nf" -entry PLOIDY_ONLY \
    -profile hpcc \
    -w "${OUT}/work" \
    --cram_dir "${PROJ}/flat_cram_dir" \
    --reference "${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna" \
    --reference_fai "${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai" \
    --metadata "${PROJ}/metadata.txt" \
    --mask_bed "${PROJ}/results/mask/mask.bed" \
    --outdir "$OUT" \
    "$@"
