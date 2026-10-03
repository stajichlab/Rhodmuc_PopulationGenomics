#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -J nf-genotype
#SBATCH -o logs/nf-genotype.%j.log

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow/26.04.6
module load singularity

export NXF_SYNTAX_PARSER=v1

PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
PIPE=/rhome/jstajich/projects/nf/nf_genotype_population

mkdir -p logs

nextflow run "${PIPE}/main.nf" \
    -profile hpcc \
    -w "${PROJ}/work_genotype_full" \
    --cram_dir "${PROJ}/flat_cram_dir" \
    --reference "${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna" \
    --reference_fai "${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai" \
    --reference_dict "${PROJ}/results/reference/dict/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.dict" \
    --metadata "${PROJ}/metadata.txt" \
    --ploidy_overrides "${PROJ}/ploidy_overrides.csv" \
    --population_sets "${PROJ}/population_sets.yaml" \
    --mask_bed "${PROJ}/results/mask/mask.bed" \
    --output_prefix RmucDH4148 \
    --snpeff_db_dir "/bigdata/stajichlab/shared/lib/snpeff_db/RmucDH4148" \
    --snpeff_genome_name "RmucDH4148" \
    "$@"
