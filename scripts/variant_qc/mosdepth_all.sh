#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 8gb -N 1 -n 1 -t 2:00:00
#SBATCH --array=1-8
#SBATCH -J mosdepth_qc --out logs/mosdepth_qc.%A_%a.log
# Read depth per strain from the current CRAMs (flat_cram_dir), for the strain-QC
# low_coverage rule. -Q 20 matches the HaplotypeCaller mapping-quality read filter.
# Each array task takes every 8th CRAM (about 40 per task).
# Output: results/variant_qc/coverage/<strain>.mosdepth.{summary,global.dist}.txt
set -euo pipefail
module load mosdepth/0.3.12
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
REF=${PROJ}/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=${PROJ}/results/variant_qc/coverage
N=${SLURM_ARRAY_TASK_COUNT:-8}
I=${SLURM_ARRAY_TASK_ID:?}
mkdir -p "$OUT"
ls "${PROJ}"/flat_cram_dir/*.cram | awk -v n="$N" -v i="$I" '(NR - 1) % n == i - 1' | while read -r cram; do
    s=$(basename "$cram" .cram)
    [[ -s "${OUT}/${s}.mosdepth.summary.txt" ]] && continue
    mosdepth -n --fast-mode -Q 20 -t "${SLURM_CPUS_ON_NODE:-4}" -f "$REF" "${OUT}/${s}" "$cram"
done
