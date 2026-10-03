#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -J snpeff_qc
#SBATCH -o logs/snpeff_qc.%j.log
# One-off: re-annotate a filtered VCF with SnpEff, same command as the
# pipeline's SNPEFF_ANNOTATE (strips stale ANN/LOF/NMD from the input first).
# Usage: sbatch scripts/variant_qc/run_snpeff.sh results/filtered/<name>.vcf.gz
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
IN=$1
OUT=${IN%.vcf.gz}.annotated.vcf.gz
DB=$(readlink -f /bigdata/stajichlab/shared/lib/snpeff_db/RmucDH4148)
SNPEFF=/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f-c51b2e46bf63786b2d9a7a7d23680791163ab39a-0.img
BCF=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
cd "${SCRATCH:?}"
singularity exec "${BCF}" bcftools annotate -x INFO/ANN,INFO/LOF,INFO/NMD -Ov "${PROJ}/${IN}" \
  | singularity exec -B "${DB}" "${SNPEFF}" snpEff -Xmx12g -XX:-UsePerfData -stats "$(basename ${OUT%.vcf.gz}).snpEff_summary.html" \
        -c "${DB}/snpEff.config" -dataDir "${DB}/data" RmucDH4148 - \
  | singularity exec "${BCF}" bgzip -@ 2 > out.vcf.gz
first=$(zcat out.vcf.gz | head -1 || true)
[[ "${first}" == "##fileformat=VCF"* ]] || { echo "bad VCF header line: ${first}" >&2; exit 1; }
singularity exec "${BCF}" tabix -p vcf out.vcf.gz
mv out.vcf.gz "${PROJ}/${OUT}"; mv out.vcf.gz.tbi "${PROJ}/${OUT}.tbi"
cp ./*.snpEff_summary.html ./*.snpEff_summary.genes.txt "${PROJ}/results/filtered/" 2>/dev/null || true
echo "wrote ${OUT}"
