#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 32G -t 2:00:00 -J rebuild_cram -o logs/rebuild_cram.%j.log
# Rebuild one strain's CRAM from a subset of its read groups (for example the
# original library only, after a mismatched library was merged in by Sarek).
# 1. samtools view -r: keep only the listed read groups (no re-alignment)
# 2. GATK MarkDuplicates on the subset: duplicate flags were set when the
#    other library was present under the same LB, so they are recomputed
# 3. write CRAM + index; report records and duplicates kept
# Writes results/preprocessing/rebuilt/<STRAIN>/<STRAIN>.rebuilt.cram
# Usage: sbatch scripts/rebuild_strain_cram.sh STRAIN in.cram RGID [RGID ...]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
SAM=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
GATK=/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-broadinstitute-gatk-4.5.0.0.img
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
S=$1; IN=$2; shift 2
OUT=results/preprocessing/rebuilt/$S
mkdir -p $OUT
W=${SCRATCH:?}/rb; mkdir -p $W
printf "%s\n" "$@" > $W/rg.txt
X="singularity exec -B $PWD,$SCRATCH"
$X $SAM bash -c "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH
  samtools view -@4 --reference $REF -R $W/rg.txt -b -o $W/sub.bam $IN
  samtools view -H $W/sub.bam | grep -v '^@RG' > $W/h.sam
  samtools view -H $W/sub.bam | grep '^@RG' | grep -F -f $W/rg.txt >> $W/h.sam
  samtools reheader $W/h.sam $W/sub.bam > $W/sub.rh.bam
  samtools index $W/sub.rh.bam"
$X $GATK gatk --java-options "-Xmx20g" MarkDuplicates -I $W/sub.rh.bam -O $W/md.bam -M $OUT/$S.rebuilt.md_metrics.txt \
    --CLEAR_DT true --TMP_DIR $W --VALIDATION_STRINGENCY SILENT
$X $SAM bash -c "export PATH=/opt/conda/envs/bcftools_samtools/bin:\$PATH
  samtools view -@4 -C --output-fmt-option version=3.0 -T $REF -o $OUT/$S.rebuilt.cram $W/md.bam
  samtools index $OUT/$S.rebuilt.cram
  echo 'read groups:'; samtools view -H $OUT/$S.rebuilt.cram | grep '^@RG' | cut -f1-5
  echo 'records in: '\$(samtools view -c -R $W/rg.txt --reference $REF $IN)'  out: '\$(samtools view -c --reference $REF $OUT/$S.rebuilt.cram)
  echo 'duplicates before (subset of old flags): '\$(samtools view -c -f 1024 -R $W/rg.txt --reference $REF $IN)'  after re-marking: '\$(samtools view -c -f 1024 --reference $REF $OUT/$S.rebuilt.cram)"
