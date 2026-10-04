#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 8G -t 2:00:00 -J merge_cram -o logs/merge_cram.%j.log
# Merge per-library markduplicates CRAMs that are the same strain into one CRAM.
# Read-group IDs and LB stay per library; SM is set to <STRAIN>_<STRAIN> (sarek form).
# Writes results/preprocessing/merged/<STRAIN>/<STRAIN>.merged.cram(.crai).
# Usage: sbatch scripts/merge_strain_crams.sh STRAIN in1.cram in2.cram [...]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
S=$1; shift
OUT=results/preprocessing/merged/$S
mkdir -p $OUT
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -s "$S" "$REF" "$OUT" "$@" <<'IN'
set -euo pipefail
export PATH=/opt/conda/envs/bcftools_samtools/bin:$PATH
S=$1; REF=$2; OUT=$3; shift 3
samtools merge --threads 4 --reference $REF -O CRAM,version=3.0 -o $SCRATCH/m.cram "$@"
samtools view -H $SCRATCH/m.cram | awk -v sm="SM:${S}_${S}" 'BEGIN{OFS="\t"} /^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/) $i=sm} {print}' > $SCRATCH/h.sam
samtools reheader $SCRATCH/h.sam $SCRATCH/m.cram > $OUT/$S.merged.cram
samtools index $OUT/$S.merged.cram
echo "samples in merged CRAM:"; samtools samples $OUT/$S.merged.cram
echo "read groups:"; samtools view -H $OUT/$S.merged.cram | grep "^@RG" | cut -f1-5
n=0; for f in "$@"; do c=$(samtools view -c --reference $REF $f); echo "input $f: $c records"; n=$((n+c)); done
m=$(samtools view -c --reference $REF $OUT/$S.merged.cram); echo "sum of inputs: $n; merged: $m"
[ "$n" -eq "$m" ] || { echo "RECORD COUNT MISMATCH" >&2; exit 1; }
echo "duplicate-flagged records (input sum vs merged):"
d=0; for f in "$@"; do d=$((d+$(samtools view -c -f 1024 --reference $REF $f))); done; echo "$d vs $(samtools view -c -f 1024 --reference $REF $OUT/$S.merged.cram)"
IN
