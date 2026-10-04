#!/usr/bin/env bash
#SBATCH -p short -c 8 --mem 16G -t 2:00:00 -J cram_v30 -o logs/cram_v30.%j.log
# Rewrite CRAMs as CRAM 3.0 in place (GATK 4.5 / htsjdk cannot read CRAM 3.1,
# which samtools 1.24 writes by default). Checks the record count and that the
# reheadered file is readable, then replaces the file and its index.
# Usage: sbatch scripts/cram_to_v30.sh file1.cram [file2.cram ...]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -s "$REF" "$@" <<'IN'
set -euo pipefail
export PATH=/opt/conda/envs/bcftools_samtools/bin:$PATH
REF=$1; shift
for f in "$@"; do
  t=$SCRATCH/$(basename $f)
  samtools view -@8 -C --output-fmt-option version=3.0 -T $REF -o $t $f
  a=$(samtools view -c -@8 --reference $REF $f); b=$(samtools view -c -@8 --reference $REF $t)
  v=$(head -c 6 $t | tail -c 2 | od -An -tu1 | tr -s ' ')
  echo "$f: records $a -> $b; CRAM major.minor bytes:$v"
  [ "$a" -eq "$b" ] || { echo "COUNT MISMATCH $f" >&2; exit 1; }
  cp $t $f.tmp && mv $f.tmp $f && samtools index $f
done
IN
