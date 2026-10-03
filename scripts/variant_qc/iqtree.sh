#!/usr/bin/env bash
#SBATCH -p epyc -c 16 --mem 32G -t 4-00:00:00 -J iqtree -o logs/iqtree.%j.log
# IQ-TREE on a SNP alignment from strain_tree.sh (same command as the pipeline IQTREE process).
# Usage: sbatch iqtree.sh <name>
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
NAME=${1:?name}
IMG=/bigdata/stajichlab/shared/singularity_cache/iqtree_3.1.2--h8471819_0.sif
cd results/strain_tree
W=${SCRATCH:?}/iq; mkdir -p $W
gzip -dc ${NAME}.mfa.gz > $W/${NAME}.mfa
run_iq() { singularity exec -B "$W" "$IMG" iqtree3 -s "$1" -st DNA --prefix "$W/${NAME}" -m GTR+ASC -B 1000 -alrt 1000 -T 16 --seed 12345 -redo; }
if ! run_iq "$W/${NAME}.mfa"; then
  [[ -s "$W/${NAME}.varsites.phy" ]] || exit 1
  run_iq "$W/${NAME}.varsites.phy"
fi
rm -f $W/${NAME}.mfa
gzip -f $W/${NAME}.varsites.phy 2>/dev/null || true
cp $W/${NAME}.* .
