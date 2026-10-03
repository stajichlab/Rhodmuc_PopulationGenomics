#!/usr/bin/env bash
#SBATCH -p short -c 2 --mem 8G -t 2:00:00 -J het_test -o logs/het_script_test.%j.log
# Reproduce custom_het_ploidy.py on one contig and test fixes.
# For each strain: the exact production command, then with -a AD, then with
# -a AD plus filters. Reports FORMAT keys, het / hom-alt counts, het fraction
# as the script computes it (het/(het+homalt)), and het per covered kb.
# Usage: sbatch het_script_test.sh CONTIG STRAIN [STRAIN ...]
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
CTG=$1; shift
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -s "$CTG" "$REF" "$@" <<'IN'
set -euo pipefail
export PATH=/opt/conda/envs/bcftools_samtools/bin:$PATH
CTG=$1; REF=$2; shift 2
summ() { # stdin VCF; label
  awk -v L="$1" -F'\t' '
    /^#/ {next}
    length($4)!=1 || length($5)!=1 {next}
    { split($9,k,":"); split($10,v,":"); for(i in k) f[k[i]]=v[i]; keys=$9
      g=f["GT"]; gsub(/\|/,"/",g); split(g,a,"/"); if (a[1]=="." ) next
      if (a[1]!=a[2]) { het++
        if ("AD" in f) { split(f["AD"],d,","); t=d[1]+d[2]; if (t>0) { ab=d[2]/t; if (ab>=0.35 && ab<=0.65) hetab++ } }
      } else homalt++
      delete f }
    END { printf "  %-34s FORMAT=%-22s het=%-6d het_AB0.35-0.65=%-6d homalt=%-6d script_frac=%.4f  het_AB/(het_AB+homalt)=%.4f\n", L, keys, het, hetab, homalt, (het+homalt? het/(het+homalt):0), (hetab+homalt? hetab/(hetab+homalt):0) }'
}
echo "bcftools: $(bcftools --version | head -1)"
for S in "$@"; do
  C=flat_cram_dir/$S.cram
  echo "== $S ($CTG)"
  bcftools mpileup -Ou -f $REF -r $CTG $C 2>/dev/null | bcftools call -mv --ploidy 2 2>/dev/null | summ "production (no -a AD)"
  bcftools mpileup -Ou -f $REF -r $CTG -a FORMAT/AD,FORMAT/DP $C 2>/dev/null | bcftools call -mv --ploidy 2 2>/dev/null | summ "with -a AD"
  bcftools mpileup -Ou -f $REF -r $CTG -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 $C 2>/dev/null | bcftools call -mv --ploidy 2 2>/dev/null \
    | bcftools view -i 'QUAL>=30 && FMT/DP>=10' 2>/dev/null | summ "-a AD, MQ20 BQ20, QUAL30 DP10"
  bcftools mpileup -Ou -f $REF -r $CTG -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 $C 2>/dev/null | bcftools call -mv --ploidy 2 2>/dev/null \
    | bcftools view -i 'QUAL>=30 && FMT/DP>=10' -T ^results/mask/mask.bed --targets-overlap 1 2>/dev/null | summ "  + repeat mask"
done
IN
