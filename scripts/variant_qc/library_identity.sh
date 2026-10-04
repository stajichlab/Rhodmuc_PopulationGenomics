#!/usr/bin/env bash
#SBATCH -p short -c 4 --mem 16G -t 2:00:00 -J lib_ident -o logs/library_identity.%j.log
# Per-read-group genotypes at population SNPs (CM179498.1) for strains with 9003
# libraries; compare each library to every strain in the joint callset.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
IMG=/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=results/variant_qc/library_identity
mkdir -p $OUT/lib
STRAINS=$(awk -F, 'NR>1 && $6~/\/9003_/{print $1}' samplesheet.csv | sort -u)
singularity exec -B "$PWD,${SCRATCH:?}" "$IMG" bash -s $STRAINS <<'IN'
set -euo pipefail
export PATH=/opt/conda/envs/bcftools_samtools/bin:$PATH
REF=refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna
OUT=results/variant_qc/library_identity
if [[ -s $OUT/pop.tsv ]]; then
  # reuse the population table; rebuild the sites file from its first 4 columns
  awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3","$4}' $OUT/pop.tsv | bgzip > $SCRATCH/sites.tsv.gz
  tabix -s1 -b2 -e2 $SCRATCH/sites.tsv.gz
else
bcftools view -r CM179498.1 -f PASS -m2 -M2 -v snps -T ^results/mask/mask.bed --targets-overlap 1 -Ou results/all.annotated.vcf.gz \
  | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GQ<20 | FMT/DP<5' > $SCRATCH/pop.bcf
bcftools index $SCRATCH/pop.bcf
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' $SCRATCH/pop.bcf | bgzip > $SCRATCH/sites.tsv.gz
tabix -s1 -b2 -e2 $SCRATCH/sites.tsv.gz
{ printf "CHROM\tPOS\tREF\tALT"; bcftools query -l $SCRATCH/pop.bcf | awk '{printf "\t%s",$1}'; echo
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $SCRATCH/pop.bcf | awk -F'\t' 'BEGIN{OFS="\t"}{for(i=5;i<=NF;i++){g=$i; gsub(/\|/,"/",g); n=split(g,a,"/"); if(a[1]=="."){$i="NA"} else {s=0; for(k=1;k<=n;k++) s+=(a[k]!="0"); $i=s/n}} print}'
} > $OUT/pop.tsv
fi
for S in "$@"; do
  C=flat_cram_dir/$S.cram
  if [[ -s $OUT/lib/$S.tsv && $OUT/lib/$S.tsv -nt $(readlink -f $C) ]]; then continue; fi
  samtools view -H $C | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^ID:/) id=substr($i,4); print id"\t"id}' > $SCRATCH/rg.$S.txt
  bcftools mpileup -f $REF -r CM179498.1 -T $SCRATCH/sites.tsv.gz -G $SCRATCH/rg.$S.txt -a FORMAT/AD -q 20 -Q 20 -d 10000 -Ob -o $SCRATCH/lib.$S.bcf $C 2>/dev/null
  { printf "CHROM\tPOS\tREF\tALT"; bcftools query -l $SCRATCH/lib.$S.bcf | awk '{printf "\t%s",$1}'; echo
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' $SCRATCH/lib.$S.bcf
  } > $OUT/lib/$S.tsv
done
IN
/usr/bin/python3.12 scripts/variant_qc/library_identity.py $OUT/pop.tsv $OUT/lib > $OUT/library_identity.tsv
cat $OUT/library_identity.tsv
