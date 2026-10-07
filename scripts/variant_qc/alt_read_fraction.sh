#!/usr/bin/bash -l
# Main-allele read fraction at each strain's ALT calls (DP >= 10), genome-wide; a pure haploid has nearly all >= 0.9.
# Used for docs/strain_identity_issues_2026-10-06.md (job 29534407). Edit the strain list below.
#SBATCH -p short -c 4 --mem 8gb -t 1:00:00 -J altfrac --out logs/altfrac.%j.log
module load bcftools
cd /bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
for s in EXF_7934 EXF_8891 EXF_9051 EXF_12768 EXF_1695 DBVPG_10619 EXF_8006; do
  bcftools view --threads 4 -f PASS -v snps -T ^results/mask/mask.bed --targets-overlap 1 -s ${s}_${s} -Ou results/RmucDH4148.all.hardfiltered.vcf.gz \
   | bcftools view -i 'GT="alt" & FMT/DP>=10' -Ou | bcftools query -f '[%AD]\n' \
   | awk -F, -v s=$s '{t=0; m=0; for(i=1;i<=NF;i++){t+=$i; if(i>1 && $i>m) m=$i} if(t>0){f=m/t; n++; b=int(f*10); if(b>9)b=9; h[b]++}} END{printf "%s n_alt=%d", s, n; for(i=0;i<10;i++) printf " %.1f-%.1f:%d", i/10,(i+1)/10,h[i]; print ""}'
done
