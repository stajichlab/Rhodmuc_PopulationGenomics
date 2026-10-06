#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 32gb -N 1 -n 1 -t 2:00:00
#SBATCH -J nearref_diff --out logs/nearref_diff.%j.log
# Differences among strains nearly identical to the reference DH4148.
# Near-ref set: strains with < MAX_ALT non-ref calls in genome.samples.tsv
# (diagnose_calls.py on the hard-filtered `all` VCF).
# Usage: sbatch scripts/variant_qc/nearref_differences.sh [MAX_ALT]
set -euo pipefail
module load bcftools
PROJ=/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
cd "$PROJ"
MAX_ALT=${1:-5000}
VCF=results/RmucDH4148.all.hardfiltered.vcf.gz
MASK=results/mask/mask.bed
OUT=results/variant_qc/nearref
mkdir -p "$OUT"
T=${SCRATCH:?}/nearref.$$
mkdir -p "$T"

awk -F'\t' -v m="$MAX_ALT" 'NR > 1 && $5 < m {print $1}' results/variant_qc/genome.samples.tsv > "$T/nr_samples.txt"
awk '{h = int(length($1) / 2); s = $1; if (substr(s, h + 1, 1) == "_" && substr(s, 1, h) == substr(s, h + 2)) s = substr(s, 1, h); print s}' \
    "$T/nr_samples.txt" > "${OUT}/nearref_strains.txt"
echo "near-ref strains: $(wc -l < "$T/nr_samples.txt")"

# Candidate sites: PASS, outside the mask, a non-ref GT in at least one near-ref strain.
# Two passes: all variant types (nearref.*) and biallelic-or-multiallelic SNPs only (nearref_snps.*).
for TYPE in all snps; do
    VT=(); SUF=""
    [[ $TYPE == snps ]] && VT=(-v snps) && SUF="_snps"
    bcftools view --threads 4 -f PASS "${VT[@]}" -T "^${MASK}" --targets-overlap 1 -S "$T/nr_samples.txt" -Ou "$VCF" \
      | bcftools view -i 'GT="alt"' -Ou \
      | bcftools query -f '%CHROM\t%POS\n' > "$T/sites_${TYPE}.tsv"
    echo "${TYPE}: candidate sites $(wc -l < "$T/sites_${TYPE}.tsv")"
    { bcftools query -l "$VCF" | paste -sd '\t'
      bcftools view -T "$T/sites_${TYPE}.tsv" -f PASS "${VT[@]}" -Ou "$VCF" \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP:%GQ:%AD]\n'
    } | /usr/bin/python3.12 scripts/variant_qc/nearref_differences.py "${OUT}/nearref_strains.txt" metadata.txt \
          results/variant_qc/divergence/clone_groups_1e-3.tsv "${OUT}/nearref${SUF}"
done
