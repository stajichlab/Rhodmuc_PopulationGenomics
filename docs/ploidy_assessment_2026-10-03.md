# Ploidy assessment and reclassification, 2026-10-03

## Summary

- Four methods were compared: two read-based ploidy tools and two checks on the joint-called VCF.
- No strain called haploid shows a diploid heterozygosity signal.
- Three strains called diploid show no heterozygosity: DBVPG_3239, DBVPG_6094 and TFCN_25-332D-2.
- These three strains were changed from diploid to haploid in `ploidy_overrides.csv`.
- The GATK population run was started again. All outputs made from the earlier call set are now out of date (see "Invalidated outputs").
- Limit: SNP data cannot separate a haploid from a fully homozygous diploid. The three changes are consistent with the data, but the data do not prove them.

## Methods

### 1. custom_het_script (pipeline module CUSTOM_HET_PLOIDY)

- **Script:** `nf_genotype_population/bin/custom_het_ploidy.py`.
- **Input:** the CRAM for each strain.
- **Variant calls:** `bcftools mpileup -Ou -f REF CRAM | bcftools call -mv --ploidy 2`. This forces diploid calls over the whole genome (no `--region`).
- **Het fraction:** heterozygous biallelic SNPs / all biallelic SNPs. A het is meant to count only if its allele balance is from 0.35 to 0.65.
- **Call:** diploid if the het fraction is above 0.01 (`--threshold`).
- **Output:** `results/ploidy_inference_all.csv` (strain, inferred_ploidy, het_fraction).
- **Defect (tested; see "Method 1 defect" below):** the denominator counts only variant sites (het + hom-alt). In a strain close to the reference this count is very small, so a few noise hets in repeats give a large fraction. An earlier version of this report said the allele-balance check was skipped because AD was missing. That was wrong: bcftools 1.24 mpileup writes FORMAT/AD by default, and the check does apply.
- **Agreement:** this method agreed with the metadata for only 136 of 319 strains (`results/ploidy_review.csv`).

### 2. nQuire (pipeline module NQUIRE)

- **Script:** `nf_genotype_population/bin/nquire_ploidy.py`.
- **Container:** nquire-a990a88.sif.
- **Steps:** convert CRAM to BAM, then `nQuire create -b BAM`, `nQuire denoise`, `nQuire lrdmodel`.
- **Call:** taken from the best-fit model only. Diploid → `diploid`. Triploid or tetraploid → `non_diploid`. nQuire has no haploid model.
- **Run:** on all 319 strains (commit 86f43b2). Logs: `logs/nquire_*.log`.
- **Output:** `results/ploidy_review_nquire.csv`. 35 strains were called diploid and 284 non_diploid.

### 3. Genotype QC on the joint VCF

- **Scripts:** `scripts/variant_qc/diagnose.sh`, `diagnose_calls.py` and `strain_qc_table.py`.
- **Input:** `results/all.annotated.vcf.gz`, which is the earlier call set.
- **Per-strain measures:** het rate, allele balance of het calls, and the share of haploid ALT calls with mixed reads.
- **Thresholds:**
  - `ploidy_mixed`: a haploid strain with at least 200 ALT calls, where more than 25% of those calls have mixed reads (the called allele has less than 80% of reads).
  - `ploidy_skewed`: a diploid strain with het rate above 5%, where more than 40% of het calls have allele balance outside 0.2–0.8.
  - `no_het_diploid` (flag only): a diploid strain with het rate below 1% and no low-coverage drop.
- **Output:** `results/variant_qc/strain_qc.tsv`, `excluded_strains.tsv`, `genome.samples.tsv`.

### 4. Read allele balance at SNPs

- **Scripts:** `scripts/variant_qc/ploidy_evidence.sh` and `allele_balance.py` (SLURM job 29354684).
- **Sites:** sites from `results/all.annotated.vcf.gz` that are FILTER=PASS, biallelic SNPs, and outside `results/mask/mask.bed`. Only the 5 largest contigs were used: CM179485.1, CM179486.1, JBZGVR010000032.1, CM179489.1 and CM179487.1. That gives 1,082,195 sites.
- **Per-strain counts:** sites with read depth (ref+alt) of at least 10 are "covered". The minor-allele fraction f = min(ref,alt)/total, counted only when the minor allele has at least 2 reads. "Balanced" sites have 0.25 ≤ f ≤ 0.5. The histogram of f uses 0.05 bins from 0.10 to 0.50.
- **Output:** `results/variant_qc/allele_balance.tsv`.
- **Combined table:** `results/variant_qc/ploidy_evidence.tsv` joins this with methods 1–3 and the metadata.

### 5. Hybrid origin of the high-het diploids (related check)

- **Script:** `scripts/variant_qc/species_check.sh` + `hybrid_check.py` (SLURM job 29354324).
- **Sites:** CM179498.1, FILTER=PASS, biallelic SNPs. Genotypes with GQ < 20 or DP < 5 were set to missing. That gives 134,131 SNPs.
- **Measure:** for each het site, check whether the ALT allele is present in the 232 kept haploids, and whether it is the major allele in 10 R. frigidialcoholis haploids or 7 R. aff. mucilaginosa haploids.
- **Output:** `results/variant_qc/hybrid_check.CM179498.tsv` and `hybrid_diploids.tsv`.

## Results

Balanced sites per 1000 covered sites (method 4):

| Group | n | Median | Max |
|---|---|---|---|
| Called haploid, kept | 232 | 0.02 | 0.85 |
| Called haploid, dropped | 33 | 0.47 | 1.68 |
| Called diploid, het rate ≥ 0.1 | 39 | 196.51 | 225.75 (min 15.39) |
| Called diploid, het rate < 0.1 | 15 | 0.15 | 26.75 |

- No called-haploid strain has a diploid-level balanced-site rate.
- nQuire calls no called-haploid strain diploid.
- EXF_1695 has the highest mixed-read rate of the kept haploids: 22.51 per 1000 covered sites. Its f distribution peaks at 0.10, not 0.5. This fits a minor mixed culture better than diploidy. This is an inference and was not tested.

### Evidence for the three reclassified strains

| Strain | GATK het rate | Balanced / 1000 | nQuire | custom_het (fraction) | Metadata |
|---|---|---|---|---|---|
| DBVPG_3239 | 0.001 | 0.14 | non_diploid | diploid (0.0132) | too low |
| DBVPG_6094 | 0.001 | 0.17 | non_diploid | diploid (0.0155) | too low |
| TFCN_25-332D-2 | 0.000 | 0.02 | non_diploid | diploid (0.2003) | too low |

- In all three, the GATK het rate and the balanced-site rate are in the haploid range.
- nQuire did not fit the diploid model for any of them.
- The only diploid call came from custom_het_script, which has the open issue described in method 1.
- All three had been listed in `results/ploidy_overrides_NEEDS_REVIEW.csv`. In `ploidy_overrides_DRAFT.csv` their draft value was haploid.

## Changes made

1. **`ploidy_overrides.csv`:** DBVPG_3239, DBVPG_6094 and TFCN_25-332D-2 changed from diploid to haploid. The counts are now 51 diploid and 268 haploid. The file before the change is `ploidy_overrides.before_2026-10-03.csv`. The file uses CRLF line endings, and those were kept.
2. **`population_sets.yaml`:** regenerated with `scripts/variant_qc/make_population_sets.py`, with a header note about this change. Group membership did not change: rmuc_core has 235 strains, rmuc_core_outgroup 252, rmuc_with_hybrids 267 and hybrid_diploids 32. The three strains are still in rmuc_core.
3. **GATK population run:** started again with `sbatch -t 3-00:00:00 run_genotype_full.sh -resume` (job 29359109). Pipeline: nf_genotype_population, branch `feature/strain-tree` (PR #4).
   - HaplotypeCaller runs again only for the three strains. The GVCFs of the other 316 strains come from cache.
   - GenomicsDBImport, GenotypeGVCFs, the hard filter, VARIANT_QC_FILTER, SnpEff, SNP_ALIGNMENT and IQTREE run again for `all` and the four YAML groups.
   - Settings: `--mask_bed results/mask/mask.bed`, `--output_prefix RmucDH4148`, `--population_mode subset` (default) and `--tree_model GTR+ASC` (default; the pipeline now refuses a model without `+ASC`).
   - The earlier run (job 29356679) was cancelled before its filter jobs started.

## Invalidated outputs

These files come from the call set made before the reclassification. Do not use them for analysis.

- `results/all.annotated.vcf.gz` (+ `.tbi`)
- `results/filtered/rmucilaginosa_qc.*`
- `results/pca/*`
- `results/strain_tree/rmuc_kept267.*`. This IQ-TREE job (29354404) failed, because the varsites rerun lacked `-st DNA`.
- `results/gvcfs/` entries for the three strains, until the rerun replaces them.

The tables in `results/variant_qc/` (`strain_qc.tsv`, `allele_balance.tsv`, `ploidy_evidence.tsv`, `hybrid_check.CM179498.tsv`, `hybrid_diploids.tsv`) are the evidence for this decision. They also come from the earlier call set. After the rerun, the per-strain QC should be recomputed from `results/RmucDH4148.all.*`.

New outputs will be:

- `results/RmucDH4148.<pop>.{qc,snps.maf}.annotated.vcf.gz`
- `results/RmucDH4148.<pop>.filter_stats.tsv`
- `results/strain_tree/RmucDH4148.<pop>.snps.maf.*`

`<pop>` is one of: all, rmuc_core, rmuc_core_outgroup, rmuc_with_hybrids, hybrid_diploids.

## Method 1 defect

Test: `scripts/variant_qc/het_script_test.sh` (SLURM job 29359874, log `logs/het_script_test.29359874.log`). It runs the production command on contig CM179498.1 (1,038,828 bp), then reruns with mapping and base quality ≥ 20, QUAL ≥ 30, DP ≥ 10, and the repeat mask.

The script computes: het_fraction = hets with allele balance 0.35–0.65 / all biallelic variant sites (het + hom-alt). It calls diploid if het_fraction > 0.01.

| Strain | Production call (genome) | AB-passing hets | het + hom-alt | Value on CM179498.1 | Hets after filters + mask |
|---|---|---|---|---|---|
| TFCN_25-332D-2 | diploid (0.200) | 3 | 126 | 0.024 | 0 |
| EXF_8984 | diploid (0.161) | 0 | 57 | 0.000 | 1 (outside AB window) |
| CCFEE_5036 | haploid (0.009) | 39 | 7,722 | 0.005 | 10 AB-passing of 6,698 sites |
| DBVPG_3857 (diploid) | diploid (0.701) | 34,407 | 47,112 | 0.730 | 32,625 AB-passing of 43,977 sites |

Cause:
- `bcftools call -v` writes only variant sites. So the denominator is mostly the number of hom-alt sites, which measures distance from the reference, not genome size.
- A strain close to the reference has very few hom-alt sites: 37 for TFCN_25-332D-2 and 24 for EXF_8984 on this 1 Mb contig.
- A handful of balanced noise hets then pass the 0.01 threshold.
- Without the mask and quality filters, those noise hets come mostly from repeats. On this contig, the mask removes all 82 filtered hets in TFCN_25-332D-2 and 15 of 16 in EXF_8984.

Population-wide effect (`results/variant_qc/ploidy_evidence.tsv` and `strain_qc.tsv`, called haploids without a low-coverage drop, n = 256):
- 115 called haploids differ from the reference at fewer than 0.05% of variant sites. Their median het_fraction is 0.172, and method 1 calls 114 of them diploid.
- For divergence from 0.5% to 5%, the median is 0.008 to 0.010.
- The Spearman correlation between divergence and het_fraction is −0.684.
- This explains most of the 162 DISAGREE rows in `results/ploidy_review.csv`.

Proposed fix:
- Count heterozygous sites per covered kb, not per variant site. This needs a denominator of callable bases (for example `samtools depth` or `mosdepth` at DP ≥ 10).
- Add mpileup `-q 20 -Q 20`, a QUAL ≥ 30 and DP ≥ 10 site filter, and the repeat mask.
- On CM179498.1, AB-passing hets per Mb after filters are: diploid DBVPG_3857 about 31,000; divergent haploid CCFEE_5036 about 10; near-reference haploids 0. The genome-wide threshold must still be set from all strains.

## Open items

- **custom_het_script:** fix the denominator and add filters (see "Method 1 defect").
- **Duplicate metadata:** `metadata.txt` has two rows for TFCN_25-332D-2.
- **Homozygous diploids:** a fully homozygous diploid cannot be detected from SNP data. Flow cytometry, or a read-depth method with an internal reference, would be needed.
