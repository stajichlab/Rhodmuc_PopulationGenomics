# Clean R. mucilaginosa population set: divergence, clonality and design notes (2026-10-03)

## Data and methods

- **Source callset:** `results/all.annotated.vcf.gz`, made before the 2026-10-03 ploidy changes. Divergence and species placement do not depend on those changes.
- **Sites:** FILTER=PASS, biallelic SNPs, outside `results/mask/mask.bed`, on the 5 largest contigs. Every 4th site was kept, which gives 270,548 sites. Genotypes with GQ < 20 or DP < 5 were set to missing.
- **Pairwise distance:** `scripts/variant_qc/divergence.sh` and `pairwise_distance.py` (SLURM job 29381744). For strains i and j: d_ij = Σ [x_i(1−x_j) + x_j(1−x_i)] / (number of sites called in both). Here x is the ALT-allele fraction (haploid 0 or 1; diploid 0, 0.5 or 1).
- **Divergence from the reference:** (variant sites − balanced hets) / callable bp, from `results/ploidy_v2/ploidy/ploidy_inference_all.csv`.
- **Plots and tables:** `scripts/variant_qc/plot_divergence.py 0.001` writes:
  - `results/variant_qc/divergence/divergence.{pdf,png}`
  - `strain_placement.tsv`
  - `clone_groups_1e-3.tsv` (clusters by single linkage at d < 0.001)

## Results

Pairwise distance (panel B):
- Core haploid vs core haploid: 0–0.045.
- Core vs *R. frigidialcoholis*: about 0.24–0.26.
- Core vs *R. aff. mucilaginosa*: about 0.26–0.28.
- No pair falls between 0.045 and 0.24, so the species boundary is clear.

Divergence from the reference (panel A): the core strains form three modes, at about 15–50, 2,000–4,000 and 5,000–10,000 SNPs per callable Mb. The other species are at about 50,000–100,000.

### The four no-metadata candidates

| Strain | Ploidy | Median d to core | Nearest other group | Decision |
|---|---|---|---|---|
| EXF_17335 | haploid (418 het/Mb) | 0.249 | *R. frigidialcoholis* EXF_7386, d = 0.0009 | *R. frigidialcoholis*; exclude |
| EXF_10854 | haploid (229 het/Mb) | 0.251 | *R. frigidialcoholis* TFCN_17-332M-1, d = 0.137 | separate lineage, close to neither group; exclude |
| TFCN_17-333M-1 | diploid (34,930 het/Mb) | 0.146 | *R. aff. mucilaginosa*, d = 0.145 | hybrid-diploid pattern; exclude |
| TFCN_186CL-2 | diploid (35,435 het/Mb) | 0.149 | *R. aff. mucilaginosa*, d = 0.145 | hybrid-diploid pattern; exclude |

None of the four is in any group. They have no metadata, so they were already excluded.

### NRRL_Y-2510

- Median d to the core is 0.039, and its nearest core strain is DBVPG_5757 (d = 0.0003).
- It is haploid by GATK, nQuire and het/Mb (56).
- The phenotype table lists it as *R. mucilaginosa*.
- It was added to `rmuc_core` and related groups as a reviewed no-metadata exception (`INCLUDE_NO_METADATA` in `make_population_sets.py`).

### Clonality (panel D, `clone_groups_1e-3.tsv`)

Nearest-neighbour distances within the core haploids have a gap between about 10^-3.4 and 10^-2.3.

At d < 0.001, the 236 `rmuc_core` strains form **21 clone groups**:
- The largest have 113, 38, 23, 17 and 10 strains.
- 9 strains are singletons.
- The largest group, CG001, spans China (marsh_tidalflat, 63 strains), the Arctic (8) and Italy (5). Its largest within-group distance is 0.00033.

## Groups (population_sets.yaml)

| Group | Strains | Contents |
|---|---|---|
| rmuc_core | 236 | kept R. mucilaginosa, haploid, no hybrids; includes NRRL_Y-2510 |
| rmuc_core_outgroup | 253 | rmuc_core + 17 haploid R. frigidialcoholis / R. aff. mucilaginosa |
| rmuc_with_hybrids | 268 | rmuc_core + 32 hybrid diploids |
| hybrid_diploids | 32 | high-het diploids |
| rmuc_pheno | 141 | rmuc_core with a phenotype row whose species is R. mucilaginosa |

- **rmuc_pheno exclusions** (`results/variant_qc/rmuc_pheno_exclusions.tsv`):
  - 80 strains have no phenotype row.
  - 15 strains have a phenotype-table species other than *R. mucilaginosa*, although their genotype is *R. mucilaginosa*: *R. sphaerocarpa* 7, *diobovata* 2, *toruloides* 2, and one each of *glutinis*, *paludigena*, *graminis*, *nothofagi*.
- **rmuc_pheno structure:**
  - The 141 strains fall in 14 clone groups, with sizes 55, 35, 16, 7, 6, 6, 6, 3, 2 and five singletons.
  - The phenotype table has three plates. The two largest clone groups are on all three plates. CG003 is on plates 1 and 2 only.
- **Phenotype analysis:** this is planned in a separate project, using the clean VCFs from here. That project should recompute MAF and missingness on its exact strain set.

## Design notes for association studies

These notes come from the numbers above. They are advice, not tested results.

1. **Effective sample size is the number of lineages, not the number of strains.**
   - 141 phenotyped strains fall in 14 clone groups, and 2 groups hold 90 of them.
   - Most variants that differ between lineages will be in near-complete linkage with lineage membership.
   - A naive per-SNP test will report lineage markers across the whole genome.
2. **Use a mixed model with a kinship matrix**, for example GEMMA, FaST-LMM or pyseer LMM.
   - Build the kinship from an LD-pruned SNP set.
   - Report the genomic inflation (λ) and QQ plots.
   - Also try fitting clone group as a random effect.
   - Compare against a permutation test that shuffles phenotypes between clone groups, not between strains.
3. **Test the lineage effect first.** Does the phenotype differ among clone groups (ANOVA or a mixed model)? A strong lineage effect with few lineages cannot be mapped to single genes by GWAS.
4. **Within-lineage tests.** CG001 (55 phenotyped strains) and CG002 (35) each have little internal variation. Within-group differences, mostly rare variants and gene presence/absence, can be tested with burden or k-mer/unitig methods (pyseer). Expect low power; replicate phenotype measurements help.
5. **Covariates.**
   - Include the phenotyping plate (PLATENUM): clone groups are not spread evenly over plates.
   - Include origin and environment only if they are not the question; they track the clone groups (CG001 is mostly China marsh_tidalflat).
6. **Culture identity.** Before use, verify the 15 strains whose phenotype-table species disagrees with the genotype, for example by ITS-sequencing the phenotyped stock. Keep EXF_1695 (possible mixed culture) out of association tests.
7. **Variant sets.**
   - Haploid genotypes only. `rmuc_core` is all haploid after the 2026-10-03 reclassification.
   - Recompute MAF on the phenotyped set and require a minor-allele count of at least 5.
   - Missingness ≤ 10%.
   - Use the SnpEff HIGH/MODERATE annotations for gene-level burden tests.
   - Run SNP-based and k-mer/presence-absence tests side by side, because accessory genes are common in yeasts.
8. **Multiple testing.** Base the threshold on the number of unique genotype patterns, or use permutation, not the raw SNP count. With few lineages, many SNPs share one pattern.
9. **Power.** 83 strains in the phenotype table are not genotyped here (other species or not sequenced). Sequencing more strains from under-sampled lineages would add more power than more strains from CG001 or CG002.
