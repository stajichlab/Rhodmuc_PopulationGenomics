# Strains close to the reference DH4148: are they different strains? (2026-10-05)

**Question (user, 2026-10-05):** TFCN_137D-4, TFCN_17-333D-2 and similar strains have very few ALT calls. Are they really different strains? Show a full table of differences.

## Answer in short

1. **No strain is a copy of DH4148.**
   - Each strain in the near-reference set has 188–810 SNPs that differ from the reference.
   - 188 of those sites are shared by the whole set. Almost every other strain in the callset also carries the ALT allele there, so the reference has an allele that no other strain has. These sites are DH4148-specific mutations or errors in the assembly; the data cannot tell which.
   - After those 188 sites are excluded, the strain closest to DH4148 is **EXF_7934** (Sweden, kitchen sink): 2 SNPs. The next closest is EXF_9051 with 83, and the median is 313.
   - **Correction (2026-10-06):** EXF_7934 and EXF_9051 are likely mixed cultures (`docs/strain_identity_issues_2026-10-06.md`, section D). Their mixed ALT calls were masked here, so their low counts are not reliable. Excluding them, the next closest is EXF_2251 with 89.
   - My earlier note ("TFCN_137D-4 has 3 ALT calls; maybe a clone of DH4148") was wrong. That count used only the `rmuc_core` MAF ≥ 0.05 SNP sites on one chromosome. Genome-wide, TFCN_137D-4 has 536 ALT SNPs, and 357 of them are not shared by the whole set.
2. **Most strains in the set are distinct from each other.**
   - The median pairwise difference is about 470 SNPs.
   - 7,096 of the 7,503 strain pairs differ at ≥ 50 SNPs.
3. **52 strains fall into 11 tight clusters (≤ 5 SNPs to the nearest member).** These are candidates for repeat isolates of one clone. Most clusters come from one collection, and many from one sampling site (for example TFCN_209-6-*, TFCN_4M-1-*, TFCN_25-332*). This analysis cannot separate repeat isolates of one strain from closely related isolates of one clonal lineage. Collection records (same sample, same plate?) are needed for that.
4. **One cluster crosses collections: NR06.**
   - TFCN_86C-3 (China, 9003_027 library, no metadata) differs from DBVPG_6742 (Italy) at 0 SNPs, and from DBVPG_4304 (Italy) at ≤ 1 SNP.
   - It differs from every other strain at ≥ 307 SNPs.
   - The most likely explanation is a sample or library mix-up. TFCN_86C-3 is one of the 4 rescued no-metadata haploids. Check the tube or library before it is used as a Chinese isolate.

## Method

Script: `scripts/variant_qc/nearref_differences.sh`, which calls `nearref_differences.py`. Ran as SLURM job 29431120.

- **Near-reference set:** strains with < 5,000 non-ref calls genome-wide in `results/variant_qc/genome.samples.tsv`. That gives 123 strains: 70 TFCN, 42 EXF and 11 DBVPG; 120 of them are in clone group CG001.
- **Input:** the hard-filtered `all` VCF, PASS sites, outside the repeat and low-complexity mask.
- **Candidate sites:** sites with a non-ref call in at least one strain of the set. There are 8,242 SNP sites (8,912 with indels).
- **Genotype rules** (the VARIANT_QC_FILTER rules of 2026-10-05):

  | Genotype | Rule |
  |---|---|
  | missing | GQ < 20 |
  | hom-ref | no DP threshold |
  | non-ref | DP ≥ 5, and the main allele has ≥ 0.8 of the reads for a haploid call |

- **Difference between two strains:** the number of candidate sites where both strains are called and the genotypes differ. Every pair was compared at 7,987–8,238 sites.
- The numbers above use SNPs only (`nearref_snps.*`). The all-variant tables (`nearref.*`) give the same picture.

## Output files (`results/variant_qc/nearref/`)

| File | Content |
|---|---|
| `nearref_snps.strains.tsv` | One row per strain: origin, environment, clone group, ALT calls split into lineage-wide / shared only within the set / shared with strains outside the set / private, nearest strain and its difference, number of strains at 0, ≤ 5 and ≤ 20 differences, median difference to the set |
| `nearref_snps.pairwise.tsv.gz` | All 7,503 pairs: differences and sites compared |
| `nearref_snps.clusters_le5.tsv` | The 11 clusters at ≤ 5 SNPs (single linkage), with the largest difference inside each cluster and the smallest difference to any strain outside it |
| `nearref_snps.sites.tsv.gz` | Per candidate site: ALT carriers and called strains inside and outside the set; lineage-wide flag |
| `nearref.*` | The same tables with indels included |

## Summary over the 123 strains (SNPs)

| Column | min | Q1 | median | Q3 | max |
|---|---|---|---|---|---|
| ALT calls vs DH4148 | 188 | 346 | 495 | 567 | 810 |
| of which lineage-wide | 152 | 182 | 183 | 183 | 187 |
| private (only this strain among all 316) | 0 | 1 | 3 | 14 | 370 |
| differences to nearest strain | 0 | 2 | 7 | 23 | 628 |
| median difference to the set | 314 | 447 | 468 | 549 | 802 |

## Clusters at ≤ 5 SNPs

"max within" is the largest difference between two members. A single-linkage cluster can hold members that differ by more than 5 SNPs. "min to outside" is the smallest difference from a member to any strain outside the cluster.

| Cluster | n | max within | min to outside | Members |
|---|---|---|---|---|
| NR01 | 11 | 10 | 6 | TFCN_152A-5, TFCN_17-332C-1, TFCN_17-332D-1, TFCN_17-332P-2, TFCN_25-332C-1, TFCN_25-332D-2, TFCN_25-332M-1, TFCN_25-332M-2, TFCN_25-332Y-1, TFCN_25-335Y-1, **TFCN_25-337M-3** |
| NR02 | 8 | 9 | 174 | TFCN_102C-1, TFCN_102D-2, TFCN_17-325D-2, TFCN_17-325D-4, TFCN_17-325P-1, TFCN_17-338D-3, TFCN_25-0-2E333-7, TFCN_25-325Y-1 |
| NR03 | 8 | 5 | 14 | TFCN_209-6-1, -2, -4, -5, -11, -12, -14, -17 |
| NR04 | 6 | 7 | 12 | **TFCN_137D-4**, TFCN_17-334Y-2, TFCN_212C-2, TFCN_213-6-2, TFCN_213-6-4, TFCN_86A-12 |
| NR05 | 4 | 2 | 72 | EXF_3417, EXF_3544, EXF_3569, EXF_3612 (all Arctic glacier ice or water) |
| NR06 | 3 | 1 | 307 | DBVPG_4304 (Italy), DBVPG_6742 (Italy), **TFCN_86C-3** (no metadata) |
| NR07 | 3 | 1 | 13 | TFCN_152A-12, **TFCN_17-333D-2**, TFCN_54D-2 |
| NR08 | 3 | 6 | 6 | TFCN_4M-1-3, TFCN_4M-1-4, TFCN_86A-3 |
| NR09 | 2 | 3 | 377 | DBVPG_6741 (Italy, soil), DBVPG_7019 (Italy) |
| NR10 | 2 | 1 | 235 | EXF_13260, EXF_13261 (Mexico, volcano) |
| NR11 | 2 | 1 | 321 | TFCN_17-334C-1, TFCN_342-3 |

All TFCN strains in these clusters are listed as China, marsh/tidal flat, except TFCN_25-332M-2 and TFCN_86C-3, which have no metadata row. NR01 and NR08 are 6 SNPs from their nearest outside strain, so the 5-SNP cut splits them from neighbours only narrowly.

## The strains asked about

| Strain | ALT vs DH4148 | lineage-wide | private | nearest strain (SNPs) | strains ≤ 5 SNPs | median to set |
|---|---|---|---|---|---|---|
| TFCN_137D-4 | 536 | 179 | 0 | TFCN_86A-12 (1) | 5 | 530 |
| TFCN_17-333D-2 | 580 | 182 | 1 | TFCN_152A-12 (1) | 2 | 559 |
| TFCN_25-337M-3 | 496 | 182 | 3 | TFCN_25-332C-1 (4) | 4 | 465 |
| TFCN_25-0-2E332-1 | 498 | 182 | 3 | TFCN_25-332C-1 (7) | 0 | 466 |
| TFCN_86C-3 | 489 | 181 | 0 | DBVPG_6742 (0) | 2 | 463 |

## TFCN_86C-3: what to check in the lab

**Evidence**

| Strain | Library / run | Metadata |
|---|---|---|
| TFCN_86C-3 | one library only: 9003_027, run 222L7FLNX (the 9003 batch) | none in the strain DB; template filled from sibling isolates: China, Guangxi, tidal flat |
| DBVPG_6742 | 23263Sta_DBVPG6742_S189, run 23263Sta | Italy |
| DBVPG_4304 | 23263Sta_DBVPG4304_S97, run 23263Sta | Italy |

- TFCN_86C-3 differs from DBVPG_6742 at 0 SNPs and from DBVPG_4304 at ≤ 1 SNP. The two DBVPG libraries come from a different run, and they match each other as closely.
- The 3 strains differ from every other strain in the callset at ≥ 307 SNPs.
- The library-identity check of 2026-10-03 (`results/variant_qc/library_identity_2026-10-03.tsv`) already listed DBVPG_4304 as the next-best match of the 9003_027 library, at the same distance as its own strain (0.0001).
- The same 9003 batch has other libraries that match a strain other than their label: DBVPG_3445, DBVPG_3538, DBVPG_6649, TFCN_2M-1-3 and TFCN_102D-1 (see `docs/straindb_fixes_proposal_2026-10-04.md`). A mix-up in that batch is therefore possible.
- The sequence data cannot separate "9003_027 contains DBVPG_6742/4304 DNA" from "TFCN_86C-3 is a true member of the same clone found in China and Italy". There is no second, independent TFCN_86C-3 library.

**Suggested check**

1. Take TFCN_86C-3 from the original TFCN stock, not the 9003 DNA.
2. Genotype 2–3 of the 198 diagnostic SNPs in `results/variant_qc/nearref/NR06_diagnostic_snps.tsv`. At these sites the 3 strains share an allele that none of the other 244 called rmuc_core strains has.
   - Original stock carries the NR06 alleles → it is the same clone as the two DBVPG strains, and the library is correct.
   - It does not → the 9003_027 library is not TFCN_86C-3. Exclude it, as for the other mismatched 9003 libraries.
3. Until then, TFCN_86C-3 stays in rmuc_core (rescued, `variant_qc_reviewed_keep.tsv`). In the de-cloned set it falls in the same near-identical group as DBVPG_6742 and DBVPG_4304, so only one of the 3 is used there.

## De-cloned group `rmuc_core_declone` (user decision 2026-10-06: build it)

`scripts/variant_qc/declone.sh` (SLURM job 29533133) compares all 247 rmuc_core strains, not only the near-reference set.
- **Input:** the rmuc_core `.qc` VCF from run 29429262, SNPs only: 435,261 sites. There is no MAF filter, so private SNPs count.
- **Comparison:** every pair was compared at ≥ 421,677 sites.
- **Groups:** single linkage at ≤ 5 SNP differences gives 12 near-identical groups holding 56 strains (`results/variant_qc/declone/rmuc_core.groups_le5.tsv`). All of them are inside CG001. Outside CG001, no pair is ≤ 5 SNPs apart.
- **Representative per group:** the strain with the fewest missing genotypes; ties go to the higher mosdepth depth.
- **Result:** 203 strains (`rmuc_core.representatives_le5.txt`). `make_population_sets.py` writes them as `rmuc_core_declone`. The other 4 groups are unchanged.

| Group | n | max within | min to outside | Strains (* = representative) |
|---|---|---|---|---|
| NI001 | 13 | 12 | 6 | TFCN_152A-5, 17-332C-1, 17-332C-2, 17-332D-1, 17-332P-1, 17-332P-2, 25-332C-1, 25-332D-2, 25-332M-1*, 25-332M-2, 25-332Y-1, 25-335Y-1, 25-337M-3 |
| NI002 | 8 | 9 | 239 | TFCN_102C-1, 102D-2, 17-325D-2, 17-325D-4, 17-325P-1, 17-338D-3, 25-0-2E333-7, 25-325Y-1* |
| NI003 | 8 | 4 | 11 | TFCN_209-6-1, -2, -4*, -5, -11, -12, -14, -17 |
| NI004 | 6 | 7 | 10 | TFCN_137D-4, 17-334Y-2*, 212C-2, 213-6-2, 213-6-4, 86A-12 |
| NI005 | 4 | 2 | 68 | EXF_3417*, EXF_3544, EXF_3569, EXF_3612 |
| NI006 | 3 | 1 | 309 | DBVPG_4304*, DBVPG_6742, TFCN_86C-3 |
| NI007 | 3 | 1 | 7 | TFCN_152A-12, 17-333D-2, 54D-2* |
| NI008 | 3 | 6 | 6 | TFCN_4M-1-3*, 4M-1-4, 86A-3 |
| NI009 | 2 | 2 | 375 | DBVPG_6741, DBVPG_7019* |
| NI010 | 2 | 1 | 244 | EXF_13260, EXF_13261* |
| NI011 | 2 | 5 | 249 | EXF_1565, EXF_1676* |
| NI012 | 2 | 1 | 308 | TFCN_17-334C-1, TFCN_342-3* |

**The cutoff is a choice, not a gap in the data.** Representatives at other cutoffs:

| Cutoff (SNPs) | Groups | Strains grouped | Representatives |
|---|---|---|---|
| ≤ 2 | 13 | 42 | 218 |
| ≤ 5 | 12 | 56 | 203 |
| ≤ 10 | 17 | 74 | 190 |
| ≤ 20 | 18 | 91 | 174 |
| ≤ 50 | 23 | 110 | 160 |

The pair counts rise smoothly: 7 pairs at 0 SNPs, 52 at ≤ 2, 113 at ≤ 5, 209 at ≤ 10. NI001 and NI008 lie 6 SNPs from their nearest outside strain. To use another cutoff, run `sbatch scripts/variant_qc/declone.sh N`, change `DECLONE` in `make_population_sets.py`, and rerun the pipeline.

**Trees and CG001.** The strain trees use the MAF ≥ 0.05 SNP set. Private and rare SNPs, which separate members of CG001, are not in it. In the tree figures CG001 is therefore a near-flat comb, and the members of one near-identical group are not always placed next to each other. The trees do not resolve structure inside CG001.

## Open decisions (not acted on)

- **TFCN_86C-3:** check the tube or library against DBVPG_6742 and DBVPG_4304.
- **EXF_7934:** withdrawn. It is a likely mixed culture (see the correction above).
