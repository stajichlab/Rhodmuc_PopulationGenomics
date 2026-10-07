# R. mucilaginosa DH4148 population genomics: summary of findings (2026-10-06)

This document summarises the work of 2026-10-04 to 2026-10-06. The detail is in the linked reports. Open work is tracked as GitHub issues (table at the end).

## 1. State of the data

- **Callset:** 316 strains, GATK 4.5 joint genotyping on the DH4148 assembly (GCA_058775505.1).
- **Group VCFs and trees:** pipeline runs 29429262 and 29533184, pipeline branch `fix/homref-dp-mask`.

| Group | Strains | Tree alignment (SNP columns) |
|---|---|---|
| rmuc_core | 247 | 235,622 |
| rmuc_core_declone | 203 | 242,736 |
| rmuc_core_outgroup | 264 | 293,720 |
| rmuc_with_hybrids | 280 | 166,617 |
| hybrid_diploids | 33 | 176,368 |
| all | 316 | 2,469 |

- **Files:**
  - group VCFs: `results/RmucDH4148.<group>.{qc,snps.maf}.annotated.vcf.gz`
  - trees: `results/strain_tree/RmucDH4148.<group>.snps.maf.treefile`
  - tree figures: the same name ending `.tree.png` / `.tree.pdf`
- **Group definitions:** `scripts/variant_qc/make_population_sets.py`. Reviewed exceptions: `variant_qc_reviewed_keep.tsv` (7 strains).

## 2. Findings

### 2.1 A genotype-filter artifact removed most hom-ref calls in strains near the reference

- **Cause:** GenotypeGVCFs gives a hom-ref genotype the depth of its gVCF reference block (block MIN_DP), not the read depth at the site. The rule `DP < 5 → missing` therefore removed 40–89 % of hom-ref calls in 6 strains with 16–28× coverage.
- **Effect:**
  - The tree alignments fell from about 17,000 to about 700 sites.
  - The same artifact made the strain-QC rule drop those 6 strains as `low_coverage`.
- **Fixes** (user decisions, 2026-10-05):
  - no DP mask on hom-ref calls; GQ < 20 still applies
  - `low_coverage` from mosdepth read depth
- **Who is affected:** haploid strains only. For a diploid, GQ < 20 already removes every hom-ref call with fewer than 5 reads.
- Report: `docs/decision_homref_dp_mask_2026-10-05.md`.

### 2.2 A tree-step bug dropped identical strains

- When IQ-TREE's +ASC model stopped on invariant columns, the pipeline reran on `varsites.phy`. IQ-TREE writes that file after it collapses identical sequences.
- So the `rmuc_with_hybrids` tree had 74 of 280 strains. The run of 2026-10-04 had the same bug: 96 of 267.
- Fix: `-keep-ident`. All trees now hold all of their strains.

### 2.3 Clone group CG001 and near-identical strains

- **CG001 is a group, not a single clone.** It has 120 strains (70 TFCN, 39 EXF, 11 DBVPG). They differ from each other at a median of about 470 SNPs.
- **No strain is a copy of DH4148.** 188 sites are fixed for ALT in the whole group and in almost all other strains. These are DH4148-specific alleles or assembly errors.
- **Near-identical groups:** 56 strains form 12 groups at ≤ 5 SNPs. Most come from one collection and often one sampling site. `rmuc_core_declone` keeps one strain per group.
  - The 5-SNP cutoff is a choice; the data show no gap.

    | Cutoff (SNPs) | Strains kept |
    |---|---|
    | ≤ 2 | 218 |
    | ≤ 5 | 203 |
    | ≤ 10 | 190 |
    | ≤ 20 | 174 |
    | ≤ 50 | 160 |
- **The trees cannot resolve CG001.** They use SNPs at MAF ≥ 0.05, and the rare SNPs that separate CG001 strains are not in that set.
- Report: `docs/nearref_strain_differences_2026-10-05.md`.

### 2.4 Strain identity: mix-ups and contamination

Full list with evidence: `docs/strain_identity_issues_2026-10-06.md` and `results/variant_qc/strain_identity_issues_2026-10-06.tsv` (66 rows).

- **SeqCoast 9003 batch.**
  - 5 libraries carry another organism than their label: DBVPG_3445, 3538, 4379, 4952 and 6649. That is 5 of the 11 strains in the batch that could be cross-checked. These libraries are excluded.
  - 3 more libraries are *R. mucilaginosa* where the StrainDB says *R. glutinis* or *R. toruloides*: DBVPG_3380, DBVPG_3985 and TFCN_152C-2.
- **TFCN_86C-3 (China, 9003 only) is identical to DBVPG_6742 and DBVPG_4304 (Italy).** It differs from them at 0–1 SNPs, and from every other strain at ≥ 307. 198 diagnostic SNPs are listed for a lab check.
- **Possible mixed cultures (8):**
  - EXF_1695 (still in rmuc_core and 3 other groups)
  - EXF_12768, EXF_7934, EXF_8891, EXF_9051, EXF_5668, TFCN_3M-1-1, EXF_14606
- **Phenotype-table species conflicts.**
  - 16 rmuc_core strains are another species in `ExRhodotorula_Phenotypes/strains.csv`.
  - They are enriched in the near-identical groups: 11/42 against 5/123 among strains with a phenotype row (Fisher OR 8.4, p = 1.5 × 10⁻⁴).
  - Either the phenotyped and sequenced cultures differ, or the phenotype table species is wrong. These data cannot tell which.
- **Labels elsewhere.**
  - StrainDB: 20 strains have the wrong species.
  - `popgen_strain_metadata.csv`: TFCN_363-1-2 is not *R. taiwanensis*.
  - Reference genomes: 5 mislabels, relabelled in Rodeo and in Species_ID_db.

### 2.5 Species ID database (sourmash)

- **Location:** `/bigdata/stajichlab/shared/projects/Rhodotorula/Species_ID_db`. Report: `REPORT_2026-10-04.md`.
- **Validation:** leave-one-out on 514 reference genomes; every label with ≥ 2 correctly named references agrees 100 %.
- **Use:** species calls for the 15 new 9003 strains, and detection of the hybrids.

## 3. Decisions made (user)

| Date | Decision |
|---|---|
| 2026-10-04 | `rmuc_pheno` dropped; association pruning moves to the phenotype project |
| 2026-10-04 | StrainDB: 12 strains → *R. frigidialcoholis*; *R. aff. mucilaginosa* stays *R. mucilaginosa* for now; DBVPG:5227 primary (CCFEE 5036 alias); manual alias file; hybrid naming `Rhodotorula <p1> x Rhodotorula <p2> hybrid` |
| 2026-10-05 | Rescue the no-metadata haploids and hybrids and the low-coverage strains; TFCN_17-333M-1 and TFCN_186CL-2 → hybrid_diploids |
| 2026-10-05 | No DP mask on hom-ref calls; IQ-TREE `-keep-ident`; `low_coverage` from mosdepth |
| 2026-10-06 | Build `rmuc_core_declone` (one strain per ≤ 5-SNP group) |

## 4. Open work (GitHub issues)

Tracking issue: [Rhodmuc_PopulationGenomics#7](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/7)

| Issue | Topic | Needs |
|---|---|---|
| [#1](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/1) | EXF_1695 (likely mixed culture) is still in rmuc_core | decision; re-isolation |
| [#2](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/2) | 16 phenotype-table species conflicts, enriched in near-identical groups | ITS on the phenotyped cultures |
| [#3](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/3) | SeqCoast 9003 batch: 5 mismatched libraries, 3 species conflicts, DBVPG_6094 | lab tube check, plate map |
| [#4](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/4) | TFCN_86C-3 = DBVPG_6742 / DBVPG_4304 | diagnostic SNPs in the original stock |
| [#5](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/5) | `metadata.txt` rows for the rescued strains | data entry |
| [#6](https://github.com/stajichlab/Rhodmuc_PopulationGenomics/issues/6) | CG001-only tree from all SNPs | analysis |
| [nf_genotype_population#7](https://github.com/stajichlab/nf_genotype_population/issues/7) | merge `fix/homref-dp-mask` (pushed) | review and merge |
| [Rhodotorula_StrainDB#1](https://github.com/stajichlab/Rhodotorula_StrainDB/issues/1) | apply the StrainDB fixes; rename the 43 hybrids? | edits to the import files |
| [ExtremeRhodotorula_Rodeyo#1](https://github.com/stajichlab/ExtremeRhodotorula_Rodeyo/issues/1) | reload the species table after 4 relabels | DB reload |
