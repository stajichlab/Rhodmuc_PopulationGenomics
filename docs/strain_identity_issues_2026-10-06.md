# Strain identity: possible mix-ups and contamination (2026-10-06)

This report collects every strain-identity problem found in the DH4148 popgen dataset up to 2026-10-06, with its evidence. Each row in the table below has a matching row in **`results/variant_qc/strain_identity_issues_2026-10-06.tsv`** (66 rows; one per strain or record). The TSV columns are: id, strain, library_or_record, category, finding, evidence, confidence, status, recommended_action and source.

**Confidence terms**

| Term | Meaning |
|---|---|
| data-supported | the measurement shows it directly |
| inferred | an explanation that fits the data but was not tested |
| unresolved | the data cannot separate the explanations |

**Data used**

- 316 genotyped strains, joint callset of 2026-10-04 (GATK, DH4148 reference)
- per-library identity check on CM179498.1 (2026-10-03)
- sourmash species ID (`Species_ID_db`, 2026-10-04)
- near-identical strain groups (2026-10-05/06)
- allele balance per strain, and the ALT-read fraction at ALT calls (SLURM job 29534407)
- the phenotype table `ExRhodotorula_Phenotypes/strains.csv`

## Summary

| Category | n | Main point | Status |
|---|---|---|---|
| A. 9003 library is another organism than its label | 5 | DBVPG_3445, 3538, 4379, 4952, 6649 | Libraries excluded, CRAMs rebuilt. Tubes not checked |
| B. 9003 library species ≠ StrainDB species | 3 | DBVPG_3380, DBVPG_3985, TFCN_152C-2 are *R. mucilaginosa* in the library | Not used here. Tube or DB wrong: unknown |
| C. Strains from two collections are identical | 1 | TFCN_86C-3 (China) = DBVPG_6742 / DBVPG_4304 (Italy) | Kept. Lab check designed |
| D. Possible mixed cultures | 8 | EXF_1695 (**still in rmuc_core**), EXF_12768, EXF_7934, EXF_8891, EXF_9051, EXF_5668, TFCN_3M-1-1, EXF_14606 | 1 kept, 7 dropped |
| E. Species label wrong or in conflict | 37 rows | 16 phenotype-table conflicts, enriched in near-identical groups; 20 StrainDB species; 1 metadata file | Fixes proposed, not applied |
| F. Reference genome mislabels | 5 | C2-MEA2-4, P1-MEA2-1, P1-RM0-1, TFCN_152C-6, public "*R. mucilaginosa*" genome | Relabelled; Rodeo DB reload pending |
| G. Name and ID errors | 5 | plate positions used as names, one isolate with two IDs, typo, date artifact, leading zeros | Resolved in project files |
| H. Identity cannot be verified | 2 rows | DBVPG_6094; strains with only a 9003 library | Provisional |

## A. 9003 libraries that are not the labelled strain

Each library was genotyped alone and compared with all strains in the callset. Sourmash agrees for all 5.

| Strain | Library | What the library is | Original library |
|---|---|---|---|
| DBVPG_3445 | 9003_005 | *R. mucilaginosa* haploid ≈ TFCN_98C-7 | hybrid diploid |
| DBVPG_3538 | 9003_007 | *R. aff. mucilaginosa* ≈ TFCN_25-395P-1 | *R. muc.* haploid (DBVPG_3382 clone) |
| DBVPG_4379 | 9003_009 | *R. muc.* × *R. aff. muc.* hybrid diploid | *R. muc.* haploid ≈ DBVPG_10842 |
| DBVPG_4952 | 9003_010 | near *R. frigidialcoholis*; sourmash also scores *R. muc.* 0.977–0.997, so it may be a hybrid or a mixture | hybrid diploid |
| DBVPG_6649 | 9003_012 | *R. mucilaginosa* haploid ≈ TFCN_17-325D-4 | hybrid diploid |

All 5 libraries carry DBVPG labels. In 4 of them, the best match is a TFCN strain: TFCN_98C-7, TFCN_25-395P-1, TFCN_17-390M-1 and TFCN_17-325D-4. The fifth matches DBVPG_6660 (*R. frigidialcoholis*). In the 9003 batch, 5 of the 11 strains with a usable second library failed this check. A mix-up during preparation of the 9003 batch fits this pattern. That is an inference; the batch records were not examined.

## B. 9003 libraries whose species disagrees with the StrainDB

DBVPG_3380 (9003_004) and DBVPG_3985 (9003_008) are *R. glutinis* in the DB. TFCN_152C-2 (9003_014) is *R. toruloides*. All 3 libraries score 1.000 against *R. mucilaginosa*. There is no genome of the original culture, so the data cannot show whether the tube or the DB record is wrong.

## C. TFCN_86C-3 = DBVPG_6742 / DBVPG_4304

- TFCN_86C-3 has a single library, 9003_027. It is 0 SNPs from DBVPG_6742 and ≤ 1 SNP from DBVPG_4304; both of those are from Italy and were sequenced on run 23263Sta.
- The 3 strains are ≥ 307 SNPs from every other strain.
- The data cannot separate a mix-up from a true clone found in both countries.
- **Check:** genotype 2–3 of the 198 diagnostic SNPs (`results/variant_qc/nearref/NR06_diagnostic_snps.tsv`) in the original TFCN stock.
- Details: `docs/nearref_strain_differences_2026-10-05.md`.

## D. Possible mixed cultures

**Measures**
- **Mixed sites per 1,000 covered:** the rate of sites where a haploid call has a second allele in the reads. The haploid median is 0.04 and the 95th percentile is 0.72. *R. frigidialcoholis* and *R. aff. mucilaginosa* sit at 0.6–0.8 because their reads map to this reference with more mismatches.
- **Main-allele fraction at ALT calls (DP ≥ 10):** in a pure haploid, nearly all ALT calls have the main allele at ≥ 0.9. Clean CG001 control EXF_8006: 313 of 322.

| Strain | Kept? | Evidence | Reading |
|---|---|---|---|
| **EXF_1695** | **yes, rmuc_core** | 22.3 mixed sites/1,000; minor allele peaks at 0.10; 85,508 of 163,435 ALT calls have the main allele at 0.8–0.9 | haploid plus a minor second genotype at about 10–20 % (inferred) |
| EXF_12768 | no | 3,787 of 4,469 ALT calls at 0.5–0.6; nQuire diploid; override haploid | 50:50 mixture of two strains, or a low-het diploid (unresolved) |
| EXF_7934 | no | 142 of 329 ALT calls at 0.5–0.8 | CG001 strain plus a second, REF-like genotype at about 20–50 % (inferred) |
| EXF_8891 | no | 166 of 440 ALT calls at 0.7–0.9 | minor second genotype at about 10–30 % (inferred) |
| EXF_9051 | no | 89 of 326 ALT calls at 0.7–0.9 | as EXF_8891 |
| EXF_5668 | no | diploid call; 209.5 mixed sites/1,000; peak at 0.20, not 0.5 | mixture, not a balanced diploid (inferred) |
| TFCN_3M-1-1 | no | *R. frigidialcoholis*; 108 mixed sites/1,000 (species median 0.72); peak at 0.10 | mixture (inferred) |
| EXF_14606 | no | sourmash: *R. glutinis* 0.995 and sp_clade_I 0.989, which are 0.924 apart | hybrid or mixture (unresolved) |

**Correction to `docs/nearref_strain_differences_2026-10-05.md`.** That report named EXF_7934 as the strain closest to DH4148, 2 SNPs apart outside the lineage-wide sites. EXF_7934 is a likely mixture, and its mixed ALT calls were masked in that analysis, so the 2-SNP figure is not reliable. The minor component carries the REF allele at those sites. A contaminating culture of DH4148 itself would give this signal, but so would any CG001 strain that lacks those ALT alleles. Not tested.

**EXF_1695 is the only likely mixed culture still in an analysis group** (`rmuc_core`, `rmuc_core_declone`, `rmuc_core_outgroup`, `rmuc_with_hybrids`). It is already listed in `rmuc_pheno_exclusions.tsv`. Removing it from the popgen groups is a decision for you; I have not changed the groups.

## E. Species labels

**E1. Phenotype table vs genotype (16 rmuc_core strains).** `ExRhodotorula_Phenotypes/strains.csv` gives another species for each of them: *R. sphaerocarpa* 7, *R. toruloides* 2, *R. diobovata* 2, and *R. nothofagi*, *R. glutinis*, *R. graminis*, *R. paludigena* and *R. taiwanensis* 1 each. All 16 genotype as *R. mucilaginosa*.

- **The conflicts cluster in the near-identical groups** (≤ 5 SNPs; `results/variant_qc/declone/`). Among the 165 rmuc_core strains with a phenotype row:

  | | Conflict | No conflict |
  |---|---|---|
  | in a near-identical group | 11 | 31 |
  | not in a group | 5 | 118 |

  Fisher exact test: OR 8.4, p = 1.5 × 10⁻⁴.
- Affected groups:
  - NI004: 4 strains
  - NI002 and NI007: 2 each
  - NI001, NI008 and NI012: 1 each
- **Possible explanations (not tested):**
  1. Some sequenced cultures were the same *R. mucilaginosa* clone, from a contaminant or a mislabelled stock, while the phenotyped cultures were the intended species.
  2. The species in the phenotype table is wrong for these isolates.
- The data cannot tell 1 from 2. An ITS check of the phenotyped cultures for these 16 strains would.
- If explanation 1 holds, the phenotype values for these strains do not belong to the sequenced genomes.

**E2/E3. StrainDB species.** The StrainDB lists 20 strains as *R. mucilaginosa*. 12 genotype as *R. frigidialcoholis*; the fix is proposed but not applied. 8 genotype as *R. aff. mucilaginosa*; they stay *R. mucilaginosa* in the DB for now, by your decision. See `docs/straindb_fixes_proposal_2026-10-04.md`, Fix A.

**E4.** `import/popgen_strain_metadata.csv` calls TFCN_363-1-2 *R. taiwanensis*. Its reads are *R. mucilaginosa* (mapping 96.4 %; sourmash 1.000). That file is not yet corrected.

## F. Reference genome mislabels (Rodeo / Species_ID_db)

| Genome | Label | k-mer ANI places it in |
|---|---|---|
| C2-MEA2-4 | *R. mucilaginosa* | *R. frigidialcoholis* |
| P1-MEA2-1 | *R. mucilaginosa* | sp.1 |
| P1-RM0-1 | *R. mucilaginosa* | sp.1 |
| TFCN_152C-6 | *R. toruloides* | sp.2 |
| `Public_genomes/Rhodotorula_mucilaginosa.fasta` | *R. mucilaginosa* | *R. aff. mucilaginosa* |

- The first 4 are relabelled in `Rhodotorula_Rodeo/genomes.csv`. The Rodeo DB species table still needs a reload; see `Rhodotorula_Rodeo/NOTE_species_relabel_2026-10-04.md`.
- In `samples.csv`, the TFCN_152C-6 ASMID still contains "toruloides".

## G. Name and ID errors (resolved in project files)

- 9003 libraries labelled TFCN_BY120-C1 and TFCN_BY120-C7 are plate positions. Both were merged into TFCN_25-0-2E333-9.
- CCFEE_5036 and DBVPG_5227 are one isolate. Merged into DBVPG_5227.
- The SeqCoast sheet has the typo `TFCN_186CL④-2`.
- The alternative ID "7-6-2003" for TFCN_7-6-3 looks like a spreadsheet date conversion (inferred).
- 11 TFCN IDs carry a leading zero in the strain metadata csv.

## H. Identity cannot be verified

- **DBVPG_6094:** its original libraries gave no usable reads, so its genotype rests on 9003_011 alone.
- **Strains with only a 9003 library:** these cannot be cross-checked: TFCN_1A-1-5, TFCN_2M-1-3, TFCN_86C-3, TFCN_25-332M-2, TFCN_17-333M-1, TFCN_186CL-2, EXF_10854 and EXF_17335, plus the 8 non-*R. mucilaginosa* strains. Given the failure rate in the 9003 batch (5 of 11), treat them as provisional.

## Suggested lab checks, by expected value

1. **ITS (or diagnostic SNPs) on the phenotyped cultures of the 16 E1 strains**, starting with group NI004 (TFCN_212C-2, 213-6-2, 213-6-4, 86A-12). This decides whether phenotype and genome data belong to the same culture.
2. **9003 batch tubes:** DBVPG_3445, 3538, 4379, 4952, 6649, 3380, 3985, TFCN_152C-2 and DBVPG_6094. Compare tube contents with the collection stock.
3. **TFCN_86C-3:** diagnostic SNPs in the original TFCN stock.
4. **Single-colony re-isolation** of EXF_1695, then EXF_7934, EXF_8891, EXF_9051 and EXF_12768, if they are needed.

## Limits

- The library-identity check used one chromosome (CM179498.1, about 129,000 sites). Distances inside clone groups (0.0001–0.0004) tie, so a "best match" other than the strain's own name inside a clone group is not evidence of a mix-up.
- Mixed-read measures cannot tell a mixed culture from a diploid when the two genotypes are present at about 50:50 (EXF_12768).
- The near-identical groups use a 5-SNP cutoff chosen by me, not by a gap in the data.
