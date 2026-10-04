# Sample ID and library fixes (2026-10-03)

This document records:
- every change to strain IDs or library assignment in this project
- the evidence for each change
- the files that were changed

## Sources compared

- **Our files:** `samplesheet.csv` (Sarek input; the FASTQs behind each CRAM), `flat_cram_dir/`, `metadata.txt`, `ploidy_overrides.csv`.
- **Strain database:** `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_StrainDB/db/strains.duckdb`, tables `strains`, `strain_aliases` and `sequencing_runs`. Also its MS2 extended metadata and `strain_summary_YPD2`.
- **SeqCoast run 9003:** `/bigdata/stajichlab/shared/projects/SeqData/SeqCoast/9003_20260831__Rhodotorula35/sample_sheet.csv`. It maps 35 libraries (9003_001–035) to strain names.
- **Phenotypes:**
  - `ExRhodotorula_Phenotypes/strains.csv`
  - the zinc screen, `Rhodotorula_Phenotyping/Heavy_Metals/ArrayedHeavyMetalScreen/`

Cross-check script: `scripts/variant_qc/sample_id_crosscheck.py`, which writes `results/variant_qc/sample_id_crosscheck.tsv`. It runs with the StrainDB pixi python, which has duckdb.
- 254 of 319 strains pass all checks.
- The flags are listed below.

## Fixes made

| Change | Evidence | Files changed |
|---|---|---|
| Libraries 9003_002 (labelled TFCN_BY120-C1) and 9003_028 (TFCN_BY120-C7) merged into **TFCN_25-0-2E333-9** | BY120-C1/C7 are plate positions, not strain IDs (confirmed by the user). The strain DB entry "BY120-C1 (BY120-C7) 25-0-2E333-9" maps to `TFCN:25-0-2E333-9`. Genetic distance C1 vs C7 = 0.00009 (identical). Each is 0.00012–0.00013 from the original library, inside its clone group. | merged CRAM `results/preprocessing/merged/TFCN_25-0-2E333-9/`; `flat_cram_dir` link; BY120 links removed; `samplesheet.csv` (lanes L002, L003); `ploidy_overrides.csv` (BY120 rows removed) |
| CCFEE_5036 merged into **DBVPG_5227** | The strain DB lists "CCFEE 5036" as one of DBVPG_5227's other collection numbers. The user confirmed they are the same isolate. Genetic distance 7×10⁻⁵, in clone group CG005. DBVPG_5227 has a phenotype row; CCFEE_5036 has none. | merged CRAM `results/preprocessing/merged/DBVPG_5227/`; `flat_cram_dir` link; CCFEE_5036 link removed; `samplesheet.csv` (lanes L003, L004); `ploidy_overrides.csv` (CCFEE_5036 row removed) |
| `TFCN_186CL④-2` → `TFCN_186CL-2` | entry error in the SeqCoast sheet (user) | `archive/sample_sheet.csv`, `archive/samples.csv` (backups `*.before_2026-10-03`). SeqCoast's own `sample_sheet.csv` was not changed, at the user's request. |

Notes on both merges:
- `scripts/merge_strain_crams.sh` merges the markduplicates CRAMs.
- Each read group keeps its own ID and LB. SM is set to `<STRAIN>_<STRAIN>`.
- The script checks that the merged record count equals the sum of the inputs, and that the duplicate-flag counts match.
- For TFCN_25-0-2E333-9: 21,615,627 records and 8,547,969 duplicate-flagged.

`scripts/variant_qc/make_population_sets.py` (`MERGED_INTO`) keeps the merged names out of every group.

## Findings not changed here

- **9003 batch not in the strain DB.** The DB `sequencing_runs` table lists none of the 29 9003 libraries that we used, and 17 of those strains have no DB entry at all. Our library-to-strain mapping matches SeqCoast's `sample_sheet.csv` exactly. Six 9003 libraries were never aligned: 004 DBVPG_3380, 008 DBVPG_3985, 014 TFCN_152C-2, 024 TFCN_7-9-2, 025 TFCN_7-9-3, 033 EXF_13474. **These six are species other than R. mucilaginosa (confirmed by the user), so they stay out of this project by design.**
- **Library identity.** Ten strains merge a 9003 library with an older library. `scripts/variant_qc/library_identity.sh` tests each library against the joint callset; see its results when complete.
- **Wrong species in the DB.** The strain DB lists 20 strains as *R. mucilaginosa*, but they genotype as *R. frigidialcoholis* (12) or *R. aff. mucilaginosa* (8), as `metadata.txt` already says. This needs fixing in the DB.
- **Phenotype-table species.** 15 `rmuc_core` strains are labelled as other species in `ExRhodotorula_Phenotypes/strains.csv` (see `docs/population_set_design_2026-10-03.md`).
- **DBVPG_Y3853 (zinc file).** This is DBVPG_3853, *R. dairenensis*, so it is not in this project.
- **TFCN_17-332D-2, TFCN_223D-8 and TFCN_17-333P-8 (zinc file).** These are different strains from the similar genotyped names (user). They have no genotype data.
- **Zinc phenotype files.** `Data/Interm/ZnArrayRun/zincmeasurementmeta.csv` holds 2 of the 5 runs: 159 strains, 92 of them in `rmuc_core`. `0.15.1_Analysis/Results/Zinc/` holds all 5 runs: 321 strains, 161 in `rmuc_core`, including TFCN_25-0-2E333-9. Use the 0.15.1 results.
- **Duplicate metadata rows (fixed).** `metadata.txt` had exact duplicate rows for six strains: TFCN_25-332D-2, DBVPG_3775, TFCN_17-332D-1, TFCN_25-332D-1, TFCN_25-332M-1 and TFCN_270H-2. One copy of each was removed (607 → 601 lines). The backup is `metadata.before_2026-10-03.txt`, and CRLF line endings were kept. Four strains (DBVPG_3045, DBVPG_3855, TFCN_17-0-2E334-4, TFCN_17Y-278-1) keep three rows each: one diploid row plus `_hap1`/`_hap2` haplotype rows. These are intended, and the pipeline excludes `_hap` rows.

## Mismatched 9003 libraries and CRAM rebuilds

The library identity check (`scripts/variant_qc/library_identity.sh`; results in `results/variant_qc/library_identity_2026-10-03.tsv`) genotyped each read group separately on CM179498.1 and compared it with all strains in the joint callset.

**Five strains had a 9003 library from a different organism.** Sarek had merged each one into the strain's original library, because `samplesheet.csv` listed it as an extra lane of the same sample:

| Strain | Original library | 9003 library (excluded) | Ploidy after rebuild |
|---|---|---|---|
| DBVPG_3445 | hybrid diploid | 9003_005: pure *R. muc.* haploid ≈ TFCN_98C-7 (China) | diploid (unchanged) |
| DBVPG_3538 | pure *R. muc.* haploid ≈ DBVPG_3382 clone | 9003_007: *R. aff. mucilaginosa* ≈ TFCN_25-395P-1 (China) | **haploid** (was diploid) |
| DBVPG_4379 | pure *R. muc.* haploid ≈ DBVPG_10842 | 9003_009: hybrid diploid | **haploid** (was diploid) |
| DBVPG_4952 | hybrid diploid | 9003_010: pure haploid near *R. frigidialcoholis* | diploid (unchanged) |
| DBVPG_6649 | hybrid diploid | 9003_012: pure *R. muc.* haploid ≈ TFCN_17-325D-4 (China) | diploid (unchanged) |

What was done:
- `scripts/rebuild_strain_cram.sh` rebuilt each CRAM from its original read groups and re-marked duplicates. The output is in `results/preprocessing/rebuilt/`.
- The five 9003 rows moved to `samplesheet_excluded_9003.csv`.
- On 2026-10-03, DBVPG_3538 and DBVPG_4379 were changed to haploid in `ploidy_overrides.csv`.

**Hybrid status.** Hybrid-like libraries peak at a minor-allele fraction of 0.40–0.50, so they are real diploid hybrids, not mixed cultures.
- DBVPG_3446 and TFCN_270H-1 are hybrids in two independent libraries.
- DBVPG_3538 and DBVPG_4379 looked hybrid only because of the merge.

**Consistent merges, kept as they are:** DBVPG_3239, DBVPG_3446, TFCN_102D-1, TFCN_25-332D-2, TFCN_270H-1 and TFCN_25-0-2E333-9.

**DBVPG_6094 is unverified.** Its original libraries gave no usable sites (the database notes say "Too Low"), so its genotype rests on 9003_011 alone.

**Single-library 9003 strains cannot be cross-checked.**
