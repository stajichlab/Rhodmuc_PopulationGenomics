# Proposed Rhodotorula_StrainDB fixes (2026-10-04) — FOR REVIEW

Status: proposal only. Nothing in `Rhodotorula_StrainDB` has been changed. User decisions of 2026-10-04 are in "Decisions" at the end.

## Where the fixes must go

`run_all_imports.sh` → `metadb.build_db` deletes `db/strains.duckdb` and rebuilds it from `import/*.csv`. A direct edit of the `.duckdb` file is lost at the next rebuild. So each fix goes into a source file:

- Species and sequencing runs come from `import/Master_List_Rhodotorula_Sequencing.csv` (columns `Species`, `FileBase`, `Notes`).
- Aliases come only from the master list `ID`/`Strain` columns, the DBVPG catalog `Other collection numbers`, and BioSample strain fields. The build has no input for manual aliases or for "same isolate" links. `enrich_metadata.py` can write `manual_correction` aliases, but `build_db` does not call it.

Evidence for every item below: `docs/sample_id_fixes_2026-10-03.md` and `results/variant_qc/sample_id_crosscheck.tsv`. The genetic species come from the 2026-10-03 callset. None of these 20 strains had a CRAM rebuilt, so the new callset should not change them. I will confirm this in the step-1 QC.

## Fix A — species for 20 strains

The DB says *R. mucilaginosa* for all 20. `metadata.txt` and the genotypes agree on another species.

| Strain | DB strain_id | Proposed `Species` |
|---|---|---|
| CHIFNET_14NJ214 | CHIFNET:14NJ214 | Rhodotorula frigidialcoholis |
| DBVPG_6660 | DBVPG:6660 | Rhodotorula frigidialcoholis |
| EXF_13250 | EXF:13250 | Rhodotorula frigidialcoholis |
| EXF_13255 | EXF:13255 | Rhodotorula frigidialcoholis |
| EXF_13601 | EXF:13601 | Rhodotorula frigidialcoholis |
| EXF_3629 | EXF:3629 | Rhodotorula frigidialcoholis |
| EXF_3633 | EXF:3633 | Rhodotorula frigidialcoholis |
| EXF_7386 | EXF:7386 | Rhodotorula frigidialcoholis |
| TFCN_134A-3 | TFCN:134A-3 | Rhodotorula frigidialcoholis |
| TFCN_17-332M-1 | TFCN:17-332M-1 | Rhodotorula frigidialcoholis |
| TFCN_1A-14 | TFCN:1A-14 | Rhodotorula frigidialcoholis |
| TFCN_3M-1-1 | TFCN:3M-1-1 | Rhodotorula frigidialcoholis |
| DBVPG_4380 | DBVPG:4380 | no change (decision 1) |
| DBVPG_4534 | DBVPG:4534 | no change (decision 1) |
| DBVPG_8043 | DBVPG:8043 | no change (decision 1) |
| EXF_5666 | EXF:5666 | no change (decision 1) |
| TFCN_25-333Y-10 | TFCN:25-333Y-10 | no change (decision 1) |
| TFCN_25-334Y-6 | TFCN:25-334Y-6 | no change (decision 1) |
| TFCN_25-395P-1 | TFCN:25-395P-1 | no change (decision 1) |
| TFCN_33A-4 | TFCN:33A-4 | no change (decision 1) |

Edit: the `Species` cell of each strain's rows in the master list. Add to `Notes`: "species from genotype, DH4148 popgen 2026-10".

For DBVPG strains, `taxon_name_full` still comes from the DBVPG catalog and will still say *R. mucilaginosa*. That is the collection's record, so I propose to leave it.

## Fix B — add the SeqCoast 9003 libraries

The DB has 0 of the 35 libraries. FASTQs: `/bigdata/stajichlab/shared/projects/SeqData/SeqCoast/9003_20260831__Rhodotorula35/9003_Illumina/`. `FileBase` pattern: `9003_0NN_SNN_R[12]_001.fastq.gz`.

### B1. Strain is in the DB, library is consistent with it (7) — add one run row each

| Library | FileBase | Strain |
|---|---|---|
| 9003_001 | 9003_001_S29_R[12]_001.fastq.gz | TFCN:102D-1 |
| 9003_003 | 9003_003_S31_R[12]_001.fastq.gz | DBVPG:3239 |
| 9003_006 | 9003_006_S34_R[12]_001.fastq.gz | DBVPG:3446 |
| 9003_011 | 9003_011_S39_R[12]_001.fastq.gz | DBVPG:6094 (unverified: no usable original library to compare) |
| 9003_019 | 9003_019_S47_R[12]_001.fastq.gz | TFCN:25-332D-2 |
| 9003_021 | 9003_021_S49_R[12]_001.fastq.gz | TFCN:270H-1 |

Also add 9003_020 (`9003_020_S48_R[12]_001.fastq.gz`) to TFCN:25-332M-2. It is the only library this project used for that strain. It genotypes as *R. mucilaginosa*. The DB's own run (ISQ0044552_S287) is noted "Too Low", so the 9003 library cannot be compared with it.

### B2. Plate-position libraries (2) — add as runs of TFCN:25-0-2E333-9

| Library | FileBase | Strain |
|---|---|---|
| 9003_002 | 9003_002_S30_R[12]_001.fastq.gz | TFCN:25-0-2E333-9 |
| 9003_028 | 9003_028_S56_R[12]_001.fastq.gz | TFCN:25-0-2E333-9 |

### B3. Library is a different organism from the strain (5)

| Library | FileBase | Labelled strain | What the library is |
|---|---|---|---|
| 9003_005 | 9003_005_S33_R[12]_001.fastq.gz | DBVPG:3445 | pure *R. muc.* haploid ≈ TFCN_98C-7 |
| 9003_007 | 9003_007_S35_R[12]_001.fastq.gz | DBVPG:3538 | *R. aff. mucilaginosa* ≈ TFCN_25-395P-1 |
| 9003_009 | 9003_009_S37_R[12]_001.fastq.gz | DBVPG:4379 | hybrid diploid |
| 9003_010 | 9003_010_S38_R[12]_001.fastq.gz | DBVPG:4952 | pure haploid near *R. frigidialcoholis* |
| 9003_012 | 9003_012_S40_R[12]_001.fastq.gz | DBVPG:6649 | pure *R. muc.* haploid ≈ TFCN_17-325D-4 |

Proposal: add the rows to the labelled strain with `Notes` = "LIBRARY MISMATCH: not this strain (library identity 2026-10-03); do not use". See Q2.

### B4. Strain is in the DB, species is not *R. mucilaginosa* (6)

DBVPG:3380 (9003_004, *R. glutinis*), DBVPG:3985 (9003_008, *R. glutinis*), TFCN:152C-2 (9003_014, *R. toruloides*), TFCN:7-9-2 (9003_024, *R. diobovata*), TFCN:7-9-3 (9003_025, *R. diobovata*), EXF:13474 (9003_033, *R. toruloides*). Add one run row each. These libraries were never aligned in this project, so I have no check of their identity.

### B5. Strain not in the DB (15) — needs new strain rows

| Library | Strain | Genotype result in this project |
|---|---|---|
| 9003_017 | TFCN_1A-1-5 | R. mucilaginosa |
| 9003_022 | TFCN_2M-1-3 | R. mucilaginosa |
| 9003_027 | TFCN_86C-3 | R. mucilaginosa |
| 9003_015 | TFCN_17-333M-1 | candidate, no metadata |
| 9003_016 | TFCN_186CL-2 (SeqCoast sheet: `186CL④-2`) | candidate, no metadata |
| 9003_029 | EXF_10854 | candidate, no metadata |
| 9003_034 | EXF_17335 | candidate, no metadata |
| 9003_013 | TFCN_152A-3 | dropped by QC (species not set) |
| 9003_018 | TFCN_211C-2 | dropped by QC |
| 9003_023 | TFCN_7-6-3 | dropped by QC |
| 9003_026 | TFCN_86A-5 | dropped by QC |
| 9003_030 | EXF_10630 | dropped by QC |
| 9003_031 | EXF_14606 | dropped by QC |
| 9003_032 | EXF_13500 | dropped by QC |
| 9003_035 | EXF_12265 | dropped by QC |

I have no isolation metadata (source, locality, date) for these 15. The master list needs `Species`; I have a genotype species for 3 only. See Q3.

## Fix C — BY120-C1 / BY120-C7 as aliases of TFCN:25-0-2E333-9

No current input makes these aliases. "BY120-x" is a plate position: `import/strain_summary_YPD2.csv` has 16 rows with BY120 positions for different strains. Its row 116 reads "BY120-C1 (BY120-C7) 25-0-2E333-9".

Two ways:
- C-i (minimal): the B2 run rows carry the link. Put "SeqCoast label TFCN_BY120-C1 / -C7 = plate position" in their `Notes`. No alias row.
- C-ii: add a small `import/manual_aliases.csv` (`alias,strain_id,note`) and a few lines in `build_db` to load it with `source='manual_correction'`, `confidence='manual'`. This also serves Fix D.

I recommend C-ii. A plate position is not a strain name, so the alias needs a note that says that.

## Fix D — CCFEE 5036 = DBVPG 5227

Current DB state:
- Two separate strain rows: `CCFEE:5036` and `DBVPG:5227`.
- Alias "CCFEE 5036" points to both: from the DBVPG catalog to DBVPG:5227, and from BioSample to CCFEE:5036.
- Runs on CCFEE:5036: SRR5223778 (PRJNA342238, SAMN06285187) and ISQ0033050_S25, both noted "Too Low".
- Runs on DBVPG:5227: 23263Sta_DBVPG5227_S109 L001/L002.

Proposal:
- Keep `DBVPG:5227` as the strain_id. It has the phenotype row, and the merged CRAM uses this name.
- In the master list, change the two CCFEE_5036 rows to `ID`/`Strain` = DBVPG_5227. Put "deposited as CCFEE 5036; same isolate" in `Notes`.
- `CCFEE:5036` then is no longer a strain row. Its name stays an alias through the DBVPG catalog and BioSample records.
- Check: SAMN06285187 then matches DBVPG:5227 through the alias.

See Q4. The CCFEE rows carry the isolation data (collection date 1996-10-08, Edmonson Point, Antarctica) and the public SRA run. Moving the rows keeps that data.

## Questions for you

1. **Name for *R. aff. mucilaginosa*.** The DB has no such species. Its informal names use tokens such as `sp_clade_I`. Options: `aff_mucilaginosa`, or a clade token.
2. **Mismatched 9003 libraries (B3).** Add them with a warning note (proposed), or leave them out of the DB?
3. **15 new strains (B5).** Do you have a source for species and isolation data, e.g. the SeqCoast manifest xlsx or the Google Sheet? Or add them with species left blank and `needs_review`?
4. **CCFEE 5036 / DBVPG 5227.** Use DBVPG:5227 as the primary ID (proposed), or CCFEE:5036?
5. **Manual alias file (C-ii).** OK to add `import/manual_aliases.csv` and the loader in `build_db`?
6. **TFCN_363-1-2.** `popgen_strain_metadata.csv` says *R. taiwanensis*. The DB says *R. mucilaginosa*. Our callset dropped it in QC, so I have no genotype species. Leave it?

## Decisions (user, 2026-10-04)

1. *R. frigidialcoholis*: use "Rhodotorula frigidialcoholis" (user: "probably okay"). The 8 *R. aff. mucilaginosa* strains stay *Rhodotorula mucilaginosa* in the DB for now (user, 2026-10-04). Only the 12 *R. frigidialcoholis* strains change species.
2. Mismatched 9003 libraries (B3): add them, with a warning in `Notes`.
3. 15 new strains (B5): the user fills in `docs/straindb_new_strains_template_2026-10-04.csv`.
4. CCFEE 5036 / DBVPG 5227: `DBVPG:5227` is primary. Add "CCFEE 5036" / "CCFEE_5036" as aliases of DBVPG:5227 through the manual alias file.
5. Manual alias file (C-ii): approved.
6. TFCN_363-1-2: the user believes it is *R. mucilaginosa*. The read data agree (below). Keep *R. mucilaginosa* in the DB. Correct `popgen_strain_metadata.csv`, which says *R. taiwanensis*.

### TFCN_363-1-2 read evidence

From `results/reports/samtools/<strain>/<strain>.md.cram.stats` (reads mapped to the DH4148 reference):

All 319 strains, grouped by genotype species and called ploidy (haploids shown; TFCN_363-1-2 is called haploid):

| Group | n | Reads mapped | Mismatch rate |
|---|---|---|---|
| **TFCN_363-1-2** | 1 | **0.964** | **0.0056** |
| R. mucilaginosa haploids | 242 | min 0.933, median 0.972 | median 0.0053, p90 0.0125 |
| R. frigidialcoholis haploids | 10 | min 0.765, median 0.832 | median 0.0663 |
| R. aff. mucilaginosa haploids | 7 | min 0.771, median 0.793 | median 0.0676 |

In the 2026-10-03 callset, its non-reference fraction at called sites is 0.013. The *R. mucilaginosa* median is 0.013, and the *R. frigidialcoholis* / *R. aff. mucilaginosa* range is 0.26–0.34. This project has no *R. taiwanensis* sample to compare against. The data show that TFCN_363-1-2 is as close to the reference as the *R. mucilaginosa* strains are, and much closer than the two sister species are.

Why it was dropped: `low_coverage`, with 55.7% of called sites at DP < 5. Its mean depth is 18.5×, against about 70× for the other strains. But mosdepth reports 98% of the genome at ≥ 5×. These two numbers disagree, and I have not found the cause. Check this strain in the step-1 QC of the new callset.

## B5 read evidence (2026-10-04)

Reads mapped and mismatch rate come from `results/reports/samtools/<strain>/<strain>.md.cram.stats`. The non-ref fraction and the QC reasons come from `results/variant_qc/strain_qc.tsv`.

| Strain | Reads mapped | Mismatch | Non-ref frac | Closest profile |
|---|---|---|---|---|
| TFCN_1A-1-5 | 0.959 | 0.0019 | 0.000 | R. mucilaginosa haploid |
| TFCN_2M-1-3 | 0.959 | 0.0049 | 0.014 | R. mucilaginosa haploid |
| TFCN_86C-3 | 0.970 | 0.0020 | 0.000 | R. mucilaginosa haploid |
| TFCN_17-333M-1 | 0.883 | 0.0324 | 0.271 | R. mucilaginosa hybrid diploid (median 0.871 / 0.0323) |
| TFCN_186CL-2 | 0.877 | 0.0342 | 0.284 | R. mucilaginosa hybrid diploid |
| EXF_10854 | 0.833 | 0.0686 | 0.258 | R. frigidialcoholis or R. aff. mucilaginosa (these numbers cannot separate the two) |
| EXF_17335 | 0.831 | 0.0685 | 0.261 | same as EXF_10854 |
| TFCN_152A-3 | 0.093 | 0.0681 | — | outside the R. mucilaginosa complex |
| TFCN_211C-2 | 0.256 | 0.0622 | — | outside the complex |
| TFCN_7-6-3 | 0.231 | 0.0633 | — | outside the complex (user: R. diobovata) |
| TFCN_86A-5 | 0.286 | 0.0798 | — | outside the complex |
| EXF_10630 | 0.271 | 0.0647 | — | outside the complex |
| EXF_14606 | 0.252 | 0.0630 | — | outside the complex |
| EXF_13500 | 0.441 | 0.0745 | — | outside the complex |
| EXF_12265 | 0.258 | 0.0593 | — | outside the complex |

For the strains outside the complex, the non-ref fraction is not useful. Most of their sites have DP < 5.

User-supplied for TFCN_7-6-3 (2026-10-04): species *Rhodotorula diobovata*, alternative ID "7-6-2003", location "China: Zhoushan City, Zhejiang Province".

Note: "7-6-2003" has the same form as the master-list IDs "7-9-2002" and "7-9-2003" for TFCN_7-9-2 and TFCN_7-9-3. This looks like a spreadsheet converting "7-6-3" to a date. Keep it as an alias so lookups still match, but mark it as a spreadsheet artifact.

`/rhome/jstajich/Strain Metadata - metadata about strains.csv` matches `import/popgen_strain_metadata.csv`:
- Same 20 columns and the same 536 strains.
- It has an extra header row and a blank row.
- 11 TFCN IDs have a leading zero, e.g. `0209-6-1`.
- It holds the CCFEE_5036 rows in a different order.
It holds none of the 15 B5 strains.

## EXF strains: catalogue search and species evidence (2026-10-04)

### Public EXF catalogue

The catalogue is at `https://catalogue.ex-genebank.com/bt_exf_index.html`, linked from `https://www.ex-genebank.com/index.php/en/fungi-2`. It has 699 species pages.
- I downloaded all 699 pages on 2026-10-04.
- They list 6,938 EXF strain IDs. The highest is **EXF-11119**.
- **None of the 6 strains is listed:** EXF-10630, EXF-10854, EXF-12265, EXF-13500, EXF-14606, EXF-17335.
- Five of them have higher numbers than any listed strain, so the public catalogue probably has not added them yet.
- EXF-10630 is inside the listed range, but it is not listed either.
- A web search for the six IDs found nothing.
- Isolation data for these strains must come from the EXF collection (Ljubljana) or from lab records.

### EXF-10630 → *Rhodotorula* sp_clade_I (user approved 2026-10-04)

Source: `ExtremeRhodotorula_DraftGenomes/EXF_external/EXF-10630.ani.txt` (fastANI, 2023-10-26) and `EXF-10630.CBS7608.txt`.
- **98.3–98.7%** ANI to the 13 strains that `Rhodotorula_Rodeo/genomes.csv` labels *R. aff. babjevae*. Examples: EXF_7045, TFCN_222A-2, TFCN_17-328C-3, TFCN_86C-7, DBVPG_4629, DBVPG_8058, EXF_14539, EXF_6481.
- 90.3% to *R. babjevae* CBS 7808.
- 88.8–89.2% to *R. graminis* (NRRL Y-2474, DBVPG_6083, DBVPG_7021, EXF_13753) and to *R. glutinis* DBVPG_6081.
- In the StrainDB, the *R. aff. babjevae* clade is named `Rhodotorula sp_clade_I` (DBVPG_4629, DBVPG_8058). EXF_10630 gets that name.

### Other EXF strains

- **EXF_10854 and EXF_17335.** Reads mapped 83%, mismatch 6.9%. This matches *R. frigidialcoholis* or *R. aff. mucilaginosa*, which the mapping numbers cannot tell apart. The sourmash database will be used to separate them.
- **EXF_12265, EXF_13500 and EXF_14606.** Reads mapped 25–44%. These strains are outside the *R. mucilaginosa* complex, and their species is not known.

## Species-ID with sourmash (2026-10-04)

The full report is `/bigdata/stajichlab/shared/projects/Rhodotorula/Species_ID_db/REPORT_2026-10-04.md`. The species calls for the 15 new strains are in `docs/straindb_new_strains_filled_2026-10-04.csv`, column `sourmash_call_2026-10-04`. Findings that change this proposal:
- TFCN_363-1-2 is confirmed as *R. mucilaginosa*. It shares almost no k-mers with *R. taiwanensis*.
- The 9003 libraries of DBVPG_3380 and DBVPG_3985 (both *R. glutinis* in the DB) and TFCN_152C-2 (*R. toruloides* in the DB) are *R. mucilaginosa*. Flag them in B4 like the B3 libraries, and add them to the lab tube check.
- The 9003 library labelled DBVPG_4952 (9003_010) matches both *R. frigidialcoholis* and *R. mucilaginosa*. It may be a hybrid or a mixture.
- Hybrids: TFCN_17-333M-1 and TFCN_186CL-2 are *R. mucilaginosa* × *R. aff. mucilaginosa*. EXF_14606 is likely *R. glutinis* × *R. aff. babjevae*.

## Decisions (user, 2026-10-04, second round)

7. **Hybrid species names.** The species name says the strain is a hybrid and names at least one parent. Name the second parent if it is known. Format: `Rhodotorula <parent 1> x Rhodotorula <parent 2> hybrid`. Applied to TFCN_17-333M-1, TFCN_186CL-2 and EXF_14606 (EXF_14606 is unconfirmed; it could be a mixed culture).
   - Open: should the 43 hybrid diploids already in the DB, now listed as *R. mucilaginosa*, be renamed the same way? Their second parent is known for some of them only.
8. **EXF_10854** is `Rhodotorula sp_clade_XIV`, a new informal token. The existing tokens are I, XI, XII and XIII.
9. **Rodeo `genomes.csv` relabels** are done. See `Rhodotorula_Rodeo/NOTE_species_relabel_2026-10-04.md`.
