# Decision: no depth mask on hom-ref genotypes; IQ-TREE keeps identical strains (2026-10-05)

**Status:** approved by the user on 2026-10-05. Applied in the pipeline on branch `fix/homref-dp-mask`, commit `64f189c` (`/rhome/jstajich/projects/nf/nf_genotype_population`).

## Decision

1. **VARIANT_QC_FILTER depth mask.**

   | Rule | Before | Now |
   |---|---|---|
   | GQ mask | `GQ < 20` → missing, all genotypes | unchanged |
   | Depth mask, genotypes with an ALT allele | `DP < 5` → missing | unchanged (`qc_min_dp = 5`) |
   | Depth mask, hom-ref genotypes | `DP < 5` → missing | none: `qc_min_dp_homref = 0` (off) |

   Set `--qc_min_dp_homref 5` to get the old rule back.

2. **IQTREE:** add `-keep-ident` to both IQ-TREE calls.

## Why: hom-ref depth mask

GenotypeGVCFs takes the `FORMAT/DP` of a hom-ref genotype from the gVCF reference block. That value is the block `MIN_DP`, not the depth at the site. A long reference block can have a low `MIN_DP` even when the strain has good coverage.

Evidence from run 29409304 (`rmuc_core`, 247 strains):

- **Tree alignments shrank.**
  - `rmuc_core` fell from 17,626 to 702 sites after the 6 rescued `low_coverage` strains were added.
  - `rmuc_core_outgroup` fell from 19,157 to 646 sites; `rmuc_with_hybrids` from 16,437 to 692.
  - The tree alignment keeps only sites with no missing genotype (`tree_max_missing = 0`).
- **Missing-genotype fraction in `rmuc_core.snps.maf`, after masking:**

  | Strain | Missing |
  |---|---|
  | TFCN_137D-4 | 0.85 |
  | TFCN_25-337M-3 | 0.81 |
  | TFCN_25-0-2E332-1 | 0.53 |
  | TFCN_363-1-2 | 0.52 |
  | TFCN_89D-3 | 0.52 |
  | TFCN_17-333D-2 | 0.48 |
  | next highest (TFCN_137C-3) | 0.40 |

- **The masked genotypes are hom-ref calls.** Counts are from the hard-filtered `all` VCF, before masking, on the largest chromosome, at the `rmuc_core` SNP sites:

  | Strain | Hom-ref calls | Hom-ref with DP<5 | ALT calls | ALT with DP<5 | Missing before mask |
  |---|---|---|---|---|---|
  | TFCN_137D-4 | 28,757 | 0.89 | 3 | 0.00 | 0 |
  | TFCN_363-1-2 | 20,047 | 0.46 | 8,713 | 0.01 | 0 |
  | TFCN_17-333D-2 | 28,745 | 0.40 | 15 | 0.00 | 0 |
  | DBVPG_10619 (control) | 20,074 | 0.00 | 8,686 | 0.00 | 0 |

- **These strains have good coverage.** mosdepth gives 16.5–27.7×, with 98–99 % of the genome at ≥ 5×; see `docs/strain_qc_2026-10-05.md`.
- **What still screens hom-ref calls:** `GQ < 20`, the allele-fraction mask (`qc_hap_min_af`) and the site-level filters.

## Why: `-keep-ident`

- `GTR+ASC` stops when a column is invariant. A column that varies only by IUPAC het codes counts as invariant.
- When that happens, the IQTREE task reruns on `<pop>.varsites.phy`.
- Without `-keep-ident`, IQ-TREE collapses identical sequences before it writes that file. The rerun tree then loses those strains:
  - the `rmuc_with_hybrids` tree had 74 of 280 strains in run 29409304
  - it had 96 of 267 strains in the run of 2026-10-04
- Test on the 281-sequence `rmuc_with_hybrids` alignment with `-keep-ident`: `varsites.phy` had 281 sequences × 398 sites, and the rerun built a tree from all 281.

## Effects and open points

- **All group VCFs change.** VARIANT_QC_FILTER, SnpEff, SNP alignment and IQ-TREE rerun for every group. The `all` callset (HaplotypeCaller → HARDFILTER) is reused.
- **Expected effect:** fewer missing genotypes and more tree sites. I have not measured the new numbers yet. Check them after the rerun.
- **Not changed:** the strain-QC `low_coverage` rule (`strain_qc.tsv`) uses the same VCF DP and has the same artifact. The 6 strains stay rescued in `variant_qc_reviewed_keep.tsv`. A change to that rule is still open.
- **New, untested observation:** TFCN_137D-4 (3 ALT calls) and TFCN_17-333D-2 (15 ALT calls) are almost identical to the reference DH4148 on the chromosome tested. They may be clones of DH4148, or mislabelled reference DNA. Check this before they are used as independent isolates. TFCN_25-337M-3 was not checked.
