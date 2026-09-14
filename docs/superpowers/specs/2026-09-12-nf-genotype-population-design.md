# nf_genotype_population: Design Spec

Date: 2026-09-12
Status: Approved for implementation planning
Author: Jason Stajich (design facilitated with Claude Code)

## Context

This project (`Rhodotorula_mucilaginosa_DH4148_ref`) runs nf-core/sarek 3.9.0
for mapping, QC, and per-sample germline variant calling against the custom
*Rhodotorula mucilaginosa* DH4148 reference (GCA_058775505.1). Sarek covers
alignment/dedup/CRAM well, but the population-genomics-specific logic this
study needs - per-strain ploidy-aware calling, population-sliced joint
genotyping, hard-filtering without a truth/known-sites resource, coverage/CNV
visualization - either isn't supported by Sarek at all, or needs different
defaults than Sarek ships.

A prior pipeline, `Rmuc_popgen_NRRLY2510` (SLURM sbatch scripts + environment
modules, ~20 scripts), implements this logic already but is not
containerized/reproducible and is hard to extend. `nf_genotype_population` is
a new, standalone Nextflow DSL2 pipeline that ports that logic into a
container-based, SLURM-native pipeline, using Sarek only for the front-end
alignment stage.

### Findings that shaped this design (verified this session)

- **Sarek has zero per-sample ploidy support** for germline HaplotypeCaller
  anywhere in its codebase (confirmed by a full-repo search of the pinned
  `3.9.0` tag). It always calls at GATK's default diploid ploidy. Its only
  `ext.args` for `GATK4_HAPLOTYPECALLER` is `-ERC GVCF` (+ optional
  `--pcr-indel-model`).
- **Sarek's only germline joint-filtering path is VQSR**
  (`VariantRecalibrator`/`ApplyVQSR`), which requires a known-sites
  truth/training resource to run at all. Rhodotorula has no dbSNP or
  gold-standard truth set, so VQSR is not viable. Sarek ships **no hard-filter
  (`VariantFiltration`) fallback** for germline joint calls - there is nothing
  to mirror here; `nf_genotype_population` must implement its own hard-filter
  step using GATK's published best-practices thresholds.
- Sarek's `GenomicsDBImport` args (`--genomicsdb-shared-posixfs-optimizations
  true --bypass-feature-reader`) are worth reusing as-is; `GenotypeGVCFs` uses
  plain GATK defaults with nothing else to copy.
- A metadata/samplesheet reconciliation (`scripts/reconcile_metadata.py`,
  this project) found: 278 of 319 sequenced strains are confirmed
  `R. mucilaginosa`, 8 are borderline `R. aff. mucilaginosa`, 12 are clearly
  other species, 21 have no metadata row, and metadata.txt's `ploidy` column
  is blank for most rows - hence this pipeline computes ploidy itself rather
  than trusting that column.
- `nQuack` (Gaynor et al. 2024) is a real, more recent alternative to nQuire
  for ploidy inference, but ships with **no Bioconda/Docker/Singularity
  container** - it is an R package installed via
  `devtools::install_github`, requiring `samtools` on PATH. A container has
  to be built in-house.

## Scope boundary with Sarek

- **Sarek** (this project) is re-configured to stop after alignment/dedup and
  CRAM output - no variant calling (`--tools ""` / `--step markduplicates`
  scope). Its `results/preprocessing/*/*.cram(.crai)` is the sole input
  contract `nf_genotype_population` depends on.
- **`nf_genotype_population`** is a separate pipeline repo. It owns
  everything downstream: ploidy inference, per-strain HaplotypeCaller,
  population-sliced joint genotyping, hard-filtering/annotation, coverage/CNV.
- This boundary isolates Sarek version upgrades from the study-specific logic;
  Sarek's CRAM output shape has been stable for years.

## Components

### 1. Ploidy inference (`PLOIDY_INFERENCE` subworkflow)

- Pluggable via `params.ploidy_method`: `custom_het_script` (Phase 1 default),
  `nquire`, `nquack` (both added Phase 2). Each method module emits a common
  contract: `[strain, inferred_ploidy, confidence_metric]`.
- `custom_het_script`: in-house het-rate/coverage-uniformity check
  (bcftools/mosdepth-based) - fast, transparent, no exotic dependency.
- `nquire`: existing bioconda/container, low integration effort.
- `nquack`: in-house container (R base + `samtools` + `devtools::install_github
  ("mgaynor1/nQuack")`) - one-time build cost, same output contract.
- **Cross-check stage**: joins inferred ploidy against `metadata.txt`'s
  existing `ploidy` column where present. Disagreements or blanks are written
  to `ploidy_review.csv` for manual sign-off. A human-edited
  `ploidy_overrides.csv` (strain -> final ploidy) is what the pipeline actually
  reads for calling - inference informs it, never silently overrides review.

### 2. Per-strain HaplotypeCaller + joint genotyping

- `GATK_HAPLOTYPECALLER`: `-ERC GVCF`, per-strain `--sample-ploidy` from
  `ploidy_overrides.csv` (1 for haploid, 2 for diploid/hybrid).
  Computationally-derived `_hap1`/`_hap2` rows are excluded from calling
  entirely - they are phased outputs of a diploid call, not independent
  calling inputs.
- `population_sets.yaml` (ported as-is from `Rmuc_popgen_NRRLY2510`, current
  format unchanged) drives `GenomicsDBImport`
  (`--genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader`)
  + `GenotypeGVCFs`, sliced per population group. **A built-in `all` group
  containing every included strain always runs alongside any sliced
  sub-population groups**, so a single all-strains callset is guaranteed to
  exist regardless of population-file completeness. Strains not yet assigned
  to a population fall into a `default_unassigned` group.
- Future (Phase 3) `PROPOSE_POPULATIONS` subworkflow: PCA/DAPC
  (`adegenet`/`poppr` in R) on the joint-genotyped `all` VCF, run at a fixed
  cadence, emitting `population_sets.proposed.yaml` as a diff against the
  current manual file for human review/merge. Never auto-overwrites the
  manual file. (Metadata Origin/Environment has shown little geographic/
  ecological structure so far, hence deriving populations empirically rather
  than from metadata.)

### 3. Filtering and annotation

- `GATK4_VARIANTFILTRATION` hard-filter (no VQSR - not viable without a
  truth/known-sites resource). Starting thresholds (GATK best-practices,
  since Sarek offers nothing better to adopt):
  - SNPs: `QD<2.0, FS>60.0, MQ<40.0, MQRankSum<-12.5, ReadPosRankSum<-8.0, SOR>3.0`
  - Indels: `QD<2.0, FS>200.0, ReadPosRankSum<-20.0, SOR>10.0`
- SnpEff annotation against the DH4148 gene models, per population group then
  merged.

### 4. Coverage / CNV interface

- `mosdepth` runs per-strain directly on Sarek's CRAM, windowed (5/10/50kb,
  matching `Rmuc_popgen_NRRLY2510`'s scheme). Feeds a ported R plotting
  script for per-strain coverage-ratio visualization.
- Interface only, in Phase 1/2: mosdepth's per-window output and the CRAMs are
  already in the shape a later CNV caller (GATK gCNV, or a simpler
  mosdepth-ratio approach) needs - no rework required when Phase 4 starts.
  Indel calling itself is already covered by GATK4's default SNP+indel joint
  calling in step 2/3 above.

### 5. Storage and reproducibility

- CRAM-only end-to-end; `nf_genotype_population` never materializes BAM.
- Per-strain GVCFs kept (compressed/indexed `.g.vcf.gz`) for re-genotyping as
  new strains are added; GenomicsDB workspaces are disposable `work/` scratch,
  not archived outputs.
- Every tool pinned to an exact container digest in `conf/modules.config`
  (mirroring Sarek's own practice). A `-profile hpcc` mirrors this project's
  SLURM/Singularity/`$SCRATCH` conventions, including the containerOptions
  quoting lesson learned this session (Groovy closures referencing
  `System.getenv()` resolve in the Nextflow *head* job's environment, not each
  per-task SLURM sub-job's environment - environment-variable references
  destined for per-task shell expansion must stay as literal, single-quoted
  strings, e.g. `containerOptions = '-B $SCRATCH:/tmp'`, the same execution-
  context class of bug as never resolving paths via `BASH_SOURCE[0]` in
  scripts destined for SLURM).

## Phased delivery plan

- **Phase 0** (done): Sarek FASTQC `/tmp` fix (JVM temp-dir env vars +
  per-task `$SCRATCH:/tmp` container bind); metadata/samplesheet
  reconciliation table and report.
- **Phase 1**: Scaffold `nf_genotype_population` repo. Sarek stopped at CRAM.
  Ploidy inference using only `custom_het_script`. Per-strain HaplotypeCaller
  with ploidy branching. `all`-group joint genotyping. Hard-filter + SnpEff.
  This is the minimum end-to-end path producing a usable SNP/indel VCF.
- **Phase 2**: Add `nquire` and `nquack` as additional pluggable
  ploidy-inference methods, cross-validated against each other, the custom
  script, and metadata. Wire in `population_sets.yaml` sliced joint genotyping
  alongside the `all` group.
- **Phase 3**: mosdepth + CNV visualization. `PROPOSE_POPULATIONS`
  (PCA/DAPC) subworkflow.
- **Phase 4**: Full CNV calling (GATK gCNV or mosdepth-ratio) as a
  first-class module.

## Out of scope for this spec

- The actual CNV-calling algorithm choice (Phase 4) - deferred until Phase 1-3
  are working and real coverage data is in hand to evaluate against.
- Automated population structure derivation beyond proposing a diff for human
  review (Phase 3) - DAPC/PCA parameters and validation are their own design
  question when that phase starts.
- Rhotodoula_phenotypes integration (copper-stress phenotyping data) - noted
  as a possible future phenotype-genotype join by strain ID, not part of this
  pipeline's scope.
