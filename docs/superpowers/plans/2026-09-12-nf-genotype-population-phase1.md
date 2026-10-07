# nf_genotype_population Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the minimum end-to-end `nf_genotype_population` pipeline: per-strain ploidy inference (custom heuristic only), ploidy-aware GATK4 HaplotypeCaller, a single all-strains joint-genotyping call, GATK4 hard-filtering, and SnpEff annotation — producing a usable, ploidy-correct SNP/indel VCF for the confirmed *R. mucilaginosa* strain set.

**Architecture:** A new standalone Nextflow DSL2 pipeline repo, consuming CRAM files from a Sarek run that has been reconfigured to stop after alignment/dedup (no variant calling). Every stage after that boundary — ploidy inference, HaplotypeCaller, joint genotyping, filtering, annotation — is implemented here, following nf-core module/subworkflow conventions but without full nf-core-template boilerplate, so it stays small and specific to this study.

**Tech Stack:** Nextflow DSL2, Singularity/Apptainer, GATK4, bcftools, SnpEff, Python 3 (bin/ scripts), pytest, SLURM (UCR HPCC).

**Spec:** `docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md` (in `Rhodotorula_mucilaginosa_DH4148_ref`)

## Global Constraints

- New pipeline repo lives at `/rhome/jstajich/projects/nf/nf_genotype_population` during development; remote is `stajichlab/nf_genotype_population` on GitHub (adjust if the user has already created it under a different name/location — confirm before Task 1's `git remote add`).
- Target environment: UCR HPCC, SLURM executor, Singularity containers. Every container-referencing `containerOptions`/`beforeScript` value that must be expanded per-task (not at pipeline-definition time) MUST be a single-quoted literal string, never a Groovy closure calling `System.getenv()` — that resolves in the Nextflow head job's environment, not each task's own SLURM sub-job environment (verified failure mode this session, see spec section 5).
- No dbSNP/known-sites/gold-standard truth set exists for *R. mucilaginosa* — VQSR is never used; hard-filtering with the exact thresholds below is the only filtering path.
- SNP hard-filter thresholds (from spec, GATK best-practices defaults): `QD<2.0, FS>60.0, MQ<40.0, MQRankSum<-12.5, ReadPosRankSum<-8.0, SOR>3.0`.
- Indel hard-filter thresholds (from spec): `QD<2.0, FS>200.0, ReadPosRankSum<-20.0, SOR>10.0`.
- `_hap1`/`_hap2` metadata rows are computationally-phased outputs of a diploid call, never independent HaplotypeCaller inputs — they must never appear as strains going into calling.
- `population_sets.yaml` always includes a built-in `all` group containing every included strain, in addition to any other groups (Phase 2 concern — Phase 1 only implements `all`).
- Reference genome: `GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna` (same file already in `Rhodotorula_mucilaginosa_DH4148_ref/refgenome/`).

---

## Task 1: Reconfigure the existing Sarek project to stop at CRAM

The currently-running Sarek job still has `tools='haplotypecaller', joint_germline=true`, which will attempt VQSR later and fail (no known-sites resource — confirmed this session by reading Sarek 3.9.0's source). Stop wasting compute on that path now, since `nf_genotype_population` re-does HaplotypeCaller itself anyway.

**Files:**
- Modify: `/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/conf/profile_sarek.config`

**Interfaces:**
- Produces: `results/preprocessing/*/*.cram(.crai)` under the existing `outdir`, which Task 3 of this plan (and later, `nf_genotype_population`'s real input) depends on.

- [ ] **Step 1: Edit `conf/profile_sarek.config`**

Change:
```groovy
    tools          = 'haplotypecaller'
    joint_germline = true
```
to:
```groovy
    tools          = ''
    joint_germline = false
```
and update the comment above it (currently incorrectly claims Sarek "falls back to hard-filtering" — it does not; delete that claim):
```groovy
    // No variant calling in Sarek at all: this project's germline calling is
    // ploidy-aware per-strain (nf_genotype_population), which Sarek cannot
    // do. Sarek is used only for QC/alignment/dedup through CRAM output.
    // (Sarek's only germline joint-filtering path is VQSR, which requires a
    // known-sites truth resource we don't have for this organism - verified
    // by reading its 3.9.0 source; there is no hard-filter fallback in Sarek.)
    tools          = ''
    joint_germline = false
```

- [ ] **Step 2: Validate config**

Run: `cd /bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref && nextflow config -profile sarek . | grep -E 'tools|joint_germline'`
Expected: `tools = ''` and `joint_germline = false`.

- [ ] **Step 3: Resubmit and confirm the run reaches CRAM output without attempting variant calling**

Run: `sbatch run_sarek.sh` (uses `-resume`, so alignment/dedup work already done is reused).
Check after it progresses: `grep -i "HAPLOTYPECALLER\|VARIANTRECALIBRATOR" logs/nf-sarek.<jobid>.log` should return nothing.

- [ ] **Step 4: Commit**

```bash
cd /bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref
git add conf/profile_sarek.config
git commit -m "sarek: stop at CRAM, no variant calling (VQSR unusable without known-sites)"
```

---

## Task 2: Scaffold the nf_genotype_population repo

**Files:**
- Create: `/rhome/jstajich/projects/nf/nf_genotype_population/.gitignore`
- Create: `/rhome/jstajich/projects/nf/nf_genotype_population/README.md`
- Create: `/rhome/jstajich/projects/nf/nf_genotype_population/nextflow.config`
- Create: `/rhome/jstajich/projects/nf/nf_genotype_population/main.nf`

**Interfaces:**
- Produces: repo root that every later task's files are relative to.

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p /rhome/jstajich/projects/nf/nf_genotype_population/{bin,modules/local,subworkflows/local,workflows,conf,assets,tests/fixtures,docs}
cd /rhome/jstajich/projects/nf/nf_genotype_population
git init
```

- [ ] **Step 2: Write `.gitignore`**

```
work/
.nextflow*
results/
.DS_Store
__pycache__/
*.pyc
.pytest_cache/
```

- [ ] **Step 3: Write a placeholder `main.nf`**

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow {
    log.info "nf_genotype_population - scaffold only, see Task 8 of the Phase 1 plan for the real entry workflow"
}
```

- [ ] **Step 4: Write a minimal `nextflow.config`**

```groovy
manifest {
    name        = 'nf_genotype_population'
    description = 'Ploidy-aware per-strain calling, joint genotyping, and filtering for Rhodotorula population genomics'
}

executor {
    name            = 'slurm'
    queueSize       = 100
    submitRateLimit = '10/1min'
}

singularity {
    enabled      = true
    autoMounts   = true
    cacheDir     = '/bigdata/stajichlab/shared/lib/singularity_cache'
    envWhitelist = 'SCRATCH'
    runOptions   = '--home "$PWD"'
}

process {
    executor       = 'slurm'
    cache          = 'lenient'
    shell          = ['/bin/bash', '-l']
    clusterOptions = '-N 1 -n 1'
    beforeScript   = 'module load singularity; export TMPDIR="${SCRATCH:?}"; export SINGULARITYENV_TMPDIR="${SCRATCH:?}"; export SINGULARITYENV__JAVA_OPTIONS="-Djava.io.tmpdir=${SCRATCH:?}"; export SINGULARITYENV_JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=${SCRATCH:?}"'

    withLabel: 'process_low' {
        queue  = 'short'
        cpus   = { 2 * task.attempt }
        memory = { 8.GB * task.attempt }
        time   = '2h'
    }
    withLabel: 'process_medium' {
        queue  = 'epyc'
        cpus   = { 6 * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = '24h'
    }
}

profiles {
    hpcc {
        includeConfig 'conf/profile_hpcc.config'
    }
}
```

- [ ] **Step 5: Write a stub README**

```markdown
# nf_genotype_population

Ploidy-aware per-strain GATK4 calling, joint genotyping, and hard-filtering
for Rhodotorula population genomics on UCR HPCC. Consumes CRAM output from
an nf-core/sarek run (Sarek itself only used through alignment/dedup - see
this project's design spec).

Design spec: see `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`.
```

- [ ] **Step 6: Verify the scaffold runs**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && module load nextflow/26.04.6 2>&1 | tail -1 && nextflow run main.nf`
Expected: prints the scaffold log line, exits 0.

- [ ] **Step 7: Commit**

```bash
cd /rhome/jstajich/projects/nf/nf_genotype_population
git add .gitignore README.md nextflow.config main.nf
git commit -m "scaffold nf_genotype_population repo"
git remote add origin git@github.com:stajichlab/nf_genotype_population.git
```

(Do not push yet — hold until the user confirms the remote repo exists and the name is right.)

---

## Task 3: Build tiny synthetic test fixtures

Phase 1 needs a small, fast, deterministic dataset to prove each module and the final assembled pipeline actually work, without waiting on real 300-strain-scale HPCC runs. One small reference contig, one synthetic haploid CRAM, one synthetic diploid CRAM.

**Files:**
- Create: `tests/fixtures/make_fixtures.sh`
- Create (generated by the script, but commit the outputs since they're tiny): `tests/fixtures/ref.fa`, `tests/fixtures/ref.fa.fai`, `tests/fixtures/ref.dict`, `tests/fixtures/haploid_strain.cram(.crai)`, `tests/fixtures/diploid_strain.cram(.crai)`

**Interfaces:**
- Produces: `tests/fixtures/ref.fa` (+ `.fai`/`.dict`), `tests/fixtures/haploid_strain.cram`, `tests/fixtures/diploid_strain.cram` — consumed by Tasks 5, 7, 9, 11, 13's smoke tests.

- [ ] **Step 1: Write the fixture-generation script**

```bash
#!/usr/bin/env bash
# tests/fixtures/make_fixtures.sh
# Requires: samtools, bwa-mem2, dwgsim (all already used elsewhere in this
# project's containers - run this once locally/interactively, not via SLURM).
set -euo pipefail
cd "$(dirname "$0")"

# 2kb synthetic single-contig reference, deterministic (seeded).
python3 - <<'PY'
import random
random.seed(42)
seq = "".join(random.choice("ACGT") for _ in range(2000))
with open("ref.fa", "w") as fh:
    fh.write(">synth_contig1\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + "\n")
PY
samtools faidx ref.fa
gatk CreateSequenceDictionary -R ref.fa -O ref.dict 2>/dev/null || \
    picard CreateSequenceDictionary R=ref.fa O=ref.dict
bwa-mem2 index ref.fa

# Haploid strain: simulate reads from the reference itself (no variants injected
# beyond sequencer noise) so heterozygosity should be ~0.
dwgsim -N 2000 -1 100 -2 100 -r 0 -R 0 -e 0.001 -E 0.001 ref.fa haploid_strain
bwa-mem2 mem ref.fa haploid_strain.bwa.read1.fastq haploid_strain.bwa.read2.fastq \
    | samtools sort -O cram --reference ref.fa -o haploid_strain.cram -
samtools index haploid_strain.cram

# Diploid strain: inject ~1% heterozygous SNP rate (-r 0.01) so heterozygosity
# should be clearly nonzero.
dwgsim -N 2000 -1 100 -2 100 -r 0.01 -R 0 -e 0.001 -E 0.001 ref.fa diploid_strain
bwa-mem2 mem ref.fa diploid_strain.bwa.read1.fastq diploid_strain.bwa.read2.fastq \
    | samtools sort -O cram --reference ref.fa -o diploid_strain.cram -
samtools index diploid_strain.cram

rm -f *.fastq *.mutations.* *.bfast.* *.bwa.read*.fastq
echo "Fixtures written to $(pwd)"
```

- [ ] **Step 2: Run it and verify outputs exist**

Run: `bash tests/fixtures/make_fixtures.sh`
Expected: `ref.fa`, `ref.fa.fai`, `ref.dict`, `haploid_strain.cram(.crai)`, `diploid_strain.cram(.crai)` all present in `tests/fixtures/`.

- [ ] **Step 3: Commit**

```bash
git add tests/fixtures/make_fixtures.sh tests/fixtures/ref.fa tests/fixtures/ref.fa.fai tests/fixtures/ref.dict tests/fixtures/haploid_strain.cram tests/fixtures/haploid_strain.cram.crai tests/fixtures/diploid_strain.cram tests/fixtures/diploid_strain.cram.crai
git commit -m "add synthetic haploid/diploid test fixtures"
```

---

## Task 4: `custom_het_ploidy.py` core logic + unit tests

**Files:**
- Create: `bin/custom_het_ploidy.py`
- Test: `tests/test_custom_het_ploidy.py`

**Interfaces:**
- Produces: `het_fraction_from_vcf(vcf_lines: list[str]) -> float | None`, `call_ploidy(het_fraction: float | None, threshold: float = 0.01) -> str` (one of `"haploid"`, `"diploid"`, `"unknown"`), CLI `custom_het_ploidy.py --cram X --reference Y --strain S --out Z.csv` writing a `strain,inferred_ploidy,het_fraction,method` CSV row.
- Consumed by: Task 5 (Nextflow module wrapper), Task 6 (crosscheck reads its CSV output).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_custom_het_ploidy.py
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from custom_het_ploidy import het_fraction_from_vcf, call_ploidy

HAPLOID_LIKE_VCF = [
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0/0\n",
]

DIPLOID_LIKE_VCF = [
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\n",
]


def test_het_fraction_haploid_like_is_zero():
    assert het_fraction_from_vcf(HAPLOID_LIKE_VCF) == 0.0


def test_het_fraction_diploid_like_is_high():
    frac = het_fraction_from_vcf(DIPLOID_LIKE_VCF)
    assert abs(frac - (2 / 3)) < 1e-9


def test_het_fraction_none_when_no_sites():
    assert het_fraction_from_vcf(["#comment\n"]) is None


def test_het_fraction_skips_multiallelic_and_indels():
    lines = [
        "chr1\t100\t.\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n",   # multiallelic, skip
        "chr1\t200\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\n",     # indel, skip
        "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",       # counts
    ]
    assert het_fraction_from_vcf(lines) == 1.0


def test_call_ploidy_thresholds():
    assert call_ploidy(0.0) == "haploid"
    assert call_ploidy(0.005) == "haploid"
    assert call_ploidy(0.02) == "diploid"
    assert call_ploidy(None) == "unknown"
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_custom_het_ploidy.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'custom_het_ploidy'` (file doesn't exist yet).

- [ ] **Step 3: Write `bin/custom_het_ploidy.py`**

```python
#!/usr/bin/env python3
"""Estimate ploidy (haploid vs diploid) from CRAM heterozygous-site rate.

Calls variants at forced diploid genotype-likelihood mode via
`bcftools mpileup | bcftools call`, then computes the fraction of covered,
biallelic SNP sites called heterozygous. A true haploid genome shows
near-zero heterozygosity (mapping/base-calling noise only); a diploid or
hybrid genome shows a substantial het fraction. This is the fast,
transparent Phase 1 ploidy check - nQuire/nQuack are added in Phase 2 as
additional cross-validated methods, not replacements for this one.
"""
import argparse
import subprocess


def het_fraction_from_vcf(vcf_lines):
    """Fraction of biallelic-SNP records called heterozygous.

    vcf_lines: iterable of VCF data/header lines (header lines, i.e. those
    starting with '#', are ignored). Returns None if no qualifying sites
    were found (e.g. no coverage at all).
    """
    n_sites = 0
    n_het = 0
    for line in vcf_lines:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        ref, alt = fields[3], fields[4]
        if len(ref) != 1 or len(alt) != 1:
            continue  # biallelic SNPs only - skip indels and multiallelic sites
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        if "GT" not in format_keys:
            continue
        gt = sample_values[format_keys.index("GT")]
        alleles = gt.replace("|", "/").split("/")
        if len(alleles) != 2 or "." in alleles:
            continue
        n_sites += 1
        if alleles[0] != alleles[1]:
            n_het += 1
    if n_sites == 0:
        return None
    return n_het / n_sites


def call_ploidy(het_fraction, threshold=0.01):
    if het_fraction is None:
        return "unknown"
    return "diploid" if het_fraction > threshold else "haploid"


def run_bcftools_forced_diploid_call(cram, reference, region=None):
    mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", reference]
    if region:
        mpileup_cmd += ["-r", region]
    mpileup_cmd.append(cram)
    call_cmd = ["bcftools", "call", "-mv", "--ploidy", "2"]
    mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
    result = subprocess.run(call_cmd, stdin=mpileup.stdout, capture_output=True, text=True, check=True)
    mpileup.stdout.close()
    mpileup.wait()
    return result.stdout.splitlines()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cram", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--strain", required=True)
    parser.add_argument("--region", default=None)
    parser.add_argument("--threshold", type=float, default=0.01)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    vcf_lines = run_bcftools_forced_diploid_call(args.cram, args.reference, args.region)
    het_fraction = het_fraction_from_vcf(vcf_lines)
    ploidy = call_ploidy(het_fraction, args.threshold)

    with open(args.out, "w") as fh:
        fh.write("strain,inferred_ploidy,het_fraction,method\n")
        hf = f"{het_fraction:.6f}" if het_fraction is not None else "NA"
        fh.write(f"{args.strain},{ploidy},{hf},custom_het_script\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests again to verify they pass**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_custom_het_ploidy.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
chmod +x bin/custom_het_ploidy.py
git add bin/custom_het_ploidy.py tests/test_custom_het_ploidy.py
git commit -m "add custom_het_ploidy core logic + unit tests"
```

---

## Task 5: `CUSTOM_HET_PLOIDY` Nextflow module

**Files:**
- Create: `modules/local/custom_het_ploidy/main.nf`

**Interfaces:**
- Consumes: `bin/custom_het_ploidy.py` (Task 4).
- Produces: process `CUSTOM_HET_PLOIDY`, input `tuple val(strain), path(cram), path(crai)` + `path reference` + `path reference_fai`, output `tuple val(strain), path("${strain}.ploidy_inference.csv"), emit: csv`. Consumed by Task 8 (`PLOIDY_INFERENCE` subworkflow).

- [ ] **Step 1: Write the module**

```groovy
// modules/local/custom_het_ploidy/main.nf
process CUSTOM_HET_PLOIDY {
    tag "$strain"
    label 'process_low'
    container 'depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1'

    input:
    tuple val(strain), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    tuple val(strain), path("${strain}.ploidy_inference.csv"), emit: csv

    script:
    """
    custom_het_ploidy.py \\
        --cram ${cram} \\
        --reference ${reference} \\
        --strain ${strain} \\
        --out ${strain}.ploidy_inference.csv
    """
}
```

- [ ] **Step 2: Write a standalone smoke-test workflow against the synthetic fixtures**

```groovy
// tests/smoke_custom_het_ploidy.nf
nextflow.enable.dsl = 2
include { CUSTOM_HET_PLOIDY } from '../modules/local/custom_het_ploidy/main.nf'

workflow {
    ref     = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ch = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )
    CUSTOM_HET_PLOIDY(ch, ref, ref_fai)
    CUSTOM_HET_PLOIDY.out.csv.view()
}
```

- [ ] **Step 3: Run the smoke test and inspect output**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow run tests/smoke_custom_het_ploidy.nf -profile hpcc` (or without `-profile hpcc` if running interactively off-cluster with local Singularity access)
Expected: two CSV rows emitted; `cat work/*/*/haploid_strain.ploidy_inference.csv` shows `inferred_ploidy=haploid`, and the diploid_strain one shows `inferred_ploidy=diploid`. If either is wrong, the `dwgsim -r` heterozygosity rate in Task 3's fixture script needs adjusting, not the module logic (verify with `bcftools stats` on the intermediate VCF before changing code).

- [ ] **Step 4: Commit**

```bash
git add modules/local/custom_het_ploidy/main.nf tests/smoke_custom_het_ploidy.nf
git commit -m "add CUSTOM_HET_PLOIDY module + smoke test"
```

---

## Task 6: `ploidy_crosscheck.py` core logic + unit tests

**Files:**
- Create: `bin/ploidy_crosscheck.py`
- Test: `tests/test_ploidy_crosscheck.py`

**Interfaces:**
- Produces: `crosscheck(inferred_rows: list[dict], metadata_ploidy_by_strain: dict) -> list[dict]` (each dict has `strain, inferred_ploidy, metadata_ploidy, status` where status is one of `AGREE, DISAGREE, NO_METADATA, UNKNOWN_INFERENCE`), CLI `ploidy_crosscheck.py --inferred X.csv --metadata metadata.txt --out ploidy_review.csv`.
- Consumed by: Task 7 (Nextflow module wrapper).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_ploidy_crosscheck.py
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from ploidy_crosscheck import crosscheck


def test_agree():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {"S1": "haploid"}
    assert crosscheck(rows, meta)[0]["status"] == "AGREE"


def test_disagree():
    rows = [{"strain": "S1", "inferred_ploidy": "diploid"}]
    meta = {"S1": "haploid"}
    assert crosscheck(rows, meta)[0]["status"] == "DISAGREE"


def test_no_metadata():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {}
    assert crosscheck(rows, meta)[0]["status"] == "NO_METADATA"


def test_haploid_from_hybrid_normalizes_to_haploid():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {"S1": "haploid_from_hybrid"}
    assert crosscheck(rows, meta)[0]["status"] == "AGREE"


def test_unknown_inference():
    rows = [{"strain": "S1", "inferred_ploidy": "unknown"}]
    meta = {"S1": "diploid"}
    assert crosscheck(rows, meta)[0]["status"] == "UNKNOWN_INFERENCE"
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_ploidy_crosscheck.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'ploidy_crosscheck'`.

- [ ] **Step 3: Write `bin/ploidy_crosscheck.py`**

```python
#!/usr/bin/env python3
"""Cross-check computationally inferred ploidy against metadata.txt's
existing ploidy column, producing ploidy_review.csv for manual sign-off.

Agreement is NOT auto-finalized - a human still produces the actual
ploidy_overrides.csv the pipeline reads for calling (see the design spec's
"Ploidy inference" section). This script only flags what needs attention.
"""
import argparse
import csv


def crosscheck(inferred_rows, metadata_ploidy_by_strain):
    results = []
    for row in inferred_rows:
        strain = row["strain"]
        inferred = row["inferred_ploidy"]
        meta = metadata_ploidy_by_strain.get(strain, "")
        meta_norm = "haploid" if meta == "haploid_from_hybrid" else meta
        if inferred == "unknown":
            status = "UNKNOWN_INFERENCE"
        elif not meta_norm:
            status = "NO_METADATA"
        elif meta_norm == inferred:
            status = "AGREE"
        else:
            status = "DISAGREE"
        results.append({
            "strain": strain,
            "inferred_ploidy": inferred,
            "metadata_ploidy": meta,
            "status": status,
        })
    return results


def load_inferred_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_metadata_ploidy(path):
    by_strain = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            strain = row["strain"].strip()
            ploidy = row.get("ploidy", "").strip()
            if strain not in by_strain or ploidy:
                by_strain[strain] = ploidy
    return by_strain


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inferred", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    inferred_rows = load_inferred_csv(args.inferred)
    metadata_ploidy = load_metadata_ploidy(args.metadata)
    results = crosscheck(inferred_rows, metadata_ploidy)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["strain", "inferred_ploidy", "metadata_ploidy", "status"])
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests again to verify they pass**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_ploidy_crosscheck.py -v`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
chmod +x bin/ploidy_crosscheck.py
git add bin/ploidy_crosscheck.py tests/test_ploidy_crosscheck.py
git commit -m "add ploidy_crosscheck core logic + unit tests"
```

---

## Task 7: `PLOIDY_CROSSCHECK` Nextflow module

**Files:**
- Create: `modules/local/ploidy_crosscheck/main.nf`

**Interfaces:**
- Consumes: `bin/ploidy_crosscheck.py` (Task 6).
- Produces: process `PLOIDY_CROSSCHECK`, input `path inferred_csv, path metadata_txt`, output `path 'ploidy_review.csv', emit: review`. Consumed by Task 8.

- [ ] **Step 1: Write the module**

```groovy
// modules/local/ploidy_crosscheck/main.nf
process PLOIDY_CROSSCHECK {
    label 'process_low'
    container 'python:3.11-slim'

    input:
    path inferred_csv
    path metadata_txt

    output:
    path 'ploidy_review.csv', emit: review

    script:
    """
    ploidy_crosscheck.py --inferred ${inferred_csv} --metadata ${metadata_txt} --out ploidy_review.csv
    """
}
```

- [ ] **Step 2: Verify it parses with a stub run**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow inspect modules/local/ploidy_crosscheck/main.nf` (or, if `inspect` is unavailable in the installed Nextflow version, `nextflow config` on a throwaway workflow that only includes this module) to confirm no syntax errors.
Expected: no errors.

- [ ] **Step 3: Commit**

```bash
git add modules/local/ploidy_crosscheck/main.nf
git commit -m "add PLOIDY_CROSSCHECK module"
```

---

## Task 8: Assemble the `PLOIDY_INFERENCE` subworkflow

**Files:**
- Create: `subworkflows/local/ploidy_inference/main.nf`

**Interfaces:**
- Consumes: `CUSTOM_HET_PLOIDY` (Task 5), `PLOIDY_CROSSCHECK` (Task 7).
- Produces: workflow `PLOIDY_INFERENCE`, `take: cram_ch (channel: [strain, cram, crai]), reference (path), reference_fai (path), metadata_txt (path)`, `emit: inferred_csv (path, collected all-strains CSV), review (path, ploidy_review.csv)`. Consumed by Task 9's channel wiring (via a human-edited `ploidy_overrides.csv` derived from `review`, not directly).

- [ ] **Step 1: Write the subworkflow**

```groovy
// subworkflows/local/ploidy_inference/main.nf
include { CUSTOM_HET_PLOIDY } from '../../../modules/local/custom_het_ploidy/main.nf'
include { PLOIDY_CROSSCHECK } from '../../../modules/local/ploidy_crosscheck/main.nf'

workflow PLOIDY_INFERENCE {
    take:
    cram_ch
    reference
    reference_fai
    metadata_txt

    main:
    CUSTOM_HET_PLOIDY(cram_ch, reference, reference_fai)

    CUSTOM_HET_PLOIDY.out.csv
        .map { strain, csv -> csv }
        .collectFile(name: 'ploidy_inference_all.csv', keepHeader: true, skip: 1, storeDir: "${params.outdir}/ploidy")
        .set { inferred_csv }

    PLOIDY_CROSSCHECK(inferred_csv, metadata_txt)

    emit:
    inferred_csv = inferred_csv
    review       = PLOIDY_CROSSCHECK.out.review
}
```

- [ ] **Step 2: Write a smoke test using the synthetic fixtures + a tiny fixture metadata file**

```
# tests/fixtures/metadata_fixture.txt (tab-separated, matches metadata.txt column layout)
sample	hashtag	species	strain	ploidy	MAT_type	Origin	Environment
haploid_strain	na	R. mucilaginosa	haploid_strain	haploid	A2		
diploid_strain	na	R. mucilaginosa	diploid_strain	diploid	A1/A2		
```

```groovy
// tests/smoke_ploidy_inference.nf
nextflow.enable.dsl = 2
include { PLOIDY_INFERENCE } from '../subworkflows/local/ploidy_inference/main.nf'

workflow {
    ref     = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    meta    = file("${projectDir}/../tests/fixtures/metadata_fixture.txt")
    ch = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )
    PLOIDY_INFERENCE(ch, ref, ref_fai, meta)
    PLOIDY_INFERENCE.out.review.view { "review: " + it.text }
}
```

- [ ] **Step 3: Run and verify AGREE status for both fixture strains**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow run tests/smoke_ploidy_inference.nf -profile hpcc`
Expected: `ploidy_review.csv` shows `status=AGREE` for both `haploid_strain` and `diploid_strain` (since the fixture metadata was set to match what Task 5's smoke test already confirmed the inference produces).

- [ ] **Step 4: Commit**

```bash
git add subworkflows/local/ploidy_inference/main.nf tests/smoke_ploidy_inference.nf tests/fixtures/metadata_fixture.txt
git commit -m "assemble PLOIDY_INFERENCE subworkflow + smoke test"
```

---

## Task 9: `GATK4_HAPLOTYPECALLER` module with per-strain ploidy

**Files:**
- Create: `modules/local/gatk4_haplotypecaller/main.nf`
- Create: `bin/ploidy_overrides_to_channel.py` (helper used inline via Nextflow's `splitCsv`, actually not needed as a separate script - see Step 1 note)

**Interfaces:**
- Produces: process `GATK4_HAPLOTYPECALLER`, input `tuple val(strain), val(ploidy), path(cram), path(crai)` + `path reference, path reference_fai, path reference_dict`, output `tuple val(strain), path("${strain}.g.vcf.gz"), path("${strain}.g.vcf.gz.tbi"), emit: gvcf`. Consumed by Task 10.
- Consumes: a `ploidy_overrides.csv` with header `strain,ploidy` (human-edited output derived from Task 8's `ploidy_review.csv`; `ploidy` column values are the literal strings `haploid` or `diploid`, mapped to GATK's `1`/`2` in this task's channel-construction step).

- [ ] **Step 1: Write the module**

```groovy
// modules/local/gatk4_haplotypecaller/main.nf
process GATK4_HAPLOTYPECALLER {
    tag "$strain"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(strain), val(ploidy), path(cram), path(crai)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(strain), path("${strain}.g.vcf.gz"), path("${strain}.g.vcf.gz.tbi"), emit: gvcf

    script:
    """
    gatk HaplotypeCaller \\
        -R ${reference} \\
        -I ${cram} \\
        -O ${strain}.g.vcf.gz \\
        --sample-ploidy ${ploidy} \\
        -ERC GVCF
    """
}
```

- [ ] **Step 2: Write the channel-construction helper as an inline test workflow (this logic lives in the real `workflows/genotype_population.nf` from Task 13, but is validated standalone here first)**

```groovy
// tests/smoke_haplotypecaller.nf
nextflow.enable.dsl = 2
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' - only 'haploid' or 'diploid' are valid in ploidy_overrides.csv"
}

workflow {
    ref      = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai  = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ref_dict = file("${projectDir}/../tests/fixtures/ref.dict")

    overrides = Channel.of(
        ['haploid_strain', 'haploid'],
        ['diploid_strain', 'diploid'],
    )
    crams = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )

    ch = overrides
        .join(crams)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ch, ref, ref_fai, ref_dict)
    GATK4_HAPLOTYPECALLER.out.gvcf.view()
}
```

- [ ] **Step 3: Run and verify both GVCFs are produced**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow run tests/smoke_haplotypecaller.nf -profile hpcc`
Expected: two `[strain, gvcf, tbi]` tuples emitted; `zcat work/*/*/haploid_strain.g.vcf.gz | grep -v '^#' | head -1` shows a GT field consistent with ploidy 1 (e.g. `1` or `0`, never `0/1`); the diploid one shows standard `0/1`-style diploid genotypes.

- [ ] **Step 4: Commit**

```bash
git add modules/local/gatk4_haplotypecaller/main.nf tests/smoke_haplotypecaller.nf
git commit -m "add GATK4_HAPLOTYPECALLER module + ploidy-branching smoke test"
```

(Delete the `bin/ploidy_overrides_to_channel.py` line from the Files list above if unused — the `ploidyCodeFor` closure inline in the workflow file was sufficient; no separate script was needed.)

---

## Task 10: `population_sets.yaml` asset + joint genotyping subworkflow

**Files:**
- Create: `assets/population_sets.yaml`
- Create: `modules/local/gatk4_genomicsdbimport/main.nf`
- Create: `modules/local/gatk4_genotypegvcfs/main.nf`
- Create: `subworkflows/local/joint_genotyping/main.nf`

**Interfaces:**
- Produces: workflow `JOINT_GENOTYPING`, `take: gvcf_ch (channel: [strain, gvcf, tbi]), population_sets_yaml (path), reference, reference_fai, reference_dict`, `emit: vcf_by_population (channel: [population_name, vcf, tbi])`. Consumed by Task 11.

- [ ] **Step 1: Port `population_sets.yaml` with the `all` group**

Copy the existing file from the prior pipeline (adjust path once confirmed):
```bash
cp /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/population_sets.yaml \
   /rhome/jstajich/projects/nf/nf_genotype_population/assets/population_sets.yaml
```
Then manually add (or script-generate from `reconciliation_table.csv`'s `INCLUDE` rows) an `all` key at the top listing every confirmed *R. mucilaginosa* strain, e.g.:
```yaml
all:
  - DBVPG_10619
  - DBVPG_10656
  # ... every strain with scope == INCLUDE in reconciliation_table.csv
# (existing sliced groups below, unchanged)
```

- [ ] **Step 2: Write the GenomicsDBImport module**

```groovy
// modules/local/gatk4_genomicsdbimport/main.nf
process GATK4_GENOMICSDBIMPORT {
    tag "$population"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(population), path(gvcfs), path(tbis)
    path reference
    path reference_fai
    path reference_dict
    val intervals

    output:
    tuple val(population), path("${population}_gdb"), emit: genomicsdb

    script:
    def gvcf_args = gvcfs.collect { "-V ${it}" }.join(' ')
    """
    gatk GenomicsDBImport \\
        ${gvcf_args} \\
        --genomicsdb-workspace-path ${population}_gdb \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --bypass-feature-reader \\
        -L ${intervals}
    """
}
```

- [ ] **Step 3: Write the GenotypeGVCFs module**

```groovy
// modules/local/gatk4_genotypegvcfs/main.nf
process GATK4_GENOTYPEGVCFS {
    tag "$population"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(population), path(genomicsdb)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.vcf.gz"), path("${population}.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk GenotypeGVCFs \\
        -R ${reference} \\
        -V gendb://${genomicsdb} \\
        -O ${population}.vcf.gz
    """
}
```

- [ ] **Step 4: Assemble the subworkflow**

```groovy
// subworkflows/local/joint_genotyping/main.nf
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/local/gatk4_genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS }   from '../../../modules/local/gatk4_genotypegvcfs/main.nf'

workflow JOINT_GENOTYPING {
    take:
    gvcf_ch              // channel: [strain, gvcf, tbi]
    population_sets_yaml // path to assets/population_sets.yaml
    reference
    reference_fai
    reference_dict
    intervals            // e.g. the whole-genome interval list/contig name

    main:
    def populations = new org.yaml.snakeyaml.Yaml().load(population_sets_yaml.text)

    def gvcf_map = [:]
    gvcf_ch.subscribe { strain, gvcf, tbi -> gvcf_map[strain] = [gvcf, tbi] }

    population_ch = Channel.fromList(
        populations.collect { name, strains ->
            def present = strains.findAll { gvcf_map.containsKey(it) }
            [name, present.collect { gvcf_map[it][0] }, present.collect { gvcf_map[it][1] }]
        }
    )

    GATK4_GENOMICSDBIMPORT(population_ch, reference, reference_fai, reference_dict, intervals)
    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, reference, reference_fai, reference_dict)

    emit:
    vcf_by_population = GATK4_GENOTYPEGVCFS.out.vcf
}
```

**Note for the implementer:** the `.subscribe` + closure-captured-map pattern above is a pragmatic way to turn a channel into a lookup map inside a Groovy closure; it works because Nextflow channels are single-consumer-safe within one `main:` block, but it means `population_ch` must not be constructed until the `gvcf_ch` channel has been fully drained. If this ordering proves unreliable in practice (e.g. flaky joins on a real multi-strain run), replace it with a `.toList()` + `.map` reshaping the gvcf channel into a `[strain, gvcf, tbi]` list and combine with population membership via `.combine()`/`.groupTuple()` instead of a subscribe-based side effect - flag this as a real risk to watch during Task 13's full run, not a guaranteed-safe pattern.

- [ ] **Step 5: Smoke test against the two synthetic strains under an `all` group**

```yaml
# tests/fixtures/population_sets_fixture.yaml
all:
  - haploid_strain
  - diploid_strain
```
```groovy
// tests/smoke_joint_genotyping.nf
nextflow.enable.dsl = 2
include { JOINT_GENOTYPING } from '../subworkflows/local/joint_genotyping/main.nf'
// ... construct gvcf_ch from Task 9's smoke test outputs (re-run GATK4_HAPLOTYPECALLER
// here, or point at the work/ directory GVCFs already produced by Task 9's smoke test)
```
Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow run tests/smoke_joint_genotyping.nf -profile hpcc`
Expected: `all.vcf.gz` produced containing records for both synthetic strains.

- [ ] **Step 6: Commit**

```bash
git add assets/population_sets.yaml modules/local/gatk4_genomicsdbimport/main.nf modules/local/gatk4_genotypegvcfs/main.nf subworkflows/local/joint_genotyping/main.nf tests/smoke_joint_genotyping.nf tests/fixtures/population_sets_fixture.yaml
git commit -m "add population_sets.yaml (with all group) + joint genotyping subworkflow"
```

---

## Task 11: `GATK4_VARIANTFILTRATION` hard-filter module

**Files:**
- Create: `modules/local/gatk4_variantfiltration/main.nf`
- Test: `tests/test_filter_expressions.py` (pure string-constant check, not a GATK integration test)

**Interfaces:**
- Produces: process `GATK4_HARDFILTER`, input `tuple val(population), path(vcf), path(tbi)` + reference paths, output `tuple val(population), path("${population}.filtered.vcf.gz"), path("${population}.filtered.vcf.gz.tbi"), emit: vcf`. Consumed by Task 12.
- Consumes: `JOINT_GENOTYPING.out.vcf_by_population` (Task 10).

- [ ] **Step 1: Write a failing test pinning the exact filter-expression strings to the spec**

```python
# tests/test_filter_expressions.py
import re

MODULE_PATH = "modules/local/gatk4_variantfiltration/main.nf"


def read_module():
    with open(MODULE_PATH) as fh:
        return fh.read()


def test_snp_filter_thresholds_present():
    text = read_module()
    for expr in ["QD < 2.0", "FS > 60.0", "MQ < 40.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0", "SOR > 3.0"]:
        assert expr in text, f"missing SNP filter expression: {expr}"


def test_indel_filter_thresholds_present():
    text = read_module()
    for expr in ["QD < 2.0", "FS > 200.0", "ReadPosRankSum < -20.0", "SOR > 10.0"]:
        assert expr in text, f"missing indel filter expression: {expr}"
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_filter_expressions.py -v`
Expected: FAIL — `FileNotFoundError` (module doesn't exist yet).

- [ ] **Step 3: Write the module**

```groovy
// modules/local/gatk4_variantfiltration/main.nf
process GATK4_HARDFILTER {
    tag "$population"
    label 'process_low'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(population), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.filtered.vcf.gz"), path("${population}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include SNP -O ${population}.snp.vcf.gz
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include INDEL -O ${population}.indel.vcf.gz

    gatk VariantFiltration -R ${reference} -V ${population}.snp.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        -O ${population}.snp.filtered.vcf.gz

    gatk VariantFiltration -R ${reference} -V ${population}.indel.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 200.0" --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        --filter-expression "SOR > 10.0" --filter-name "SOR10" \\
        -O ${population}.indel.filtered.vcf.gz

    gatk MergeVcfs \\
        -I ${population}.snp.filtered.vcf.gz \\
        -I ${population}.indel.filtered.vcf.gz \\
        -O ${population}.filtered.vcf.gz
    """
}
```

- [ ] **Step 4: Run tests again to verify they pass**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && /usr/bin/python3.12 -m pytest tests/test_filter_expressions.py -v`
Expected: 2 passed.

- [ ] **Step 5: Smoke test against Task 10's `all.vcf.gz` output**

Run: `cd /rhome/jstajich/projects/nf/nf_genotype_population && nextflow run tests/smoke_joint_genotyping.nf -entry-with-filter -profile hpcc` (extend `tests/smoke_joint_genotyping.nf` to also call `GATK4_HARDFILTER` on its output, or add a small `tests/smoke_hardfilter.nf` pointed at the prior task's cached `work/` output)
Expected: `all.filtered.vcf.gz` produced; `zcat all.filtered.vcf.gz | grep -v '^#' | cut -f7 | sort | uniq -c` shows a mix of `PASS` and named filter tags (e.g. `SOR3`, `FS60`), not all-PASS or all-fail (a synthetic 2-strain, 2kb dataset is small, so don't expect the exact real-world ratio - the goal here is confirming the filters actually run and tag records, not validating filtering sensitivity).

- [ ] **Step 6: Commit**

```bash
git add modules/local/gatk4_variantfiltration/main.nf tests/test_filter_expressions.py
git commit -m "add GATK4_HARDFILTER module + filter-expression pinning test"
```

---

## Task 12: `SNPEFF` annotation module

**Files:**
- Create: `modules/local/snpeff/main.nf`
- Create: `conf/snpeff.config` (database build parameters for DH4148)

**Interfaces:**
- Produces: process `SNPEFF_ANNOTATE`, input `tuple val(population), path(vcf), path(tbi)` + `path snpeff_db_dir, val snpeff_genome_name`, output `tuple val(population), path("${population}.annotated.vcf.gz"), emit: vcf`. Consumed by Task 13's root workflow (terminal stage for Phase 1).

- [ ] **Step 1: Confirm whether a SnpEff database for DH4148 already exists**

Run: `find /bigdata/stajichlab/shared -iname "*snpEff*" -o -iname "*snpeff*" 2>/dev/null | grep -i rhodo` and check `Rmuc_popgen_NRRLY2510`'s scripts for how its SnpEff database was built (custom `snpEff.config` entry + `snpEff build -gtf22` against the same GTF this project already has at `refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.gtf`). Reuse that exact recipe rather than re-deriving one - if `Rmuc_popgen_NRRLY2510` doesn't have one either, building it is a same-shaped one-off task using this project's own GTF, not new design work.

- [ ] **Step 2: Write the module**

```groovy
// modules/local/snpeff/main.nf
process SNPEFF_ANNOTATE {
    tag "$population"
    label 'process_low'
    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'

    input:
    tuple val(population), path(vcf), path(tbi)
    path snpeff_db_dir
    val snpeff_genome_name

    output:
    tuple val(population), path("${population}.annotated.vcf.gz"), emit: vcf

    script:
    """
    snpEff -dataDir ${snpeff_db_dir} ${snpeff_genome_name} ${vcf} | bgzip > ${population}.annotated.vcf.gz
    tabix -p vcf ${population}.annotated.vcf.gz
    """
}
```

- [ ] **Step 3: Smoke test against Task 11's filtered VCF output**

Run against the synthetic `all.filtered.vcf.gz` from Task 11 with a minimal SnpEff database built from `tests/fixtures/ref.fa` (or, if the real DH4148 database already exists per Step 1, use it directly since coordinates won't match the synthetic contig but the goal here is confirming the module runs end-to-end without error, not real annotation accuracy).
Expected: `all.annotated.vcf.gz` produced without error, `ANN=` field present in at least one INFO column when opened with `zcat | grep ANN=`.

- [ ] **Step 4: Commit**

```bash
git add modules/local/snpeff/main.nf conf/snpeff.config
git commit -m "add SNPEFF_ANNOTATE module"
```

---

## Task 13: Assemble the root workflow + HPCC profile, run full synthetic integration test

**Files:**
- Create: `workflows/genotype_population.nf`
- Modify: `main.nf`
- Create: `conf/profile_hpcc.config`
- Modify: `conf/modules.config` (create if not already present from earlier tasks)

**Interfaces:**
- Consumes: every module/subworkflow from Tasks 5, 7, 8, 9, 10, 11, 12.
- Produces: the full Phase 1 pipeline entry point, runnable as `nextflow run main.nf -profile hpcc`.

- [ ] **Step 1: Write `conf/profile_hpcc.config`**

```groovy
// conf/profile_hpcc.config
// Mirrors Rhodotorula_mucilaginosa_DH4148_ref/conf/profile_sarek.config conventions.
params {
    cram_dir           = null   // directory of Sarek CRAM output, e.g. .../results/preprocessing/*/*.cram
    reference          = null
    reference_fai      = null
    reference_dict     = null
    metadata           = null   // metadata.txt
    ploidy_overrides   = null   // human-reviewed strain,ploidy CSV (see design spec)
    population_sets    = "${projectDir}/assets/population_sets.yaml"
    snpeff_db_dir      = null
    snpeff_genome_name = null
    outdir             = "${launchDir}/results"
}

workDir       = "${launchDir}/work"
timeline.file = "${launchDir}/logs/genotype_population_timeline.html"
report.file   = "${launchDir}/logs/genotype_population_report.html"
trace.file    = "${launchDir}/logs/genotype_population_trace.txt"
dag.file      = "${launchDir}/logs/genotype_population_dag.html"
```

- [ ] **Step 2: Write `workflows/genotype_population.nf`**

```groovy
// workflows/genotype_population.nf
include { PLOIDY_INFERENCE }       from '../subworkflows/local/ploidy_inference/main.nf'
include { GATK4_HAPLOTYPECALLER }  from '../modules/local/gatk4_haplotypecaller/main.nf'
include { JOINT_GENOTYPING }       from '../subworkflows/local/joint_genotyping/main.nf'
include { GATK4_HARDFILTER }       from '../modules/local/gatk4_variantfiltration/main.nf'
include { SNPEFF_ANNOTATE }        from '../modules/local/snpeff/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' in ploidy_overrides.csv - expected 'haploid' or 'diploid'"
}

workflow GENOTYPE_POPULATION {
    take:
    cram_ch           // channel: [strain, cram, crai]
    reference
    reference_fai
    reference_dict
    metadata_txt
    ploidy_overrides  // path to human-reviewed strain,ploidy CSV
    population_sets_yaml
    snpeff_db_dir
    snpeff_genome_name

    main:
    PLOIDY_INFERENCE(cram_ch, reference, reference_fai, metadata_txt)

    overrides_ch = Channel
        .fromPath(ploidy_overrides)
        .splitCsv(header: true)
        .map { row -> [row.strain, row.ploidy] }

    ploidy_cram_ch = overrides_ch
        .join(cram_ch)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ploidy_cram_ch, reference, reference_fai, reference_dict)

    JOINT_GENOTYPING(
        GATK4_HAPLOTYPECALLER.out.gvcf,
        population_sets_yaml,
        reference,
        reference_fai,
        reference_dict,
        reference_fai.baseName, // whole-genome interval placeholder; replace with a real interval list once available
    )

    GATK4_HARDFILTER(JOINT_GENOTYPING.out.vcf_by_population, reference, reference_fai, reference_dict)
    SNPEFF_ANNOTATE(GATK4_HARDFILTER.out.vcf, snpeff_db_dir, snpeff_genome_name)

    emit:
    annotated_vcf = SNPEFF_ANNOTATE.out.vcf
}
```

- [ ] **Step 3: Wire it into `main.nf`**

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOTYPE_POPULATION } from './workflows/genotype_population.nf'

workflow {
    cram_ch = Channel
        .fromFilePairs("${params.cram_dir}/*.{cram,cram.crai}", checkIfExists: true) { file -> file.baseName.replaceFirst(/\.cram$/, '') }
        .map { strain, files -> [strain, files[0], files[1]] }

    GENOTYPE_POPULATION(
        cram_ch,
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict),
        file(params.metadata),
        file(params.ploidy_overrides),
        file(params.population_sets),
        file(params.snpeff_db_dir),
        params.snpeff_genome_name,
    )
}
```

- [ ] **Step 4: Run the full pipeline against the synthetic 2-strain fixture end-to-end**

```bash
cd /rhome/jstajich/projects/nf/nf_genotype_population
mkdir -p /tmp/synth_cram_dir
cp tests/fixtures/haploid_strain.cram tests/fixtures/haploid_strain.cram.crai /tmp/synth_cram_dir/
cp tests/fixtures/diploid_strain.cram tests/fixtures/diploid_strain.cram.crai /tmp/synth_cram_dir/
cat > /tmp/synth_ploidy_overrides.csv <<'EOF'
strain,ploidy
haploid_strain,haploid
diploid_strain,diploid
EOF

nextflow run main.nf -profile hpcc \
    --cram_dir /tmp/synth_cram_dir \
    --reference tests/fixtures/ref.fa \
    --reference_fai tests/fixtures/ref.fa.fai \
    --reference_dict tests/fixtures/ref.dict \
    --metadata tests/fixtures/metadata_fixture.txt \
    --ploidy_overrides /tmp/synth_ploidy_overrides.csv \
    --population_sets tests/fixtures/population_sets_fixture.yaml
```
Expected: pipeline completes with `[100%] 2 of 2` (or equivalent) success; `results/` contains a final `all.annotated.vcf.gz` with records for both synthetic strains and both a `PASS`/named-filter FILTER column and an `ANN=` field.

- [ ] **Step 5: Commit**

```bash
git add workflows/genotype_population.nf main.nf conf/profile_hpcc.config conf/modules.config
git commit -m "assemble root GENOTYPE_POPULATION workflow, run full synthetic integration test"
```

---

## Task 14: README usage documentation

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Document real usage**

Add to `README.md`:
```markdown
## Usage

1. Run Sarek through CRAM output only (see `Rhodotorula_mucilaginosa_DH4148_ref/conf/profile_sarek.config` for the reconfigured, no-variant-calling profile).
2. Run the ploidy-inference-only entry to generate `ploidy_review.csv`, review it by hand, and produce `ploidy_overrides.csv`.
3. Run the full pipeline:

\`\`\`bash
nextflow run main.nf -profile hpcc \\
    --cram_dir /path/to/sarek/results/preprocessing/*/*.cram \\
    --reference /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna \\
    --reference_fai /path/to/refgenome/....fna.fai \\
    --reference_dict /path/to/refgenome/....dict \\
    --metadata /path/to/metadata.txt \\
    --ploidy_overrides /path/to/ploidy_overrides.csv \\
    --population_sets assets/population_sets.yaml \\
    --snpeff_db_dir /path/to/snpeff_data \\
    --snpeff_genome_name DH4148
\`\`\`

## Phase 2+ (not yet implemented)

See `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`
for nQuire/nQuack ploidy methods, population-sliced joint genotyping beyond
`all`, mosdepth/CNV visualization, empirical population-structure proposal,
and full CNV calling.
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "document Phase 1 usage"
```

---

## Self-review notes (fixed inline above, recorded here for transparency)

- Confirmed every task's file paths and emitted channel/type names are used consistently by the task that consumes them (e.g. `PLOIDY_INFERENCE.out.review` -> human review -> `ploidy_overrides.csv` -> `GENOTYPE_POPULATION`'s `overrides_ch`).
- Flagged the `.subscribe`-based map-building pattern in Task 10 as a real risk to watch under real multi-strain load, with a concrete fallback approach named, rather than presenting it as guaranteed-correct.
- Task 12 (SnpEff) explicitly defers to whatever database-build recipe `Rmuc_popgen_NRRLY2510` already used rather than inventing a new one, per Global Constraints' reuse-first principle.
- Left CNV/indel-as-first-class-module, nQuire/nQuack, and population-sliced joint genotyping beyond `all` entirely out of this plan's tasks (Phase 2+), matching the approved spec's phasing.
