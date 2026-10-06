#!/usr/bin/env python3
"""Per-strain QC table and keep/drop lists from diagnose_calls.py output.

Drop reasons (any one drops the strain):
  other_species   metadata species is not "R. mucilaginosa"
  no_metadata     strain absent from metadata.txt
  low_coverage    < 50% of the genome at >= 5x read depth (mosdepth, MAPQ >= 20),
                  or > 50% of called genotypes have GQ<20.
                  Before 2026-10-05 the rule was "> 50% of called genotypes have
                  GQ<20 or DP<5". FORMAT/DP of a hom-ref genotype is the gVCF
                  block MIN_DP, so that rule dropped strains with good coverage
                  (docs/decision_homref_dp_mask_2026-10-05.md).
  ploidy_mixed    haploid call set, but >25% of >=200 ALT calls have mixed reads
  ploidy_skewed   diploid with het rate > 5% but >40% of hets outside AB 0.2-0.8
Flag only (kept):
  no_het_diploid  diploid call set with het rate < 1% and adequate coverage:
                  haploid or homozygous diploid; consider haploid in overrides
Usage: strain_qc_table.py genome.samples.tsv metadata.txt ploidy_overrides.csv outdir coverage.tsv
coverage.tsv: strain, mean_depth, frac_ge5x (mosdepth_summary.py)
"""
import csv, sys

samples_tsv, metadata, overrides, outdir, coverage = sys.argv[1:6]
cov = {r["strain"]: r for r in csv.DictReader(open(coverage), delimiter="\t")}
species = {}
for row in open(metadata):
    f = row.rstrip("\n").split("\t")
    if len(f) > 4 and f[3] != "strain":
        species[f[3]] = f[2]
override = {r["strain"]: r["ploidy"] for r in csv.DictReader(open(overrides))}

rows = []
DEC, RSN = 15, 16  # column index of decision, reasons
for x in csv.DictReader(open(samples_tsv), delimiter="\t"):
    name = x["sample"]
    h = len(name) // 2
    strain = name[:h] if name[h] == "_" and name[:h] == name[h + 1:] else name
    n = {k: int(v) for k, v in x.items() if k not in ("sample", "ploidy")}
    called = max(1, n["called"])
    lowq = max(n["gq_lt20"], n["dp_lt5"]) / called
    het = n["het"] / called
    abext = n["dip_het_ABext"] / max(1, n["het"])
    mixed = n["hap_alt_mixedAD"] / max(1, n["hap_alt_called"])
    sp = species.get(strain)
    c = cov[strain]  # KeyError = no mosdepth result for this strain: fail loudly
    ge5 = float(c["frac_ge5x"])
    drop, flag = [], []
    if sp is None:
        drop.append("no_metadata")
    elif sp != "R. mucilaginosa":
        drop.append("other_species")
    if ge5 < 0.5 or n["gq_lt20"] / called > 0.5:
        drop.append("low_coverage")
    if x["ploidy"] == "1" and n["hap_alt_called"] >= 200 and mixed > 0.25:
        drop.append("ploidy_mixed")
    if x["ploidy"] == "2" and het > 0.05 and abext > 0.4:
        drop.append("ploidy_skewed")
    if x["ploidy"] == "2" and het < 0.01 and "low_coverage" not in drop:
        flag.append("no_het_diploid")
    rows.append([strain, name, sp or "NA", override.get(strain, "NA"), x["ploidy"],
                 f"{lowq:.3f}", f"{n['gq_lt20']/called:.3f}", f"{n['dp_lt5']/called:.3f}",
                 f"{n['nonref']/called:.3f}", f"{het:.3f}", f"{abext:.3f}", f"{mixed:.3f}",
                 str(n["hap_alt_called"]), c["mean_depth"], f"{ge5:.3f}", "drop" if drop else "keep",
                 ",".join(drop + flag) or "-"])

hdr = ["strain", "vcf_sample", "species", "override_ploidy", "called_ploidy",
       "frac_lowq_gt", "frac_gq_lt20", "frac_dp_lt5", "frac_nonref", "het_rate",
       "het_AB_extreme", "hap_alt_mixed_reads", "hap_alt_calls", "mosdepth_mean", "genome_ge5x",
       "decision", "reasons"]
with open(f"{outdir}/strain_qc.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in sorted(rows, key=lambda r: (r[DEC], r[RSN], r[0])):
        o.write("\t".join(r) + "\n")
with open(f"{outdir}/excluded_strains.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in sorted(rows, key=lambda r: (r[RSN], r[0])):
        if r[DEC] == "drop":
            o.write("\t".join(r) + "\n")
with open(f"{outdir}/rmucilaginosa_qc.samples.txt", "w") as o:
    for r in rows:
        if r[DEC] == "keep":
            o.write(r[1] + "\n")
from collections import Counter
print("keep:", sum(r[DEC] == "keep" for r in rows), "drop:", sum(r[DEC] == "drop" for r in rows))
print(Counter(r[RSN] for r in rows if r[DEC] == "drop"))
print("flagged kept:", [r[0] for r in rows if r[DEC] == "keep" and r[RSN] != "-"])
