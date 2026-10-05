#!/usr/bin/env python3
"""Sample lists used by species_check.sh / hybrid_check.py and pca.sh, from strain_qc.tsv.

  kept_haploid      decision keep, called haploid
  kept_dip_highhet  decision keep, called diploid, het_rate > 0.01
  frig_haploid      species R. frigidialcoholis, called haploid (any decision)
  affmuc_haploid    species R. aff. mucilaginosa, called haploid (any decision)

These rules reproduce the hand-made lists of 2026-10-02 exactly.
Strains in variant_qc_reviewed_keep.tsv count as decision keep (reasons must match).
Usage: make_qc_groups.py strain_qc.tsv outdir [reviewed_keep.tsv]
"""
import csv
import sys

qc_tsv, outdir = sys.argv[1:3]
q = list(csv.DictReader(open(qc_tsv), delimiter="\t"))
if len(sys.argv) > 3:
    reviewed = {r["strain"]: r["drop_reasons"] for r in csv.DictReader(
        (l for l in open(sys.argv[3]) if not l.startswith("#")), delimiter="\t")}
    for r in q:
        if r["strain"] in reviewed:
            assert r["reasons"] == reviewed[r["strain"]], (r["strain"], r["reasons"])
            r["decision"] = "keep"
groups = {
    "kept_haploid": lambda r: r["decision"] == "keep" and r["called_ploidy"] == "1",
    "kept_dip_highhet": lambda r: r["decision"] == "keep" and r["called_ploidy"] == "2"
    and float(r["het_rate"]) > 0.01,
    "frig_haploid": lambda r: r["species"] == "R. frigidialcoholis" and r["called_ploidy"] == "1",
    "affmuc_haploid": lambda r: r["species"] == "R. aff. mucilaginosa" and r["called_ploidy"] == "1",
}
for g, keep in groups.items():
    s = sorted(r["vcf_sample"] for r in q if keep(r))
    with open(f"{outdir}/{g}.txt", "w") as o:
        o.write("".join(x + "\n" for x in s))
    print(g, len(s))
