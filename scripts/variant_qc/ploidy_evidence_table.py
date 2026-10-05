#!/usr/bin/env python3
"""Per-strain ploidy evidence: strain QC + read allele balance + per-Mb het call + nQuire + metadata.

Replaces the hand-joined table of 2026-10-02. The old 'custom'/'custom_het' columns
(het-ratio method, found defective) are replaced by the per-Mb call from CUSTOM_HET_PLOIDY.

Usage: ploidy_evidence_table.py strain_qc.tsv allele_balance.tsv ploidy_inference_all.csv \
           ploidy_review_nquire.csv metadata.txt ploidy_overrides.csv > ploidy_evidence.tsv
"""
import csv
import sys

qc_tsv, ab_tsv, permb_csv, nquire_csv, metadata, overrides = sys.argv[1:7]
ab = {r["sample"]: r for r in csv.DictReader(open(ab_tsv), delimiter="\t")}
permb = {r["strain"]: r for r in csv.DictReader(open(permb_csv))}
nq = {r["strain"]: r["inferred_ploidy"] for r in csv.DictReader(open(nquire_csv))}
override = {r["strain"]: r["ploidy"] for r in csv.DictReader(open(overrides))}
meta = {}
for row in open(metadata):
    f = row.rstrip("\r\n").split("\t")
    if len(f) > 4 and f[3] != "strain":
        meta.setdefault(f[3], f[4])

hdr = ["strain", "sample", "species", "called_ploidy", "decision", "reasons", "het",
       "covered_sites", "balanced_per_kcov", "mixed_per_kcov", "peak", "f45", "f30",
       "nquire", "meta", "override", "per_mb_call", "het_per_mb", "callable_bp"]
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(hdr)
for r in csv.DictReader(open(qc_tsv), delimiter="\t"):
    s, a, p = r["strain"], ab.get(r["vcf_sample"], {}), permb.get(r["strain"], {})
    w.writerow([s, r["vcf_sample"], r["species"], r["called_ploidy"], r["decision"], r["reasons"],
                r["het_rate"], a.get("covered", "NA"), a.get("balanced_per_kcov", "NA"),
                a.get("mixed_per_kcov", "NA"), a.get("peak_bin", "NA"),
                a.get("frac_f_045_050", "NA"), a.get("frac_f_030_037", "NA"),
                nq.get(s, "NA"), meta.get(s, "NA"), override.get(s, "NA"),
                p.get("inferred_ploidy", "NA"), p.get("het_per_mb", "NA"), p.get("callable_bp", "NA")])
