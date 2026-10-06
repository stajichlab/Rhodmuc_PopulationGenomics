#!/usr/bin/env python3
"""Differences among the strains that are nearly identical to the reference DH4148.

stdin: first line = tab-separated VCF sample names; then one line per site:
  %CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP:%GQ:%AD]
The sites are the candidate sites: PASS, outside the repeat mask, and a non-ref
call in at least one near-ref strain (see nearref_differences.sh).

Genotype rules (the VARIANT_QC_FILTER rules of 2026-10-05):
  missing  GT has '.', or GQ < 20
  hom-ref  all alleles 0 (no DP threshold: gVCF block MIN_DP)
  non-ref  DP >= 5, and for a haploid/hom call the main allele has >= 0.8 of AD
Two strains differ at a site if both are called and their genotypes differ.
All non-candidate sites are hom-ref or missing in every near-ref strain, so they
add no differences among them.

Usage: nearref_differences.py nearref_samples.txt metadata.txt clone_groups.tsv OUTPREFIX < query.tsv
Writes:
  OUTPREFIX.strains.tsv        one row per near-ref strain
  OUTPREFIX.pairwise.tsv.gz    strain_a, strain_b, diffs, compared (near-ref pairs)
  OUTPREFIX.sites.tsv.gz       per candidate site: ALT carriers in and outside the set
"""
import gzip
import sys

import numpy as np

nr_file, metadata, clones, out = sys.argv[1:5]


def strain_of(name):
    h = len(name) // 2
    return name[:h] if len(name) > 2 and name[h] == "_" and name[:h] == name[h + 1:] else name


nearref = {l.strip() for l in open(nr_file) if l.strip()}
meta = {}
for row in open(metadata):
    f = row.rstrip("\r\n").split("\t")
    if len(f) > 7 and f[3] != "strain":
        meta.setdefault(f[3], (f[6], f[7]))
cg = {}
for r in open(clones):
    f = r.rstrip("\n").split("\t")
    if f[0] != "clone_group":
        cg[f[2]] = f[0]

samples = sys.stdin.readline().rstrip("\n").split("\t")
strains = [strain_of(s) for s in samples]
nr_idx = [i for i, s in enumerate(strains) if s in nearref]
missing_nr = nearref - {strains[i] for i in nr_idx}
if missing_nr:
    sys.stderr.write(f"near-ref strains not in VCF: {sorted(missing_nr)}\n")
is_nr = np.zeros(len(samples), bool)
is_nr[nr_idx] = True


def call(field):
    """-1 missing, 0 hom-ref, k>0 an id for the genotype (sorted allele string)."""
    p = field.split(":")
    gt = p[0].replace("|", "/")
    if "." in gt:
        return None
    gq = p[2] if len(p) > 2 else "."
    if gq == "." or int(gq) < 20:
        return None
    al = gt.split("/")
    if all(a == "0" for a in al):
        return "0"
    dp = p[1] if len(p) > 1 else "."
    if dp == "." or int(dp) < 5:
        return None
    if len(set(al)) == 1 and len(p) > 3 and p[3] != ".":
        ad = [int(x) for x in p[3].split(",") if x != "."]
        tot = sum(ad)
        if tot == 0 or max(ad) / tot < 0.8:
            return None
    return "/".join(sorted(al))


cols, site_rows = [], []
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    g = [call(x) for x in f[4:]]
    ids, v = {"0": 0}, np.full(len(g), -1, np.int16)
    for i, x in enumerate(g):
        if x is not None:
            v[i] = ids.setdefault(x, len(ids))
    alt = v > 0
    called = v >= 0
    n_nr_called = int(called[is_nr].sum())
    n_nr_alt = int(alt[is_nr].sum())
    if n_nr_alt == 0:
        continue
    n_out_alt = int(alt[~is_nr].sum())
    site_rows.append((f[0], f[1], f[2], f[3], n_nr_alt, n_nr_called, n_out_alt,
                      int(called[~is_nr].sum())))
    cols.append(v)

G = np.vstack(cols).T  # samples x sites
S = np.array([[r[4], r[5], r[6]] for r in site_rows])  # nr_alt, nr_called, out_alt
lineage_wide = S[:, 0] >= 0.9 * S[:, 1]
alt_all = (G > 0)
alt_count_all = alt_all.sum(axis=0)

GN = G[nr_idx]
CN = GN >= 0
names = [strains[i] for i in nr_idx]
n = len(names)
D = np.zeros((n, n), np.int32)
N = np.zeros((n, n), np.int32)
for i in range(n):
    both = CN[i] & CN
    D[i] = ((GN[i] != GN) & both).sum(axis=1)
    N[i] = both.sum(axis=1)

hdr = ["strain", "Origin", "Environment", "clone_group", "candidate_sites_called",
       "alt_calls", "alt_lineage_wide", "alt_shared_nearref_only", "alt_shared_outside_set",
       "alt_private", "nearest_strain", "nearest_diffs", "nearest_compared",
       "n_strains_0_diffs", "n_strains_le5_diffs", "n_strains_le20_diffs", "median_diffs_to_set"]
with open(out + ".strains.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for k, s in enumerate(names):
        a = GN[k] > 0
        lw = int((a & lineage_wide).sum())
        private = int((a & (alt_count_all == 1)).sum())
        outside = int((a & ~lineage_wide & (S[:, 2] > 0)).sum())
        nr_only = int(a.sum()) - lw - private - outside
        d = D[k].astype(float)
        d[k] = np.inf
        j = int(np.argmin(d))
        others = np.delete(D[k], k)
        o.write("\t".join(map(str, [
            s, *meta.get(s, ("NA", "NA")), cg.get(s, "NA"), int(CN[k].sum()), int(a.sum()), lw,
            nr_only, outside, private, names[j], int(D[k, j]), int(N[k, j]),
            int((others == 0).sum()), int((others <= 5).sum()), int((others <= 20).sum()),
            int(np.median(others))])) + "\n")

with gzip.open(out + ".pairwise.tsv.gz", "wt") as o:
    o.write("strain_a\tstrain_b\tdiffs\tcompared\n")
    for i in range(n):
        for j in range(i + 1, n):
            o.write(f"{names[i]}\t{names[j]}\t{D[i, j]}\t{N[i, j]}\n")

with gzip.open(out + ".sites.tsv.gz", "wt") as o:
    o.write("chrom\tpos\tref\talt\tnearref_alt\tnearref_called\toutside_alt\toutside_called\tlineage_wide\n")
    for r, lw in zip(site_rows, lineage_wide):
        o.write("\t".join(map(str, r)) + f"\t{int(lw)}\n")

sys.stderr.write(f"{n} near-ref strains, {len(site_rows)} candidate sites, "
                 f"{int(lineage_wide.sum())} lineage-wide\n")
