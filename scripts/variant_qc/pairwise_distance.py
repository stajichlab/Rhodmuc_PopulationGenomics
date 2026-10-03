#!/usr/bin/env python3
"""Pairwise genetic distance between strains (mixed ploidy).

stdin: header line of VCF sample names, then one line per site with
`[%GT\t]` (GQ<20 / DP<5 genotypes already set to missing).
Each genotype becomes an ALT-allele fraction x in {0, 0.5, 1} (haploid 0/1,
diploid 0/0.5/1). Distance between strains i and j over sites called in both:
  d_ij = sum x_i(1-x_j) + x_j(1-x_i) / n_ij
i.e. the chance that one allele drawn from each strain differs. For a diploid
het strain d_ii > 0. Divergence from the reference is mean(x_i).

Writes <out>.dist.tsv.gz (matrix), <out>.n.tsv.gz (sites called in both),
<out>.ref_div.tsv (per-strain mean ALT fraction and called sites).
Usage: pairwise_distance.py OUTPREFIX < gt.tsv
"""
import gzip, sys
import numpy as np
out = sys.argv[1]
names = sys.stdin.readline().rstrip("\n").rstrip("\t").split("\t")
code = {}
def enc(g):
    v = code.get(g)
    if v is None:
        a = g.replace("|", "/").split("/")
        if "." in a:
            v = -1.0
        else:
            v = sum(1 for x in a if x != "0") / len(a)
        code[g] = v
    return v
rows = []
for line in sys.stdin:
    rows.append(np.fromiter((enc(g) for g in line.rstrip("\n").rstrip("\t").split("\t")), dtype=np.float32, count=len(names)))
X = np.vstack(rows).T            # strains x sites
del rows
M = (X >= 0).astype(np.float32)
X = np.where(M > 0, X, 0).astype(np.float32)
Y = M - X                        # REF-allele fraction where called
N = M @ M.T
D = (X @ Y.T + Y @ X.T) / np.maximum(N, 1)
with gzip.open(out + ".dist.tsv.gz", "wt") as o:
    o.write("strain\t" + "\t".join(names) + "\n")
    for i, n in enumerate(names):
        o.write(n + "\t" + "\t".join(f"{v:.6f}" for v in D[i]) + "\n")
with gzip.open(out + ".n.tsv.gz", "wt") as o:
    o.write("strain\t" + "\t".join(names) + "\n")
    for i, n in enumerate(names):
        o.write(n + "\t" + "\t".join(str(int(v)) for v in N[i]) + "\n")
with open(out + ".ref_div.tsv", "w") as o:
    o.write("sample\tcalled_sites\tmean_alt_fraction\n")
    for i, n in enumerate(names):
        c = M[i].sum()
        o.write(f"{n}\t{int(c)}\t{X[i].sum() / max(c, 1):.6f}\n")
print(f"strains={len(names)} sites={X.shape[1]}")
