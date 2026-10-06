#!/usr/bin/env python3
"""Pairwise SNP differences among all strains of a group VCF, near-identical groups, one representative each.

stdin: first line = tab-separated VCF sample names; then `[%GT\t]` per SNP site
(from the group .qc VCF: genotypes already masked by VARIANT_QC_FILTER).
Two strains differ at a site if both are called and the genotypes differ
(allele order ignored, so 0|1 = 0/1). Single linkage at <= CUTOFF differences.
Representative per group: the strain with the fewest missing genotypes; ties ->
higher mosdepth mean depth (coverage.tsv), then name.

Usage: declone.py CUTOFF coverage.tsv metadata.txt OUTPREFIX < gt.tsv
Writes:
  OUTPREFIX.pairwise.tsv.gz      strain_a, strain_b, diffs, compared (all pairs)
  OUTPREFIX.groups_leCUTOFF.tsv  near_identical_group, size, strain, representative (1/0),
                                 missing_frac, mosdepth_mean, Origin, Environment,
                                 max_diffs_within, min_diffs_to_outside
  OUTPREFIX.representatives_leCUTOFF.txt  one strain per group + every ungrouped strain
"""
import gzip
import sys

import numpy as np

cutoff = int(sys.argv[1])
coverage, metadata, out = sys.argv[2:5]


def strain_of(name):
    h = len(name) // 2
    return name[:h] if len(name) > 2 and name[h] == "_" and name[:h] == name[h + 1:] else name


cov = {}
for l in open(coverage):
    f = l.rstrip("\n").split("\t")
    if f[0] != "strain":
        cov[f[0]] = float(f[1])
meta = {}
for row in open(metadata):
    f = row.rstrip("\r\n").split("\t")
    if len(f) > 7 and f[3] != "strain":
        meta.setdefault(f[3], (f[6], f[7]))

names = [strain_of(s) for s in sys.stdin.readline().rstrip("\n").rstrip("\t").split("\t")]
ns = len(names)
codes = {}


def enc(g):
    v = codes.get(g)
    if v is None:
        a = g.replace("|", "/").split("/")
        v = -1 if "." in a else codes.setdefault("/".join(sorted(a)), len(codes))
        codes[g] = v
    return v


cols = []
for line in sys.stdin:
    cols.append(np.fromiter((enc(g) for g in line.rstrip("\n").rstrip("\t").split("\t")),
                            dtype=np.int16, count=ns))
G = np.vstack(cols).T  # strains x sites
del cols
M = (G >= 0).astype(np.float32)
N = M @ M.T
same = np.zeros((ns, ns), np.float64)
for k in sorted(set(v for v in codes.values() if v >= 0)):
    I = (G == k).astype(np.float32)
    if I.any():
        same += I @ I.T
D = np.rint(N - same).astype(np.int64)
N = np.rint(N).astype(np.int64)
miss = 1 - M.mean(axis=1)

with gzip.open(f"{out}.pairwise.tsv.gz", "wt") as o:
    o.write("strain_a\tstrain_b\tdiffs\tcompared\n")
    for i in range(ns):
        for j in range(i + 1, ns):
            o.write(f"{names[i]}\t{names[j]}\t{D[i, j]}\t{N[i, j]}\n")

par = list(range(ns))


def find(x):
    while par[x] != x:
        par[x] = par[par[x]]
        x = par[x]
    return x


for i in range(ns):
    for j in range(i + 1, ns):
        if D[i, j] <= cutoff:
            par[find(i)] = find(j)
groups = {}
for i in range(ns):
    groups.setdefault(find(i), []).append(i)
multi = sorted((sorted(g, key=lambda i: names[i]) for g in groups.values() if len(g) > 1),
               key=lambda g: (-len(g), names[g[0]]))

reps = []
with open(f"{out}.groups_le{cutoff}.tsv", "w") as o:
    o.write("near_identical_group\tsize\tstrain\trepresentative\tmissing_frac\tmosdepth_mean\t"
            "Origin\tEnvironment\tmax_diffs_within\tmin_diffs_to_outside\n")
    for k, g in enumerate(multi, 1):
        rep = min(g, key=lambda i: (round(miss[i], 4), -cov.get(names[i], 0), names[i]))
        reps.append(names[rep])
        inside = set(g)
        mx = max(D[a, b] for a in g for b in g if a != b)
        mn = min((D[a, b] for a in g for b in range(ns) if b not in inside), default=-1)
        for i in g:
            o.write(f"NI{k:03d}\t{len(g)}\t{names[i]}\t{int(i == rep)}\t{miss[i]:.4f}\t"
                    f"{cov.get(names[i], 'NA')}\t{meta.get(names[i], ('NA', 'NA'))[0]}\t"
                    f"{meta.get(names[i], ('NA', 'NA'))[1]}\t{mx}\t{mn}\n")
grouped = {i for g in multi for i in g}
reps += [names[i] for i in range(ns) if i not in grouped]
with open(f"{out}.representatives_le{cutoff}.txt", "w") as o:
    o.write("".join(s + "\n" for s in sorted(reps)))

off = D[np.triu_indices(ns, 1)]
sys.stderr.write(f"{ns} strains, {G.shape[1]} sites; pairs compared min {N[np.triu_indices(ns, 1)].min()}\n")
for c in (0, 2, 5, 10, 20, 50):
    sys.stderr.write(f"pairs <= {c}: {(off <= c).sum()}\n")
sys.stderr.write(f"<= {cutoff}: {len(multi)} groups holding {len(grouped)} strains; "
                 f"{len(reps)} representatives (of {ns})\n")
