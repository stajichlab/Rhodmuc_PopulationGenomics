#!/usr/bin/env python3
"""Clone groups: single-linkage clusters at pairwise distance < CUTOFF among one population.

Groups are numbered CG001.. by size (largest first), ties by first strain name.
Usage: clone_groups.py dist.tsv.gz population_sets.yaml POP metadata.txt CUTOFF > clone_groups.tsv
dist.tsv.gz: square matrix from pairwise_distance.py (row/col names = VCF sample names).
"""
import gzip
import sys

dist_gz, yaml_path, pop, metadata, cutoff = sys.argv[1:6]
cutoff = float(cutoff)


def strain_of(name):
    h = len(name) // 2
    return name[:h] if name[h] == "_" and name[:h] == name[h + 1:] else name


members, cur = [], None
for line in open(yaml_path):
    s = line.strip()
    if s.endswith(":") and not s.startswith("-"):
        cur = s[:-1]
    elif s.startswith("- ") and cur == pop:
        members.append(s[2:].strip())
members = set(members)

meta = {}
for row in open(metadata):
    f = row.rstrip("\r\n").split("\t")
    if len(f) > 7 and f[3] != "strain":
        meta.setdefault(f[3], (f[6], f[7]))

with gzip.open(dist_gz, "rt") as fh:
    cols = [strain_of(c) for c in fh.readline().rstrip("\n").split("\t")[1:]]
    keep = [i for i, c in enumerate(cols) if c in members]
    parent = {cols[i]: cols[i] for i in keep}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for line in fh:
        f = line.rstrip("\n").split("\t")
        a = strain_of(f[0])
        if a not in parent:
            continue
        v = f[1:]
        for i in keep:
            b = cols[i]
            if b != a and v[i] not in ("", "NA", "nan") and float(v[i]) < cutoff:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb

missing = members - set(parent)
if missing:
    sys.stderr.write(f"not in distance matrix: {sorted(missing)}\n")
groups = {}
for s in parent:
    groups.setdefault(find(s), []).append(s)
ordered = sorted((sorted(g) for g in groups.values()), key=lambda g: (-len(g), g[0]))
print("clone_group\tsize\tstrain\tOrigin\tEnvironment")
for k, g in enumerate(ordered, 1):
    for s in g:
        o, e = meta.get(s, ("NA", "NA"))
        print(f"CG{k:03d}\t{len(g)}\t{s}\t{o}\t{e}")
sys.stderr.write(f"{len(parent)} strains, {len(ordered)} groups, sizes {[len(g) for g in ordered[:8]]}, "
                 f"singletons {sum(len(g) == 1 for g in ordered)}\n")
