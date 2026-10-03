#!/usr/bin/env python3
"""PCA of a biallelic SNP set with mixed ploidy.

Input: `bcftools +dosage -- -t GT` output (.gz or plain). Dosage per sample is
divided by the sample ploidy (strain_qc.tsv called_ploidy), so haploid and
diploid strains share one 0..1 alt-allele-fraction scale. Missing (-1) is
mean-imputed per site. Sites are centred and scaled by sqrt(p(1-p))
(Patterson et al. 2006). PCs come from the sample x sample Gram matrix.

Usage: pca_snps.py dosage.tsv.gz strain_qc.tsv metadata.txt groups.tsv outprefix [exclude.txt]
groups.tsv: vcf_sample<TAB>group label used for colour.
"""
import csv, gzip, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dos_f, qc_f, meta_f, grp_f, out = sys.argv[1:6]
excl = set(open(sys.argv[6]).read().split()) if len(sys.argv) > 6 else set()

qc = {r["vcf_sample"]: r for r in csv.DictReader(open(qc_f), delimiter="\t")}
meta = {r["strain"]: r for r in csv.DictReader(open(meta_f), delimiter="\t")}
grp = dict(l.rstrip("\n").split("\t") for l in open(grp_f))

op = gzip.open if dos_f.endswith(".gz") else open
with op(dos_f, "rt") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    names = [h.split("]", 1)[1] for h in head[4:]]
    keep = [i for i, n in enumerate(names) if n not in excl]
    names = [names[i] for i in keep]
    ploidy = np.array([float(qc[n]["called_ploidy"]) for n in names])
    rows = []
    for line in fh:
        v = np.array(line.rstrip("\n").split("\t")[4:], dtype=np.float32)[keep]
        rows.append(v)
X = np.vstack(rows)                      # sites x samples
del rows
miss = X < 0
X = X / ploidy[None, :]
X[miss] = np.nan
p = np.nanmean(X, axis=1)
ok = np.isfinite(p) & (p > 0) & (p < 1)  # polymorphic within this sample set
X, p, miss = X[ok], p[ok], miss[ok]
X = np.where(np.isnan(X), p[:, None], X)
X = (X - p[:, None]) / np.sqrt(p * (1 - p))[:, None]
G = X.T @ X / X.shape[0]
w, V = np.linalg.eigh(G)
order = np.argsort(w)[::-1]
w, V = w[order], V[:, order]
var = w / w.sum()
npc = 10
pcs = V[:, :npc] * np.sqrt(w[:npc])

with open(out + ".eigenvec.tsv", "w") as o:
    o.write("sample\tstrain\tgroup\tploidy\tOrigin\tEnvironment\t" + "\t".join(f"PC{i+1}" for i in range(npc)) + "\n")
    for j, n in enumerate(names):
        s = qc[n]["strain"]; m = meta.get(s, {})
        o.write(f"{n}\t{s}\t{grp.get(n,'NA')}\t{int(ploidy[j])}\t{m.get('Origin','')}\t{m.get('Environment','')}\t"
                + "\t".join(f"{x:.5f}" for x in pcs[j]) + "\n")
with open(out + ".eigenval.tsv", "w") as o:
    o.write(f"# samples={len(names)} sites={X.shape[0]}\nPC\teigenvalue\tvar_fraction\n")
    for i in range(npc):
        o.write(f"PC{i+1}\t{w[i]:.5f}\t{var[i]:.4f}\n")

def scatter(ax, a, b, labels, title):
    cats = sorted(set(labels), key=lambda c: -labels.count(c))
    cmap = plt.get_cmap("tab20")
    for k, c in enumerate(cats):
        idx = [j for j, l in enumerate(labels) if l == c]
        ax.scatter(pcs[idx, a], pcs[idx, b], s=14, color=cmap(k % 20), label=f"{c} ({len(idx)})",
                   alpha=0.8, edgecolors="none")
    ax.set_xlabel(f"PC{a+1} ({var[a]*100:.1f}%)"); ax.set_ylabel(f"PC{b+1} ({var[b]*100:.1f}%)")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, loc="best", markerscale=1.2)

def top(vals, n=12):
    from collections import Counter
    c = Counter(vals); t = {k for k, _ in c.most_common(n)}
    return [v if v in t else "other" for v in vals]

lab_g = [grp.get(n, "NA") for n in names]
lab_env = top([meta.get(qc[n]["strain"], {}).get("Environment", "") or "NA" for n in names])
lab_ori = top([meta.get(qc[n]["strain"], {}).get("Origin", "") or "NA" for n in names])
fig, axs = plt.subplots(2, 2, figsize=(14, 12))
scatter(axs[0, 0], 0, 1, lab_g, "group")
scatter(axs[0, 1], 0, 2, lab_g, "group")
scatter(axs[1, 0], 0, 1, lab_env, "Environment (top 12)")
scatter(axs[1, 1], 0, 1, lab_ori, "Origin (top 12)")
fig.suptitle(f"{out.split('/')[-1]}: {len(names)} strains, {X.shape[0]} SNPs")
fig.tight_layout()
fig.savefig(out + ".png", dpi=130)
print(f"samples={len(names)} sites={X.shape[0]} var(PC1..5)=" + ",".join(f"{v:.3f}" for v in var[:5]))
