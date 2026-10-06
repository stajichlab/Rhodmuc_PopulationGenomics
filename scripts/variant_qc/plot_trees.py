#!/usr/bin/env python3
"""Draw IQ-TREE strain trees as PNG and PDF next to the treefile.

Midpoint-rooted phylogram, ladderized. Tip dot colour = strain collection
(TFCN, EXF, DBVPG, CHIFNET, other); the label gives strain, origin, species if
not R. mucilaginosa, "hybrid" for hybrid_diploids members, and the near-identical
group (NIxxx, <= 5 SNPs; declone.sh) when a group file is given. Internal nodes with
SH-aLRT >= 80 and UFBoot >= 95 get a small dark dot.

Usage: plot_trees.py metadata.txt population_sets.yaml [clusters.tsv] -- tree1.treefile [tree2 ...]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

args = sys.argv[1:]
sep = args.index("--")
opts, trees = args[:sep], args[sep + 1:]
metadata, yaml_path = opts[0], opts[1]
clusters = opts[2] if len(opts) > 2 else None

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
COLL = [("TFCN", "#2a78d6"), ("EXF", "#eb6834"), ("DBVPG", "#1baf7a"),
        ("CHIFNET", "#eda100"), ("other", "#4a3aa7")]
COLOR = dict(COLL)


def strain_of(name):
    h = len(name) // 2
    return name[:h] if len(name) > 2 and name[h] == "_" and name[:h] == name[h + 1:] else name


meta = {}
for row in open(metadata):
    f = row.rstrip("\r\n").split("\t")
    if len(f) > 7 and f[3] != "strain":
        meta.setdefault(f[3], (f[2], f[6]))
hybrids, cur = set(), None
for line in open(yaml_path):
    s = line.strip()
    if s.endswith(":") and not s.startswith("-") and not s.startswith("#"):
        cur = s[:-1]
    elif s.startswith("- ") and cur == "hybrid_diploids":
        hybrids.add(s[2:].strip())
cl = {}
if clusters:
    for r in open(clusters):
        f = r.rstrip("\n").split("\t")
        if f[0] not in ("cluster", "near_identical_group"):
            cl[f[2]] = f[0]


def label(tip):
    s = strain_of(tip)
    if s == "reference":
        return "reference DH4148", "other"
    sp, orig = meta.get(s, ("no metadata", ""))
    parts = [s]
    if orig and orig != "NA":
        parts.append(orig)
    if sp != "R. mucilaginosa":
        parts.append(sp)
    if s in hybrids:
        parts.append("hybrid")
    if s in cl:
        parts.append(cl[s])
    coll = s.split("_")[0] if s.split("_")[0] in COLOR else "other"
    return "  ·  ".join(parts), coll


def support_ok(conf_name):
    if not conf_name:
        return False
    try:
        sh, uf = conf_name.split("/")
        return float(sh) >= 80 and float(uf) >= 95
    except ValueError:
        return False


for path in trees:
    t = Phylo.read(path, "newick")
    # IQ-TREE writes "SH/UFBoot" as the internal node name; keep it before rooting moves things
    for c in t.get_nonterminals():
        if c.confidence is not None and c.name is None:
            c.name = str(c.confidence)
    t.root_at_midpoint()
    t.ladderize()
    tips = t.get_terminals()
    n = len(tips)
    depth = t.depths()
    if max(depth.values()) == 0:
        depth = t.depths(unit_branch_lengths=True)
    y = {tip: i for i, tip in enumerate(tips)}

    def ypos(c):
        if c in y:
            return y[c]
        y[c] = sum(ypos(k) for k in c.clades) / len(c.clades)
        return y[c]
    ypos(t.root)

    height = max(6, 0.13 * n + 1.5)
    fig, ax = plt.subplots(figsize=(14, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    xmax = max(depth.values())
    for c in t.find_clades():
        if c.clades:
            ys = [y[k] for k in c.clades]
            ax.plot([depth.get(c, 0)] * 2, [min(ys), max(ys)], color=INK2, lw=0.8)
            for k in c.clades:
                ax.plot([depth.get(c, 0), depth[k]], [y[k]] * 2, color=INK2, lw=0.8)
            if c is not t.root and support_ok(c.name):
                ax.plot(depth[c], y[c], "o", ms=2.5, color=INK, zorder=3)
    used = []
    for tip in tips:
        text, coll = label(tip.name)
        used.append(coll)
        ax.plot(depth[tip], y[tip], "o", ms=5, color=COLOR[coll], mec=SURFACE, mew=0.6, zorder=4)
        ax.text(depth[tip] + xmax * 0.008, y[tip], text, va="center", fontsize=6.5, color=INK)
    ax.set_ylim(n, -1)
    ax.set_xlim(-xmax * 0.01, xmax * 1.55)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK2)
    ax.tick_params(axis="x", colors=INK2, labelsize=8)
    ax.set_xlabel("substitutions per SNP site (midpoint-rooted)", color=INK2, fontsize=9)
    name = path.split("/")[-1].replace(".snps.maf.treefile", "")
    ax.set_title(f"{name}: {n} tips  ·  dot on node = SH-aLRT ≥ 80 and UFBoot ≥ 95",
                 loc="left", fontsize=11, color=INK)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=6, color=col, label=k)
               for k, col in COLL if k in used]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8, title="collection",
              title_fontsize=8, labelcolor=INK)
    out = path[: -len(".treefile")]
    fig.tight_layout()
    fig.savefig(out + ".tree.pdf")
    fig.savefig(out + ".tree.png", dpi=150)
    plt.close(fig)
    print(out + ".tree.{pdf,png}", n, "tips")
