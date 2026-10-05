#!/usr/bin/env python3
"""Divergence plots and per-strain placement table.

Inputs (run from the project root):
  results/variant_qc/divergence/all_5contigs.dist.tsv.gz   pairwise distances
  results/variant_qc/strain_qc.tsv                         species, QC decision
  results/ploidy_inference_all.csv        callable_bp, SNP sites
  population_sets.yaml                                     rmuc_core, hybrid_diploids
  ExRhodotorula_Phenotypes/strains.csv                     phenotype-table species
Outputs:
  results/variant_qc/divergence/divergence.{pdf,png}
  results/variant_qc/divergence/strain_placement.tsv
Usage: plot_divergence.py [CLONE_DISTANCE]
"""
import csv, gzip, sys
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

D_DIR = "results/variant_qc/divergence"
CANDIDATES = ["EXF_10854", "EXF_17335", "TFCN_17-333M-1", "TFCN_186CL-2"]
NO_META_HAP = ["TFCN_1A-1-5", "TFCN_25-332M-2", "TFCN_2M-1-3", "TFCN_86C-3", "TFCN_BY120-C1", "TFCN_BY120-C7"]
EXTRA = ["NRRL_Y-2510"]
clone_d = float(sys.argv[1]) if len(sys.argv) > 1 else None

with gzip.open(f"{D_DIR}/all_5contigs.dist.tsv.gz", "rt") as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
    D = np.array([[float(x) for x in l.rstrip("\n").split("\t")[1:]] for l in fh])
qc = {r["vcf_sample"]: r for r in csv.DictReader(open("results/variant_qc/strain_qc.tsv"), delimiter="\t")}
strain = [qc[s]["strain"] for s in samples]
idx = {s: i for i, s in enumerate(strain)}
v2 = {r["strain"]: r for r in csv.DictReader(open("results/ploidy_inference_all.csv"))}
pops = yaml.safe_load(open("population_sets.yaml"))["Populations"]
core, hyb = set(pops["rmuc_core"]), set(pops["hybrid_diploids"])
ph = {r["STRAIN"]: r["ASSIGNEDSPECIES"] for r in
      csv.DictReader(open("/bigdata/stajichlab/shared/projects/Rhodotorula/ExRhodotorula_Phenotypes/strains.csv"))}


def group(s):
    q = qc[samples[idx[s]]]
    if s in CANDIDATES: return "candidate (no metadata)"
    if s in NO_META_HAP: return "no-metadata haploid (6)"
    if s in EXTRA: return "NRRL_Y-2510"
    if s in core: return "rmuc_core"
    if s in hyb: return "hybrid diploid"
    if q["species"] == "R. frigidialcoholis": return "R. frigidialcoholis"
    if q["species"] == "R. aff. mucilaginosa": return "R. aff. mucilaginosa"
    return "other dropped"


grp = [group(s) for s in strain]
COL = {"rmuc_core": "#4C72B0", "hybrid diploid": "#DD8452", "R. frigidialcoholis": "#55A868",
       "R. aff. mucilaginosa": "#8172B3", "candidate (no metadata)": "#C44E52",
       "no-metadata haploid (6)": "#DA8BC3", "NRRL_Y-2510": "#000000", "other dropped": "#BBBBBB"}
ORDER = list(COL)

# divergence from reference: hom-alt-ish SNPs per callable Mb (v2 counts, forced diploid
# call: n_snp_sites includes unbalanced hets, so this is an upper bound for diploids)
div = np.array([(int(v2[s]["n_snp_sites"]) - int(v2[s]["n_het"])) / int(v2[s]["callable_bp"]) * 1e6
                if s in v2 and int(v2[s]["callable_bp"]) > 0 else np.nan for s in strain])

core_hap = [i for i, s in enumerate(strain) if grp[i] == "rmuc_core" and qc[samples[i]]["called_ploidy"] == "1"]
frig = [i for i, g in enumerate(grp) if g == "R. frigidialcoholis" and qc[samples[i]]["called_ploidy"] == "1"]
aff = [i for i, g in enumerate(grp) if g == "R. aff. mucilaginosa" and qc[samples[i]]["called_ploidy"] == "1"]


def nearest(i, pool):
    pool = [j for j in pool if j != i]
    j = min(pool, key=lambda j: D[i, j])
    return strain[j], D[i, j]


with open(f"{D_DIR}/strain_placement.tsv", "w") as o:
    o.write("strain\tgroup\tqc_decision\tphenotype_table_species\tsnps_per_Mb_vs_ref\t"
            "median_dist_core_haploids\tnearest_core\tnearest_core_dist\t"
            "nearest_frig\tnearest_frig_dist\tnearest_affmuc\tnearest_affmuc_dist\n")
    for i, s in enumerate(strain):
        q = qc[samples[i]]
        nc = nearest(i, core_hap); nf = nearest(i, frig); na = nearest(i, aff)
        md = np.median([D[i, j] for j in core_hap if j != i])
        o.write(f"{s}\t{grp[i]}\t{q['decision']}:{q['reasons']}\t{ph.get(s, 'NA')}\t{div[i]:.1f}\t{md:.4f}\t"
                f"{nc[0]}\t{nc[1]:.4f}\t{nf[0]}\t{nf[1]:.4f}\t{na[0]}\t{na[1]:.4f}\n")

fig, ax = plt.subplots(2, 2, figsize=(15, 11))
# A: divergence from reference
a = ax[0, 0]
bins = np.logspace(0, np.log10(np.nanmax(div)) + 0.1, 45)
stack = [div[[i for i, g in enumerate(grp) if g == k]] for k in ORDER]
a.hist([x[np.isfinite(x)] for x in stack], bins=bins, stacked=True, color=[COL[k] for k in ORDER],
       label=[f"{k} ({len(x)})" for k, x in zip(ORDER, stack)])
a.set_xscale("log")
for s in CANDIDATES + EXTRA:
    a.annotate(s, (div[idx[s]], 0), xytext=(0, 40 + 12 * (CANDIDATES + EXTRA).index(s)), textcoords="offset points",
               fontsize=7, rotation=0, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5))
a.set_xlabel("SNPs vs DH4148 reference per callable Mb (log)")
a.set_ylabel("strains")
a.set_title("A. Divergence from the reference, all 319 strains", fontsize=10)
a.legend(fontsize=7)
# B: pairwise distances
b = ax[0, 1]
cc = [D[i, j] for k, i in enumerate(core_hap) for j in core_hap[k + 1:]]
cf = [D[i, j] for i in core_hap for j in frig]
caf = [D[i, j] for i in core_hap for j in aff]
fa = [D[i, j] for i in frig for j in aff]
bb = np.linspace(0, max(max(cc), max(cf), max(caf), max(fa)) * 1.02, 80)
for vals, lab, c in [(cc, "core haploid vs core haploid", COL["rmuc_core"]),
                     (cf, "core vs R. frigidialcoholis", COL["R. frigidialcoholis"]),
                     (caf, "core vs R. aff. mucilaginosa", COL["R. aff. mucilaginosa"]),
                     (fa, "R. frigidialcoholis vs R. aff. mucilaginosa", "#999999")]:
    b.hist(vals, bins=bb, alpha=0.6, color=c, label=f"{lab} (n={len(vals)})", density=True)
for s in CANDIDATES + EXTRA:
    b.axvline(np.median([D[idx[s], j] for j in core_hap]), color=COL["candidate (no metadata)"] if s in CANDIDATES else "k", lw=0.8, ls="--")
b.set_xlabel("pairwise distance (share of variant sites that differ)")
b.set_ylabel("density")
b.set_title("B. Pairwise distances; dashed = median distance of each candidate / NRRL_Y-2510 to core", fontsize=10)
b.legend(fontsize=7)
# C: placement
c = ax[1, 0]
xm = np.array([np.median([D[i, j] for j in core_hap if j != i]) for i in range(len(strain))])
ym = np.array([min(D[i, j] for j in frig + aff if j != i) for i in range(len(strain))])
for k in ORDER:
    ii = [i for i, g in enumerate(grp) if g == k]
    c.scatter(xm[ii], ym[ii], s=18, color=COL[k], label=k, alpha=0.85, edgecolors="none")
for s in CANDIDATES + EXTRA:
    i = idx[s]; c.annotate(s, (xm[i], ym[i]), fontsize=7, xytext=(4, 4), textcoords="offset points")
c.set_xlabel("median distance to rmuc_core haploids")
c.set_ylabel("distance to nearest R. frigidialcoholis / R. aff. mucilaginosa haploid")
c.set_title("C. Placement of each strain", fontsize=10)
c.legend(fontsize=7)
# D: nearest neighbour within core haploids
d = ax[1, 1]
nn = np.array([min(D[i, j] for j in core_hap if j != i) for i in core_hap])
d.hist(np.log10(np.maximum(nn, 1e-5)), bins=60, color=COL["rmuc_core"])
if clone_d:
    d.axvline(np.log10(clone_d), color="r", ls="--", lw=1, label=f"clone threshold {clone_d}")
    d.legend(fontsize=8)
d.set_xlabel("log10 distance to nearest other core haploid (floor 1e-5)")
d.set_ylabel("core haploid strains")
d.set_title(f"D. Nearest-neighbour distance within rmuc_core haploids (n={len(core_hap)})", fontsize=10)
fig.suptitle("Divergence of DH4148-mapped strains (5 largest contigs, PASS unmasked biallelic SNPs, every 4th site)")
fig.tight_layout()
fig.savefig(f"{D_DIR}/divergence.png", dpi=130)
fig.savefig(f"{D_DIR}/divergence.pdf")
print("written", D_DIR)
