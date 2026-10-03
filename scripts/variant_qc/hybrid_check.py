#!/usr/bin/env python3
"""For each high-het diploid: at its het sites (biallelic SNP, PASS), what
fraction of ALT alleles are the majority allele in each reference group?
stdin: bcftools query -f '[%GT\t]\n' with header line of sample names.
argv: groups dir. Writes TSV to stdout."""
import sys
gd = sys.argv[1]
grp = {g: set(open(f"{gd}/{g}.txt").read().split()) for g in
       ["kept_haploid", "kept_dip_highhet", "frig_haploid", "affmuc_haploid"]}
names = sys.stdin.readline().rstrip("\n").rstrip("\t").split("\t")
idx = {g: [i for i, n in enumerate(names) if n in s] for g, s in grp.items()}
dips = idx["kept_dip_highhet"]
from collections import Counter
c = {i: Counter() for i in dips}
def altfreq(gts, ii):
    a = n = 0
    for i in ii:
        for x in gts[i].replace("|", "/").split("/"):
            if x in ("0", "1"):
                n += 1; a += x == "1"
    return a / n if n else None
for line in sys.stdin:
    g = line.rstrip("\n").rstrip("\t").split("\t")
    f = {k: altfreq(g, idx[k]) for k in ("kept_haploid", "frig_haploid", "affmuc_haploid")}
    if None in f.values():
        continue
    for i in dips:
        if g[i].replace("|", "/") not in ("0/1", "1/0"):
            continue
        cc = c[i]; cc["het"] += 1
        cc["muc_hap_alt>0.5"] += f["kept_haploid"] > 0.5
        cc["muc_hap_alt=0"] += f["kept_haploid"] == 0
        cc["frig_alt>0.5"] += f["frig_haploid"] > 0.5
        cc["affmuc_alt>0.5"] += f["affmuc_haploid"] > 0.5
        cc["muc_hap_alt=0&frig>0.5"] += f["kept_haploid"] == 0 and f["frig_haploid"] > 0.5
        cc["muc_hap_alt=0&affmuc>0.5"] += f["kept_haploid"] == 0 and f["affmuc_haploid"] > 0.5
keys = ["het", "muc_hap_alt>0.5", "muc_hap_alt=0", "frig_alt>0.5", "affmuc_alt>0.5",
        "muc_hap_alt=0&frig>0.5", "muc_hap_alt=0&affmuc>0.5"]
print("sample\t" + "\t".join(keys))
for i in dips:
    h = c[i]["het"] or 1
    print(names[i] + "\t" + str(c[i]["het"]) + "\t" +
          "\t".join(f"{c[i][k]/h:.3f}" for k in keys[1:]))
