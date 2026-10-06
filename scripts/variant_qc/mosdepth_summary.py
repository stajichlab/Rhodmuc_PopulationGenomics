#!/usr/bin/env python3
"""One row per strain from mosdepth_all.sh output: mean depth and genome fractions at >= 5x / 10x.

Usage: mosdepth_summary.py results/variant_qc/coverage > results/variant_qc/coverage.tsv
"""
import glob
import os
import sys

d = sys.argv[1]
print("strain\tmean_depth\tfrac_ge5x\tfrac_ge10x")
for summ in sorted(glob.glob(os.path.join(d, "*.mosdepth.summary.txt"))):
    s = os.path.basename(summ)[: -len(".mosdepth.summary.txt")]
    mean = next(l.split("\t")[3] for l in open(summ) if l.startswith("total\t"))
    frac = {}
    for l in open(os.path.join(d, f"{s}.mosdepth.global.dist.txt")):
        f = l.split("\t")
        if f[0] == "total" and f[1] in ("5", "10"):
            frac[f[1]] = f[2].strip()
    print(f"{s}\t{mean}\t{frac['5']}\t{frac['10']}")
