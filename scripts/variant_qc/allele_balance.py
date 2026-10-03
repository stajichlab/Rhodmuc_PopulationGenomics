#!/usr/bin/env python3
"""Per-sample read allele balance at biallelic SNPs (ploidy evidence).

stdin: header line of sample names (tab separated), then one line per site:
  [%AD\t]   (ref,alt read counts per sample)
For each sample, at sites with total depth >= MIN_DP and minor-allele reads >= 2:
  minor fraction f = min(ref,alt)/total.
Reports, per sample:
  covered       sites with depth >= MIN_DP
  mixed         sites with f >= 0.10
  balanced      sites with 0.25 <= f <= 0.5  (expected for diploid hets)
  mixed_per_kcov, balanced_per_kcov  per 1000 covered sites
  peak_bin      modal bin of f among mixed sites (bins 0.10-0.15 ... 0.45-0.50)
  frac_f_045_050, frac_f_030_037  share of mixed sites with f in the
                diploid (0.45-0.50) or triploid (0.30-0.37) window
"""
import sys
MIN_DP = 10
names = sys.stdin.readline().rstrip("\n").rstrip("\t").split("\t")
n = len(names)
cov = [0]*n; mixed = [0]*n; bal = [0]*n
bins = [[0]*8 for _ in range(n)]
for line in sys.stdin:
    f = line.rstrip("\n").rstrip("\t").split("\t")
    for i in range(n):
        ad = f[i]
        if ad == "." or "," not in ad:
            continue
        a, b = ad.split(",")[:2]
        if a == "." or b == ".":
            continue
        a = int(a); b = int(b); t = a + b
        if t < MIN_DP:
            continue
        cov[i] += 1
        m = a if a < b else b
        if m < 2:
            continue
        fr = m / t
        if fr >= 0.10:
            mixed[i] += 1
            bins[i][min(7, int((fr - 0.10) / 0.05))] += 1
            if fr >= 0.25:
                bal[i] += 1
print("sample\tcovered\tmixed\tbalanced\tmixed_per_kcov\tbalanced_per_kcov\tpeak_bin\tfrac_f_045_050\tfrac_f_030_037\t"
      + "\t".join(f"bin_{0.10+0.05*k:.2f}" for k in range(8)))
for i in range(n):
    c = cov[i] or 1; mx = mixed[i] or 1
    pk = max(range(8), key=lambda k: bins[i][k])
    # 0.30-0.37 ~ bins 0.30-0.35 plus 0.35-0.40 (approximate window)
    tri = (bins[i][4] + bins[i][5]) / mx
    print(f"{names[i]}\t{cov[i]}\t{mixed[i]}\t{bal[i]}\t{1000*mixed[i]/c:.2f}\t{1000*bal[i]/c:.2f}\t"
          f"{0.10+0.05*pk:.2f}\t{bins[i][7]/mx:.3f}\t{tri:.3f}\t" + "\t".join(str(x) for x in bins[i]))
