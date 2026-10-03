#!/usr/bin/env python3
"""Per-sample genotype QC and multiallelic-site diagnostics.

Reads `bcftools query` output on stdin with this format:
  %CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO/DP[\t%GT:%AD:%DP:%GQ]\n
The first line must be the sample list (tab separated), written by the caller.

Writes:
  <prefix>.samples.tsv     per-sample counts
  <prefix>.multiallelic.tsv  summary table of multiallelic-site classes
  <prefix>.site_dp.tsv     INFO/DP histogram (PASS sites), 100-read bins
"""
import sys
from collections import Counter, defaultdict

prefix = sys.argv[1]
mask = None
if len(sys.argv) > 2:
    # optional merged BED: in-memory interval lookup per contig
    import bisect
    mask = defaultdict(list)
    for line in open(sys.argv[2]):
        c, s, e = line.split()[:3]
        mask[c].append((int(s), int(e)))
    starts = {c: [s for s, _ in v] for c, v in mask.items()}

    def in_mask(c, pos0):
        iv = mask.get(c)
        if not iv:
            return False
        i = bisect.bisect_right(starts[c], pos0) - 1
        return i >= 0 and iv[i][0] <= pos0 < iv[i][1]
else:
    def in_mask(c, pos0):
        return False

samples = sys.stdin.readline().rstrip("\n").split("\t")
ns = len(samples)
ploidy = [0] * ns
keys = ["called", "missing", "nonref", "het", "gq_lt20", "dp_lt5",
        "hap_alt_mixedAD", "hap_alt_called", "dip_het_ABext", "het_multi_alt",
        "uses_allele_ge2"]
st = [Counter() for _ in range(ns)]
multi = Counter()
dp_hist = Counter()
n_sites = Counter()

for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    chrom, pos, ref, alt, filt, qual, sdp = f[:7]
    alts = alt.split(",")
    nalt = len(alts)
    is_pass = filt == "PASS"
    masked = in_mask(chrom, int(pos) - 1)
    n_sites["all"] += 1
    n_sites["pass"] += is_pass
    n_sites["masked"] += masked
    if is_pass and sdp != ".":
        dp_hist[int(sdp) // 100 * 100] += 1
    allele_ac = Counter()
    carriers_hap_ge2 = 0
    carriers_dip_12 = 0
    for i, g in enumerate(f[7:]):
        gt, ad, dp, gq = g.split(":")
        sep = "/" if "/" in gt else ("|" if "|" in gt else None)
        al = gt.split(sep) if sep else [gt]
        if not ploidy[i]:
            ploidy[i] = len(al)
        s = st[i]
        if "." in al:
            s["missing"] += 1
            continue
        s["called"] += 1
        if gq != "." and int(gq) < 20:
            s["gq_lt20"] += 1
        if dp != "." and int(dp) < 5:
            s["dp_lt5"] += 1
        ai = [int(a) for a in al]
        for a in ai:
            allele_ac[a] += 1
        if any(a > 0 for a in ai):
            s["nonref"] += 1
        if max(ai) >= 2:
            s["uses_allele_ge2"] += 1
        ads = [int(x) for x in ad.split(",")] if ad != "." else None
        tot = sum(ads) if ads else 0
        if len(ai) == 1:
            if ai[0] > 0 and tot > 0:
                s["hap_alt_called"] += 1
                frac = ads[ai[0]] / tot
                if frac < 0.8:  # >20% of reads support some other allele
                    s["hap_alt_mixedAD"] += 1
            if ai[0] >= 2:
                carriers_hap_ge2 += 1
        else:
            if ai[0] != ai[1]:
                s["het"] += 1
                if 0 not in ai:
                    s["het_multi_alt"] += 1
                    carriers_dip_12 += 1
                if tot > 0:
                    ab = ads[ai[1]] / tot
                    if ab < 0.2 or ab > 0.8:
                        s["dip_het_ABext"] += 1
    if nalt > 1:
        has_star = "*" in alts
        kind = "snp" if all(len(a) == 1 for a in alts if a != "*") and len(ref) == 1 else "indel/mixed"
        ac_extra = sorted((allele_ac[k] for k in range(2, nalt + 1)))
        min_extra = ac_extra[0] if ac_extra else 0
        rare = "minor_alt_AC<=2" if min_extra <= 2 else "minor_alt_AC>2"
        multi[(filt == "PASS", kind, "star" if has_star else "nostar", rare,
               "masked" if masked else "unmasked")] += 1
        multi[("any_haploid_carries_allele>=2", carriers_hap_ge2 > 0)] += 1
        multi[("any_diploid_het_without_ref", carriers_dip_12 > 0)] += 1
        n_sites["multiallelic"] += 1

with open(prefix + ".samples.tsv", "w") as o:
    o.write("sample\tploidy\t" + "\t".join(keys) + "\n")
    for i, name in enumerate(samples):
        o.write(f"{name}\t{ploidy[i]}\t" + "\t".join(str(st[i][k]) for k in keys) + "\n")
with open(prefix + ".multiallelic.tsv", "w") as o:
    o.write("# " + " ".join(f"{k}={v}" for k, v in n_sites.items()) + "\n")
    for k, v in sorted(multi.items(), key=lambda kv: -kv[1]):
        o.write("\t".join(map(str, k)) + f"\t{v}\n")
with open(prefix + ".site_dp.tsv", "w") as o:
    for k in sorted(dp_hist):
        o.write(f"{k}\t{dp_hist[k]}\n")
