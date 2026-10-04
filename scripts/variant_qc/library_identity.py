#!/usr/bin/env python3
"""Identity of each sequencing library (read group) against the population callset.

Inputs (from library_identity.sh):
  pop.tsv : header 'CHROM POS REF ALT <samples...>', then per site the
            population genotypes as ALT fraction (0, 0.5, 1, or NA)
  lib/*.tsv : per strain, 'CHROM POS REF ALT [AD per read group]' from
            bcftools mpileup -G (one pseudo-sample per read group); header
            line lists read-group names.
For each library, at sites with >= MIN_DP reads: x = ALT read fraction.
  mixed_frac = share of sites with 0.2 < x < 0.8 (a pure haploid library ~0)
  distance to every population sample = mean |x - g| over shared sites
Writes TSV to stdout: library, strain, n_sites, mixed_frac, own_dist,
best_match, best_dist, second_match, second_dist.
"""
import glob, os, sys
import numpy as np
MIN_DP = 5
pop_f, lib_dir = sys.argv[1], sys.argv[2]
with open(pop_f) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    samples = h[4:]
    keys, G = {}, []
    for i, l in enumerate(fh):
        f = l.rstrip("\n").split("\t")
        keys[(f[0], f[1])] = (i, f[3])
        G.append([np.nan if v == "NA" else float(v) for v in f[4:]])
G = np.array(G, dtype=np.float32)
strain_of = {s: s.rsplit("_", 1)[0] if s.count("_") >= 3 else s for s in samples}
print("library\tstrain\tn_sites\tmixed_frac\town_dist\tbest_match\tbest_dist\tsecond_match\tsecond_dist")
for fn in sorted(glob.glob(os.path.join(lib_dir, "*.tsv"))):
    strain = os.path.basename(fn)[:-4]
    with open(fn) as fh:
        libs = fh.readline().rstrip("\n").split("\t")[4:]
        X = np.full((len(libs), G.shape[0]), np.nan, dtype=np.float32)
        for l in fh:
            f = l.rstrip("\n").split("\t")
            k = keys.get((f[0], f[1]))
            if k is None:
                continue
            i, palt = k
            alts = f[3].split(",")
            if palt not in alts:
                ai = None
            else:
                ai = alts.index(palt) + 1
            for j, ad in enumerate(f[4:]):
                if ad in (".", ""):
                    continue
                c = [int(v) if v != "." else 0 for v in ad.split(",")]
                ref = c[0]; alt = c[ai] if ai is not None and ai < len(c) else 0
                if ref + alt >= MIN_DP:
                    X[j, i] = alt / (ref + alt)
    for j, lib in enumerate(libs):
        x = X[j]; ok = ~np.isnan(x)
        mixed = float(np.mean((x[ok] > 0.2) & (x[ok] < 0.8))) if ok.any() else float("nan")
        both = ok[:, None] & ~np.isnan(G)
        d = np.where(both, np.abs(np.nan_to_num(x)[:, None] - np.nan_to_num(G)), 0).sum(0) / np.maximum(both.sum(0), 1)
        d[both.sum(0) < 100] = np.nan
        order = np.argsort(np.where(np.isnan(d), 9, d))
        own = [d[i] for i, s in enumerate(samples) if s.startswith(strain + "_") or s == strain]
        print(f"{lib}\t{strain}\t{int(ok.sum())}\t{mixed:.3f}\t{own[0] if own else float('nan'):.4f}\t"
              f"{samples[order[0]]}\t{d[order[0]]:.4f}\t{samples[order[1]]}\t{d[order[1]]:.4f}")
