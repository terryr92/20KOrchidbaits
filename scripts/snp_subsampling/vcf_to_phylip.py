#!/usr/bin/env python3
"""
Convert a (filtered, biallelic) VCF to a simple PHYLIP alignment.
Encoding:
- 0/0 -> '0'
- 1/1 -> '1'
- missing -> '?'
This is a minimal representation suitable for consistent tree comparisons.
"""
import sys

vcf = sys.argv[1]
out_phy = sys.argv[2]

samples = []
seqs = None
n_sites = 0

with open(vcf, "r") as f:
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            parts = line.strip().split("\t")
            samples = parts[9:]
            seqs = {s: [] for s in samples}
            continue
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        fmt = parts[8].split(":")
        gt_i = fmt.index("GT") if "GT" in fmt else None
        if gt_i is None:
            raise RuntimeError("GT not found in FORMAT column")

        gts = parts[9:]
        for s, cell in zip(samples, gts):
            fields = cell.split(":")
            gt = fields[gt_i]
            if gt in ("0/0","0|0"):
                seqs[s].append("0")
            elif gt in ("1/1","1|1"):
                seqs[s].append("1")
            else:
                seqs[s].append("?")

        n_sites += 1

# Write PHYLIP
with open(out_phy, "w") as w:
    w.write(f"{len(samples)} {n_sites}\n")
    for s in samples:
        w.write(f"{s[:10].ljust(10)} {''.join(seqs[s])}\n")

print(f"DONE: PHYLIP -> {out_phy} (sites={n_sites})")
