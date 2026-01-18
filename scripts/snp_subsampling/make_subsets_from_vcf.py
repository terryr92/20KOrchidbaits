#!/usr/bin/env python3
import gzip
import random
from pathlib import Path
import sys

vcf_gz = Path(sys.argv[1])  # Diagnostic VCF (subset positions already)
outdir = Path(sys.argv[2])
reps = int(sys.argv[3]) if len(sys.argv) > 3 else 100
fractions = [1.0, 0.75, 0.50, 0.25, 0.10]
seed = 123

outdir.mkdir(parents=True, exist_ok=True)
random.seed(seed)

header, records = [], []
with gzip.open(vcf_gz, "rt") as f:
    for line in f:
        if line.startswith("#"):
            header.append(line)
        else:
            records.append(line)

n = len(records)
for frac in fractions:
    k = max(1, int(round(n * frac)))
    for r in range(reps):
        sub = random.sample(records, k)
        tag = f"{int(frac*100):03d}_rep{r:03d}"
        out = outdir / f"subset_{tag}.vcf"
        out.write_text("".join(header + sub))

print(f"DONE: subsets written to {outdir} (n_sites={n})")
